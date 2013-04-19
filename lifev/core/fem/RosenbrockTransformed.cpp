//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Rosenbrock method with variable change

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 01-04-2013
 */



template<UInt s>
RosenbrockTransformed<s>::RosenbrockTransformed(Real g, const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
												const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat,
												UInt order, boost::shared_ptr<Epetra_Comm>& Comm, const string& solvParam)
:M_g(g), M_A(A), M_C(C), M_gammai(gammai), M_a(a), M_m(m), M_mhat(mhat), M_p(order)
{
	initMembers();
	setSolver(Comm, solvParam);
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::solve(RightHandSide& Fun, vector_Type& y0, Real t0, Real TF, Real dt_init)
{
	boost::shared_ptr<RightHandSide> FunPtr(new RightHandSide(Fun));

	solve(FunPtr, y0, t0, TF, dt_init);
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::solve( boost::shared_ptr<RightHandSide> Fun, vector_Type& y0, Real t0, Real TF, Real dt_init)
{
	ofstream output("test_ros3p.txt");

	Real t(t0);						//time t_k
	Real dt(dt_init);
	Real dt_old(dt_init);

	Int n= y0.size();
	vector_Type y(y0);					//y contains the solution y_k
	vector<vector_Type> U(s,y0);
	matrix_Type I( getIdentity(y0.map(), n) );
	I.matrixPtr()->FillComplete();
	matrixPtr_Type B(new matrix_Type(y0.map(), n, false));			//Linear system matrix
	vector_Type ytmp(y0.map());				//temporary variable
	vectorPtr_Type Utmp(new vector_Type(y0.map()));
	vectorPtr_Type rhs(new vector_Type(y0.map()));					//rhs will be the right side
	Real err_n;							//error at step n
	Real err_n_1;						//error at step n-1
	Real fac_max = 5.0;					//maximal value for this factor, dt(k+1) < dt(k)*fac_max
	Int k = 1;							//iteration counter
	bool rejected = false;					//used to know if a step is rejected two times consecutively

	output << t << " " << y[0] << " " << y[1] << " " << dt << " " <<rejected<<"\n";

	//First step, to set err_n_1
	cout<<"Begin of iteration k = 0\n";

	zero(U);
	updateMatrix(B, Fun->getJac(y), I*(1.0/(dt*M_g)));
	M_solver.setOperator(B);

	computeStages<RightHandSide>(U, y, ytmp, rhs, Utmp, Fun, dt);

	y = y + combLin(U, M_m);
	*Utmp = combLin(U, M_m - M_mhat );
	err_n_1 = Utmp->norm2();
	t = t + dt;

	cout<<"t(k) = "<<t-dt<<"\n";
	cout<<"dt(k) = "<<dt<<"\n";
	cout<<"err_n_1 = "<<err_n_1<<"\n";
	cout<<"dt(k+1) = "<<dt<<"\n";
	cout<<"Iteration 0 finished."<<"\n\n";

	output << t << " " << y[0] << " " << y[1] << " " << dt << " " <<rejected<<"\n";

	while (t < TF)
	{
		cout<<"Begin of iteration k = "<<k<<"\n";
		cout<<"t("<<k<<") = "<<t<<"\n";
		cout<<"dt("<<k<<") = "<<dt<<"\n";

		zero(U);
		updateMatrix(B, Fun->getJac(y), I*(1.0/(dt*M_g)));
		M_solver.setOperator(B);								//Computing the inverse, which will be used s times
		computeStages<RightHandSide>(U, y, ytmp, rhs, Utmp, Fun, dt);

		if( computeError(U, Utmp, err_n, err_n_1, fac_max, dt, dt_old, TF-t, y.norm2(), rejected) )
		{
			rejected = true;
			continue;
		}
		else
		{
			y = y + combLin(U, M_m);							//upgrading the solution
			t = t + dt;									//upgrading the time
			k++;

			cout<<"dt("<<k<<") = "<<dt<<"\n";
			cout<<"Iteration "<<k-1<<" finished."<<"\n\n";

			output << t << " " << y[0] << " " << y[1] << " " <<dt<< " " << rejected <<"\n";

			rejected = false;
		}

	}

	y0 = y;

	output.close();

}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::computeStages(vector<vector_Type>& U, const vector_Type& y, vector_Type& ytmp, vectorPtr_Type& rhs, vectorPtr_Type& Utmp,
											 boost::shared_ptr<RightHandSide> Fun, Real dt)
{
	for (int i = 0; i<s; i++)
	{
		ytmp = y + combLin(U, M_A, i);				//ytmp = y0 + sum_{j=1}^{i-1} A(i,j)*U(:,j)
		*Utmp = combLin(U, M_C, i)/dt;				//Utmp = sum_{j=1}^{i-1} C(i,j)*U(:,j)/dt
		Fun->computeRhs( ytmp, 0.0, *rhs);
		*rhs = *rhs + *Utmp;
		M_solver.setRightHandSide(rhs);
		M_solver.solve(Utmp);
		U[i] = *Utmp;
	}
}

template<UInt s>
void RosenbrockTransformed<s>::updateMatrix(matrixPtr_Type& B, const matrix_Type& J, const matrix_Type& aI)
{
	*B = J*(-1.0);
	*B += aI;
}

template<UInt s>
bool RosenbrockTransformed<s>::computeError(const vector<vector_Type>& U, vectorPtr_Type& Utmp, Real& err_n, Real& err_n_1,
				  	  	  	  	  	  	  	  Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected)
{
	Real Tol = M_absTol + M_relTol * ynorm;			//Tol = atol + rtol*|y_k|
	Real fac;
	*Utmp =  combLin(U, M_m - M_mhat );				//difference with the embedded method
	err_n = Utmp->norm2();						//norm of the error
	if (err_n == 0.0)							//here we set fac, where dt(k+1) = fac*dt(k)
		fac = fac_max;							//if the actual error is zero we set fac to its maximal value
	else if (err_n_1 == 0.0)					//if the previous error was zero and the actual is not then fac~1
		fac = M_S;
	else										//formula to compute fac, takes in account Tol, errors and time steps
		fac = M_S * std::pow( (Tol*err_n_1)/(err_n*err_n) , M_p_1 ) * ( dt / dt_old ) ;

	if (err_n > Tol)							//the step is rejected
	{
		if (rejected)						//if it is the second time that it is rejected we
			dt /= M_D;							//divide the time step by
		else									//else dt = dt*fac, if the previous step has not grown too much then
			dt *= fac;							//fac<1, if fac>1 it will be rejected one more time and dt will be
													//divided by 10.
		return true;							//this timestep has been rejected
	}
	else
	{
		err_n_1 = err_n;
		dt_old = dt;
		dt = min<Real>( Trem, min<Real>( fac, fac_max )*dt );	// dt(k+1) = min( TF-t, fac*dt(k), fac_max*dt(k) )
		return false;
	}
}

template<UInt s>
matrix_Type RosenbrockTransformed<s>::getIdentity(const map_Type& map, Int n)
{
	matrix_Type I(map, n, false);
	int* Indices = new int[n];
	double* Values = new double[n];

	for( int i=0; i<n; i++)
		Indices[i] = i;

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
			Values[j] = (double)(i==j);

		I.matrixPtr()->InsertGlobalValues (i, n, Values, Indices);
	}

	I.globalAssemble();

	delete Indices;
	delete Values;

	return I;
}

template<UInt s>
vector_Type RosenbrockTransformed<s>::combLin(const vector<vector_Type>& U, const MatrixSmall<s,s>& M, const Int line)
{
	vector_Type sum = M(line,0)*U[0];

	for(int i=1; i<s; i++)
		sum = sum + M(line,i)*U[i];

	return sum;
}

template<UInt s>
vector_Type RosenbrockTransformed<s>::combLin(const vector<vector_Type>& U, const VectorSmall<s>& v)
{
	vector_Type sum = v(0)*U[0];

	for(int i=1; i<s; i++)
		sum = sum + v(i)*U[i];

	return sum;
}

template<UInt s>
void RosenbrockTransformed<s>::zero(vector<vector_Type>& U)
{
	for(int i=0; i<s; i++)
		U[i] = U[i]*0.0;
}

template<UInt s>
void RosenbrockTransformed<s>::initMembers()
{
	M_D = 1.5;
	M_S = 0.95;
	M_absTol = 0.0000001;
	M_relTol = 0.0000001;
	M_p_1 = 1.0/M_p;
}

template<UInt s>
void RosenbrockTransformed<s>::setSolver(boost::shared_ptr<Epetra_Comm>& Comm, const string& solvParam)
{
	Teuchos::RCP< Teuchos::ParameterList > List = Teuchos::rcp ( new Teuchos::ParameterList );
	List = Teuchos::getParametersFromXmlFile ( solvParam );
	M_solver.setCommunicator(Comm);
	M_solver.setParameters(*List);
}

template<UInt s>
void RosenbrockTransformed<s>::setMethod(Real g, const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
		   	   	   	   	   	   	   	   	   const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat, UInt order)
{
	M_g = g;
	M_A = A;
	M_C = C;
	M_gammai = gammai;
	M_a = a;
	M_m = m;
	M_mhat = mhat;
	M_p = (double)(order);
}





