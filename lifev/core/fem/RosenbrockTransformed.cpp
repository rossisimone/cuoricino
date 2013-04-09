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
	//ofstream output("test_ros3p_vectorEpetra.txt");

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

	//output << t << " " << y[0] << " " << y[1] << " " << dt << " " <<rejected<<"\n";

	//First step, to set err_n_1
	//cout<<"Begin of iteration k = 0\n";

	zero(U);
	updateMatrix(B, Fun->getJac(y), I*(1.0/(dt*M_g)));
	M_solver.setOperator(B);

	computeStages<RightHandSide>(U, y, ytmp, rhs, Utmp, Fun, dt);

	y = y + combLin(U, M_m);
	*Utmp = combLin(U, M_m - M_mhat );
	err_n_1 = Utmp->norm2();
	t = t + dt;

	/*
	cout<<"t(k) = "<<t-dt<<"\n";
	cout<<"dt(k) = "<<dt<<"\n";
	cout<<"err_n_1 = "<<err_n_1<<"\n";
	cout<<"dt(k+1) = "<<dt<<"\n";
	cout<<"Iteration 0 finished."<<"\n\n";
	*/

	//output << t << " " << y[0] << " " << y[1] << " " << dt << " " <<rejected<<"\n";

	while (t < TF)
	{
		//cout<<"Begin of iteration k = "<<k<<"\n";
		//cout<<"t("<<k<<") = "<<t<<"\n";
		//cout<<"dt("<<k<<") = "<<dt<<"\n";

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

			//cout<<"dt("<<k<<") = "<<dt<<"\n";
			//cout<<"Iteration "<<k-1<<" finished."<<"\n\n";

			//output << t << " " << y[0] << " " << y[1] << " " <<dt<< " " << rejected <<"\n";

			rejected = false;
		}

	}

	y0 = y;

	//output.close();

	//cout<<"Computed with VECTOR EPETRA\n";
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::solve( boost::shared_ptr<RightHandSide> Fun, vector<Real>& y0, Real t0, Real TF, Real dt_init)
{
	//ofstream output("test_ros3p_vectorReal.txt");

	Real t(t0);						//time t_k
	Real dt(dt_init);
	Real dt_old(dt_init);

	UInt n= y0.size();
	vector<Real> y(y0);					//y contains the solution y_k
	vector<vector<Real> > U(s,vector<Real>(n,0.0));
	vector<vector<Real> > I( getIdentity(n) );
	vector<vector<Real> > Psys(I),Qsys(I),Usys(I),Lsys(I);
	vector<vector<Real> > B(n,vector<Real>(n,0.0));			//Linear system matrix
	vector<Real> ytmp(y0);				//temporary variable
	vector<Real> Utmp(y0);
	vector<Real> rhs(y0);					//rhs will be the right side
	Real err_n;							//error at step n
	Real err_n_1;						//error at step n-1
	Real fac_max = 5.0;					//maximal value for this factor, dt(k+1) < dt(k)*fac_max
	Int k = 1;							//iteration counter
	bool rejected = false;					//used to know if a step is rejected two times consecutively

	//output << t << " " << y[0] << " " << y[1] << " " << dt << " " <<rejected<<"\n";

	//First step, to set err_n_1
	//cout<<"Begin of iteration k = 0\n";

	zero(U, n);
	B = I;
	mult(B, 1.0/(dt*M_g), n);
	minusequal(B, Fun->getJac(y), n);
	LU(B, Psys, Qsys, Lsys, Usys, n);

	computeStages<RightHandSide>(U, y, ytmp, rhs, Utmp, Fun, dt, Lsys, Usys, Psys, Qsys, n);

	y = sum(y, combLin(U, M_m, n), n);
	Utmp = combLin(U, M_m - M_mhat, n);

	err_n_1 = norm2(Utmp, n);
	t = t + dt;
	/*
	cout<<"t(k) = "<<t-dt<<"\n";
	cout<<"dt(k) = "<<dt<<"\n";
	cout<<"err_n_1 = "<<err_n_1<<"\n";
	cout<<"dt(k+1) = "<<dt<<"\n";
	cout<<"Iteration 0 finished."<<"\n\n";
	*/

	//output << t << " " << y[0] << " " << y[1] << " " << dt << " " <<rejected<<"\n";

	while (t < TF)
	{
		/*
		cout<<"Begin of iteration k = "<<k<<"\n";
		cout<<"t("<<k<<") = "<<t<<"\n";
		cout<<"dt("<<k<<") = "<<dt<<"\n";
		*/

		zero(U, n);
		B = I;
		mult(B, 1.0/(dt*M_g), n);
		minusequal(B, Fun->getJac(y), n);
		LU(B, Psys, Qsys, Lsys, Usys, n);								//Computing the inverse, which will be used s times
		computeStages<RightHandSide>(U, y, ytmp, rhs, Utmp, Fun, dt, Lsys, Usys, Psys, Qsys, n);

		if( computeError(U, Utmp, err_n, err_n_1, fac_max, dt, dt_old, TF-t, norm2(y,n), rejected, n) )
		{
			rejected = true;
			continue;
		}
		else
		{
			y = sum(y, combLin(U, M_m, n), n);							//upgrading the solution
			t = t + dt;									//upgrading the time
			k++;

			/*
			cout<<"dt("<<k<<") = "<<dt<<"\n";
			cout<<"Iteration "<<k-1<<" finished."<<"\n\n";
			*/

			//output << t << " " << y[0] << " " << y[1] << " " <<dt<< " " << rejected <<"\n";

			rejected = false;
		}

	}

	y0 = y;

	//output.close();

	//cout<<"Computed with VECTOR<REAL>\n";

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

template<UInt s>
void RosenbrockTransformed<s>::LU(vector<vector<Real> >& A, vector<vector<Real> >& P, vector<vector<Real> >& Q,
								  vector<vector<Real> >& L, vector<vector<Real> >& U, UInt n)
{

	vector< vector<Real> > I(getIdentity(n));
	vector< vector<Real> > Minv(I), Pk(I), Qk(I), Mk(I), Mkinv(I);
	P = I;
	Q = I;

	for(int k = 1; k<n; k++)
	{
		Pk = I;
		Qk = I;
		pivot(A, k, n, Pk, Qk);
		A = mult(Pk,A,n);
		A = mult(A,Qk,n);
		Mk = I;
		MGauss(A, k, n, Mk);
		Mkinv = I;
		mult(Mkinv, 2.0, n);
		minusequal(Mkinv, Mk, n);
		A = mult(Mk, A, n);
		P = mult(Pk, P, n);
		Q = mult(Q, Qk, n);
		Minv = mult(Minv, Pk, n);
		Minv = mult(Minv, Mkinv, n);
	}

	U = triu(A, n);
	L = mult(P, Minv, n);

}

template<UInt s>
vector< vector<Real> > RosenbrockTransformed<s>::getIdentity(UInt n)
{
	vector< vector<Real> > I(n, vector<Real>(n,0.0) );

	for(int i=0; i<n; i++)
		I[i][i] = 1.0;

	return I;
}

template<UInt s>
void RosenbrockTransformed<s>::pivot(const vector<vector<Real> >& A, UInt k, UInt n, vector<vector<Real> >& Pk, vector<vector<Real> >& Qk)
{
	int i = k-1;
	int j = k-1;
	Real max = abs(A[i][j]);

	for(int l = k-1; l<n; l++)
	{
		for(int c = k-1; c<n; c++)
		{
			if( abs(A[l][c]) > max)
			{
				max = abs(A[l][c]);
				i = l;
				j = c;
			}
		}
	}

	Pk[i][i] = 0.;	Pk[k-1][k-1] = 0.;	Pk[k-1][i] = 1.;	Pk[i][k-1] = 1.;
	Qk[j][j] = 0.;	Qk[k-1][k-1] = 0.;	Qk[k-1][j] = 1.;	Qk[j][k-1] = 1.;
}

template<UInt s>
void RosenbrockTransformed<s>::MGauss(const vector<vector<Real> >& A, UInt k, UInt n, vector<vector<Real> >& Mk)
{
	for(int i = k; i<n; i++)
		Mk[i][k-1] = -A[i][k-1]/A[k-1][k-1];
}

template<UInt s>
void RosenbrockTransformed<s>::mult(vector<vector<Real> >& B, Real a, UInt n)
{
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			B[i][j] *= a;
}

template<UInt s>
void RosenbrockTransformed<s>::minusequal(vector<vector<Real> >& A, const vector<vector<Real> >& B, UInt n)
{
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			A[i][j] -= B[i][j];
}

template<UInt s>
vector< vector<Real> > RosenbrockTransformed<s>::mult(const vector<vector<Real> >& A, const vector<vector<Real> >& B, UInt n)
{
	vector<vector<Real> > C(n, vector<Real>(n,0.0));

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			for(int k=0; k<n; k++)
				C[i][j] += A[i][k]*B[k][j];
		}
	}

	return C;
}

template<UInt s>
vector<Real> RosenbrockTransformed<s>::mult(const vector<vector<Real> >& A, const vector<Real>& b, UInt n)
{
	vector<Real> c(n,0.0);

	for(int i=0; i<n; i++)
	{
		for(int k=0; k<n; k++)
			c[i] += A[i][k]*b[k];
	}

	return c;
}

template<UInt s>
vector<vector<Real> > RosenbrockTransformed<s>::triu(const vector<vector<Real> >& B, UInt n)
{
	vector< vector<Real> > U(n, vector<Real>(n,0.0));

	for(int i=0; i<n; i++)
		for(int j=i; j<n; j++)
			U[i][j] = B[i][j];

	return U;
}

template<UInt s>
void RosenbrockTransformed<s>::zero(vector<vector<Real> >& U, UInt n)
{
	for(int i=0; i<s; i++)
		for(int j=0; j<n; j++)
			U[i][j] = 0.0;
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::computeStages(vector< vector<Real> >& U, const vector<Real>& y, vector<Real>& ytmp,
											 vector<Real>& rhs, vector<Real>& Utmp, boost::shared_ptr<RightHandSide> Fun,
											 Real dt, vector< vector<Real> >& Lsys, vector< vector<Real> >& Usys,
											 vector< vector<Real> >& Psys, vector< vector<Real> >& Qsys, UInt n)
{
	for (int i = 0; i<s; i++)
	{
		ytmp = sum(y, combLin(U, M_A, i, n), n);				//ytmp = y0 + sum_{j=1}^{i-1} A(i,j)*U(:,j)
		Utmp = combLin(U, M_C, i, n);				//Utmp = sum_{j=1}^{i-1} C(i,j)*U(:,j)/dt
		multequal(Utmp, 1.0/dt, n);
		Fun->computeRhs( ytmp, 0.0, rhs);
		rhs = sum(rhs, Utmp, n);
		rhs = mult(Psys, rhs, n);
		rhs = solvel(Lsys,rhs, n);
		rhs = solveu(Usys,rhs, n);
		rhs = mult(Qsys, rhs, n);
		U[i] = rhs;
	}
}

template<UInt s>
vector<Real> RosenbrockTransformed<s>::combLin(const vector<vector<Real> >& U, const MatrixSmall<s,s>& M, const Int line, UInt n)
{
	vector<Real> sum(n,0.0);

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<s; j++)
			sum[i] += M(line,j)*U[j][i];
	}

	return sum;
}

template<UInt s>
vector<Real> RosenbrockTransformed<s>::combLin(const vector< vector<Real> >& U, const VectorSmall<s>& v, UInt n)
{
	vector<Real> sum(n,0.0);

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<s; j++)
			sum[i] += v(j)*U[j][i];
	}

	return sum;
}

template<UInt s>
vector<Real> RosenbrockTransformed<s>::sum(const vector<Real>& a, const vector<Real>& b, UInt n)
{
	vector<Real> sum(n,0.0);

	for(int i=0; i<n; i++)
		sum[i] = a[i] + b[i];

	return sum;
}

template<UInt s>
void RosenbrockTransformed<s>::multequal(vector<Real>& v, Real a, UInt n)
{
	for(int i=0; i<n; i++)
		v[i] *= a;
}

template<UInt s>
vector<Real> RosenbrockTransformed<s>::solveu(const vector<vector<Real> >& U, const vector<Real>& b, UInt n)
{
	vector<Real> x(n,0.0);
	Real tmp;

	x[n-1] = b[n-1]/U[n-1][n-1];

	for(int i=2; i<=n; i++)
	{
		tmp = 0.;
		for(int k = n-i+1; k<=n-1; k++)
			tmp += U[n-i][k]*x[k];
		x[n-i] = (b[n-i]-tmp)/U[n-i][n-i];
	}

	return x;
}

template<UInt s>
vector<Real> RosenbrockTransformed<s>::solvel(const vector<vector<Real> >& L, const vector<Real>& b, UInt n)
{
	vector<Real> x(n,0.0);
	Real tmp;

	x[0] = b[0]/L[0][0];

	for(int i=1; i<n; i++)
	{
		tmp = 0.;
		for(int k = 0; k<=i-1; k++)
			tmp += L[i][k]*x[k];
		x[i] = (b[i]-tmp)/L[i][i];
	}

	return x;
}

template<UInt s>
Real RosenbrockTransformed<s>::norm2(const vector<Real>& v, UInt n)
{
	Real norm = 0.0;

	for(int i=0; i<n; i++)
		norm += v[i]*v[i];

	return std::sqrt(norm);
}

template<UInt s>
bool RosenbrockTransformed<s>::computeError(const vector<vector<Real> >& U, vector<Real>& Utmp, Real& err_n, Real& err_n_1,
				  	  	  	  	  	  	  	  Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected, UInt n)
{
	Real Tol = M_absTol + M_relTol * ynorm;			//Tol = atol + rtol*|y_k|
	Real fac;
	Utmp =  combLin(U, M_m - M_mhat, n );				//difference with the embedded method
	err_n = norm2(Utmp, n);						//norm of the error
	if (err_n == 0.0)							//here we set fac, where dt(k+1) = fac*dt(k)
		fac = fac_max;							//if the actual error is zero we set fac to its maximal value
	else if (err_n_1 == 0.0)					//if the previous error was zero and the actual is not then fac~1
		fac = M_S;
	else										//formula to compute fac, takes in account Tol, errors and time steps
		fac = M_S * pow( (Tol*err_n_1)/(err_n*err_n) , M_p_1 ) * ( dt / dt_old ) ;

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
void RosenbrockTransformed<s>::disp(vector<vector<Real> >& M, string name, UInt n)
{
	cout<<name<<" = \n";
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
			cout<<"   "<<M[i][j];
		cout<<"\n";
	}
	cout<<"\n\n";
}





