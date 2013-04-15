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
												UInt order)
:M_g(g), M_A(A), M_C(C), M_gammai(gammai), M_a(a), M_m(m), M_mhat(mhat), M_p(order)
{
	initMembers();
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::solve(RightHandSide& Fun, VectorLU& y0, Real t0, Real TF, Real& dt_init)
{
	boost::shared_ptr<RightHandSide> FunPtr(new RightHandSide(Fun));

	solve(FunPtr, y0, t0, TF, dt_init);
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::solve(RightHandSide& Fun, vector<Real>& y0, Real t0, Real TF, Real& dt_init)
{
	boost::shared_ptr<RightHandSide> FunPtr(new RightHandSide(Fun));

	VectorLU y0LU(y0);
	solve(FunPtr, y0LU, t0, TF, dt_init);
	y0 = y0LU.getVector();
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::solve( boost::shared_ptr<RightHandSide> Fun, vector<Real>& y0, Real t0, Real TF, Real& dt_init)
{
	VectorLU y0LU(y0);
	solve(Fun, y0LU, t0, TF, dt_init);
	y0 = y0LU.getVector();
}

template<UInt s>
template<typename RightHandSide>
void RosenbrockTransformed<s>::solve( boost::shared_ptr<RightHandSide> Fun, VectorLU& y0, Real t0, Real TF, Real& dt_init)
{
	//ofstream output("test_ros3p_vectorLU.txt");

	Real t(t0);						//time t_k
	Real dt(dt_init);
	Real dt_old(dt_init);

	UInt n= y0.size();
	VectorLU y(y0);					//y contains the solution y_k
	MatrixLU U(n, s);
	MatrixLU I(n);
	MatrixLU Psys(n),Qsys(n),Usys(n),Lsys(n);
	MatrixLU B(n);			//Linear system matrix
	VectorLU ytmp(y0);				//temporary variable
	VectorLU Utmp(y0);
	VectorLU rhs(y0);					//rhs will be the right hand side
	Real err_n;							//error at step n
	Real err_n_1;						//error at step n-1
	Real fac_max = 5.0;					//maximal value for this factor, dt(k+1) < dt(k)*fac_max
	Int k = 1;							//iteration counter
	bool rejected = false;					//used to know if a step is rejected two times consecutively

	//output << t << " " << y[0] << " " << y[1] << " " << dt << " " <<rejected<<"\n";

	//First step, to set err_n_1
	//cout<<"Begin of iteration k = 0\n";

	B = I/(dt*M_g) - MatrixLU(Fun->getJac(y.getVector()));
	B.LU(Psys, Qsys, Lsys, Usys);

	computeStages<RightHandSide>(U, y, ytmp, rhs, Utmp, Fun, dt, Lsys, Usys, Psys, Qsys);

	y = y + U.timesVectorSmall<s>(M_m);
	Utmp = U.timesVectorSmall<s>(M_m-M_mhat);

	err_n_1 = Utmp.norm2();
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


		U.setZero();
		B = I/(dt*M_g) - MatrixLU(Fun->getJac(y.getVector()));
		B.LU(Psys, Qsys, Lsys, Usys);								//Computing the inverse, which will be used s times
		computeStages<RightHandSide>(U, y, ytmp, rhs, Utmp, Fun, dt, Lsys, Usys, Psys, Qsys);

		if( computeError(U, Utmp, err_n, err_n_1, fac_max, dt, dt_old, TF-t, y.norm2(), rejected) )
		{
			rejected = true;
			continue;
		}
		else
		{
			y = y + U.timesVectorSmall<s>(M_m);
			t = t + dt;									//upgrading the time
			k++;


			//cout<<"dt("<<k<<") = "<<dt<<"\n";
			//cout<<"Iteration "<<k-1<<" finished."<<"\n\n";


			//output << t << " " << y[0] << " " << y[1] << " " <<dt<< " " << rejected <<"\n";

			rejected = false;
		}

	}

	y0 = y;
	dt_init = dt;

	//output.close();
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
template<typename RightHandSide>
void RosenbrockTransformed<s>::computeStages(MatrixLU& U, const VectorLU& y, VectorLU& ytmp, VectorLU& rhs, VectorLU& Utmp,
		   boost::shared_ptr<RightHandSide> Fun, Real dt, MatrixLU& Lsys, MatrixLU& Usys,
		   MatrixLU& Psys, MatrixLU& Qsys)
{
	for (UInt i = 0; i<s; i++)
	{
		ytmp = y + U.timesVectorSmall<s>(M_A.extract(i));				//ytmp = y0 + sum_{j=1}^{i-1} A(i,j)*U(:,j)
		Utmp = U.timesVectorSmall<s>(M_C.extract(i));				//Utmp = sum_{j=1}^{i-1} C(i,j)*U(:,j)/dt
		Utmp /= dt;
		Fun->computeRhs( ytmp.getVector(), 0.0, rhs.getVector());
		rhs = Psys*(rhs+Utmp);
		rhs = Lsys.solveL(rhs);
		rhs = Usys.solveU(rhs);
		rhs = Qsys*rhs;
		U.setCol(i,rhs);
	}
}

template<UInt s>
bool RosenbrockTransformed<s>::computeError(const MatrixLU& U, VectorLU& Utmp, Real& err_n, Real& err_n_1,
					  Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected)
{
	Real Tol = M_absTol + M_relTol * ynorm;			//Tol = atol + rtol*|y_k|
	Real fac;
	Utmp =  U.timesVectorSmall<s>(M_m-M_mhat);				//difference with the embedded method
	err_n = Utmp.norm2();					//norm of the error
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

