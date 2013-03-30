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
    @brief 0D test with the Fitz-Hugh Nagumo model.

    @date 01−2013
    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @contributor
    @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <fstream>
#include <string>
#include <algorithm>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/heart/solver/IonicModels/IonicFitzHughNagumo.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace std;
using namespace LifeV;


void EulerExplicit (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output);

void ROS3P (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& S, const Real& D,
			const Real& I, std::ofstream& output);

template<UInt n, UInt s>
void RosenbrockTransformed( IonicFitzHughNagumo model, const VectorSmall<n>& y0, Real t0, Real TF, Real dt_init,
							Real g,	const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
							const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat,
							const Real& S, const Real& D, const Real& Iapp, std::ofstream& output );

MatrixSmall<2,2> Invert(const MatrixSmall<2,2>& A);
template<UInt Dim1, UInt Dim2> void setCol(MatrixSmall<Dim1,Dim2>& A, const VectorSmall<Dim1>& b, const Int& j);


Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }


    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    std::cout << "Importing parameters list...";
    Teuchos::ParameterList FHNParameterList = * ( Teuchos::getParametersFromXmlFile ( "FitzHughNagumoParameters.xml" ) );
    std::cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Fitz-Hugh Nagumo. The           //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for Fitz-Hugh Nagumo Model with parameters ... ";
    IonicFitzHughNagumo  model ( FHNParameterList );
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    model.showMe();

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//

    Real TF ( FHNParameterList.get ("TF", 300.0) );
    Real dt ( FHNParameterList.get ("dt", 0.01) );
    Real S ( FHNParameterList.get ("S", 0.9) );
    Real D ( FHNParameterList.get ("D", 10.0) );

    cout << "Time parameters : " << endl;
    cout << "TF = " << TF << endl;
    cout << "dt = " << dt << endl;

    string filename;
    if(FHNParameterList.get ("Rosen", 1.0) == 0.0)
    	filename = "output_ExpEuler.txt";
    else
    	filename = "output_ROS3P.txt";

    ofstream output(filename.c_str());

    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n\n\n";

    LifeChrono chrono;
    chrono.start();

    if ( FHNParameterList.get ("Rosen", 1.0) == 0.0 )
        EulerExplicit (dt, TF, model, FHNParameterList.get ("Iapp", 2000.0), output);
    else
        ROS3P (dt, TF, model, S, D, FHNParameterList.get ("Iapp", 2000.0), output);

    chrono.stop();
    std::cout << "\n...Time loop ends.\n";
    std::cout << "\nElapsed time : " << chrono.diff() << std::endl;
    std::cout << "Solution written on file: " << filename << "\n";
    //********************************************//
    // Close exported file.                       //
    //********************************************//
    output.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}

void EulerExplicit (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output)
{
    std::vector<Real> unknowns ( model.Size(), 0.0);
    unknowns.at(0) = 98.999200000000002;
    unknowns.at(1) = 0.023763800000000;

    std::vector<Real> rhs ( model.Size(), 0.0);
    Real Iapp;

    cout << "Computing using Explicit Euler" << endl;

    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if ( t > 10.0 && t < 10.1 )
        {
            Iapp = 0.0;
        }
        else
        {
            Iapp = 0.0;
        }

        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        model.computeRhs ( unknowns, Iapp, rhs);

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time.                          //
        //********************************************//
        unknowns.at (0) = unknowns.at (0)  + dt * rhs.at (0);
        unknowns.at (1) = unknowns.at (1)  + dt * rhs.at (1);


        //********************************************//
        // Writes solution on file.                   //
        //********************************************//
        output << t << " " << unknowns.at (0) << " " << unknowns.at (1) << "\n";

        //********************************************//
        // Update the time.                           //
        //********************************************//
        t = t + dt;
    }
}

void ROS3P (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& S, const Real& D,
			const Real& I, std::ofstream& output)
{
	cout << "Computing using ROS3P." << endl;

	// Initialization of ROS3P parameters
	MatrixSmall<3,3> A;
	MatrixSmall<3,3> C;
	VectorSmall<3> gammai;
	VectorSmall<3> a;
	VectorSmall<3> m;
	VectorSmall<3> mhat;
	Real g(0.7886751345948129);

	A(1,0) = 1.267949192431123;			A(2,0) = 1.267949192431123;
	C(1,0) = -1.607695154586736;		C(2,0) = -3.464101615137755;		C(2,1) = -1.732050807568877;
	gammai(0) = 0.7886751345948129;		gammai(1) = -0.2113248654051871;	gammai(2) = -1.077350269189626;
	a(0) = 0.0;							a(1) = 1.0;							a(2) = 1.0;
	m(0) = 2.0;							m(1) = 0.5773502691896258;			m(2) = 0.4226497308103742;
	mhat(0) = 2.113248654051871;		mhat(1) = 1.0;						mhat(2) = 0.4226497308103742;

	//Problem initial values
	VectorSmall<2> y0;
	Real t0(0.0);
	y0(0) = 98.999200000000002;		y0(1) = 0.023763800000000;

	//Call of the general Rosenbrock, passing ROS3P parameters, right side and initial values.
	//Implemented with variable changes to avoid some matrix*vector multiplications
	RosenbrockTransformed < 2, 3 > ( model, y0, t0, TF, dt, g, A, C, gammai, a, m, mhat, S, D, I, output );

}

template<UInt n, UInt s>
void RosenbrockTransformed( IonicFitzHughNagumo model, const VectorSmall<n>& y0, Real t0, Real TF, Real dt_init,
							Real g,	const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
							const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat,
							const Real& S, const Real& D, const Real& Iapp, std::ofstream& output )
{
	// Constants
	Int N = (TF-t0)/(100*dt_init);	//Initial size of the vector containing the solution
	MatrixSmall<n,n> I;				//Identity matrix
	for(int i=0; i<n; i++)
		I(i,i) = 1.0;

	//Problem variables
	Real It(0.0);					//Applied current
	VectorSmall<n> y(y0);			//y contains the solution y_k
	Real t(t0);						//time t_k
	Real dt(dt_init);				//time step dt_k
	Real dt_old(dt_init);
	Real h = 1.0e-8;

	//Stepsize control variables
	//Real S = 0.9;						//Safety parameter for the new time step
	Real abs_tol = 0.0000001;			//absolute error tolerance
	Real rel_tol = 0.00001;				//relative error tolerance
	Real p = 3.0; 						//order of the method, which is also phat+1
	Real p_1 = 1/p;						//used in computations
	//Real D = 2.0;

	//The Rosenbrock method requires the multiplication J*(sum_{j=1}^{i-1} alpha(i,j)*K_j)
	//This multiplication can be avoided using U_i = sum_{j=1}^{i-1} alpha(i,j)*K_j
	//In the following the U_i variables are used

	//Auxiliary variables
	MatrixSmall<n,s> U;					//matrix which columns are the U_i
	MatrixSmall<n,n> B;					//Linear system matrix
	VectorSmall<n> ytmp;				//temporary variable
	VectorSmall<s> line;
	VectorSmall<n> Utmp;				//temporary variable
	VectorSmall<n> rhs;					//rhs will be the right side
	Real err_n;							//error at step n
	Real err_n_1;						//error at step n-1
	Real fac;							//factor used for the new time step, dt(k+1) = dt(k) * fac
	Real fac_max = 5.0;					//maximal value for this factor, dt(k+1) < dt(k)*fac_max
	Real Tol;							//tolerance, which depends on abs_tol, rel_tol and y_k
	Int k = 1;							//iteration counter
	Int rejected = 0;					//used to know if a step is rejected two times consecutively
	Int tot_rej = 0;
	Int dub_rej = 0;

	output << t << " " << y(0) << " " << y(1) << " " << dt << " " <<rejected<<"\n";

	//First step, to set err_n_1
	cout<<"Begin of iteration k = 0"<<endl;
	cout<<"Manip matrix..."<<endl;
	B = I/(dt*g);
	cout<<"Calling class member..."<<endl;
	B = B - model.computeJ( t, y);
	cout<<"Inverting B...";
	B = Invert(B);
	cout<<"Done !"<<endl;
	for (int i = 0; i<s; i++)
	{
		cout<<"In for...\nExtracting column..."<<endl;
		line = A.extract(i);
		cout<<"Multiplying by U..."<<endl;
		ytmp = U*line;
		cout<<"Summing..."<<endl;
		ytmp = y + ytmp;
		cout<<"All togheter..."<<endl;
		ytmp = y + U*A.extract(i); 				//ytmp = y0 + sum_{j=1}^{i-1} A(i,j)*U(:,j)
		cout<<"One more time..."<<endl;
		Utmp = (U*C.extract(i))/dt;				//Utmp = sum_{j=1}^{i-1} C(i,j)*U(:,j)/dt
		cout<<"Computing rhs...";
		model.computeRhs ( y, It, rhs);
		cout<<"Done!\nUpdating Utmp..."<<endl;
		Utmp = B*( rhs + Utmp );
		setCol<n,s>(U, Utmp, i);
		cout<<"End for cycle"<<endl;
	}
	cout<<"Updating solution..."<<endl;
	y = y + U*m;
	Utmp = U*m-U*mhat;
	err_n_1 = Utmp.norm();
	t = t + dt;

	cout<<"t(k) = "<<t-dt<<endl;
	cout<<"dt(k) = "<<dt<<endl;
	cout<<"err_n_1 = "<<err_n_1<<endl;
	cout<<"dt(k+1) = "<<dt<<endl;
	cout<<"Iteration 0 finished."<<endl<<endl;

	output << t << " " << y(0) << " " << y(1) << " " << dt << " " <<rejected<<"\n";

	while (t < TF)
	{
		cout<<"Begin of iteration k = "<<k<<endl;
		cout<<"t("<<k<<") = "<<t<<endl;
		cout<<"dt("<<k<<") = "<<dt<<endl;

		U *= 0.0;									//U variables initilized at 0
		B = I/(dt*g) - model.computeNumJ( t, y, h );		//Building linear system matrix
		B = Invert(B);								//Computing the inverse, which will be used s times

		for (int i = 0; i<s; i++)					//Computing the s stages
		{
			ytmp = y + U*A.extract(i);				//Point where f will be evalued
			Utmp = (U*C.extract(i))/dt;				//U_i = sum_{j=1}^{i-1} C(i,j)*U_j
			model.computeRhs( y, It, rhs);			//rhs = f(ytmp)
			Utmp = B*( rhs + Utmp );				//solving the linear system, Utmp = U_i
			setCol<n,s>(U, Utmp, i);				//Putting U_i in the ith column of U
		}

		Tol = abs_tol + rel_tol * y.norm();			//Tol = atol + rtol*|y_k|
		Utmp =  U*(m-mhat) ;						//difference with the embedded method
		err_n = Utmp.norm();						//norm of the error
		if (err_n == 0.0)							//here we set fac, where dt(k+1) = fac*dt(k)
			fac = fac_max;							//if the actual error is zero we set fac to its maximal value
		else if (err_n_1 == 0.0)					//if the previous error was zero and the actual is not then fac~1
			fac = S;
		else										//formula to compute fac, takes in account Tol, errors and time steps
			fac = S * pow( (Tol*err_n_1)/(err_n*err_n) , p_1 ) * ( dt / dt_old ) ;

		cout<<"Tol = "<<Tol<<endl;
		cout<<"err_n_1 = "<<err_n_1<<endl;
		cout<<"err_n   = "<<err_n<<endl;
		cout<<"fac = "<<fac<<endl;

		if (err_n > Tol)							//the step is rejected
		{
			cout<<"Rejected step at time "<<t<<" with dt = "<<dt<<" and fac = "<<fac<<endl<<endl;
			if (rejected == 1)						//if it is the second time that it is rejected we
			{
				dt /= D;							//divide the time step by D
				dub_rej++;							//increase the counter of double rejections
			}
			else									//else dt = dt*fac, if the previous step has not grown too much then
				dt *= fac;							//fac<1, if fac>1 it will be rejected one more time and dt will be
													//divided by 10.
			rejected = 1;							//this timestep has been rejected
			tot_rej++;								//increase the total number of rejections
			continue;								//redo the iteration with the new time step
		}

		//This part of the cycle is executed only if the time step has been accepted

		y = y + U*m;								//upgrading the solution
		t = t + dt;									//upgrading the time

		//setting the new time step and upgrading variables
		err_n_1 = err_n;
		dt_old = dt;
		dt = min<Real>( TF-t, min<Real>( fac, fac_max )*dt );	// dt(k+1) = min( TF-t, fac*dt(k), fac_max*dt(k) )
		k++;

		cout<<"dt("<<k<<") = "<<dt<<endl;
		cout<<"Iteration "<<k-1<<" finished."<<endl<<endl;

		output << t << " " << y(0) << " " << y(1) << " " <<dt<< " " << rejected <<"\n";

		rejected = 0;								//here the step has been accepted, reset rejections to 0
	}

	cout<<"Total rejections : "<<tot_rej<<endl;
	cout<<"Double rejections : "<<dub_rej<<endl<<endl;


}

MatrixSmall<2,2> Invert(const MatrixSmall<2,2>& A)
{
	Real det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
	MatrixSmall<2,2> B;
	B(0,0) = A(1,1);
	B(0,1) = -A(0,1);
	B(1,0) = -A(1,0);
	B(1,1) = A(0,0);
	B = B/det;

	return B;
}

template<UInt Dim1, UInt Dim2>
void setCol(MatrixSmall<Dim1,Dim2>& A, const VectorSmall<Dim1>& b, const Int& j)
{
	for( Int i=0; i<Dim1; i++)
		A(i,j) = b(i);
}













