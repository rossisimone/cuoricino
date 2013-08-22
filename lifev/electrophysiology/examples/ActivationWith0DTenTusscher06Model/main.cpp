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
    @brief 0D test with the minimal model

    @date 01−2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
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

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher06.hpp>
#include <lifev/electrophysiology/solver/XbModels/XbNegroniLascano96.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using std::cout;
using std::endl;
using namespace LifeV;

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
    Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( "Parameters.xml" ) );
    Teuchos::ParameterList nlList    = * ( Teuchos::getParametersFromXmlFile ( "NegroniLascano96Parameters.xml" ) );

    std::cout << " Done!" << endl;
    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for TenTusscher 2006 Model with parameters ... ";
    IonicTenTusscher06  ionicModel;
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    ionicModel.showMe();


    //********************************************//
    // Initialize the solution to 0. The model    //
    // consist of three state variables. Xe.Size()//
    // returns the number of state variables of   //
    // the model. rStates is the reference to the //
    // the vector states                          //
    //********************************************//
    std::cout << "Initializing solution vector...";
    std::vector<Real> states (ionicModel.restingConditions());
    std::vector<Real>& rStates = states;
    std::cout << " Done!" << endl;


    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//
    std::cout << "Initializing rhs..." ;
    std::vector<Real> rhs (ionicModel.Size(), 0);

    std::cout << " Done! "  << endl;




    //********************************************//
    // The model needs as external informations   //
    // the contraction velocity and the Calcium   //
    // concentration.                             //
    //********************************************//
    Real Iapp (0.0);


    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//
    Real TF (list.get("TF",100.));
    Real dt (list.get("dt",100.));


    //********************************************//
    // Xbs                                         //
    //********************************************//
	std::cout << "Building Constructor for NegrpniLascano96 Model with parameters ... ";
    XbNegroniLascano96  xb ( nlList );
	xb.showMe();

    std::cout << "Initializing Xb solution vector...";
    std::vector<Real> XbStates (xb.Size(), 0);
//    XbStates.at (0) = 7.0;
//    XbStates.at (1) = 8.0;
//    XbStates.at (2) = 1.0;
//    XbStates.at (3) = 0.021553043080281;


    std::cout << "Initializing Xb rhs..." ;
    std::vector<Real> XbRhs (xb.Size(), 0);

    Real vel  ( 0.0 );
    Real Ca   ( 0.0 );

    Real dgammaf=vel;

    std::ofstream XbOutput ("XbOutput.txt");




    cout << "Potential: " << rStates.at (0) << endl;
    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";

    Real gammaf( list.get("gfi",0.) );
    Real gamma_f_initial = gammaf;

	Real alpha( list.get("alpha",1.0) );

	Real pwr2( list.get("power2",2.0) );
	Real eta( list.get("eta",100.0) );
	Real Cai_diast =  1.1703 / 10;

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    string filename = "output.txt";
    std::ofstream output ("output.txt");
    output << 0 << "\t";
    for ( int it (0); it <= ionicModel.Size() - 1; it++)
     {
         output << rStates.at (it) << "\t";
     }
    output << 0.0 << "\t";
    output << gammaf << "\n";


    XbOutput << XbStates[0] << "\t";
	XbOutput << XbStates[1] << "\t";
	XbOutput << XbStates[2] << "\t";
	XbOutput << 0.0 << "\t";
	XbOutput << dgammaf << "\t";
	XbOutput << gammaf << "\n";



    ///********************************
    // DEBUGGING
    //////////////////////////
    Real V = states[0];
    Real m = states[1];
	Real h = states[2];
    Real j = states[3];
    Real d = states[4];
    Real f = states[5];
    Real f2 = states[6];
    Real fCass = states[7];
    Real r = states[8];
    Real s = states[9];
    Real Xr1 = states[10];
    Real Xr2 = states[11];
    Real Xs = states[12];
    Real Nai = states[13];
    Real Ki = states[14];
    Real Cai = states[15];
    Real Cass = states[16];
    Real Casr = states[17];
    Real Rprime = states[18];

    Real itot = ionicModel.Itot(V, m, h, j, d, f, f2, fCass, r, s, Xr1, Xr2, Xs, Nai, Ki, Cai, Cass );
    Real iK1 = ionicModel.IK1(V,Ki);
    Real ito = ionicModel.Ito(V,r,s,Ki);
    Real iKr = ionicModel.IKr(V,Xr1,Xr2,Ki);
    Real iKs = ionicModel.IKs(V,Xs,Ki,Nai);
    Real iCaL = ionicModel.ICaL(V,d,f,f2,fCass,Cass);
    Real iNa = ionicModel.INa(V,m,h,j,Nai);
    Real ibNa = ionicModel.IbNa(V,Nai);
    Real iNaCa = ionicModel.INaCa(V,Nai,Cai);
    Real ibCa = ionicModel.IbCa(V,Cai);
    Real ipK = ionicModel.IpK(V,Ki);
    Real ipCa = ionicModel.IpCa(Cai);

    string filename2 = "currents.txt";
    std::ofstream output2 ("currents.txt");
    output2 << 0 << "\t";
    output2 << itot << "\t";
    output2 << iK1 << "\t";
    output2 << ito << "\t";
    output2 << iKr << "\t";
    output2 << iKs << "\t";
    output2 << iCaL << "\t";
    output2 << iNa << "\t";
    output2 << ibNa << "\t";
    output2 << iNaCa << "\t";
    output2 << ibCa << "\t";
    output2 << ipK << "\t";
    output2 << ipCa << "\n";






    std::cout << "\nItot: " << ionicModel.Itot(V, m, h, j, d, f, f2, fCass, r, s, Xr1, Xr2, Xs, Nai, Ki, Cai, Cass );
    std::cout << "\nIK1: " << ionicModel.IK1(V,Ki);
    std::cout << "\nIto: " << ionicModel.Ito(V,r,s,Ki);
    std::cout << "\nIKr: " << ionicModel.IKr(V,Xr1,Xr2,Ki);
    std::cout << "\nIKs: " << ionicModel.IKs(V,Xs,Ki,Nai);
    std::cout << "\nICaL: " << ionicModel.ICaL(V,d,f,f2,fCass,Cass);
    std::cout << "\nINa: " << ionicModel.INa(V,m,h,j,Nai);
    std::cout << "\nIbNa: " << ionicModel.IbNa(V,Nai);
    std::cout << "\nINaCa: " << ionicModel.INaCa(V,Nai,Cai);
    std::cout << "\nIbCa: " << ionicModel.IbCa(V,Cai);
    std::cout << "\nIpK: " << ionicModel.IpK(V,Ki);
    std::cout << "\nIpCa: " << ionicModel.IpCa(Cai);


    int savestep( ( list.get("savestep",1.) / dt ) );
    int pacingstepstart( ( list.get("pacingstep",1000.) / dt ) );

    int iter(0);

    std::cout << "\nMembrane capacitance: " << ionicModel.membraneCapacitance() << "\n";
    std::cout.precision(16);
    for ( Real t = 0; t < TF; )
    {
        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if( iter % pacingstepstart == 0 )
		{
        	std::cout << "\nStarting stimulation: " << t << "\n";
        	Iapp = list.get("Iapp",100.);
            ionicModel.setAppliedCurrent(Iapp);
		}
        int iter2( iter - 1.0 / dt );
        if( iter2 % pacingstepstart == 0 && iter >= 0 )
		{
        	std::cout << "\nStopping stimulation: " << t << "\n";
        	Iapp = 0.0;
            ionicModel.setAppliedCurrent(Iapp);
		}

    	iter++;


//        std::cout << "\nIapp: "<<Iapp ;
//
//        std::cout << "\nIonic Iapp: "<< ionicModel.appliedCurrent() ;
//        std::cout << "\n: ";
//
//        ionicModel.computeRhs ( states, rhs);
//        std::cout << "\nrhs : "<< rhs[0] ;
//        std::cout << "\n: ";
//        rhs[0] += Iapp;
        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time.                          //
        //********************************************//
        V = states[0];
        m = states[1];
    	h = states[2];
        j = states[3];
        d = states[4];
        f = states[5];
        f2 = states[6];
        fCass = states[7];
        r = states[8];
        s = states[9];
        Xr1 = states[10];
        Xr2 = states[11];
        Xs = states[12];
        Nai = states[13];
        Ki = states[14];
        Cai = states[15];
        Cass = states[16];
        Casr = states[17];
        Rprime = states[18];

        ionicModel.solveOneStep(states,dt);

//		rStates.at (0) = rStates.at (0)  + dt * rRhs.at (0);
//        ionicModel.computeGatingVariablesWithRushLarsen( states, dt);
//
//        int offset = 1 + ionicModel.numberOfGatingVariables();
//        for ( int j (0); j < ( ionicModel.Size() - offset ); j++)
//			{
//				rStates.at (j+offset) = rStates.at (j+offset)  + dt * rRhs.at (j+offset);
//
//			}

                //std::vector<Real> BErhs   ( xb.computeBackwardEuler( XbStates, Ca, vel, dt ) );
        //XbStates.at (2) = BErhs.at (2);

        Real FLR;
        Real I4f = (gammaf+1.0) * (gammaf+1.0);
    	if(  I4f > 0.87277 && I4f < 1.334)
    	{

			Real d0 = -4.333618335582119e3;
			Real d1 = 2.570395355352195e3;
			Real e1 = -2.051827278991976e3;
			Real d2 = 1.329536116891330e3;
			Real e2 = 0.302216784558222e3;
			Real d3 = 0.104943770305116e3;
			Real e3 = 0.218375174229422e3;
			Real l0 = 1.95;

			FLR = d0/2 + d1 * std::sin(I4f * l0)
							  + e1 * std::cos(I4f * l0)
							  + d2 * std::sin(2 * I4f * l0)
							  + e2 * std::cos(2 * I4f * l0)
							  + d3 * std::sin(3 * I4f * l0)
							  + e3 * std::cos(3 * I4f * l0);
    	}
    	else
    	{
    		FLR = 0.0;
    	}


        //vel = 0;
        Ca = Cai * 1000;

    	xb.computeRhs( XbStates, Ca, vel, XbRhs, FLR );


    	std::vector<Real> XbRhs2 (xb.Size(), 0);
    	std::vector<Real> Xb2 (xb.Size(), 0);

    	Xb2[0] = XbStates.at (0)  + dt * XbRhs.at (0);
    	Xb2[1] = XbStates.at (1)  + dt * XbRhs.at (1);
    	Xb2[2] = XbStates.at (2)  + dt * XbRhs.at (2);
    	xb.computeRhs( Xb2, Ca, vel, XbRhs2, FLR );

        XbStates.at (0) = XbStates.at (0)  + 0.5 * dt * ( XbRhs.at (0) + XbRhs2.at (0) ) ;
        XbStates.at (1) = XbStates.at (1)  + 0.5 * dt * ( XbRhs.at (1) + XbRhs2.at (1) ) ;
        XbStates.at (2) = XbStates.at (2)  + 0.5 * dt * ( XbRhs.at (2) + XbRhs2.at (2) ) ;



    	//Cai = Cai / 1e-3;

    	Real Pa;
    	if(Ca < Cai_diast)
    	{
    		Pa = 0.0;
    		//Pa = alpha * ( XbStates[1] + XbStates[2] );

    	}
    	else
    	{
    		Real C50=std::pow(10, 6-5.33);
    		//Pa = alpha * FLR * std::pow( Cai - Cai_diast, pwr2);// * ( Cai - Cai_diast);
    		Pa = - alpha * alpha * alpha * alpha * FLR * std::pow( Cai , pwr2) / ( std::pow( Cai , pwr2) + std::pow( C50 , pwr2));// * ( Cai - Cai_diast);
    		    		//Pa = alpha * ( XbStates[1] + XbStates[2] );
    	}
    	//std::cout << "\n\ntime: " << t << "\n alpha: " << alpha << ", FLR: " << FLR<< ", Cai-Cai_d: " << Cai-Cai_diast  << ", Pa: " << Pa;
    	Real dW = - 2.0 / ( gammaf + 1.0 );
    	Real dW0= - 2.0 / ( gamma_f_initial + 1.0 );
    	//std::cout << "\ndW: " << dW << ", dW0: " << dW0;
    	Real activeViscosity = eta;// * std::pow(Ca, pwr);// * Cai * Cai * Cai * Cai;
    	//PARAM: power: 3, 3, eta =50000, alpha = -50
    	//Real dCa2 = std::pow( ionicModel.dCai(V, Nai, Cai, Casr, Cass), pwr );
    	//Real f = std::exp( - ( ( 1000. * ( t / 1000 - std::floor(t/1000)  ) )
    	//				* ( 1000. * ( t / 1000 - std::floor(t/1000) ) ) ) / 50 );
    	//f = f * (1 - Cai_diast) + Cai_diast;


//    	if( f < 0.5)
//    	{
//    		Real activeViscosity = eta / f / f / 4.0;  ;
//    	}
//    	else
//    	{
//    		f = 1.0;
//    	}
		//Real activeViscosity = eta /  f;

    	dgammaf = (Pa + dW0 - dW) / activeViscosity;
    	vel = dgammaf;
    	//std::cout << "\neta: " << activeViscosity << ", dCa2: " << dCa2 << ", gamma_f: " << gammaf;
    	gammaf = gammaf + dt * dgammaf;


        //********************************************//
        // Writes solution on file.                   //
        //********************************************//

        t = t + dt;
        if( iter % savestep == 0  ){
        output << t << "\t";
        for ( int it (0); it <= ionicModel.Size() - 1; it++)
        {
            output << rStates.at (it) << "\t";
        }
        output << dgammaf << "\t";
        output << gammaf << "\n";

        XbOutput << XbStates[0] << "\t";
        XbOutput << XbStates[1] << "\t";
        XbOutput << XbStates[2] << "\t";
    	XbOutput << Pa << "\t";
        XbOutput << dgammaf << "\t";
        XbOutput << gammaf << "\n";


        ///********************************
        // DEBUGGING
        //////////////////////////


        itot = ionicModel.Itot(V, m, h, j, d, f, f2, fCass, r, s, Xr1, Xr2, Xs, Nai, Ki, Cai, Cass );
        iK1 = ionicModel.IK1(V,Ki);
        ito = ionicModel.Ito(V,r,s,Ki);
        iKr = ionicModel.IKr(V,Xr1,Xr2,Ki);
        iKs = ionicModel.IKs(V,Xs,Ki,Nai);
        iCaL = ionicModel.ICaL(V,d,f,f2,fCass,Cass);
        iNa = ionicModel.INa(V,m,h,j,Nai);
        ibNa = ionicModel.IbNa(V,Nai);
        iNaCa = ionicModel.INaCa(V,Nai,Cai);
        ibCa = ionicModel.IbCa(V,Cai);
        ipK = ionicModel.IpK(V,Ki);
        ipCa = ionicModel.IpCa(Cai);

        output2 << t << "\t";
         output2 << itot << "\t";
         output2 << iK1 << "\t";
         output2 << ito << "\t";
         output2 << iKr << "\t";
         output2 << iKs << "\t";
         output2 << iCaL << "\t";
         output2 << iNa << "\t";
         output2 << ibNa << "\t";
         output2 << iNaCa << "\t";
         output2 << ibCa << "\t";
         output2 << ipK << "\t";
         output2 << ipCa << "\n";

//         std::cout << t << "\t";
//          std::cout << itot << "\t";
//          std::cout << iK1 << "\t";
//          std::cout << ito << "\t";
//          std::cout << iKr << "\t";
//          std::cout << iKs << "\t";
//          std::cout << iCaL << "\t";
//          std::cout << iNa << "\t";
//          std::cout << ibNa << "\t";
//          std::cout << iNaCa << "\t";
//          std::cout << ibCa << "\t";
//          std::cout << ipK << "\t";
//          std::cout << ipCa << "\n";
        }






        //********************************************//
        // Update the time.                           //
        //********************************************//

    }
    std::cout << "\n...Time loop ends.\n";
    std::cout << "Solution written on file: " << filename << "\n";
    //********************************************//
    // Close exported file.                       //
    //********************************************//
    output.close();
    XbOutput.close();
    output2.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
