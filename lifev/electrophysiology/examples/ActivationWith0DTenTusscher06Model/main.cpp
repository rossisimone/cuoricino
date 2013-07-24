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
    std::vector<Real>& rRhs = rhs;
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






    cout << "Potential: " << rStates.at (0) << endl;
    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";

    Real gammaf( list.get("gfi",0.) );
    Real gamma_f_initial = gammaf;
    Real l0 = 1.95;
	Real alpha( list.get("alpha",1.0) );
	Real eta( list.get("eta",100.0) );
	Real Cai_diast =  states[15]/1e-4;

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    string filename = "output.txt";
    std::ofstream output ("output.txt");
    output << 0 << ", ";
    for ( int it (0); it <= ionicModel.Size() - 1; it++)
     {
         output << rStates.at (it) << ", ";
     }
     output << gammaf << "\n";

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
    output2 << 0 << ", ";
    output2 << itot << ", ";
    output2 << iK1 << ", ";
    output2 << ito << ", ";
    output2 << iKr << ", ";
    output2 << iKs << ", ";
    output2 << iCaL << ", ";
    output2 << iNa << ", ";
    output2 << ibNa << ", ";
    output2 << iNaCa << ", ";
    output2 << ibCa << ", ";
    output2 << ipK << ", ";
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
    int iter(0);

    bool FE( list.get("FE", false ) );


    std::cout << "\nMembrane capacitance: " << ionicModel.membraneCapacitance() << "\n";
    std::cout.precision(16);
    for ( Real t = 0; t < TF; )
    {
    	iter++;
        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if( t > 50 && t < 51 ) Iapp = list.get("Iapp",100.);
        else if( t > 1050 && t < 1051 ) Iapp = list.get("Iapp",100.);
        else if( t > 2050 && t < 2051 ) Iapp = list.get("Iapp",100.);
        else Iapp = 0;
        //std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
//        std::cout << "\nIapp: "<<Iapp ;
        ionicModel.setAppliedCurrent(Iapp);

//        std::cout << "\nIonic Iapp: "<< ionicModel.appliedCurrent() ;
//        std::cout << "\n: ";

        ionicModel.computeRhs ( states, rhs);
//        std::cout << "\nrhs : "<< rhs[0] ;
//        std::cout << "\n: ";
        rhs[0] += Iapp;


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
        for(int it2(0) ; it2 < ionicModel.Size(); it2++)
        {
        	if( states[it2] != states[it2] )
           	{
        		std::cout << "\n\nI'm Nan " << it2 << "\n";
        		return 0;
        	}
        }
//        if(!FE)
//        {
//			rStates.at (0) = rStates.at (0)  + dt * rRhs.at (0);
////			rStates.at (1) = ionicModel.M_INF(V)
////						   - ( ionicModel.M_INF(V) - m)
////						   * std::exp(-dt/ ionicModel.TAU_M(V));
////			rStates.at (2) = ionicModel.H_INF(V)
////						   - ( ionicModel.H_INF(V) - h)
////						   * std::exp(-dt/ ionicModel.TAU_H(V));
////			rStates.at (3) = ionicModel.jinf(V)
////						   - ( ionicModel.jinf(V) - states[3])
////						   * std::exp(-dt/ ionicModel.TAU_J(V));
////			rStates.at (4) = ionicModel.dinf(V)
////						   - ( ionicModel.dinf(V) - states[4])
////						   * std::exp(-dt/ ionicModel.TAU_D(V));
////			rStates.at (5) = ionicModel.finf(V)
////						   - ( ionicModel.finf(V) - states[5])
////						   * std::exp(-dt/ ionicModel.TAU_F(V));
////			rStates.at (6) = ionicModel.f2inf(V)
////						   - ( ionicModel.f2inf(V) - states[6])
////						   * std::exp(-dt/ ionicModel.tf2(V));
////			rStates.at (7) = ionicModel.f2inf(states[16])
////						   - ( ionicModel.f2inf(states[16]) - states[7])
////						   * std::exp(-dt/ ionicModel.tf2(states[16]));
////			rStates.at (8) = ionicModel.rinf(V)
////						   - ( ionicModel.rinf(V) - states[8])
////						   * std::exp(-dt/ ionicModel.tr(V));
////			rStates.at (9) = ionicModel.sinf(V)
////						   - ( ionicModel.sinf(V) - states[9])
////						   * std::exp(-dt/ ionicModel.ts(V));
////			rStates.at (10) = ionicModel.Xr1inf(V)
////						   - ( ionicModel.Xr1inf(V) - states[10])
////						   * std::exp(-dt/ ionicModel.tXr1(V));
////			rStates.at (11) = ionicModel.Xr2inf(V)
////						   - ( ionicModel.Xr2inf(V) - states[11])
////						   * std::exp(-dt/ ionicModel.tXr2(V));
////			rStates.at (12) = ionicModel.Xsinf(V)
////						   - ( ionicModel.Xsinf(V) - states[12])
////						   * std::exp(-dt/ ionicModel.tXs(V));
//
//        	ionicModel.computeGatingVariablesWithRushLarsen(states, dt);
//			rStates.at (13) = rStates.at (13)  + dt * rRhs.at (13);
//			rStates.at (14) = rStates.at (14)  + dt * rRhs.at (14);
//
////			Real cabufc = ionicModel.CaBuf( states[15]  );
////			Real dCai = dt * ionicModel.dCai(V, Nai, Cai, Casr, Cass);
////			Real bc = ionicModel.getBufc() - cabufc - dCai - Cai + ionicModel.getKbufc();
////			Real cc = ionicModel.getKbufc() * ( cabufc + dCai + Cai);
////			rStates.at(15) = std::sqrt( bc * bc + 4. * cc ) / 2;
////
//			states[15] = ionicModel.solveCai(V,Nai,Cai,Casr,Cass,dt);
////			Real cassbufc = ionicModel. CaSSBuf( Cass );
////			Real dCass = dt * ionicModel.dCaSS(V, Cai, Casr, Cass, Rprime, d, f, f2, fCass);
////			Real bcss = ionicModel.getBufss() - cassbufc - dCass - Cass + ionicModel.getKbufss();
////			Real ccss = ionicModel.getKbufss() * ( cassbufc + dCass + Cass);
////			rStates.at(16) = ( std::sqrt( bcss * bcss + 4. * ccss ) - bcss )/ 2;
////
//			states[16] = ionicModel.solveCaSS(Cai,Casr,Cass,Rprime,V,d,f,f2,fCass,dt);
////			Real CaCSQN = ionicModel.CaCSQN( Casr );
////			Real dCasr = dt * ionicModel.dCaSR(V, Cai, Casr, Cass, Rprime);
////			Real bcsr = ionicModel.getBufsr() - CaCSQN - dCasr - Casr + ionicModel.getKbufsr();
////			Real ccsr = ionicModel.getKbufsr() * ( CaCSQN + dCasr + Casr);
////			rStates.at(17) = ( std::sqrt( bcsr * bcsr + 4. * ccsr ) - bcsr ) / 2;
//			states[17] = ionicModel.solveCaSR(Cai, Casr, Cass, Rprime, dt);
//
//			rStates.at (18) = rStates.at (18)  + dt * rRhs.at (18);
//    	}
//        else
//        {


//		rStates.at (0) = rStates.at (0)  + dt * rRhs.at (0);
//        ionicModel.computeGatingVariablesWithRushLarsen( states, dt);
//
//        int offset = 1 + ionicModel.numberOfGatingVariables();
//        for ( int j (0); j < ( ionicModel.Size() - offset ); j++)
//			{
//				rStates.at (j+offset) = rStates.at (j+offset)  + dt * rRhs.at (j+offset);
//
//			}
//        }

        ionicModel.solveOneStep(states,dt);
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

    	Cai = Cai / 1e-4;

    	Real Pa;
    	if(Cai < Cai_diast)
    	{
    		Pa = 0.0;
    	}
    	else
    	{
    		Pa = alpha * FLR * ( Cai - Cai_diast) * ( Cai - Cai_diast);
    	}
    	//std::cout << "\n\ntime: " << t << "\n alpha: " << alpha << ", FLR: " << FLR<< ", Cai-Cai_d: " << Cai-Cai_diast  << ", Pa: " << Pa;
    	Real dW = - 2.0 / ( gammaf + 1.0 );
    	Real dW0= - 2.0 / ( gamma_f_initial + 1.0 );
    	//std::cout << "\ndW: " << dW << ", dW0: " << dW0;
    	Real activeViscosity = eta * Cai * Cai * Cai * Cai * Cai * Cai * Cai * Cai;
    	Real dgammaf = (Pa + dW0 - dW) / activeViscosity;
    	//std::cout << "\neta: " << activeViscosity << ", dgf: " << dgammaf << ", gamma_f: " << gammaf;
    	gammaf = gammaf + dt * dgammaf;


//
////        V = states[0];
//        states[0] = ionicModel.solveV(V,m,h,j,d,f,f2,fCass,r,s,Xr1,Xr2,Xs,Nai,Ki,Cai,Cass,Iapp,dt);
//        //        m = states[1];
//        states[1] = ionicModel.solveM(V, m, dt);
////    	h = states[2];
//        states[2] = ionicModel.solveH(V, h, dt);
//		//        j = states[3];
//        states[3] = ionicModel.solveJ(V, j, dt);
////        d = states[4];
//        states[4] = ionicModel.solveD(V, d, dt);
//
//        //        f = states[5];
//        states[5] = ionicModel.solveF(V,f,dt);
////        f2 = states[6];
//        states[6] = ionicModel.solveF2(V,f2,dt);
////        fCass = states[7];
//        states[7] =  ionicModel.solveFCaSS(V,fCass,dt);
////        r = states[8];
//        states[8] =  ionicModel.solveR(V,r,dt);
////        s = states[9];
//        states[9] = ionicModel.solveS(V,s,dt);
////        Xr1 = states[10];
//        states[10] =  ionicModel.solveXr1(V,Xr1,dt);
////        Xr2 = states[11];
//        states[11] =  ionicModel.solveXr2(V,Xr2,dt);
////        Xs = states[12];
//        states[12] =  ionicModel.solveXs(V,Xs,dt);
//
////        Nai = states[13];
//        states[13] =  ionicModel.solveNai(V,m,h,j,Nai,Cai,dt);
////        Ki = states[14];
//        states[14] = ionicModel.solveKi(V,r,s,Xr1,Xr2,Xs,Nai,Ki,Iapp,dt);
//        //        Cai = states[15];
//        states[15] = ionicModel.solveCai(V,Nai,Cai,Casr,Cass,dt);
////        Cass = states[16];
//        states[16] = ionicModel.solveCaSS(Cai,Casr,Cass,Rprime,V,d,f,f2,fCass, dt);
////        Casr = states[17];
//		states[17] = ionicModel.solveCaSR(Cai,Casr,Cass,Rprime,dt);
////        Rprime = states[18];
//        states[18] = ionicModel.solveRR(Casr,Cass,Rprime,dt);
//

        //********************************************//
        // Writes solution on file.                   //
        //********************************************//

        t = t + dt;
        if( iter % savestep == 0  ){
        output << t << ", ";
        for ( int it (0); it <= ionicModel.Size() - 1; it++)
        {
            output << rStates.at (it) << ", ";
        }
        output << gammaf << "\n";



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

        output2 << t << ", ";
         output2 << itot << ", ";
         output2 << iK1 << ", ";
         output2 << ito << ", ";
         output2 << iKr << ", ";
         output2 << iKs << ", ";
         output2 << iCaL << ", ";
         output2 << iNa << ", ";
         output2 << ibNa << ", ";
         output2 << iNaCa << ", ";
         output2 << ibCa << ", ";
         output2 << ipK << ", ";
         output2 << ipCa << "\n";

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
    output2.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
