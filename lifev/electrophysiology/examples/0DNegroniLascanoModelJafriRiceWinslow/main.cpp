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
	  @brief Ionic model based on Jafri, Rice And Winslow model coupled with
	  @ Negroni Lascano crossbridge model
	  @date 04-2013
	  @author Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

	  @contributors
	  @mantainer Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
	  @last update 03-2013
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

#include <lifev/electrophysiology/solver/XbModels/XbNegroniLascano96.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicJafriRiceWinslow.hpp>
#include <lifev/core/LifeV.hpp>

#include <lifev/electrophysiology/solver/StimulationProtocol.hpp>

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
    Teuchos::ParameterList nlParameterList    = * ( Teuchos::getParametersFromXmlFile ( "NegroniLascano96Parameters.xml" ) );
    Teuchos::ParameterList ionicParameterList = * ( Teuchos::getParametersFromXmlFile ( "JafriRiceWinslowParameters.xml" ) );
    Teuchos::ParameterList pacingPParameterList = * ( Teuchos::getParametersFromXmlFile ( "StimulationParameters.xml" ) );
    std::cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//

	std::cout << "Building Constructor for NegrpniLascano96 Model with parameters ... ";
    XbNegroniLascano96  xb ( nlParameterList );
    IonicJafriRiceWinslow ionicModel ( ionicParameterList );
    StimulationProtocol   stimulation ( pacingPParameterList );
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//

	xb.showMe();
	stimulation.showMe();


	//********************************************//
	// Initialize the solution to 0. The model    //
	// consist of 31 state variables.             //
	// model.Size() returns the number of state   //
	// variables of the model.                    //
	//********************************************//

	std::cout << "Initializing solution vector...";
    std::vector<Real> XbStates (xb.Size(), 0);
    std::vector<Real> states (ionicModel.Size(), 0);
	states.at (0)  = - 84.1638;
    states.at (1)  = 0.0328302;
    states.at (2)  = 0.988354;
    states.at (3)  = 0.992540;
	states.at (4)  = 0.000928836;
    states.at (5)  = 10.2042;
    states.at (6)  = 143.727;
	states.at (7)  = 5.4;
	states.at (8)  = 9.94893e-11;
    states.at (9)  = 1.243891;
    states.at (10) = 1.36058e-4;
    states.at (11) = 1.17504;
    states.at (12) = 0.762527;
	states.at (13) = 1.19168e-3;
    states.at (14) = 6.30613e-9;
    states.at (15) = 0.236283;
    states.at (16) = 0.997208;
    states.at (17) = 6.38897e-5;
    states.at (18) = 1.535e-9;
    states.at (19) = 1.63909e-14;
    states.at (20) = 6.56337e-20;
    states.at (21) = 9.84546e-21;
    states.at (22) = 2.72826e-3;
	states.at (23) = 6.99215e-7;
    states.at (24) = 6.71989e-11;
    states.at (25) = 2.87031e-15;
    states.at (26) = 4.59752e-20;
	states.at (27) = 0.0;
    states.at (28) = 0.998983;
    states.at (29) = 0.00635;
	states.at (30) = 0.13598;
    std::cout << " Done!" << endl;


    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//

	std::cout << "Initializing rhs..." ;
    std::vector<Real> XbRhs (xb.Size(), 0);
    std::vector<Real> rhs (ionicModel.Size(), 0);
    std::cout << " Done! "  << endl;


    //********************************************//
    // The model needs as external informations   //
    // the contraction velocity and the Calcium   //
    // concentration.                             //
    //********************************************//

    Real vel  ( 0.0 );
    Real Ca   ( 0.0 );
    Real Iapp ( 0.0 );
    Real X    ( 1.05 );

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//

    Real TF     = nlParameterList.get( "endTime", 5.0 );
    Real dt     = nlParameterList.get( "timeStep", 5.77e-5 );

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//

	string filename = "output.txt";
    std::ofstream output ("output.txt");
    string XbFilename = "XbOutput.txt";
    std::ofstream XbOutput ("XbOutput.txt");


    //********************************************//
    // Time loop starts.                          //
    //********************************************//

    int iter(0);
    int savedt( ionicParameterList.get( "savedt", 1.0) / dt );
    int NbStimulus ( 0 );

	std::cout << "Time loop starts...\n";
    for ( Real t = 0; t < TF; )
    {
		// Stimuli for ionic model

    	stimulation.pacingProtocolChoice( t, dt, NbStimulus, Iapp ); // Protocol stimulation
    	// The list of protocols are described in the StimulationProtocol.hpp


		// Velocity of motion

        xb.computeVelocity( dt, X, vel );


        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//

    	Ca = states.at(8)*1000;
    	// Because the concentration is in mM in the ionic model
    	// and in the Xb model it sould be in uM.

        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        
		xb.computeRhs             ( XbStates, Ca, vel, XbRhs );
        ionicModel.computeRhs     ( states, Iapp, rhs );
        std::vector<Real> gateInf ( ionicModel.gateInf( states ) );


        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time for all the variables     //
        // different that the ionic concentrations.   //
        // This ones are treated with Euler implicit  //
        // method and Newton algorithm.               //
        //********************************************//
        
		for(int j(0); j <= 30; ++j)
        {
    		if ( ( j <= 4 ) || ( j >= 12 ) )
    			states.at (j) = states.at (j)  + dt * rhs.at (j);
        }
		states.at (5)  = ionicModel.computeNewtonNa    (states, dt, 10);
		states.at (6)  = ionicModel.computeNewtonKi    (states, dt, 10);
		states.at (7)  = ionicModel.computeNewtonKo    (states, dt, 10);
		states.at (8)  = ionicModel.computeNewtonCai   (states, dt, 10);
		states.at (9)  = ionicModel.computeNewtonCaNSR (states, dt, 10);
		states.at (10) = ionicModel.computeNewtonCaSS  (states, dt, 10);
		states.at (11) = ionicModel.computeNewtonCaJSR (states, dt, 10);
		
        XbStates.at (0) = XbStates.at (0)  + dt * XbRhs.at (0);
        XbStates.at (1) = XbStates.at (1)  + dt * XbRhs.at (1);
//        XbStates.at (2) = XbStates.at (2)  + dt * XbRhs.at (2);

		// Implicit method
        std::vector<Real> BErhs   ( xb.computeBackwardEuler( XbStates, Ca, vel, dt ) );
//        XbStates.at (0) = BErhs.at (0);
//        XbStates.at (1) = BErhs.at (1);
        XbStates.at (2) = BErhs.at (2);

        //********************************************//
        // Writes solution on file.                   //
        //********************************************//

        iter++;
        if( iter % savedt == 0)
        {
        	for ( int j (0); j < ionicModel.Size() - 1; j++)
        	{
            	output << states.at (j) << ", ";
        	}
        	output << states.at ( ionicModel.Size() - 1 ) << "\n";

        	XbOutput << t << ", " << X << ", " << XbStates.at (0) << ", " << XbStates.at (1) << ", " << XbStates.at (2) << "\n";
        }

        //********************************************//
        // Update the time.                           //
        //********************************************//

        t = t + dt;
    }
    std::cout << "\n...Time loop ends.\n";
    std::cout << "Solution written on file: " << filename << "\n";

    //********************************************//
    // Close exported file.                       //
    //********************************************//

    XbOutput.close();
    output.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}

