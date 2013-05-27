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
	  @brief Ionic model based on Ten Tusscher 2004 model coupled with
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

#include <lifev/electrophysiology/solver/XbModels/XbNegroniLascano96TTJRW.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher.hpp>
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
    Teuchos::ParameterList nlParameterList      = * ( Teuchos::getParametersFromXmlFile ( "NegroniLascano96Parameters.xml" ) );
    Teuchos::ParameterList ionicParameterList   = * ( Teuchos::getParametersFromXmlFile ( "TenTusscherParameters.xml" ) );
    Teuchos::ParameterList pacingPParameterList = * ( Teuchos::getParametersFromXmlFile ( "StimulationParameters.xml" ) );
    std::cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//

	std::cout << "Building Constructor for NegrpniLascano96 Model with parameters ... ";
    XbNegroniLascano96TTJRW  xb     ( nlParameterList );
    IonicTenTusscher    ionicModel  ( ionicParameterList );
    StimulationProtocol stimulation ( pacingPParameterList );
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//

	xb.showMe();
	stimulation.showMe();


    //********************************************//
    // Initialize the solution to 0. The model    //
    // consist of three state variables. Xe.Size()//
    // returns the number of state variables of   //
    // the model. rStates is the reference to the //
    // the vector states                          //
    //********************************************//

	std::cout << "Initializing solution vector...";
    std::vector<Real> XbStates (xb.Size(), 0);
    std::vector<Real> XbStatesTemp (xb.Size(), 0);
    std::vector<Real> states (ionicModel.Size(), 0);
    states.at (0)  = - 86.2;
    states.at (1)  = 0.0;
    states.at (2)  = 0.75;
    states.at (3)  = 0.75;
    states.at (4)  = 0.0;
    states.at (5)  = 1.0;
    states.at (6)  = 0.0;
    states.at (7)  = 0.0;
    states.at (8)  = 1.0;
    states.at (9)  = 1.0;
    states.at (10) = 0.0;
    states.at (11) = 1.0;
    states.at (12) = 1.0;
    states.at (13) = 2e-4;
    states.at (14) = 0.2;
    states.at (15) = 11.6;
    states.at (16) = 138.3;
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
    std::vector<Real> XbRhsTemp (xb.Size(), 0);
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
    Real X    ( 1.045 );
    Real L    = nlParameterList.get("inLength", 1.05);

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//

    Real TF     = nlParameterList.get("endTime", 100.0);
    Real dt     = nlParameterList.get("timeStep", 0.001);

    Real tStim  ( 0 );

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//

    string filename        = "output.txt";
    string filenameStimPro = "outputStimPro.txt";
    string XbFilename      = "XbOutput.txt";

    std::ofstream output        ("output.txt");
    std::ofstream outputStimPro ("outputStimPro.txt");
    std::ofstream XbOutput      ("XbOutput.txt");


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


		//********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//

    	Ca = states.at(13)*1000;
    	// Because the concentration is in mM in the ionic model
    	// and in the Xb model it sould be in uM.

        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//

        xb.computeRhs             ( XbStates, Ca, L, vel, XbRhs );
        ionicModel.computeRhs     ( states, Iapp, rhs );
        std::vector<Real> gateInf ( ionicModel.gateInf( states ) );

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time for the concentration     //
        // and Rush and Larsen for the gating         //
        // variables                                  //
        //********************************************//

        for ( int j (0); j < ionicModel.Size(); j++)
        {
        	if ( j < 13 && j != 0 )
        		states.at (j) = gateInf.at(j-1) + ( states.at (j) - gateInf.at(j-1) ) * exp( dt * rhs.at(j) );
        	else
        	    states.at (j) = states.at (j)   + dt * rhs.at (j);
        }

        XbStatesTemp = XbStates;
        XbRhsTemp    = XbRhs;

        XbStates.at (0) = XbStates.at (0)  + dt * XbRhs.at (0);
        XbStates.at (1) = XbStates.at (1)  + dt * XbRhs.at (1);
        XbStates.at (2) = XbStates.at (2)  + dt * XbRhs.at (2);

        xb.computeHalfSarcomereLength( t + dt, L );
        xb.computeX        ( dt, L, X );
        xb.computeVelocity ( X, L, vel );
        xb.computeRhs      ( XbStates, Ca, L, vel, XbRhs );

        XbStates.at (0) = XbStatesTemp.at (0) + 0.5 * dt * ( XbRhs.at (0) + XbRhsTemp.at (0) );
        XbStates.at (1) = XbStatesTemp.at (1) + 0.5 * dt * ( XbRhs.at (1) + XbRhsTemp.at (1) );
        XbStates.at (2) = XbStatesTemp.at (2) + 0.5 * dt * ( XbRhs.at (2) + XbRhsTemp.at (2) );


		// Implicit method
//        xb.computeBackwardEuler( XbStates, Ca, vel, dt );

        //********************************************//
        // Writes solution on file.                   //
        //********************************************//

		iter++;
		if( iter % savedt == 0 )
		{
			for ( int j (0); j < ionicModel.Size() - 1; j++)
			{
				output << states.at (j) << ", ";
			}
			output << states.at ( ionicModel.Size() - 1 ) << "\n";

			XbOutput << t << ", " << L << ", " << X << ", " << XbStates.at (0) << ", " << XbStates.at (1) << ", " << XbStates.at (2) << "\n";
		}

		tStim = stimulation.timeSt();

		if ( t >= tStim && t <= tStim + dt )
			outputStimPro << t << "," << states.at(0) << "," << NbStimulus << "\n";


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
    outputStimPro.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
