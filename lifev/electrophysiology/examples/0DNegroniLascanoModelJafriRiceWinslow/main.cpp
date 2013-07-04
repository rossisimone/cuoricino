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

#include <lifev/electrophysiology/solver/XbModels/XbNegroniLascano96TTJRW.hpp>
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
    XbNegroniLascano96TTJRW  xb ( nlParameterList );
    IonicJafriRiceWinslow ionicModel ( ionicParameterList );
    StimulationProtocol   stimulation ( pacingPParameterList );
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//

    stimulation.showMe();
	xb.showMe();

	//********************************************//
	// Initialize the solution to 0. The model    //
	// consist of 31 state variables.             //
	// model.Size() returns the number of state   //
	// variables of the model.                    //
	//********************************************//

	std::cout << "Initializing solution vector...";
	std::vector<Real> XbStates (xb.Size(), 0);
	std::vector<Real> XbStatesTemp (xb.Size(), 0);

    std::vector<Real> states (ionicModel.Size(), 0);
	states.at (0)  = -86.1638;
    states.at (1)  = 3.28302e-2;
    states.at (2)  = 0.988354;
    states.at (3)  = 0.99254;
	states.at (4)  = 9.28836e-4;
    states.at (5)  = 10.2042;
    states.at (6)  = 143.727;
	states.at (7)  = 5.4;
	states.at (8)  = 9.94893e-5;
    states.at (9)  = 1.24891;
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
    states.at (21) = 9.084546e-21;
    states.at (22) = 2.72826e-3;
	states.at (23) = 6.99215e-7;
    states.at (24) = 6.71989e-11;
    states.at (25) = 2.87031e-15;
    states.at (26) = 4.59752e-20;
	states.at (27) = 0.0;
    states.at (28) = 0.998983;
    states.at (29) = 0.0; // Computation of this value is done with XbStates.
	states.at (30) = 135.9813e-3;


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
    Real              rhsCoupling;
    std::cout << " Done! "  << endl;


    //********************************************//
    // The model needs as external informations   //
    // the contraction velocity and the Calcium   //
    // concentration.                             //
    //********************************************//

    Real vel     ( 0.0 );
    Real Ca      ( 0.0 );
    Real Bi      ( 0.0 );
    Real CmdnTot = ionicModel.cmdnTot()*1e-3;   // Term 1e-3 is to adjust in order to have consistent units
    Real KmCmdn  = ionicModel.constmCmdn()*1e-3;
    Real Iapp    ( 0.0 );
    Real Lm      = nlParameterList.get("inLength", 1.05);
    Real L       = nlParameterList.get("inLength", 1.05);
    Real X       = L - xb.Hc();
    Real F       ( 0.0 );

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//

    Real TF     = nlParameterList.get( "endTime", 5.0 );
    Real dt     = nlParameterList.get( "timeStep", 5.77e-5 );

    Real tStim  ( 0 );

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//

    string filename        = "output1350.txt";
    string filenameStimPro = "outputStimPro.txt";
    string XbFilename      = "XbOutput1350.txt";

    std::ofstream output        ("output1350.txt");
    std::ofstream outputStimPro ("outputStimPro.txt");
    std::ofstream XbOutput      ("XbOutput1350.txt");


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

    	Ca = states.at(8)*1000;
    	Bi = 1 / ( 1 + CmdnTot * KmCmdn / ( ( KmCmdn + Ca ) * ( KmCmdn + Ca ) ) );

    	// Because the concentration is in mM in the ionic model
    	// and in the Xb model it should be in uM.

        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        
        xb.computeRhs             ( XbStates, Ca, L, vel, XbRhs );
        ionicModel.computeRhs     ( states, Iapp, rhs );
        xb.computeCoupling        ( XbStates, Ca, Bi, vel, rhsCoupling );

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time for all the variables     //
        // different that the ionic concentrations.   //
        // This ones are treated with Euler implicit  //
        // method and Newton algorithm.               //
        //********************************************//
        
        for(int j(0); j <= 30; ++j)
        {
        	if ( j != 8 && j != 29)
        		states.at (j) = states.at (j) + dt * rhs.at (j);
        	else if ( j == 8 )
        		states.at (j) = states.at (j) + dt * ( rhs.at (j) + rhsCoupling );
        }


        XbStatesTemp = XbStates;
        XbRhsTemp    = XbRhs;

        XbStates.at (0) = XbStates.at (0)  + dt * XbRhs.at (0);
        XbStates.at (1) = XbStates.at (1)  + dt * XbRhs.at (1);
        XbStates.at (2) = XbStates.at (2)  + dt * XbRhs.at (2);

        Ca = states.at(8)*1000;

        xb.computeTotalMuscleLength( XbStates, t, dt, X, Lm, L, vel, F );

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
        	output << t << ", ";

        	for ( int j (0); j < ionicModel.Size() - 1; j++)
        	{
        		output << states.at (j) << ", ";
        	}

        	output << states.at ( ionicModel.Size() - 1 ) << "\n";

        	XbOutput << t << ", " << rhsCoupling << ", " << Lm << ", " << L << ", " << X << ", "<< F << ", "
        			 << XbStates.at (0) << ", " << XbStates.at (1) << ", " << XbStates.at (2) << "\n";
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
    std::cout << "Solution written on file: " << filename << ", " << XbFilename << " and " << filenameStimPro << "\n";

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

