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
	  @brief Ionic model based on Jafri, Rice And Winslow model.
	  @date 03-2013
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

#include <lifev/heart/solver/XbModels/XbNegroniLascano96.hpp>
#include <lifev/heart/solver/IonicModels/IonicTenTusscher.hpp>
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
    Teuchos::ParameterList NLParameterList    = * ( Teuchos::getParametersFromXmlFile ( "NegroniLascano96Parameters.xml" ) );
    Teuchos::ParameterList IonicParameterList = * ( Teuchos::getParametersFromXmlFile ( "TenTusscherParameters.xml" ) );
    std::cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//

	std::cout << "Building Constructor for NegrpniLascano96 Model with parameters ... ";
    XbNegroniLascano96  xb ( NLParameterList );
    IonicTenTusscher  ionicModel ( IonicParameterList );
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//

	xb.showMe();


    //********************************************//
    // Initialize the solution to 0. The model    //
    // consist of three state variables. Xe.Size()//
    // returns the number of state variables of   //
    // the model. rStates is the reference to the //
    // the vector states                          //
    //********************************************//

	std::cout << "Initializing solution vector...";
    std::vector<Real> XbStates (xb.Size(), 0);
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

    Real TF     = NLParameterList.get("endTime", 100.0);
    Real dt     = NLParameterList.get("timeStep", 0.001);
    Real timeSt = IonicParameterList.get( "stimuliTime", 10.0 );
    Real stInt  = IonicParameterList.get( "stimuliInterval", 1000.0 );


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
    int savedt( IonicParameterList.get( "savedt", 1.0) / dt );

	std::cout << "Time loop starts...\n";
    for ( Real t = 0; t < TF; )
    {
		// Stimuli for ionic model
    	if ( t >= timeSt && t <= timeSt + 1.0 )
    	{
    	   	Iapp = -52.0;
    	   	if ( t >= timeSt + 1.0 - dt && t <= timeSt + 1.0 )
    	   		timeSt = timeSt + stInt;
    	}
    	else
    	   	Iapp = 0;

		// Velocity of motion

        xb.computeVelocity( dt, X, vel );


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
        xb.computeRhs             ( XbStates, Ca, vel, XbRhs );
        ionicModel.computeRhs     ( states, Iapp, rhs );
        std::vector<Real> gateInf ( ionicModel.gateInf( states ) );

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time.                          //
        //********************************************//
        for ( int j (0); j < ionicModel.Size(); j++)
        {
        	if ( j < 13 && j != 0 )
        		states.at (j) = gateInf.at(j-1) + ( states.at (j) - gateInf.at(j-1) ) * exp( dt * rhs.at(j) );
        	else
        	    states.at (j) = states.at (j)   + dt * rhs.at (j);
        }
        XbStates.at (0) = XbStates.at (0)  + dt * XbRhs.at (0);
        XbStates.at (1) = XbStates.at (1)  + dt * XbRhs.at (1);
//        XbStates.at (2) = XbStates.at (2)  + dt * XbRhs.at (2);

		// Implicit method
        std::vector<Real> BErhs    ( xb.computeBackwardEuler( XbStates, Ca, vel, dt ) );
//		XbStates.at (0) = BErhs.at (0);
//		XbStates.at (1) = BErhs.at (1);
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

			XbOutput << X << ", " << XbStates.at (0) << ", " << XbStates.at (1) << ", " << XbStates.at (2) << "\n";
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
