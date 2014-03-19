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
    @brief 0D test with the Jafri Rice Winslow model

    @date 04−2013
    @author Luis Miguel de Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

    @contributor
    @mantainer Luis Miguel de Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
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

#include <lifev/electrophysiology/solver/IonicModels/IonicJafriRiceWinslow.hpp>
#include <lifev/core/LifeV.hpp>

#include <lifev/electrophysiology/stimulus/StimulusPacingProtocol.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using std::endl;
using namespace LifeV;

Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
    {
        std::cout << "% using MPI" << endl;
    }


    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    std::cout << "Importing parameters list...";
    Teuchos::ParameterList ionicMParameterList = * ( Teuchos::getParametersFromXmlFile ( "JafriRiceWinslowParameters.xml" ) );
    Teuchos::ParameterList pacingPParameterList = * ( Teuchos::getParametersFromXmlFile ( "StimulationParameters.xml" ) );
    std::cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Jafri Rice Winslow 1998. The    //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//

    std::cout << "Building Constructor for JafriRiceWinslow Model with parameters ... ";
    IonicJafriRiceWinslow model       ( ionicMParameterList );
//    StimulationProtocol   stimulation ( pacingPParameterList );
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//

    model.showMe();
//    stimulation.showMe();


    //********************************************//
    // Initialize the solution to 0. The model    //
    // consist of 31 state variables.             //
    // model.Size() returns the number of state   //
    // variables of the model.                    //
    //********************************************//

    std::cout << "Initializing solution vector...";
    std::vector<Real> unknowns  (model.Size(), 0 );
    model.initialize(unknowns);
    std::cout << " Done!" << endl;


    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//

    std::cout << "Initializing rhs..." ;
    std::vector<Real> rhs  (model.Size(), 0);
    std::cout << " Done! "  << endl;


    //********************************************//
    //Initialization of the applied current       //
    //********************************************//

    Real Iapp (0.0);



    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//

    Real TF     = ionicMParameterList.get ( "endTime", 5.0 );
    Real dt     = ionicMParameterList.get ( "timeStep", 5.77e-5 );

    Real tStim  ( 0 );

    //********************************************//
    // Setup pacing protocol                      //
    //********************************************//
    StimulusPacingProtocol stimulus;
    stimulus.setParameters(pacingPParameterList);
    stimulus.setTimeStep(dt);

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//

    std::string filename             = "output.txt";
    std::string filenameStimPro      = "outputStimPro.txt";

    std::ofstream output        ("output3.txt");
    std::ofstream outputStimPro ("outputStimPro.txt");

    //********************************************//
    // Time loop starts.                          //
    //********************************************//

    std::cout << "Time loop starts...\n";

    int iter (0);
    int savedt ( ionicMParameterList.get ( "savedt", 1.0) / dt );
    int NbStimulus ( 0 );

    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Set the stimulus over time                 //
        // according to different pacing protocol     //
        //********************************************//

        Iapp =  stimulus.pacingProtocolChoice ( t ); // Protocol stimulation
//        std::cout << "\ntime = " << t << " ms.       " << std::endl;
//        std::cout << "Iapp = " << Iapp << std::endl;
//        std::cout << "time step = " <<stimulus.timeStep() << std::endl;
//        std::cout << "Stimulus starting time = " <<stimulus.startingTimeStimulus() << std::endl;
//        std::cout << "Stimulus duration = " <<stimulus.stimDuration()<< std::endl;
//        std::cout << "V = " << unknowns[0]<< std::endl;


        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        model.setAppliedCurrent (Iapp);
        model.computeRhs ( unknowns, rhs );
        model.addAppliedCurrent(rhs);

        //********************************************//
        // Writes solution on file.                   //
        //********************************************//

        iter++;
        if ( iter % savedt == 0)
        {
            output << t << ", " << unknowns.at (0) << ", " << unknowns.at (1) << ", "
                   << unknowns.at (2) << ", " << unknowns.at (3) << ", "
                   << unknowns.at (4) << ", " << unknowns.at (5) << ", "
                   << unknowns.at (6) << ", " << unknowns.at (7) << ", "
                   << unknowns.at (8) << ", " << unknowns.at (9) << ", "
                   << unknowns.at (10) << ", " << unknowns.at (11) << ", "
                   << unknowns.at (12) << ", " << unknowns.at (13) << ", "
                   << unknowns.at (14) << ", " << unknowns.at (15) << ", "
                   << unknowns.at (16) << ", " << unknowns.at (17) << ", "
                   << unknowns.at (18) << ", " << unknowns.at (19) << ", "
                   << unknowns.at (20) << ", " << unknowns.at (21) << ", "
                   << unknowns.at (22) << ", " << unknowns.at (23) << ", "
                   << unknowns.at (24) << ", " << unknowns.at (25) << ", "
                   << unknowns.at (26) << ", " << unknowns.at (27) << ", "
                   << unknowns.at (28) << ", " << unknowns.at (29) << ", "
                   << unknowns.at (30) << "\n";
        }


        //tStim = stimulation.timeSt();

        if ( t >= tStim && t <= tStim + dt )
        {
            outputStimPro << t << "," << unknowns.at (0) << "," << NbStimulus << "\n";
        }


        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time for all the variables     //
        // different that the ionic concentrations.   //
        // This ones are treated with Euler implicit  //
        // method and Newton algorithm.               //
        //********************************************//
        for (int j (0); j <= 30; ++j)
        {
            unknowns.at (j) = unknowns.at (j)   + dt * rhs.at (j);
            if(unknowns[j]!=unknowns[j])
            {
            	std::cout << "\n FOUND NAN! at " << j <<"\n Stopping";
            	break;
            }
        }


        unknowns.at (5)  = model.computeNewtonNa    (unknowns, dt, 10);
        unknowns.at (6)  = model.computeNewtonKi    (unknowns, dt, 10);
        unknowns.at (7)  = model.computeNewtonKo    (unknowns, dt, 10);
        unknowns.at (8)  = model.computeNewtonCai   (unknowns, dt, 10);
        unknowns.at (9)  = model.computeNewtonCaNSR (unknowns, dt, 10);
        unknowns.at (10) = model.computeNewtonCaSS  (unknowns, dt, 10);
        unknowns.at (11) = model.computeNewtonCaJSR (unknowns, dt, 10);

        //********************************************//
        // Update the time.                           //
        //********************************************//

        t = t + dt;
    }

    std::cout << "\n...Time loop ends.\n";
    std::cout << "Solution written on file: " << filename << " and " << filenameStimPro << "\n";

    //********************************************//
    // Close exported file.                       //
    //********************************************//

    output.close();
    outputStimPro.close();


    //! Finalizing Epetra communicator
    MPI_Finalize();
    Real returnValue;

    if (std::abs (unknowns.at (5) - 10.2042) > 1e-4 )
    {
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}
