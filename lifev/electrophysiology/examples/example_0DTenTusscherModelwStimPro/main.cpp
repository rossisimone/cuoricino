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
    @brief 0D test with the Ten Tusscher model with general stimulation protocol

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
#include <sstream>

#include <lifev/core/array/MatrixEpetra.hpp>

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

    cout << "Importing parameters list...";
    Teuchos::ParameterList ionicMParameterList  = * ( Teuchos::getParametersFromXmlFile ( "TenTusscherParameters.xml" ) );
    Teuchos::ParameterList pacingPParameterList = * ( Teuchos::getParametersFromXmlFile ( "StimulationParameters.xml" ) );
    cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Ten Tusscher ionic model.       //
    // It creates also a stimulation protocol     //
    // object that defines the way to create the  //
    // stimulus excitation.                       //
    // The model input are the parameters. Pass   //
    // the parameter list in the constructor      //
    //********************************************//

    cout << "Building Constructor for TenTusscher Model with parameters ... ";
    IonicTenTusscher    model       ( ionicMParameterList );
    StimulationProtocol stimulation ( pacingPParameterList );
    cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//

    model.showMe();
    stimulation.showMe();


    //********************************************//
    // Initialize the solution to 0. The model    //
    // consist of 17 state variables.             //
    // model.Size() returns the number of state   //
    // variables of the model.                    //
    //********************************************//

    cout << "Initializing solution vector...";
    std::vector<Real> unknowns (model.Size(), 0 );
    unknowns[0]  = - 86.2;
    unknowns[1]  = 0.0;
    unknowns[2]  = 0.75;
    unknowns[3]  = 0.75;
    unknowns[4]  = 0.0;
    unknowns[5]  = 1.0;
    unknowns[6]  = 0.0;
    unknowns[7]  = 0.0;
    unknowns[8]  = 1.0;
    unknowns[9]  = 1.0;
    unknowns[10] = 0.0;
    unknowns[11] = 1.0;
    unknowns[12] = 1.0;
    unknowns[13] = 2e-4;
    unknowns[14] = 0.2;
    unknowns[15] = 11.6;
    unknowns[16] = 138.3;
    cout << " Done!" << endl;


    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//

    cout << "Initializing rhs..." ;
    std::vector<Real> rhs (model.Size(), 0);
    cout << " Done! "  << endl;

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

    Real TF     = ionicMParameterList.get ( "endTime", 5.0 );
    Real dt     = ionicMParameterList.get ( "timeStep", 1e-3 );

    Real tStim  ( 0 );

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//

    string filename             = "output.txt";
    string filenameStimPro      = "outputStimPro.txt";

    std::ofstream output        ("output.txt");
    std::ofstream outputStimPro ("outputStimPro.txt");

    //********************************************//
    // Time loop starts.                          //
    //********************************************//

    cout << "Time loop starts...\n";

    int iter (0);
    int savedt ( ionicMParameterList.get ( "savedt", 1.0) / dt );
    int NbStimulus ( 0 );

    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Gives the appropriate Iapp (stimulation .  //
        // current ) according to time variable.      //
        //********************************************//

        stimulation.pacingProtocolChoice ( t, dt, NbStimulus, Iapp ); // Protocol stimulation
        // The list of protocols are described in the StimulationProtocol.hpp

        cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//

        model.computeRhs ( unknowns, Iapp, rhs );
        std::vector<Real> gateInf ( model.gateInf ( unknowns ) );

        //********************************************//
        // Writes solution on file.                   //
        //********************************************//

        iter++;

        if ( iter % savedt == 0)
        {
            output  << t << ", " << unknowns.at (0) << ", " << unknowns.at (1) << ", "
                    << unknowns.at (2) << ", " << unknowns.at (3) << ", "
                    << unknowns.at (4) << ", " << unknowns.at (5) << ", "
                    << unknowns.at (6) << ", " << unknowns.at (7) << ", "
                    << unknowns.at (8) << ", " << unknowns.at (9) << ", "
                    << unknowns.at (10) << ", " << unknowns.at (11) << ", "
                    << unknowns.at (12) << ", " << unknowns.at (13) << ", "
                    << unknowns.at (14) << ", " << unknowns.at (15) << ", "
                    << unknowns.at (16) << "\n";
        }

        tStim = stimulation.timeSt();

        if ( t >= tStim && t <= tStim + dt )
        {
            outputStimPro << t << "," << unknowns.at (0) << "," << NbStimulus << "\n";
        }


        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time for the concentration     //
        // and Rush and Larsen for the gating         //
        // variables                                  //
        //********************************************//

        for (int j (0); j <= 16; ++j)
        {
            if ( j < 13 && j != 0 )
            {
                unknowns.at (j) = gateInf.at (j - 1) + ( unknowns.at (j) - gateInf.at (j - 1) ) * exp ( dt * rhs.at (j) );
            }
            else
            {
                unknowns.at (j) = unknowns.at (j)   + dt * rhs.at (j);
            }

        }

        //********************************************//
        // Update the time.                           //
        //********************************************//

        t = t + dt;
    }

    cout << "\n...Time loop ends.\n";
    cout << "Solution written on file: " << filename << " and " << filenameStimPro << "\n";

    //********************************************//
    // Close exported file.                       //
    //********************************************//

    output.close();
    outputStimPro.close();


    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
