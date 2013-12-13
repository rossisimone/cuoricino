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

#include <lifev/electrophysiology/solver/IonicModels/IonicNoblePurkinje.hpp>
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

    //std::cout << "Importing parameters list...";
    // Teuchos::ParameterList NLParameterList = *( Teuchos::getParametersFromXmlFile( "NegroniLascano96Parameters.xml" ) );
    //std::cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for NegrpniLascano96 Model with parameters ... ";
    IonicNoblePurkinje  ionicModel;
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
    std::vector<Real> states (ionicModel.Size(), 0);
    ionicModel.initialize (states);
    //    states.at (0) = -80.0;
    //    states.at (1) = ionicModel.mInf(-80.0);
    //    states.at (2) = ionicModel.nInf(-80.0);
    //    states.at (3) = ionicModel.hInf(-80.0);
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
    Real Iapp (0);


    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//
    Real TF (4000);
    Real dt (0.005);
    int useRushLarsen (1);

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    string filename = "output.txt";
    std::ofstream output ("output.txt");

    cout << "Potential: " << rStates.at (0) << endl;
    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";
    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if ( t > 1.5 && t < 2 )
        {
            Iapp = 240.0;
        }
        else
        {
            Iapp = 0;
        }
        //        Iapp = 0;
        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        ionicModel.setAppliedCurrent (Iapp);
        if (useRushLarsen)
        {
            ionicModel.computeGatingVariablesWithRushLarsen (states, dt);
            Real RHS = ionicModel.computeLocalPotentialRhs ( states);
            RHS += Iapp;
            //********************************************//
            // Use forward Euler method to advance the    //
            // solution in time.                          //
            //********************************************//

            for ( int j (0); j < 1; j++)
            {
                rStates.at (j) = rStates.at (j)  + dt * (RHS);
            }
        }
        else
        {
            ionicModel.computeRhs ( states, rhs);
            rhs[0] += Iapp;
            //********************************************//
            // Use forward Euler method to advance the    //
            // solution in time.                          //
            //********************************************//

            for ( int j (0); j < ionicModel.Size(); j++)
            {
                rStates.at (j) = rStates.at (j)  + dt * rRhs.at (j);
            }
        }

        //********************************************//
        // Writes solution on file.                   //
        //********************************************//
        output << t << ", ";
        for ( int j (0); j < ionicModel.Size() - 1; j++)
        {
            output << rStates.at (j) << ", ";
        }
        output << rStates.at ( ionicModel.Size() - 1 ) << "\n";

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
    output.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    Real returnValue;

    if (std::abs (rStates.at ( ionicModel.Size() - 2 ) - 0.456093) > 1e-4 )
    {
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}
