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

using namespace LifeV;

#define SolutionTestNorm  -2.848213312500001e+03

Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
    {
    	std::cout << "% using MPI" << std::endl;
    }


    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//
    std::cout << "Importing parameters list...";
    Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( "Parameters.xml" ) );
    std::cout << " Done!" << std::endl;
    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for TenTusscher 2006 Model with parameters ... ";
    IonicTenTusscher06  ionicModel;
    std::cout << " Done!" << std::endl;


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
    std::vector<Real> states (ionicModel.restingConditions() );
    std::cout << " Done!" << std::endl;


    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//
    std::cout << "Initializing rhs..." ;
    std::vector<Real> rhs (ionicModel.Size(), 0);
    std::cout << " Done! "  << std::endl;




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
    Real TF (list.get ("TF", 100.) );
    Real dt (list.get ("dt", 100.) );


    //********************************************//
    // We record the norm of the solution to      //
    // check the failure of the test              //
    //********************************************//
    Real SolutionNorm = states[0];


    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    std::string filename = "output.txt";
    std::ofstream output ("output.txt");
    output << 0 << "\t";

    for ( int j (0); j < ionicModel.Size() - 1; j++)
    {
        output << states[j] << "\t";
    }
    output << states[ ionicModel.Size() - 1 ] << "\n";

    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";


    int savestep ( ( list.get ("savestep", 1.) / dt ) );
    int pacingstepstart ( ( list.get ("pacingstep", 1000.) / dt ) );
    int pacingstepstop ( ( list.get ("pacingstep", 1001.) / dt ) );

    int iter (1);



    for ( Real t = 0; t < TF; )
    {
        if ( t > 50 && t < 52 )
        {
            Iapp = 50.;
        }
        else
        {
            Iapp = 0;
        }
        std::cout << "\r " << t << " ms.       " << std::flush;


        iter++;
        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        ionicModel.computeRhs ( states, rhs);
        rhs[0] += Iapp;



        states[0] = states[0]  + dt * rhs[0];
        ionicModel.computeGatingVariablesWithRushLarsen ( states, dt);

        int offset = 1 + ionicModel.numberOfGatingVariables();
        for ( int index (0); index < ( ionicModel.Size() - offset ); index++)
        {
            states[index + offset] += dt * rhs[index + offset];

        }


        //********************************************//
        // Update the time.                           //
        //********************************************//
        t = t + dt;


        //********************************************//
        // Writes solution on file.                   //
        //********************************************//
        if ( iter % savestep == 0  )
        {
            output << t << "\t";
            for ( int index (0); index < ionicModel.Size() - 1; index++)
            {
                output << states[index] << "\t";
            }
            output << states[ ionicModel.Size() - 1 ] << "\n";

            //********************************************//
            // Update the norm of the solution to check   //
            // test failure                               //
            //********************************************//
            SolutionNorm += states[0];
        }

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

    Real err = std::abs (SolutionNorm - SolutionTestNorm) / std::abs(SolutionTestNorm);
    if ( err > 1e-2 )
    {
    	std::cout << "\nTest Failed: " <<  err <<"\n" << "\nSolution Norm: " <<  SolutionNorm << "\n";
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}

#undef SolutionTestNorm
