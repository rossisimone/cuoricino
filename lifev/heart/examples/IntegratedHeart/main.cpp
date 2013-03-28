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
 *  @file
 *  @brief File containing the Integrated Heart example
 *
 *  @date 2013-03-27
 *  @author Toni Lassila <toni.lassila@epfl.ch>
 *  @maintainer Toni Lassila <toni.lassila@epfl.ch>
 *
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <sys/stat.h>
#include <sys/types.h>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/multiscale/solver/MultiscaleSolver.hpp>

using namespace LifeV;
using namespace Multiscale;


class HeartApplication : public MultiscaleSolver
{
public:

    bool
    solveProblem ( const Real& referenceSolution, const Real& tolerance )
    {
//        typedef RegionMesh<LinearTetra>                         mesh_Type;
//        typedef boost::function < Real (const Real& /*t*/,
//                const Real &   x,
//                const Real &   y,
//                const Real& /*z*/,
//                const ID&   /*i*/ ) >   function_Type;
//
//        typedef HeartETAMonodomainSolver< mesh_Type, IonicAlievPanfilov >       monodomainSolver_Type;
//        typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;

#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 8000 ) << "HeartApplication::solveProblem() \n";
#endif

        // Chrono definitions
        LifeChrono buildUpdateChrono;
        LifeChrono solveChrono;
        LifeChrono saveChrono;
        LifeChrono updateSolutionChrono;
        LifeChrono globalChrono;
        Real       totalSimulationTime (0);
        Real       timeStepTime (0);

        for ( ; M_globalData->dataTime()->canAdvance(); M_globalData->dataTime()->updateTime() )
        {
            // Global chrono start
            globalChrono.start();

            if ( M_comm->MyPID() == 0 )
            {
                std::cout << std::endl;
                std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
                std::cout << "                   INTEGRATED HEART SIMULATOR                " << std::endl;
                std::cout << "             time = " << M_globalData->dataTime()->time() << " s;"
                        << " time step number = " << M_globalData->dataTime()->timeStepNumber() << std::endl;
                std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
            }

            // Solve the electrophysiology problem

            // Pass the activation Gamma_f to the LV material law

            // Update the FSI models
            buildUpdateChrono.start();
            if ( M_globalData->dataTime()->isFirstTimeStep() )
            {
                M_model->buildModel();
            }
            else
            {
                M_model->updateModel();
            }
            buildUpdateChrono.stop();

            // Solve the activated FSI model
            solveChrono.start();
            M_model->solveModel();
            solveChrono.stop();

            // Update the solution
            updateSolutionChrono.start();
            M_model->updateSolution();
            updateSolutionChrono.stop();

            // Save the solution
            saveChrono.start();
            if ( M_globalData->dataTime()->timeStepNumber() % multiscaleSaveEachNTimeSteps == 0 || M_globalData->dataTime()->isLastTimeStep() )
            {
                M_model->saveSolution();
            }
            saveChrono.stop();

            // Global chrono stop
            globalChrono.stop();

            // Compute time step time
            timeStepTime = globalChrono.globalDiff ( *M_comm );

            // Updating total simulation time
            totalSimulationTime += timeStepTime;

            if ( M_comm->MyPID() == 0 )
            {
                std::cout << " MS-  Total iteration time:                    " << timeStepTime << " s" << std::endl;
                std::cout << " MS-  Total simulation time:                   " << totalSimulationTime << " s" << std::endl;
            }

            // Save CPU time
            saveCPUTime ( timeStepTime, buildUpdateChrono.globalDiff ( *M_comm ), solveChrono.globalDiff ( *M_comm ),
                    updateSolutionChrono.globalDiff ( *M_comm ), saveChrono.globalDiff ( *M_comm ) );
        }

        return multiscaleExitFlag;
    }

private:

};

int main (int argc, char** argv)
{
    //Setup main communicator
    boost::shared_ptr< Epetra_Comm > comm;

    //Setup MPI variables
    Int numberOfProcesses (1);
    Int rank (0);

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );

    MPI_Comm_size ( MPI_COMM_WORLD, &numberOfProcesses );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
#endif

    if ( rank == 0 )
    {
        std::cout << "MPI Processes: " << numberOfProcesses << std::endl;
    }

#ifdef EPETRA_MPI
    if ( rank == 0 )
    {
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    }
    comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    std::cout << "SERIAL Epetra Initialization ... " << std::endl;
    comm.reset ( new Epetra_SerialComm() );
#endif

    // Setup IntegratedHeart solvers
    bool exitFlag = EXIT_SUCCESS;
    HeartApplication IH;

    // Set the communicator
    IH.setCommunicator ( comm );

    // Command line parameters
    GetPot commandLine ( argc, argv );
    std::string dataFile      = commandLine.follow ( "./Run_IntegratedHeart.dat", 2, "-f", "--file" );
    bool verbose              = commandLine.follow ( true, 2, "-s", "--showme" );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    Real referenceSolution    = -1;
    UInt coresPerNode         = commandLine.follow (  1, 2, "-ns", "--nodesize" );
    Real tolerance            = commandLine.follow (  1e-8, 2, "-t", "--tolerance" );

    if ( coresPerNode > static_cast<UInt> ( numberOfProcesses ) )
    {
        coresPerNode = numberOfProcesses;
    }

    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    // Setup the problem
    IH.setupProblem ( dataFile, problemFolder, coresPerNode );

    // Display problem information
    if ( verbose )
    {
        IH.showMe();
    }

    // TODO: Electro-mechanical-fluid coupling algorithm
    exitFlag = IH.solveProblem ( referenceSolution, tolerance );

#ifdef HAVE_MPI
    if ( rank == 0 )
    {
        std::cout << "MPI Finalization" << std::endl;
    }
    MPI_Finalize();
#endif

    return exitFlag;
}




