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

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/multiscale/solver/MultiscaleSolver.hpp>

#include <lifev/heart/solver/HeartETAMonodomainSolver.hpp>
#include <lifev/heart/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/heart/utility/HeartUtility.hpp>

using namespace LifeV;
using namespace Multiscale;

typedef VectorEpetra 			vector_Type;
typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
typedef RegionMesh<LinearTetra> mesh_Type;
typedef boost::function < Real ( const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/ ) > function_Type;
typedef IonicMinimalModel								ionicModel_Type;
typedef boost::shared_ptr<ionicModel_Type> 				ionicModelPtr_Type;
typedef HeartETAMonodomainSolver< mesh_Type, ionicModel_Type > monodomainSolver_Type;
typedef boost::shared_ptr< monodomainSolver_Type > monodomainSolverPtr_Type;
typedef boost::function < Real (const Real& /*t*/,
                                const Real &   x,
                                const Real &   y,
                                const Real& /*z*/,
                                const ID&   /*i*/ ) >   function_Type;


using namespace LifeV;

Real smoothing (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    return ( 0.5 + 0.5 * ( std::tanh ( - ( z - 50 ) / 5.0 ) ) );

}



class HeartApplication : public MultiscaleSolver
{
public:

    bool
    solveProblem ( )
    {

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

        //********************************************//
        // Import parameters from an xml list. Use    //
        // Teuchos to create a list from a given file //
        // in the execution directory.                //
        //********************************************//

        if ( M_comm->MyPID() == 0 )
        {
            std::cout << "\nImporting electrophysiology model parameter lists...";
        }
        Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << " Done!" << endl;
        }

        //********************************************//
        // Creates a new model object representing the//
        // model from Aliev and Panfilov 1996.  The   //
        // model input are the parameters. Pass  the  //
        // parameter list in the constructor          //
        //********************************************//
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << "\nBuilding the ionic model ... ";
        }
        ionicModelPtr_Type  IonicModel ( new ionicModel_Type () );
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << " Done!" << endl;
        }

        //********************************************//
        // In the parameter list we need to specify   //
        // the mesh name and the mesh path.           //
        //********************************************//
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << "\nImporting the mesh for the elctrophysiology ... ";
        }
        std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
        std::string meshPath = monodomainList.get ("mesh_path", "./");
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << " Done!" << endl;
        }

        //********************************************//
        // We need the GetPot datafile for to setup   //
        // the preconditioner.                        //
        //********************************************//
//        GetPot command_line (argc, argv);
//        const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
        GetPot dataFile;

        //********************************************//
        // We create three solvers to solve with:     //
        // 1) Operator Splitting method               //
        // 2) Ionic Current Interpolation             //
        // 3) State Variable Interpolation            //
        //********************************************//
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << "Building Monodomain Solvers... ";
        }

        monodomainSolverPtr_Type electro ( new monodomainSolver_Type ( meshName, meshPath, dataFile, IonicModel ) );
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << " Done... ";
        }

        //********************************************//
        // Setting up the initial condition form      //
        // a given function.                          //
        //********************************************//
        if (M_comm->MyPID() == 0 )
        {
            cout << "\nSetting initial condition for the electrophysiology...:  " ;
        }

        HeartUtility::setValueOnBoundary( *(electro -> potentialPtr() ), electro -> fullMeshPtr(), 1.0, 43 );
        HeartUtility::setValueOnBoundary( *(electro -> potentialPtr() ), electro -> fullMeshPtr(), 1.0, 45 );

        function_Type f = &smoothing;
        vectorPtr_Type smoother( new vector_Type( electro -> potentialPtr() -> map() ) );
        electro -> feSpacePtr() -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *smoother , 0);
        (*smoother) *= *(electro -> potentialPtr() );
        electro -> setPotentialPtr(smoother);

        //setting up initial conditions
        * ( electro -> globalSolution().at (1) ) = 1.0;
        * ( electro -> globalSolution().at (2) ) = 1.0;
        * ( electro -> globalSolution().at (3) ) = 0.021553043080281;

        if ( M_comm->MyPID() == 0 )
        {
            cout << "Done! \n" ;
        }

        //********************************************//
        // Importing Fiber direction                  //
        //********************************************//
        if ( M_comm->MyPID() == 0 )
        {
            cout << "\nImporting fibers:  " ;
        }

        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
        ( new FESpace< mesh_Type, MapEpetra > ( electro -> localMeshPtr(), "P1", 3, electro -> commPtr() ) );

        boost::shared_ptr<VectorEpetra> fiber ( new VectorEpetra ( Space3D -> map() ) );
    	std::string nm = monodomainList.get("fiber_file","FiberDirection") ;
        HeartUtility::importFibers( fiber, nm, electro -> localMeshPtr() );
        electro -> setFiberPtr(fiber);
        //delete Space3D;
        if ( M_comm->MyPID() == 0 )
        {
            cout << "Done! \n" ;
        }

        //********************************************//
        // Setting up the time data                   //
        //********************************************//
        electro -> setParameters ( monodomainList );

        //********************************************//
        // Create the global matrix: mass + stiffness //
        //********************************************//
        if ( M_comm->MyPID() == 0 )
        {
            cout << "\nSetup operators:  " ;
        }

        electro -> setupLumpedMassMatrix();
        electro -> setupStiffnessMatrix();
        electro -> setupGlobalMatrix();
        if ( M_comm->MyPID() == 0 )
        {
            cout << "Done! \n" ;
        }


        //********************************************//
        // Setting up the SVI solver                  //
        //********************************************//
        //********************************************//
        // Creating exporters to save the solution    //
        //********************************************//
        ExporterHDF5< RegionMesh <LinearTetra> > exporterElectro;

        electro -> setupExporter ( exporterElectro, "Electrophysiology" );

        electro -> exportSolution ( exporterElectro, 0);
        //********************************************//
        // Creating exporters to save the solution    //
        //********************************************//

        if ( M_comm->MyPID() == 0 )
        {
            cout << "\nstart solving:  " ;
        }

        Real saveStep = monodomainList.get ("saveStep", 1.0);

        //********************************************//
        // //setup activation;                        //
        //********************************************//

        Real maxCalciumLikeVariable = 0.838443;
        Real minCalciumLikeVariable = 0.021553;

        ExporterHDF5< RegionMesh <LinearTetra> > exporterActivation;

        vectorPtr_Type gammaf( new  vector_Type( *( electro -> globalSolution().at(3) ) ) );

        exporterActivation.setMeshProcId ( electro -> localMeshPtr(), M_comm -> MyPID() );
        //exporterActivation.exportPID (electro -> localMeshPtr(), M_comm );
        exporterActivation.setPrefix ("gammaf");
        exporterActivation.addVariable ( ExporterData<mesh_Type>::ScalarField,  "gamma_f", electro -> feSpacePtr(), gammaf, UInt (0) );
        exporterActivation.postProcess(0.0);

        Real beta = -0.3;
        Real electroDT = electro -> timeStep();
//        for ( ; M_globalData->dataTime()->canAdvance(); M_globalData->dataTime()->updateTime() )
//        {
        Real time = 0.;
        for ( ; time < electro -> endTime();  )
        {
        	electro   -> solveOneSplittingStep( exporterElectro, M_globalData -> dataTime() -> time() );
        	( *gammaf ) = *( electro -> globalSolution().at(3) );
        	HeartUtility::rescaleVector(*gammaf,minCalciumLikeVariable,maxCalciumLikeVariable,beta);
        	time += electroDT;
//            // Global chrono start
//            globalChrono.start();
//
//            if ( M_comm->MyPID() == 0 )
//            {
//                std::cout << std::endl;
//                std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
//                std::cout << "                   INTEGRATED HEART SIMULATOR                " << std::endl;
//                std::cout << "             time = " << M_globalData->dataTime()->time() << " s;"
//                        << " time step number = " << M_globalData->dataTime()->timeStepNumber() << std::endl;
//                std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
//            }
//
//            //********************************************//
//            // Solving the electrophysiology system      //
//            //********************************************//
//
//            //********************************************//
//            // Computing activation and passing it to FSI //
//            //********************************************//
//
//            //********************************************//
//            // Updating the FSI                           //
//            //********************************************//
//
//            buildUpdateChrono.start();
//            if ( M_globalData->dataTime()->isFirstTimeStep() )
//            {
//                M_model->buildModel();
//            }
//            else
//            {
//                M_model->updateModel();
//            }
//            buildUpdateChrono.stop();
//
//            //********************************************//
//            // Solving the activated FSI in the ventricle //
//            //********************************************//
//
//            solveChrono.start();
//            M_model->solveModel();
//            solveChrono.stop();
//
//            updateSolutionChrono.start();
//            M_model->updateSolution();
//            updateSolutionChrono.stop();
//
//            //********************************************//
//            // Exporting the FSI solution                 //
//            //********************************************//
//
//            saveChrono.start();
//            if ( M_globalData->dataTime()->timeStepNumber() % multiscaleSaveEachNTimeSteps == 0 || M_globalData->dataTime()->isLastTimeStep() )
//            {
//                M_model->saveSolution();
//            }
//            saveChrono.stop();
//
//            // Global chrono stop
//            globalChrono.stop();
//
//            // Compute time step time
//            timeStepTime = globalChrono.globalDiff ( *M_comm );
//
//            // Updating total simulation time
//            totalSimulationTime += timeStepTime;
            totalSimulationTime = time;

        	//
            if ( M_comm->MyPID() == 0 )
            {
//                std::cout << " MS-  Total iteration time:                    " << timeStepTime << " s" << std::endl;
                std::cout << " MS-  Total simulation time:                   " << totalSimulationTime << " s" << std::endl;
            }
//
//            // Save CPU time
//            saveCPUTime ( timeStepTime, buildUpdateChrono.globalDiff ( *M_comm ), solveChrono.globalDiff ( *M_comm ),
//                    updateSolutionChrono.globalDiff ( *M_comm ), saveChrono.globalDiff ( *M_comm ) );
        }

        exporterElectro.closeFile();
        exporterActivation.closeFile();
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
    UInt coresPerNode         = commandLine.follow (  1, 2, "-ns", "--nodesize" );

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
    exitFlag = IH.solveProblem ( );

#ifdef HAVE_MPI
    if ( rank == 0 )
    {
        std::cout << "MPI Finalization" << std::endl;
    }
    MPI_Finalize();
#endif

    return exitFlag;
}




