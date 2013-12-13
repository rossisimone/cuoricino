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
    @brief

    @date
    @author
    @contributor
    @mantainer
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

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicLuoRudyI.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher06.hpp>
#include <lifev/electrophysiology/util/CardiacStimulusPMJ.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <sys/stat.h>

using namespace LifeV;

Int main ( Int argc, char** argv )
{

    typedef RegionMesh<LinearTetra>                                  mesh_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >            function_Type;

    typedef ElectroIonicModel    ionicModel_Type;
    typedef boost::shared_ptr<ionicModel_Type>                       ionicModelPtr_Type;
    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type > monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >               monodomainSolverPtr_Type;

    typedef VectorEpetra                                             vector_Type;
    typedef boost::shared_ptr<vector_Type>                           vectorPtr_Type;

    bool verbose = false;

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );

    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "work", 2, "-w", "--work" );

    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( Comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    //********************************************//
    // Starts the chronometer.                    //
    //********************************************//


    LifeChrono chronoinitialsettings;

    if ( Comm->MyPID() == 0 )
    {
        chronoinitialsettings.start();
    }

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    GetPot command_line (argc, argv);
    const string monodomain_datafile_name = command_line.follow ("MonodomainSolverParamList.xml", 2, "-f", "--file");
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );

    std::string ionic_model ( monodomainList.get ("ionic_model", "minimalModel") );

    ionicModelPtr_Type  model;
    if ( ionic_model == "LuoRudyI" )
    {
        model.reset ( new IonicLuoRudyI() );
    }
    else if ( ionic_model == "TenTusscher06")
    {
        model.reset (new IonicTenTusscher06() );
    }
    else
    {
        model.reset ( new IonicMinimalModel() );
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    std::string meshName = command_line.follow ("default", 2, "-m", "--model") + ".mesh";

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    const string preconditioner_datafile_name = command_line.follow ("MonodomainSolverPreconditioner", 2, "-p", "--prec");
    GetPot dataFile (preconditioner_datafile_name);

    //********************************************//
    // We create three solvers to solve with:     //
    // 1) Operator Splitting method               //
    // 2) Ionic Current Interpolation             //
    // 3) State Variable Interpolation            //
    //********************************************//
    monodomainSolverPtr_Type solver ( new monodomainSolver_Type ( meshName, problemFolder, dataFile, model ) );

    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//

    // Initial pacing
    solver -> initializePotential();
    solver -> initializeAppliedCurrent();
    solver -> setInitialConditions();

    CardiacStimulusPMJ stimulus;
    stimulus.setRadius ( monodomainList.get ("applied_current_radius", 0.2) );
    stimulus.setTotalCurrent ( monodomainList.get ("applied_total_current", 0.2) );

    std::string purkinjeFile = problemFolder + command_line.follow ("default", 2, "-m", "--model") + "_activation.txt";
    stimulus.setPMJFromFile (purkinjeFile);

    solver -> setAppliedCurrentFromCardiacStimulus (stimulus, 0.0);

    //********************************************//
    // Setting up the time data                   //
    //********************************************//
    solver -> setParameters ( monodomainList );

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( solver -> localMeshPtr(), "P1", 3, solver -> commPtr() ) );

    boost::shared_ptr<VectorEpetra> fiber ( new VectorEpetra ( Space3D -> map() ) );
    std::string nm = problemFolder + monodomainList.get ("fiber_file", "FiberDirection") ;
    HeartUtility::importFibers ( fiber, nm, solver -> localMeshPtr() );

    solver -> setFiberPtr (fiber);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    bool lumpedMass = monodomainList.get ("lumped_mass", false);
    if ( lumpedMass)
    {
        solver -> setupLumpedMassMatrix();
    }
    else
    {
        solver -> setupMassMatrix();
    }


    solver -> setupStiffnessMatrix ( solver -> diffusionTensor() );
    solver -> setupGlobalMatrix();

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporter;

    solver -> setupExporter ( exporter, monodomainList.get ("outputFile", command_line.follow ("default", 2, "-m", "--model") + "_solution") );
    exporter.setPostDir (problemFolder);
    solver -> exportSolution ( exporter, 0);

    //********************************************//
    // Activation time                            //
    //********************************************//
    vectorPtr_Type activationTimeVector ( new vector_Type ( solver -> potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;

    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId (solver -> localMeshPtr(), solver -> commPtr() ->MyPID() );
    activationTimeExporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Activation Time",
                                        solver -> feSpacePtr(), activationTimeVector, UInt (0) );
    activationTimeExporter.setPrefix (command_line.follow ("default", 2, "-m", "--model") + "_activationTime");
    activationTimeExporter.setPostDir (problemFolder);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

    Real dt = monodomainList.get ("timeStep", 0.1);
    Real TF = monodomainList.get ("endTime", 150.0);
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;
    Int k (0);

    Real timeReac = 0.0;
    Real timeDiff = 0.0;
    Real timeReacDiff = 0.0;
    LifeChrono chrono;

    std::string solutionMethod = monodomainList.get ("solutionMethod", "splitting");


    for ( Real t = 0.0; t < TF; )
    {
        solver -> setAppliedCurrentFromCardiacStimulus ( stimulus, t );

        if ( solutionMethod == "splitting" )
        {
            chrono.reset();
            chrono.start();
            solver->solveOneReactionStepFE( );
            chrono.stop();

            timeReac += chrono.globalDiff ( *Comm );

            (*solver->rhsPtrUnique() ) *= 0.0;
            solver->updateRhs();

            chrono.reset();
            chrono.start();
            solver->solveOneDiffusionStepBE();
            chrono.stop();
            timeDiff += chrono.globalDiff ( *Comm );
        }
        else if ( solutionMethod == "ICI" )
        {
            chrono.reset();
            chrono.start();
            solver -> solveOneStepGatingVariablesFE();
            solver -> solveOneICIStep();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }
        else if ( solutionMethod == "Mixed" )
        {
            chrono.reset();
            chrono.start();
            solver -> solveOneStepGatingVariablesFE();
            solver -> solveOneMixedStep();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }
        else if ( solutionMethod == "SVI" )
        {
            chrono.reset();
            chrono.start();
            solver -> solveOneStepGatingVariablesFE();
            solver -> solveOneSVIStep();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }

        //register activation time
        k++;
        t = t + dt;
        solver -> registerActivationTime (*activationTimeVector, t, 1.0);

        if ( k % iter == 0 )
        {
            solver -> exportSolution (exporter, t);
        }
    }

    exporter.closeFile();
    activationTimeExporter.postProcess (0);
    activationTimeExporter.closeFile();

    chronoinitialsettings.stop();

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();

    return ( EXIT_SUCCESS );
}
