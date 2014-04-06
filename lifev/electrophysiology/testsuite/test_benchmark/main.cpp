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
    @brief Electrophysiology benchmark of Niederer et al. 2011

    @date 01−2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

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


#include <lifev/electrophysiology/testsuite/test_benchmark/benchmarkUtility.hpp>



#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <sys/stat.h>

using namespace LifeV;

//This final time was computed using gcc on Ubuntu/Linux
//On Mac Os X 10.9 Mavercik with clang 503.0.38 the test has final time 40.
#define finalActivationTime  39.19

Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
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
    // Some typedefs                              //
    //********************************************//

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef ElectroIonicModel                                        ionicModel_Type;
    typedef boost::shared_ptr<ionicModel_Type>                       ionicModelPtr_Type;
    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type > monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >               monodomainSolverPtr_Type;

    typedef VectorEpetra                                             vector_Type;
    typedef boost::shared_ptr<vector_Type>                           vectorPtr_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    //********************************************//
    // Staring the chronometers                   //
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

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // Creates a new model object representing the//
    // ionic model.  The model parameters are the //
    // default ones.                              //
    //********************************************//
    std::string ionic_model ( monodomainList.get ("ionic_model", "minimalModel") );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nIonic_Model:" << ionic_model;
    }


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nBuilding Constructor for " << ionic_model << " Model with parameters ... ";
    }
    ionicModelPtr_Type  model;
    Real activationThreshold = BenchmarkUtility::chooseIonicModel(model, ionic_model, *Comm );



    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = monodomainList.get ("mesh_path", "./");

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot dataFile  (argc, argv);

    //********************************************//
    // We create monodomain solver.               //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    monodomainSolverPtr_Type solver ( new monodomainSolver_Type ( meshName, meshPath, dataFile , model ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " solver done... ";
    }

    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nInitializing potential and gating variables:  " ;
    }

    // Initial pacing
    solver -> initializePotential();
    solver -> initializeAppliedCurrent();
    solver -> setInitialConditions();


    //********************************************//
    // Set initial stimulus.                      //
    //********************************************//
    function_Type stimulus;
    BenchmarkUtility::setStimulus(stimulus, ionic_model);
    solver -> setAppliedCurrentFromFunction (stimulus, 0.0);

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the monodomain solver           //
    //********************************************//
    solver -> setParameters ( monodomainList );

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nSetting fibers:  " ;
    }

    VectorSmall<3> fibers;
    fibers[0] = 0.0;
    fibers[1] = 0.0;
    fibers[2] = 1.0;
    solver -> setupFibers ( fibers );

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nSetup operators:  " ;
    }

    bool lumpedMass = monodomainList.get ("LumpedMass", true);
    LifeChrono timer;
    matrixPtr_Type hlmass;
    if( lumpedMass)
    { 
        solver -> setLumpedMassMatrix(false);
	    MPI_Barrier (MPI_COMM_WORLD);
    	timer.start();
		solver -> setupMassMatrix();
		hlmass.reset(new matrix_Type( *(solver -> massMatrixPtr() ) ) );

		timer.stop();
		solver -> setLumpedMassMatrix(lumpedMass);
    }

    timer.reset();


    timer.start();
    MPI_Barrier (MPI_COMM_WORLD);
    solver -> setupMassMatrix();
    timer.stop();
    solver -> setupStiffnessMatrix ( solver -> diffusionTensor() );
    solver -> setupGlobalMatrix();
    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporter;

    solver -> setupExporter ( exporter, monodomainList.get ("OutputFile", "Solution") );
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
    activationTimeExporter.setPrefix ("ActivationTime");
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

    if ( Comm->MyPID() == 0 )
    {
    std::cout << "\nSolving the monodomain using " << solutionMethod;
    }

    for ( Real t = 0.0; t < TF; )
    {
        //********************************************//
        // Updating the applied current               //
        //********************************************//
        solver -> setAppliedCurrentFromFunction ( stimulus, t );

        //********************************************//
        // Solving using the operator splitting method//
        //********************************************//
        if ( solutionMethod == "splitting" )
        {
			chrono.reset();
			chrono.start();
			if(ionic_model != "MinimalModel")
				solver->solveOneReactionStepRL();
			else solver->solveOneReactionStepFE();
			chrono.stop();

			timeReac += chrono.globalDiff( *Comm );

			(*solver->rhsPtrUnique()) *= 0.0;
			solver->updateRhs();

			chrono.reset();
			chrono.start();
			solver->solveOneDiffusionStepBE();
			chrono.stop();
			timeDiff += chrono.globalDiff( *Comm );
		}
        //********************************************//
        // Solving using the L-ICI method             //
        //********************************************//
        else if( solutionMethod == "L-ICI" )
        {
			chrono.reset();
			chrono.start();
			if(ionic_model != "MinimalModel")
				solver -> solveOneStepGatingVariablesRL();
			else
				solver -> solveOneStepGatingVariablesFE();
			solver -> solveOneICIStep();
			chrono.stop();
			timeReacDiff += chrono.globalDiff( *Comm );
        }
        //********************************************//
        // Solving using the ICI method             //
        //********************************************//
        else if( solutionMethod == "ICI" )
        {
			chrono.reset();
			chrono.start();
			if(ionic_model != "MinimalModel")
				solver -> solveOneStepGatingVariablesRL();
			else
				solver -> solveOneStepGatingVariablesFE();
			solver -> solveOneICIStep(*hlmass);
			chrono.stop();
			timeReacDiff += chrono.globalDiff( *Comm );
        }
        //********************************************//
        // Solving using the SVI  method             //
        //********************************************//
        else if( solutionMethod == "SVI" )
        {
			chrono.reset();
			chrono.start();
			if(ionic_model != "MinimalModel" && ionic_model != "Fox")
				solver -> solveOneStepGatingVariablesRL();
			else
				solver -> solveOneStepGatingVariablesFE();
        	solver -> solveOneSVIStep();
			chrono.stop();
			timeReacDiff += chrono.globalDiff( *Comm );
        }

        //register activation time
        k++;
        t = t + dt;

        //********************************************//
        // Save the activation time                   //
        //********************************************//
        solver -> registerActivationTime (*activationTimeVector, t, activationThreshold);


        //********************************************//
        // export the solution                        //
        //********************************************//
        if ( k % iter == 0 )
        {
	        if ( Comm->MyPID() == 0 )
	        {
	            std::cout << "\nTime : " << t;
	        }
            solver -> exportSolution (exporter, t);
        }

    }

    Real normSolution = ( ( solver -> globalSolution().at (0) )->norm2() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\n2-norm of potential solution: " << normSolution;
    }

    //********************************************//
    // Closing exporters                          //
    //********************************************//

    exporter.closeFile();
    activationTimeExporter.postProcess (0);
    activationTimeExporter.closeFile();

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nExporting fibers ...  ";
    }

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    solver -> exportFiberDirection(problemFolder);


    //********************************************//
    // Destroy the solver                         //
    //********************************************//
    solver.reset();

    Real fullActivationTime = activationTimeVector -> maxValue();
    std::cout << std::setprecision(8) << fullActivationTime;
    if ( Comm->MyPID() == 0 )
    {
        chronoinitialsettings.stop();
        std::cout << "\n\n\nTotal lapsed time : " << chronoinitialsettings.diff() << std::endl;
        if ( solutionMethod == "splitting" )
        {
            std::cout << "Diffusion time : " << timeDiff << std::endl;
            std::cout << "Reaction time : " << timeReac << std::endl;
        }
        else if ( solutionMethod == "ICI" )
        {
            std::cout << "Solution time : " << timeReacDiff << std::endl;
        }
        else if ( solutionMethod == "SVI" )
        {
            std::cout << "Solution time : " << timeReacDiff << std::endl;
        }

        std::cout << "\n\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();


    //********************************************//
    // Checking if the test failed                //
    //********************************************//
    Real returnValue;

    Real err = std::abs (fullActivationTime - finalActivationTime) / std::abs(finalActivationTime);
    if ( err > 2.5e-2 )
    {
    	std::cout << "\nTest Failed: " <<  err <<"\n" << "\nSolution Norm: " <<  fullActivationTime << "\n";
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}

