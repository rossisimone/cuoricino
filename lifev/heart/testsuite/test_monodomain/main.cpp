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
    @brief 0D test with the Negroni Lascano model of 1996.

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

#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/heart/solver/HeartETAMonodomainSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/heart/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"


using namespace LifeV;

Real smoothing (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    return ( 0.5 + 0.5 * ( std::tanh ( - ( z - 50 ) / 5.0 ) ) );
}


Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    //********************************************//
    // Starts the chronometer.                    //
    //********************************************//
    //  LifeChrono chronoinitialsettings;
    //  chronoinitialsettings.start();

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef HeartETAMonodomainSolver< mesh_Type, IonicMinimalModel >        monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
    typedef VectorEpetra				vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;


    LifeChrono chronoinitialsettings;

    if ( Comm->MyPID() == 0 )
      	chronoinitialsettings.start();

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
    // model from Aliev and Panfilov 1996.  The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building Constructor for Minimal Model with parameters ... ";
    }
    boost::shared_ptr<IonicMinimalModel>  model ( new IonicMinimalModel() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


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
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // We create three solvers to solve with:     //
    // 1) Operator Splitting method               //
    // 2) Ionic Current Interpolation             //
    // 3) State Variable Interpolation            //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    monodomainSolverPtr_Type splitting ( new monodomainSolver_Type ( meshName, meshPath, dataFile, model ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }


    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nInitializing potential:  " ;
    }

    HeartUtility::setValueOnBoundary( *(splitting -> potentialPtr() ), splitting -> fullMeshPtr(), 1.0, 43 );
    HeartUtility::setValueOnBoundary( *(splitting -> potentialPtr() ), splitting -> fullMeshPtr(), 1.0, 45 );

    function_Type f = &smoothing;
    vectorPtr_Type smoother( new vector_Type( splitting -> potentialPtr() -> map() ) );
    splitting -> feSpacePtr() -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *smoother , 0);
    (*smoother) *= *(splitting -> potentialPtr() );
    splitting -> setPotentialPtr(smoother);

    //setting up initial conditions
    * ( splitting -> globalSolution().at (1) ) = 1.0;
    * ( splitting -> globalSolution().at (2) ) = 1.0;
    * ( splitting -> globalSolution().at (3) ) = 0.021553043080281;

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the time data                   //
    //********************************************//
    splitting -> setParameters ( monodomainList );

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nImporting fibers:  " ;
    }

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( splitting -> localMeshPtr(), "P1", 3, splitting -> commPtr() ) );

    boost::shared_ptr<VectorEpetra> fiber ( new VectorEpetra ( Space3D -> map() ) );
	std::string nm = monodomainList.get("fiber_file","FiberDirection") ;
    HeartUtility::importFibers( fiber, nm, splitting -> localMeshPtr() );
    splitting -> setFiberPtr(fiber);

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

    splitting -> setupLumpedMassMatrix();
    splitting -> setupStiffnessMatrix();
    splitting -> setupGlobalMatrix();
    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }
    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;

    splitting -> setupExporter ( exporterSplitting, monodomainList.get ("OutputFile", "Splitting") );

    splitting -> exportSolution ( exporterSplitting, 0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

    Real saveStep = monodomainList.get ("saveStep", 1.0);

//    splitting   -> solveSplitting ( exporterSplitting, saveStep );


    monodomainSolver_Type::vectorPtr_Type dtVec ( new VectorEpetra ( splitting->feSpacePtr() -> map(), LifeV::Unique ) );
    ExporterHDF5<mesh_Type> Exp;
    Exp.setMeshProcId ( splitting -> localMeshPtr(), splitting -> commPtr() -> MyPID() );
    Exp.setPrefix (monodomainList.get ("OutputTimeSteps", "TimeSteps"));
    Exp.addVariable ( ExporterData<mesh_Type>::ScalarField,  "dt", splitting->feSpacePtr(), dtVec, UInt (0) );

    Real dt = monodomainList.get ("timeStep", 0.1);
    Real TF = monodomainList.get ("endTime", 150.0);
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;
    Int meth = monodomainList.get ("meth", 1.0);
    Real dt_min = dt/50.0;
    Int k(0),j(0);
    Int nodes;

    Real timeReac = 0.0;
    Real timeDiff = 0.0;
    LifeChrono chrono;

    for ( Real t = 0.0; t < TF; )
    {
        t = t + dt;

        chrono.start();
        if(meth==1)
            splitting->solveOneReactionStepROS3P(dtVec, dt_min);
        else
           	splitting->solveOneReactionStepFE();
        chrono.stop();
        timeReac += chrono.diff();

        (*splitting->rhsPtrUnique()) *= 0.0;
        splitting->updateRhs();

        chrono.start();
        splitting->solveOneDiffusionStepBE();
        chrono.stop();
        timeDiff += chrono.diff();

        if( k % iter == 0 )
        {
           	splitting -> exportSolution (exporterSplitting, t);
           	Exp.postProcess (t);
        }

        nodes = dtVec->epetraVector().MyLength();
        j = dtVec->blockMap().GID (0);
        dt_min = (*dtVec)[j];
        for(int i=1; i<nodes; i++)
        {
        	j = dtVec->blockMap().GID (i);
        	if(dt_min>(*dtVec)[j])
        		dt_min = (*dtVec)[j];
        }

        k++;

        if ( Comm->MyPID() == 0 )
        	std::cout<<"\n\n\nActual time : "<<t<<std::endl<<std::endl<<std::endl;

    }

    exporterSplitting.closeFile();
    Exp.closeFile();


    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    splitting -> exportFiberDirection();


    if ( Comm->MyPID() == 0 )
    {
    	chronoinitialsettings.stop();
    	std::cout << "\n\n\nTotal lapsed time : " << chronoinitialsettings.diff() << std::endl;
    	std::cout<<"Diffusion time : "<<timeDiff<<std::endl;
    	std::cout<<"Reaction time : "<<timeReac<<std::endl;
        cout << "\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
