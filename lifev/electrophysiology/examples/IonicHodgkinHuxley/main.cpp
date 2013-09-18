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
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroIonicSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

//#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicHodgkinHuxley.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
// ---------------------------------------------------------------
// In order to use the ETA framework, a special version of the
// FESpace structure must be used. It is called ETFESpace and
// has basically the same role as the FESpace.
// ---------------------------------------------------------------

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>


// ---------------------------------------------------------------
// The most important file to include is the Integrate.hpp file
// which contains all the definitions required to perform the
// different integrations.
// ---------------------------------------------------------------

//#include <lifev/eta/expression/Integrate.hpp>
//
//#include <lifev/eta/expression/ExpressionDot.hpp>


using std::cout;
using std::endl;
using namespace LifeV;

//Initial condition: spherical from pacing site
Real Stimulus2 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    Real pacingSite_X = monodomainList.get ("pacingSite_X", 0.);
    Real pacingSite_Y = monodomainList.get ("pacingSite_Y", 0.);
    Real pacingSite_Z = monodomainList.get ("pacingSite_Z", 1. );
    Real stimulusRadius = 0.1; // monodomainList.get ("stimulusRadius", 0.1);

    if (  ( ( x - pacingSite_X ) * ( x - pacingSite_X ) +  ( y - pacingSite_Y ) * ( y - pacingSite_Y ) +  ( z - pacingSite_Z ) * ( z - pacingSite_Z )  )
                    <= ( stimulusRadius * stimulusRadius ) )
        return 1.0;
    else
        return 0.0;
}

//Initiation of a plane wave
Real PlaneWave (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if ( x<= 0.05 )
        return 1.0;
    else if( x<= 0.1)
        return 1.0*( 0.1 - x )/(0.05);
    else
        return 0.0;
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

    typedef RegionMesh<LinearTetra>            mesh_Type;
    typedef boost::shared_ptr<mesh_Type>       meshPtr_Type;
    typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;
    typedef VectorEpetra                       vector_Type;
    typedef FESpace< mesh_Type, MapEpetra >    feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>    feSpacePtr_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& z,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, IonicHodgkinHuxley >        monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    //   Teuchos::ParameterList APParameterList = *( Teuchos::getParametersFromXmlFile( "AlievPanfilovParameters.xml" ) );
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
        std::cout << "Building Constructor for Aliev Panfilov Model with parameters ... ";
    }
    boost::shared_ptr<IonicHodgkinHuxley>  model ( new IonicHodgkinHuxley() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    // +-----------------------------------------------+
   // |               Loading the mesh                |
   // +-----------------------------------------------+
   if ( Comm->MyPID() == 0 )
   {
       std::cout << std::endl << "[Loading the mesh]" << std::endl;
   }

   meshPtr_Type fullMeshPtr ( new mesh_Type ( Comm ) );

   VectorSmall<3> meshDim;
   meshDim[0] =  monodomainList.get ("meshDim_X", 10 );
   meshDim[1] =  monodomainList.get ("meshDim_Y", 10 );
   meshDim[2] =  monodomainList.get ("meshDim_Z", 10 );
   VectorSmall<3> domain;
   domain[0] =  monodomainList.get ("domain_X", 1. );
   domain[1] =  monodomainList.get ("domain_Y", 1. );
   domain[2] =  monodomainList.get ("domain_Z", 1. );

   // Building the mesh from the source
   regularMesh3D ( *fullMeshPtr,
                   1,
                   meshDim[0], meshDim[1], meshDim[2],
                   false,
                   domain[0], domain[1], domain[2],
                   0.0, 0.0, 0.0 );

   if ( Comm->MyPID() == 0 )
   {
       std::cout << "Mesh size  : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).maxH << std::endl;
   }
   if ( Comm->MyPID() == 0 )
   {
       std::cout << "Partitioning the mesh ... " << std::endl;
   }
   meshPtr_Type meshPtr;
   {
       MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, Comm );
       meshPtr = meshPart.meshPartition();
   }
   fullMeshPtr.reset(); //Freeing the global mesh to save memory

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

    monodomainSolverPtr_Type splitting ( new monodomainSolver_Type ( dataFile, model, meshPtr ) );
    const feSpacePtr_Type FESpacePtr =  splitting->feSpacePtr(); //FE Space


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

//    function_Type f = &Stimulus;
    function_Type f = &Stimulus2;
    splitting -> setPotentialFromFunction ( f );

    //setting up initial conditions
    * ( splitting -> globalSolution().at (1) ) = 0.052932485257250;
    * ( splitting -> globalSolution().at (2) ) = 0.317676914060697;
    * ( splitting -> globalSolution().at (3) ) = 0.596120753508460;

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
    VectorSmall<3> fibers;
    fibers[0] =  monodomainList.get ("fiber_X", std::sqrt (2) / 2.0 );
    fibers[1] =  monodomainList.get ("fiber_Y", std::sqrt (2) / 2.0 );
    fibers[2] =  monodomainList.get ("fiber_Z", 0.0 );

    splitting ->setupFibers (fibers);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    splitting -> setupLumpedMassMatrix();
    splitting -> setupStiffnessMatrix();
    splitting -> setupGlobalMatrix();

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;

    string filenameSplitting =  monodomainList.get ("OutputFile", "MinMod" );
    filenameSplitting += "Splitting";

    splitting -> setupPotentialExporter ( exporterSplitting, filenameSplitting );

    splitting -> exportSolution ( exporterSplitting, 0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

    Real Savedt = monodomainList.get ("saveStep", 0.1);
    Real timeStep = monodomainList.get ("timeStep", 0.01);
    Real endTime = monodomainList.get ("endTime", 10.);
    Real initialTime = monodomainList.get ("initialTime", 0.);

    vectorPtr_Type previousPotential0Ptr( new vector_Type ( FESpacePtr->map() ) );
    *(previousPotential0Ptr) = *(splitting->globalSolution().at(0));

    int iter((Savedt / timeStep)+ 1e-9);
    int nbTimeStep (1);
    int k(0);
    if (endTime > timeStep) {
        for (Real t = initialTime; t < endTime;){

            t += timeStep;
            k++;
            if (nbTimeStep==1) {
                    splitting->solveOneReactionStepFE();
                    (*(splitting->rhsPtrUnique())) *= 0;
                    splitting->updateRhs();
                    splitting->solveOneDiffusionStepBE();
                    splitting->exportSolution(exporterSplitting, t);
            }else{
                *(previousPotential0Ptr) = *(splitting->globalSolution().at(0));
                splitting->solveOneReactionStepFE(2);
                (*(splitting->rhsPtrUnique())) *= 0;
                splitting->updateRhs();
                splitting->solveOneDiffusionStepBDF2(previousPotential0Ptr);
                splitting->solveOneReactionStepFE(2);
                if (k % iter == 0) {
                    splitting->exportSolution(exporterSplitting, t);
                }
            }
            nbTimeStep++;
        }
    }

    exporterSplitting.closeFile();


    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    splitting -> exportFiberDirection();


    if ( Comm->MyPID() == 0 )
    {
        cout << "\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
