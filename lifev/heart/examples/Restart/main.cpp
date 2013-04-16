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
#include <lifev/core/filter/Exporter.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/heart/solver/HeartETAMonodomainSolver.hpp>
#include <lifev/heart/solver/HeartIonicSolver.hpp>
#include <lifev/core/util/HeartUtility.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/heart/solver/IonicModels/IonicAlievPanfilov.hpp>
#include <lifev/heart/solver/IonicModels/IonicMinimalModel.hpp>
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


Real Stimulus (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    //      return ( 0.5 + 0.5 * ( std::tanh( - 300 * ( ( x - 0.4 ) * ( x - 0.6 ) + ( y - 0.4 ) * ( y - 0.6 ) ) ) ) );
    return ( 0.5 + 0.5 * ( std::tanh ( - 10 * ( ( x - 30 ) * ( x - 50 ) ) ) ) );
    //      return ( exp( -0.001 * ( ( x - 50 ) * ( x - 50 ) + ( y - 37 ) * ( y - 37 ) + ( z - 42 ) * ( z - 42 )  ) ) );
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
    typedef boost::shared_ptr< RegionMesh<LinearTetra> >	   meshPtr_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef  IonicAlievPanfilov                             ionicModel_Type;
    typedef HeartETAMonodomainSolver< mesh_Type, ionicModel_Type >      monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
    typedef ExporterData<mesh_Type> exporterData_Type;

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
    boost::shared_ptr<ionicModel_Type>  model ( new ionicModel_Type() );
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

    function_Type f = &Stimulus;
    splitting -> setPotentialFromFunction ( f );

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the time data                   //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\n\nsetting paramters:  " ;
    }

    splitting -> setParameters ( monodomainList );
    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\n\nsetting fibers:  " ;
    }

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( splitting -> localMeshPtr(), "P1", 3, splitting -> commPtr() ) );

    boost::shared_ptr<VectorEpetra> fiber ( new VectorEpetra ( Space3D -> map() ) );
    std::string fibersDirectory = monodomainList.get ("fiber_path", "./" );
    std::string fibersFile = monodomainList.get ("fiber_file", "fibers.dat" );

    HeartUtility::importFibers(fiber,fibersFile, fibersDirectory);

    splitting -> setFiberPtr(fiber);

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\n\nexporting fibers:  \n\n" ;
    }

    splitting -> exportFiberDirection();

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//


//    ExporterHDF5<mesh_Type> exp;
//    exp.setMeshProcId ( splitting -> localMeshPtr(), splitting -> commPtr() -> MyPID() );
//    exp.setPrefix ("FiberDirection");
//    exp.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", Space3D, splitting -> fiberPtr(), UInt (0) );
//    exp.postProcess (0);
//    exp.closeFile();

    if ( Comm->MyPID() == 0 )
    {
        cout << "\n\nDone! \n" ;
    }

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        cout << "\n Preparing the restarter ..." ;
    }

    std::string name ( "FiberDirection" );
    monodomainSolver_Type::vectorPtr_Type newFiber ( new VectorEpetra ( Space3D -> map(), LifeV::Unique ) );



//    exporterData_Type impData (exporterData_Type::VectorField, "fibers.00000", Space3D,
//                               newFiber, UInt (0), exporterData_Type::UnsteadyRegime);
//
//    if ( Comm->MyPID() == 0 )
//    {
//        cout << "Done! \n" ;
//    }
//
//    if ( Comm->MyPID() == 0 )
//    {
//        cout << "\n Restarting ..." ;
//    }
//
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >  hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                  hdf5FilterPtr_Type;
//
//    //    filterPtr_Type importer( new hdf5Filter_Type(dataFile, name) );
//    filterPtr_Type importer ( new hdf5Filter_Type() );
//
//    importer -> setMeshProcId ( splitting -> localMeshPtr(), Comm -> MyPID() );
    std::string const Name = "fiberdirection";
//    importer-> setPrefix (name);
//    if ( Comm->MyPID() == 0 )
//    {
//        cout << "Done! \n" ;
//    }
//
//    importer -> readVariable (impData);
//
//    importer -> closeFile();

    HeartUtility::importFibers(newFiber, name, splitting -> localMeshPtr() );
//    HeartUtility::importFibers(newFiber, splitting -> localMeshPtr() );

    ExporterHDF5<mesh_Type> Exp;
    Exp.setMeshProcId ( splitting -> localMeshPtr(), splitting -> commPtr() -> MyPID() );
    Exp.setPrefix (Name);
//    (*newFiber) *= 100.0;
    Exp.addVariable ( ExporterData<mesh_Type>::VectorField,  "ciccia", Space3D, newFiber, UInt (0) );
    Exp.postProcess (0);
    Exp.closeFile();



    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;

    splitting -> setupExporter ( exporterSplitting, "Splitting" );

    splitting -> exportSolution ( exporterSplitting, 0);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    splitting -> setupLumpedMassMatrix();
    splitting -> setupStiffnessMatrix();
    splitting -> setupGlobalMatrix();

    splitting -> solveICI (exporterSplitting);
    exporterSplitting.closeFile();


    //================================================
    //
    //================================================
    std::string sol ( "Splitting" );
    monodomainSolver_Type::vectorPtr_Type newSol ( new VectorEpetra ( splitting -> feSpacePtr() -> map(), LifeV::Unique ) );
    exporterData_Type ImportData (exporterData_Type::ScalarField, "Variable0.00001", splitting -> feSpacePtr(),
                                  newSol, UInt (0), exporterData_Type::UnsteadyRegime);

    //================================================
    //
    //================================================
    filterPtr_Type Imp ( new hdf5Filter_Type() );
    //
    Imp -> setMeshProcId ( splitting -> localMeshPtr(), Comm -> MyPID() );
    Imp-> setPrefix (sol);
    //
    Imp -> readVariable (ImportData);
    Imp -> closeFile();
    //================================================
    //
    //================================================
    splitting -> setPotentialPtr (newSol);
    ( * ( splitting -> globalSolution().at (1) ) ) *= 0;

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporter2;

    splitting -> setupExporter ( exporter2, "exp2" );

    splitting -> exportSolution ( exporter2, 0);

    splitting -> solveICI (exporter2);
    exporter2.closeFile();

    //================================================
    //
    //================================================
    typedef VectorEpetra                           vector_Type;
    typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;

    //    boost::shared_ptr<mesh_Type> lmesh( new mesh_Type( Comm ) );
    //    lmesh ->se
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > space (new FESpace<mesh_Type, MapEpetra> (  splitting -> localMeshPtr(), "P1", 1, Comm ) );


    //    vectorPtr_Type vec( new vector_Type( splitting -> feSpacePtr() -> map(), Unique ) );
    vectorPtr_Type vec ( new vector_Type ( space -> map(), Unique ) );
//    meshPtr_Type mesh1 ( new mesh_Type ( Comm ) );
//    meshPtr_Type fullMesh1 ( new mesh_Type ( Comm ) );
//    MeshUtility::fillWithFullMesh (mesh1, fullMesh1, meshName, meshPath);

    for ( UInt j (0); j < vec -> epetraVector().MyLength() ; ++j)
    {

        if ( splitting -> fullMeshPtr() -> point ( vec -> blockMap().GID (j) ).markerID() == 15 )
        {
            //if ( vec -> blockMap().LID ( vec -> blockMap().GID (j) ) != -1 )
            //{
                (*vec) ( vec -> blockMap().GID (j) ) = 1.0;
            //}
        }

        if ( splitting -> fullMeshPtr() -> point ( vec -> blockMap().GID (j) ).markerID() == 14 )
        {
            if ( vec -> blockMap().LID ( vec -> blockMap().GID (j) ) != -1 )
            {
                (*vec) ( vec -> blockMap().GID (j) ) = 2.0;
            }
        }
        if ( splitting -> fullMeshPtr() -> point ( vec -> blockMap().GID (j) ).markerID() == 13 )
        {
            if ( vec -> blockMap().LID ( vec -> blockMap().GID (j) ) != -1 )
            {
                (*vec) ( vec -> blockMap().GID (j) ) = 3.0;
            }
        }
        if ( splitting -> fullMeshPtr() -> point ( vec -> blockMap().GID (j) ).markerID() == 37 )
        {
            if ( vec -> blockMap().LID ( vec -> blockMap().GID (j) ) != -1 )
            {
                (*vec) ( vec -> blockMap().GID (j) ) = 4.0;
            }
        }
    }



    ExporterHDF5<mesh_Type> Exp1;
    Exp1.setMeshProcId ( splitting -> localMeshPtr(), splitting -> commPtr() -> MyPID() );
    Exp1.setPrefix ("pippo");
    Exp1.addVariable ( ExporterData<mesh_Type>::ScalarField,  "ciccia", splitting -> feSpacePtr(), vec, UInt (0) );
    Exp1.postProcess (0);
    Exp1.closeFile();



    if ( Comm->MyPID() == 0 )
    {
        cout << "\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
