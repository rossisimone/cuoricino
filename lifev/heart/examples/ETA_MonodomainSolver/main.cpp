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
#include <lifev/heart/solver/HeartIonicSolver.hpp>

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


Real Stimulus (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    return ( 0.5 + 0.5 * ( std::tanh ( - 300 * ( ( x - 0.4 ) * ( x - 0.6 ) + ( y - 0.4 ) * ( y - 0.6 ) ) ) ) );
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

    typedef HeartETAMonodomainSolver< mesh_Type, IonicAlievPanfilov >       monodomainSolver_Type;
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
    Teuchos::ParameterList APParameterList = * ( Teuchos::getParametersFromXmlFile ( "AlievPanfilovParameters.xml" ) );
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
    if ( Comm->MyPID() == 0 ) std::cout << "Building Constructor for Aliev Panfilov Model with parameters ... ";
	boost::shared_ptr<IonicAlievPanfilov>  model( new IonicAlievPanfilov( APParameterList ) );
	if ( Comm->MyPID() == 0 ) std::cout << " Done!" << endl;


	//********************************************//
	// In the parameter list we need to specify   //
	// the mesh name and the mesh path.           //
	//********************************************//
	std::string meshName = monodomainList.get("mesh_name","lid16.mesh");
	std::string meshPath = monodomainList.get("mesh_path","./");

	//********************************************//
	// We need the GetPot datafile for to setup   //
	// the preconditioner.                        //
	//********************************************//
	GetPot command_line(argc, argv);
	const string data_file_name = command_line.follow("data", 2, "-f", "--file");
	GetPot dataFile(data_file_name);

	//********************************************//
	// We create three solvers to solve with:     //
	// 1) Operator Splitting method               //
	// 2) Ionic Current Interpolation             //
	// 3) State Variable Interpolation            //
	//********************************************//
	if ( Comm->MyPID() == 0 )  std::cout << "Building Monodomain Solvers... ";

	monodomainSolverPtr_Type splitting( new monodomainSolver_Type( meshName, meshPath, dataFile, model ) );
	if ( Comm->MyPID() == 0 )  std::cout << " Splitting solver done... ";
	monodomainSolverPtr_Type ici( new monodomainSolver_Type( dataFile, model, splitting ->  meshPtr() ) );
	if ( Comm->MyPID() == 0 )  std::cout << " ICI solver done!\n. ";


	//********************************************//
	// Setting up the initial condition form      //
	// a given function.                          //
	//********************************************//
	if ( Comm->MyPID() == 0 )  cout << "\nInitializing potential:  " ;

	function_Type f = &Stimulus;
	splitting -> setPotentialFromFunction( f );
	ici -> copyPotential( splitting -> potentialPtr() );
	if ( Comm->MyPID() == 0 ) cout << "Done! \n" ;

	//********************************************//
	// Setting up the time data                   //
	//********************************************//
	splitting -> setParameters( monodomainList );
	ici -> setParameters( monodomainList );

	//********************************************//
	// Create a fiber direction                   //
	//********************************************//
	VectorSmall<3> fibers;
	fibers[0]=  monodomainList.get("fiber_X", std::sqrt(2) / 2.0 );
	fibers[1]=  monodomainList.get("fiber_Y", std::sqrt(2) / 2.0 );
	fibers[2]=  monodomainList.get("fiber_Z", 0.0 );

	splitting ->setupFibers(fibers);
	ici -> setFiberPtr( splitting -> fiberPtr() );

	//********************************************//
	// Create the global matrix: mass + stiffness //
	//********************************************//
	splitting -> setupMassMatrix();
	splitting -> setupStiffnessMatrix();
	splitting -> setupGlobalMatrix();

	ici -> setupMassMatrix();
	ici -> setupStiffnessMatrix();
	ici -> setupGlobalMatrix();

	//********************************************//
	// Setting up the SVI solver                  //
	//********************************************//
	monodomainSolverPtr_Type svi( new monodomainSolver_Type( *ici ) );

	//********************************************//
	// Creating exporters to save the solution    //
	//********************************************//
	ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;
    ExporterHDF5< RegionMesh <LinearTetra> > exporterICI;
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSVI;

    splitting -> setupExporter ( exporterSplitting, "Splitting" );
    ici -> setupExporter ( exporterICI, "ICI" );
    svi -> setupExporter ( exporterSVI, "SVI" );

    splitting -> exportSolution ( exporterSplitting, 0);
    ici -> exportSolution ( exporterICI, 0);
    svi -> exportSolution ( exporterSVI, 0);


    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }
    splitting   -> solveSplitting ( exporterSplitting );
    exporterSplitting.closeFile();

    ici         -> solveICI ( exporterICI );
    exporterICI.closeFile();

    svi         -> solveSVI ( exporterSVI );
    exporterSVI.closeFile();


    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    svi -> exportFiberDirection();


    if ( Comm->MyPID() == 0 )
    {
        cout << "\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
