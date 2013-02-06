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


Real Stimulus(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
	{
		return ( 0.5 + 0.5 * ( std::tanh( - 300 * ( ( x - 0.4 ) * ( x - 0.6 ) + ( y - 0.4 ) * ( y - 0.6 ) ) ) ) );
	}



Int main( Int argc, char** argv )
{

   //! Initializing Epetra communicator
	MPI_Init(&argc, &argv);
	//if ( Comm->MyPID() == 0 ) cout << "% using MPI" << endl;

    //********************************************//
	// Starts the chronometer.                    //
	//********************************************//
//	LifeChrono chronoinitialsettings;
//	chronoinitialsettings.start();

	typedef RegionMesh<LinearTetra> 						mesh_Type;
    typedef boost::function< Real(const Real& /*t*/,
    								const Real&   x,
    								const Real&   y,
    								const Real& /*z*/,
    								const ID&   /*i*/ ) > 	function_Type;

    function_Type f = &Stimulus;

	boost::shared_ptr<Epetra_Comm>  Comm( new Epetra_MpiComm(MPI_COMM_WORLD) );

	//********************************************//
	// Import parameters from an xml list. Use    //
	// Teuchos to create a list from a given file //
	// in the execution directory.                //
	//********************************************//

	std::cout << "Importing parameters list...";
    Teuchos::ParameterList APParameterList = *( Teuchos::getParametersFromXmlFile( "AlievPanfilovParameters.xml" ) );
	std::cout << " Done!" << endl;

	//********************************************//
	// Creates a new model object representing the//
	// model from Negroni and Lascano 1996. The   //
	// model input are the parameters. Pass  the  //
	// parameter list in the constructor          //
	//********************************************//
	std::cout << "Building Constructor for Aliev Panfilov Model with parameters ... ";
	boost::shared_ptr<IonicAlievPanfilov>  model( new IonicAlievPanfilov( APParameterList ) );
	std::cout << " Done!" << endl;


	std::string meshName = APParameterList.get("mesh_name","lid16.mesh");
	std::string meshPath = APParameterList.get("mesh_path","./");
	std::cout << "Building Constructor forMonodomain Solver... ";
	typedef HeartETAMonodomainSolver< mesh_Type, IonicAlievPanfilov > 		monodomainSolver_Type;
	typedef boost::shared_ptr< monodomainSolver_Type > 	monodomainSolverPtr_Type;

	GetPot command_line(argc, argv);
	const string data_file_name = command_line.follow("data", 2, "-f", "--file");
	GetPot dataFile(data_file_name);

	monodomainSolverPtr_Type splitting( new monodomainSolver_Type( meshName, meshPath, dataFile, model ) );
	monodomainSolverPtr_Type ici( new monodomainSolver_Type( dataFile, model, splitting ->  meshPtr() ) );
	monodomainSolverPtr_Type svi( new monodomainSolver_Type( dataFile, model, splitting ->  meshPtr() ) );

	//monodomainSolverPtr_Type mono1( new monodomainSolver_Type ( new IonicAlievPanfilov() ) );
	std::cout << " Done!" << endl;


	cout << "\nInitializing potential:  " ;
	splitting -> setPotentialFromFunction( f );
	ici -> copyPotential( splitting -> potentialPtr() );
	svi -> setPotentialFromFunction( f );

	cout << "Done! \n" ;

	splitting -> setTimeStep( APParameterList.get("dt",0.01) );
	splitting -> setEndTime ( APParameterList.get("endTime",60.0) );
	splitting -> setDiffusionCoeff( APParameterList.get("diffusion",0.001) );

	cout << "\nSplitting diffusion coefficient: " << splitting -> diffusionCoeff();
	ici -> setTimeStep( splitting -> timeStep() );
	ici -> setEndTime( splitting -> endTime() );
	ici -> setDiffusionCoeff( splitting -> diffusionCoeff() );
	cout << "\nICI diffusion coefficient: " << ici -> diffusionCoeff();

	svi -> setTimeStep( splitting -> timeStep() );
	svi -> setEndTime( splitting -> endTime() );
	svi -> setDiffusionCoeff( splitting -> diffusionCoeff() );
	cout << "\nSVI diffusion coefficient: " << svi -> diffusionCoeff();


	splitting -> setupGlobalMatrix();
	ici -> setupGlobalMatrix();
	svi -> setupGlobalMatrix();

	ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;
    ExporterHDF5< RegionMesh <LinearTetra> > exporterICI;
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSVI;

    splitting -> setupExporter( exporterSplitting, "Splitting" );
    ici -> setupExporter( exporterICI, "ICI" );
    svi -> setupExporter( exporterSVI, "SVI" );

    splitting -> exportSolution( exporterSplitting, 0);
    ici -> exportSolution( exporterICI, 0);
    svi -> exportSolution( exporterSVI, 0);


    cout << "\nstart solving:  " ;
    splitting 	-> solveSplitting( exporterSplitting );
    ici 		-> solveICI( exporterICI );
    svi 		-> solveSVI( exporterSVI );

    exporterSplitting.closeFile();
    exporterICI.closeFile();
    exporterSVI.closeFile();

    //********************************************//
	// Import mesh.				                  //
	//********************************************//
//	if (!Comm->MyPID()) std::cout << "My PID = " << Comm->MyPID() << std::endl;
//	//Create the mesh data and read and partitioned the mesh
//	boost::shared_ptr< RegionMesh <LinearTetra> > meshPtr( new mesh_Type( Comm ) );
//	std::string meshName = APParameterList.get("mesh_name","lid16.mesh");
//	std::string meshPath = APParameterList.get("mesh_path","./");
//
//	MeshUtility::fillWithFullMesh( meshPtr, meshName, meshPath );



//	uFESpace->interpolate( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type >( Stimulus ), *Sol , 0);



	//********************************************//
	// The model needs as external informations   //
	// the contraction velocity and the Calcium   //
	// concentration.                             //
	//********************************************//
//	boost::shared_ptr<VectorEpetra>  Iapp( new VectorEpetra(uSpace->map()) );
//	*Iapp = 0;

//	typedef HeartETAMonodomainSolver< RegionMesh<LinearTetra> > monodomainSolver_Type;
//	typedef boost::shared_ptr<HeartETAMonodomainSolver< RegionMesh<LinearTetra> > > monodomainSolverPtr_Type;
	//monodomainSolverPtr_Type mono1( new monodomainSolver_Type ( new IonicAlievPanfilov() ) );
	//monodomainSolverPtr_Type mono2( new monodomainSolver_Type ( new IonicMinimalModel() ) );
	//mono->M_ionicModel->showMe();
	//mono1->M_ionicModel->showMe();
	//mono2->M_ionicModel->showMe();

//	 exporter.closeFile();
////    std::cout << std::endl;
//    timer.stop();
//    std::cout << "Time: "<< timer.diff() << ", on rank: " << Comm->MyPID() <<std::endl;
    //! Finalizing Epetra communicator
//    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return( EXIT_SUCCESS );
}
