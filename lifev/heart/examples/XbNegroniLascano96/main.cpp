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



#include <string>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/heart/solver/HeartMonodomainSolver.hpp>
#include <lifev/heart/solver/HeartIonicSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/heart/solver/XbModels/XbNegroniLascano96.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using std::cout;
using std::endl;
using namespace LifeV;


Real Calcium(const Real& t, const Real& x, const Real& y, const Real&/*z*/, const ID& /*i*/)
	{
		return ( std::exp(- 0.01 * ( t - 30.0 ) * ( t - 30.0 ) ) * (x * y) );
	}


Int main( Int argc, char** argv )
{
   //! Initializing Epetra communicator
	MPI_Init(&argc, &argv);
	boost::shared_ptr<Epetra_Comm>  Comm( new Epetra_MpiComm(MPI_COMM_WORLD) );

	if ( Comm->MyPID() == 0 ) cout << "% using MPI" << endl;

	//********************************************//
	// Type definitions                           //
	//********************************************//
	typedef RegionMesh<LinearTetra> mesh_Type;
	typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
	typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    //********************************************//
	// Starts the chronometer.                    //
	//********************************************//
	LifeChrono chronoinitialsettings;
	chronoinitialsettings.start();


	//********************************************//
	// Import mesh.				                  //
	//********************************************//
	// Using the datafile and GetPot

	//first create the GetPot
	GetPot command_line(argc, argv);
	const string data_file_name = command_line.follow("data", 2, "-f", "--file");
	GetPot dataFile(data_file_name);


	if (!Comm->MyPID()) std::cout << "My PID = " << Comm->MyPID() << std::endl;


	//Create the mesh data and read the mesh
	MeshData meshData;
	meshData.setup( dataFile,"discretization/space" );
	boost::shared_ptr<mesh_Type > fullMeshPtr(new mesh_Type( Comm ) );
	readMesh(*fullMeshPtr, meshData);
	bool verbose = 1;

	//! Construction of the partitioned mesh
	MeshPartitioner< mesh_Type >   meshPart(fullMeshPtr, Comm);


	//********************************************//
	// Building map Epetra  to define distributed //
	// vectors. We need to first construct the FE //
	// space, which needs the quadrature rules.   //
	//********************************************//

	//! Define the quadrature rules
	const ReferenceFE*    refFE_u = &feTetraP1;
	const QuadratureRule* qR_u = &quadRuleTetra15pt;
	const QuadratureRule* bdQr_u = &quadRuleTria3pt;

    //! Construction of the FE spaces
    if (verbose) std::cout << "Building the FE space ... " << std::flush;
    feSpacePtr_Type uFESpacePtr( new feSpace_Type( meshPart, *refFE_u, *qR_u, *bdQr_u, 1, Comm) );
    std::cout << " Done!" << endl;

    //Create the Map
	MapEpetra localMap(uFESpacePtr->map());




	//********************************************//
	// Import parameters from an xml list. Use    //
	// Teuchos to create a list from a given file //
	// in the execution directory.                //
	//********************************************//

	std::cout << "Importing parameters list...";
    Teuchos::ParameterList NLParameterList = *( Teuchos::getParametersFromXmlFile( "NegroniLascano96Parameters.xml" ) );
	std::cout << " Done!" << endl;


	//********************************************//
	// Creates a new model object representing the//
	// model from Negroni and Lascano 1996. The   //
	// model input are the parameters. Pass  the  //
	// parameter list in the constructor          //
	//********************************************//
	std::cout << "Building Constructor for NegrpniLascano96 Model with parameters ... ";
	XbNegroniLascano96  xb( NLParameterList );
	std::cout << " Done!" << endl;


	//********************************************//
	// Show the parameters of the model as well as//
	// other informations  about the object.      //
	//********************************************//
	xb.showMe();



	//********************************************//
	// Initialize the solution to 0. The model    //
	// consist of three state variables. Xe.Size()//
	// returns the number of state variables of   //
	// the model. rStates is the reference to the //
	// the vector states                          //
	//********************************************//

	std::cout << "Initializing solution vector...";
	//std::vector<vectorPtr_Type> states(xb.Size(), *( new vectorPtr_Type( new VectorEpetra(localMap) ) ) );
	std::vector<vectorPtr_Type> states;
	for(int k = 0; k < xb.Size(); ++k ){
		states.push_back( *(new vectorPtr_Type( new VectorEpetra(localMap) ) ) );
	}
	std::cout << " Done!" << endl;



	//********************************************//
	// Initialize the rhs to 0. The rhs is the    //
	// vector containing the numerical values of  //
	// the time derivatives of the state          //
	// variables, that is, the right hand side of //
	// the differential equation.                 //
	//********************************************//
	std::cout << "Initializing rhs..." ;
	std::vector<vectorPtr_Type> rhs;
	for(int k = 0; k < xb.Size(); ++k ){
		rhs.push_back( *(new vectorPtr_Type( new VectorEpetra(localMap) ) ) );
	}
	std::cout << " Done! "  << endl;


	//********************************************//
	// The model needs as external informations   //
	// the contraction velocity and the Calcium   //
	// concentration.                             //
	//********************************************//
	boost::shared_ptr<VectorEpetra>  vel( new VectorEpetra(localMap) );
	*vel = 0;
	boost::shared_ptr<VectorEpetra>  Ca( new VectorEpetra(localMap) );
	*Ca = 0.0;



	//********************************************//
	// Simulation starts on t=0 and ends on t=TF. //
	// The timestep is given by dt                //
	//********************************************//
	Real TF = dataFile("discretization/time/endtime",100);
	Real dt = dataFile("discretization/time/timestep",0.1);

	//********************************************//
	// Creating the exporter to save the solution //
	//********************************************//
    //! Setting generic Exporter postprocessing
    ExporterHDF5< RegionMesh <LinearTetra> > exporter( dataFile, meshPart.meshPartition(), "solution", Comm->MyPID() );

    boost::shared_ptr<VectorEpetra> state1( new VectorEpetra( ( *( states.at(0).get() ) ), Repeated ) );
    boost::shared_ptr<VectorEpetra> state2( new VectorEpetra( ( *( states.at(1).get() ) ), Repeated ) );
    boost::shared_ptr<VectorEpetra> state3( new VectorEpetra( ( *( states.at(2).get() ) ), Repeated ) );
    boost::shared_ptr<VectorEpetra> stateCa( new VectorEpetra( *Ca, Repeated ) );

    exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "TCa", uFESpacePtr,
                           state1, UInt(0) );
    exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "TCas", uFESpacePtr,
                           state2, UInt(0) );
    exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "Ts", uFESpacePtr,
                           state3, UInt(0) );
    exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "Ca", uFESpacePtr,
                           stateCa, UInt(0) );


    exporter.postProcess( 0 );

    MPI_Barrier(MPI_COMM_WORLD);
    chronoinitialsettings.stop();


    //********************************************//
	// Time loop starts.                          //
	//********************************************//
    short int iter=0;
	std::cout << "Time loop starts...\n";
	for( Real t = 0; t < TF; ){
		iter++;
		//********************************************//
		// Compute Calcium concentration. Here it is  //
		// given as a function of time and space.     //
		//********************************************//
		uFESpacePtr->interpolate( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type >( Calcium ), *Ca , t);
		std::cout << "\r " << t << " ms.       "<< std::flush;


		MPI_Barrier(MPI_COMM_WORLD);
		//********************************************//
		// Compute the rhs using the model equations  //
		//********************************************//
		xb.computeRhs( states, *Ca, *vel, rhs);

		//********************************************//
		// Use forward Euler method to advance the    //
		// solution in time.                          //
		//********************************************//

		*( states.at(0) ) = *( states.at(0) ) + dt * ( *( rhs.at(0) ) );
		*( states.at(1) ) = *( states.at(1) ) + dt * ( *( rhs.at(1) ) );
		*( states.at(2) ) = *( states.at(2) ) + dt * ( *( rhs.at(2) ) );

		//********************************************//
		// Update the time.                           //
		//********************************************//
		t = t + dt;


		//********************************************//
		// Writes solution on file.                   //
		//********************************************//
			if( iter%100 == 0 ){
				*state1 = *( states.at(0) );
				*state2 = *( states.at(1) );
				*state3 = *( states.at(2) );
				*stateCa = *Ca;

				exporter.postProcess( t );
			}

			MPI_Barrier(MPI_COMM_WORLD);

		}

	//********************************************//
	// Close exported file.                       //
	//********************************************//
    exporter.closeFile();
    std::cout << std::endl;

    //! Finalizing Epetra communicator
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return( EXIT_SUCCESS );
}
