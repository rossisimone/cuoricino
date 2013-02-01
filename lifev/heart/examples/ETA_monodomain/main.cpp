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
#include <lifev/heart/solver/HeartMonodomainSolver.hpp>
#include <lifev/heart/solver/HeartIonicSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/heart/solver/IonicModels/IonicAlievPanfilov.hpp>
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


Real Stimulus(const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
	{
		return ( 0.5 + 0.5 * ( std::tanh( - 300 * ( ( x - 0.4 ) * ( x - 0.6 ) + ( y - 0.4 ) * ( y - 0.6 ) ) ) ) );
	}
Real fZero( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/ )
{
    return 0.0;
}


Int main( Int argc, char** argv )
{

   //! Initializing Epetra communicator
	MPI_Init(&argc, &argv);
	LifeChrono timer;
	timer.start();
	boost::shared_ptr<Epetra_Comm>  Comm( new Epetra_MpiComm(MPI_COMM_WORLD) );

	if ( Comm->MyPID() == 0 ) cout << "% using MPI" << endl;

	//********************************************//
	// Type definitions                           //
	//********************************************//
	typedef RegionMesh<LinearTetra> mesh_Type;
	typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
	typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    //********************************************//
	// Starts the chronometer.                    //
	//********************************************//
	LifeChrono chronoinitialsettings;
	chronoinitialsettings.start();

	//********************************************//
	// Import parameters from an xml list. Use    //
	// Teuchos to create a list from a given file //
	// in the execution directory.                //
	//********************************************//

	std::cout << "Importing parameters list...";
    Teuchos::ParameterList APParameterList = *( Teuchos::getParametersFromXmlFile( "AlievPanfilovParameters.xml" ) );
	std::cout << " Done!" << endl;

	//********************************************//
	// Import mesh.				                  //
	//********************************************//
	if (!Comm->MyPID()) std::cout << "My PID = " << Comm->MyPID() << std::endl;
	//Create the mesh data and read and partitioned the mesh
	boost::shared_ptr< RegionMesh <LinearTetra> > meshPtr( new mesh_Type( Comm ) );
	std::string meshName = APParameterList.get("mesh_name","lid16.mesh");
	std::string meshPath = APParameterList.get("mesh_path","./");

	MeshUtility::fillWithFullMesh( meshPtr, meshName, meshPath );


	std::cout << " -- Building ETFESpaces ... " << std::flush;

	boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace
		( new ETFESpace< mesh_Type, MapEpetra, 3, 1 >(meshPtr,&feTetraP1, Comm));
 	boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace( new FESpace< mesh_Type, MapEpetra >(meshPtr, "P1", 1, Comm) );


	std::cout << " done ! " << std::endl;
	std::cout << " ---> Dofs: " << uSpace->dof().numTotalDof() << std::endl;

    std::cout << " -- Defining the matrices ... " << std::flush;

    boost::shared_ptr<matrix_Type> stiffnessMatrix (new matrix_Type( uSpace->map() ));
    boost::shared_ptr<matrix_Type> massMatrix (new matrix_Type( uSpace->map() ));

    *stiffnessMatrix *=0.0;
    *massMatrix *= 0.0;

    std::cout << " done! " << std::endl;

    std::cout << " -- Assembling the Laplace matrix ... " << std::flush;

    {
	   using namespace ExpressionAssembly;


	   integrate(  elements(uSpace->mesh()),
				   quadRuleTetra4pt,
				   uSpace,
				   uSpace,
				   dot( grad(phi_i) , grad(phi_j) )
		   )
		   >> stiffnessMatrix;

	   integrate(  elements(uSpace->mesh()),
				   quadRuleTetra4pt,
				   uSpace,
				   uSpace,
				   phi_i*phi_j
		   )
		   >> massMatrix;
    }
    std::cout << " done! " << std::endl;

    std::cout << " -- Closing the matrix ... " << std::flush;

	//********************************************//
	// Simulation starts on t=0 and ends on t=TF. //
	// The timestep is given by dt                //
	//********************************************//
	Real TF = APParameterList.get("endTime",1.0);
	Real dt = APParameterList.get("dt",0.1);
	Real D  = APParameterList.get("diffusion",0.001);
	int n  = APParameterList.get("nf" , 1);
    stiffnessMatrix->globalAssemble();
    massMatrix->globalAssemble();
     std::cout << " done ! " << std::endl;
     *stiffnessMatrix *= D;
     *massMatrix *= 1.0 / dt;
     *stiffnessMatrix += (*massMatrix);
     VectorEpetra Grhs(uSpace->map(), Repeated);
     Grhs = 0.0;
     Grhs = (*massMatrix) * Grhs;


     // Boundary conditions
     //BCHandler bcH;
     //BCFunctionBase uZero(fZero);
     //bcH.addBC( "perimeter",   	20,		Essential,	Full,	uZero,  1 );
     //bcH.addBC( "bottom",   2, 	Essential, 	Full,   uZero, 	1 );
     //bcH.addBC( "top",    	1,  	Essential, 	Full,   uZero, 	1 );



     //bcH.bcUpdate( *meshPtr, uFESpace->feBd(), uFESpace->dof() );

     //VectorEpetra rhsbc(rhs, Unique);
     boost::shared_ptr<VectorEpetra> rhsptr( new VectorEpetra( Grhs, Unique ) );

     //bcManage(*stiffnessMatrix,  *rhsbcptr, *uFESpace->mesh(), uFESpace->dof(), bcH, uFESpace->feBd(), 1.0, 0.0);


     boost::shared_ptr<VectorEpetra> Sol( new VectorEpetra( uSpace->map() ) );


  	GetPot command_line(argc, argv);
  	const string data_file_name = command_line.follow("data", 2, "-f", "--file");
  	GetPot dataFile(data_file_name);


    std::cout << "Setting up SolverAztecOO... " << std::flush;
    SolverAztecOO linearSolver1;
    linearSolver1.setCommunicator( Comm );
    linearSolver1.setDataFromGetPot( dataFile, "solver" );
    linearSolver1.setTolerance( 1e-10 );
    linearSolver1.setupPreconditioner(dataFile, "prec");
    std::cout << "done" << std::endl;
    linearSolver1.setMatrix( *stiffnessMatrix );





     ExporterHDF5< RegionMesh <LinearTetra> > exporter;
	 exporter.setMeshProcId( meshPtr, Comm->MyPID() );
	 exporter.setPrefix("solution");

	 exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "u", uFESpace, Sol, UInt(0) );



	uFESpace->interpolate( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type >( Stimulus ), *Sol , 0);


		exporter.postProcess( 0 );


//	//********************************************//
//	// Building map Epetra  to define distributed //
//	// vectors. We need to first construct the FE //
//	// space, which needs the quadrature rules.   //
//	//********************************************//

	//********************************************//
	// Creates a new model object representing the//
	// model from Negroni and Lascano 1996. The   //
	// model input are the parameters. Pass  the  //
	// parameter list in the constructor          //
	//********************************************//
	std::cout << "Building Constructor for Aliev Panfilov Model with parameters ... ";
	IonicAlievPanfilov  model( APParameterList );
	std::cout << " Done!" << endl;
//
//
//	//********************************************//
//	// Show the parameters of the model as well as//
//	// other informations  about the object.      //
//	//********************************************//
//	model.showMe();
//
//
//
	//********************************************//
	// Initialize the solution to 0. The model    //
	// consist of three state variables. Xe.Size()//
	// returns the number of state variables of   //
	// the model. rStates is the reference to the //
	// the vector states                          //
	//********************************************//

	std::cout << "Initializing solution vector...";
	std::vector<vectorPtr_Type> unknowns;
	unknowns.push_back( Sol );
	for(int k = 1; k < model.Size(); ++k ){
		unknowns.push_back( *(new vectorPtr_Type( new VectorEpetra(uSpace->map()) ) ) );
	}
	std::cout << " Done!" << endl;
//
//
//
	//********************************************//
	// Initialize the rhs to 0. The rhs is the    //
	// vector containing the numerical values of  //
	// the time derivatives of the state          //
	// variables, that is, the right hand side of //
	// the differential equation.                 //
	//********************************************//
	std::cout << "Initializing rhs..." ;
	std::vector<vectorPtr_Type> rhs;
	rhs.push_back( rhsptr );
	for(int k = 1; k < model.Size(); ++k ){
		rhs.push_back( *(new vectorPtr_Type( new VectorEpetra(uSpace->map()) ) ) );
	}
	std::cout << " Done! "  << endl;
//
//
	//********************************************//
	// The model needs as external informations   //
	// the contraction velocity and the Calcium   //
	// concentration.                             //
	//********************************************//
	boost::shared_ptr<VectorEpetra>  Iapp( new VectorEpetra(uSpace->map()) );
	*Iapp = 0;
//
//

//
//	//********************************************//
//	// Creating the exporter to save the solution //
//	//********************************************//
//    //! Setting generic Exporter postprocessing
//    ExporterHDF5< RegionMesh <LinearTetra> > exporter;
//    exporter.setMeshProcId( meshPtr, Comm->MyPID() );
//    exporter.setPrefix("solution");
//
//    boost::shared_ptr<VectorEpetra> V( new VectorEpetra( ( *( unknowns.at(0).get() ) ), Repeated ) );
//    boost::shared_ptr<VectorEpetra> r( new VectorEpetra( ( *( unknowns.at(1).get() ) ), Repeated ) );
//
//
//    exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "V", uFESpacePtr,
//                           V, UInt(0) );
//    exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "r", uFESpacePtr,
//                           r, UInt(0) );
//
//
//
//    exporter.postProcess( 0 );
//

//    chronoinitialsettings.stop();
//
//
//    cout << "nbdof: " << uFESpacePtr.get()->fe().nbFEDof() << ", nbquad: " << uFESpacePtr.get()->fe().nbQuadPt() << endl ;
//    //********************************************//
//	// Time loop starts.                          //
//	//********************************************//
//    short int iter=0;
//	std::cout << "\nTime loop starts...\n";
//	for( Real t = 0; t < TF; ){
//		iter++;
//		//********************************************//
//		// Compute Calcium concentration. Here it is  //
//		// given as a function of time and space.     //
//		//********************************************//
//		uFESpacePtr->interpolate( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type >( Stimulus ), *Iapp , t);
//		std::cout << "\r " << t << " ms.       "<< std::flush;
//
//
//		MPI_Barrier(MPI_COMM_WORLD);
//		//********************************************//
//		// Compute the rhs using the model equations  //
//		//********************************************//
//		model.computeRhs( unknowns, *Iapp, rhs);
//
//		//********************************************//
//		// Use forward Euler method to advance the    //
//		// solution in time.                          //
//		//********************************************//
//
//		*( unknowns.at(0) ) = *( unknowns.at(0) ) + dt * ( *( rhs.at(0) ) );
//		*( unknowns.at(1) ) = *( unknowns.at(1) ) + dt * ( *( rhs.at(1) ) );
//
//		//********************************************//
//		// Update the time.                           //
//		//********************************************//
//		t = t + dt;
//
//
//		//********************************************//
//		// Writes solution on file.                   //
//		//********************************************//
//			if( iter%100 == 0 ){
//				*V = *( unknowns.at(0) );
//				*r = *( unknowns.at(1) );
//
//				exporter.postProcess( t );
//			}
//
//			MPI_Barrier(MPI_COMM_WORLD);
//
//		}


	for( Real t = 0; t <  TF; ){

			//MPI_Barrier(MPI_COMM_WORLD);

				model.computeRhs( unknowns, *Iapp, rhs);

				*( unknowns.at(0) ) = *( unknowns.at(0) ) + dt * ( *( rhs.at(0) ) );
				*( unknowns.at(1) ) = *( unknowns.at(1) ) + dt * ( *( rhs.at(1) ) );

			//MPI_Barrier(MPI_COMM_WORLD);

			*rhsptr = (*massMatrix) * (*Sol);
			linearSolver1.solveSystem( *rhsptr, *Sol, stiffnessMatrix );

			t = t + dt;
			//MPI_Barrier(MPI_COMM_WORLD);
			exporter.postProcess( t );

	}


	//********************************************//
	// Close exported file.                       //
	//********************************************//
    //exporter.closeFile();
	 exporter.closeFile();
    std::cout << std::endl;
    timer.stop();
    std::cout << "Time: "<< timer.diff() << ", on rank: " << Comm->MyPID() <<std::endl;
    //! Finalizing Epetra communicator
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return( EXIT_SUCCESS );
}
