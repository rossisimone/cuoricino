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
  @brief Class for solving the Monodomain equations in electrophysiology.

  @date 11-2007
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>

  @contributors J.Castelneau (INRIA), Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @mantainer Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @last update 11-2010
 */

#ifndef _HEARTETAMONODOMAINSOLVER_H_
#define _HEARTETAMONODOMAINSOLVER_H_

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/heart/solver/IonicModels/HeartIonicModel.hpp>

#include <lifev/core/util/LifeChrono.hpp>
#include <boost/shared_ptr.hpp>
#include <lifev/core/fem/FESpace.hpp>


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>

namespace LifeV
{

//! monodomainSolver - Class featuring the usual solver for monodomain equations

template< typename Mesh, typename IonicModel >
class HeartETAMonodomainSolver
{

public:

    //! @name Type definitions
    //@{

	 //stypedef typename Mesh				meth_Type;
	 typedef Mesh							mesh_Type;
	 typedef boost::shared_ptr< mesh_Type > 	 	meshPtr_Type;

	 typedef VectorEpetra											vector_Type;
	 typedef boost::shared_ptr<VectorEpetra> 				    vectorPtr_Type;
	 typedef std::vector<vectorPtr_Type>							vectorOfPtr_Type;

	 typedef MatrixEpetra<Real> 									matrix_Type;
	 typedef boost::shared_ptr<matrix_Type>						matrixPtr_Type;
	 typedef boost::shared_ptr<Epetra_Comm>						commPtr_Type;

	 typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 > 						 	ETFESpace_Type;
	 typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > 	ETFESpacePtr_Type;
	 typedef FESpace< mesh_Type, MapEpetra >					 		  	feSpace_Type;
	 typedef boost::shared_ptr<feSpace_Type> 						 			feSpacePtr_Type;

	 typedef boost::shared_ptr<LinearSolver>						linearSolverPtr_Type;

	 typedef ExporterHDF5< RegionMesh <LinearTetra> >			exporter_Type;
	 typedef boost::shared_ptr<exporter_Type>						exporterPtr_Type;

	typedef LifeV::Preconditioner             basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack       prec_Type;
	typedef boost::shared_ptr<prec_Type>      precPtr_Type;


	typedef HeartIonicModel 										superIonicModel;
	typedef boost::shared_ptr<IonicModel>		ionicModelPtr_Type;

	typedef Teuchos::ParameterList					list_Type;
	 //    typedef HeartMonodomainData data_type;
//    typedef Real ( *Function ) ( const Real& t,
//                                 const Real& x,
//                                 const Real& y,
//                                 const Real& z,
//                                 const ID& id );
//    typedef boost::function<Real ( Real const&, Real const&, Real const&,
//                                   Real const&, ID const& )> source_Type;
//
//    typedef Mesh mesh_Type;
//    typedef BCHandler                               bcHandlerRaw_Type;
//    typedef boost::shared_ptr<bcHandlerRaw_Type>    bcHandler_Type;
//
//    typedef typename SolverType::matrix_type        matrix_Type;
//    typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;
//    typedef typename SolverType::vector_type        vector_Type;
//
//    typedef typename SolverType::prec_raw_type      preconditionerRaw_Type;
//    typedef typename SolverType::prec_type          preconditioner_Type;


     //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param dataType
     * @param potential FE space
     * @param bcHandler boundary conditions for potential
     * @param Epetra communicator
     */
	HeartETAMonodomainSolver(){}

	HeartETAMonodomainSolver( commPtr_Type 		comm	);

	HeartETAMonodomainSolver( ionicModelPtr_Type model	);

	HeartETAMonodomainSolver( ionicModelPtr_Type model,
								 commPtr_Type 		comm	);

	HeartETAMonodomainSolver( list_Type 			list,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model	);

	HeartETAMonodomainSolver( list_Type 			list,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model,
								 commPtr_Type 		comm	);

	HeartETAMonodomainSolver( std::string		meshName,
								 std::string 		meshPath,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model	);

	HeartETAMonodomainSolver( std::string 		meshName,
								 std::string 		meshPath,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model,
								 commPtr_Type 		comm	);


//    HeartETAMonodomainSolver(  const data_type& dataType,
//                      	  	  	  FESpace<Mesh, MapEpetra>& uFESpace,
//                      	  	  	  BCHandler& bcHandler,
//                      	  	  	  boost::shared_ptr<Epetra_Comm>& comm);

    //! Destructor
    virtual ~HeartETAMonodomainSolver() {}



    //@}

    //! @name Methods
    //@{

    void setup(std::string meshName, std::string meshPath, GetPot& dataFile, short int ionicSize);

    inline const Real surfaceVolumeRatio()	const { return M_surfaceVolumeRatio; 	}
    inline const Real membraneCapacitance() 	const { return M_membraneCapacitance; 	}
    inline const Real initialTime() 			const { return M_initialTime; 			}
    inline const Real timeStep() 				const { return M_timeStep; 			}
    inline const Real endTime() 				const { return M_endTime; 				}
    inline const Real diffusionCoeff()		const { return M_diffusionCoeff; 		}

    inline const ionicModelPtr_Type 	ionicModelPtr		()	const { return M_ionicModelPtr; 			}
    inline const commPtr_Type 			commPtr				()	const { return M_commPtr; 				}
    inline const meshPtr_Type			meshPtr				()	const { return M_meshPtr; 				}
    inline const ETFESpacePtr_Type 	ETFESpacePtr		()	const { return M_ETFESpacePtr; 		}
    inline const feSpacePtr_Type		feSpacePtr			()	const { return M_feSpacePtr; 			}
    inline const matrixPtr_Type		massMatrixPtr		()	const { return M_massMatrixPtr; 		}
    inline const matrixPtr_Type		stiffnessMatrixPtr()	const { return M_stiffnessMatrixPtr;	}
    inline const matrixPtr_Type		globalMatrixPtr	()	const { return M_globalMatrixPtr; 		}
    inline const vectorPtr_Type		rhsPtr				()	const { return M_rhsPtr; 				}
    inline const vectorPtr_Type		rhsPtrUnique		()	const { return M_rhsPtrUnique; 		}
    inline const vectorPtr_Type		potentialPtr		()	const { return M_potentialPtr; 		}
    inline const linearSolverPtr_Type	linearSolverPtr	()	const { return M_linearSolverPtr; 		}
    inline const vectorOfPtr_Type&		globalSolution		()	const { return M_globalSolution; 		}
    inline const vectorOfPtr_Type&		globalRhs			()	const { return M_globalRhs; 			}


    inline  void setSurfaceVolumeRatio	( const Real & p ) { this -> M_surfaceVolumeRatio = p; }
    inline  void setMembraneCapacitance	( const Real & p ) { this -> M_membraneCapacitance = p;}
    inline  void setInitialTime 			( const Real & p ) { this -> M_initialTime = p;  		}
    inline  void setTimeStep 				( const Real & p ) { this -> M_timeStep = p;  			}
    inline  void setEndTime 				( const Real & p ) { this -> M_endTime = p;  			}
    inline  void setDiffusionCoeff		( const Real & p ) { this -> M_diffusionCoeff = p;  	}


    //@}

private:
    Real 				M_surfaceVolumeRatio;
    Real				M_membraneCapacitance;

    ionicModelPtr_Type	M_ionicModelPtr;

    commPtr_Type 		M_commPtr;
	meshPtr_Type		M_meshPtr;
	ETFESpacePtr_Type   M_ETFESpacePtr;
	feSpacePtr_Type		M_feSpacePtr;
	matrixPtr_Type		M_massMatrixPtr;
	matrixPtr_Type		M_stiffnessMatrixPtr;
	matrixPtr_Type		M_globalMatrixPtr;

	Real				M_initialTime;
	Real				M_endTime;
	Real				M_timeStep;
	Real				M_diffusionCoeff;

	vectorPtr_Type			M_rhsPtr;
	vectorPtr_Type			M_rhsPtrUnique;
	vectorPtr_Type			M_potentialPtr;

	linearSolverPtr_Type	M_linearSolverPtr;

	//exporterPtr_Type		M_exporter;

	vectorOfPtr_Type		M_globalSolution;
	vectorOfPtr_Type		M_globalRhs;








}; // class MonodomainSolver



//
// IMPLEMENTATION
//
// ===================================================
//! Constructors
// ===================================================
template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::HeartETAMonodomainSolver( commPtr_Type 		comm	):
M_surfaceVolumeRatio(2400.0),//cm^-1
M_membraneCapacitance(1.0),//uF
M_diffusionCoeff(1e-4),
M_initialTime(0.0),
M_timeStep(0.01),
M_endTime(100.0),
M_commPtr( comm ),
M_meshPtr( new mesh_Type( M_commPtr ) ),
M_linearSolverPtr ( new LinearSolver() ),
M_globalSolution ( *(new vectorOfPtr_Type() ) ),
M_globalRhs ( *(new vectorOfPtr_Type() ) )
{
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::HeartETAMonodomainSolver(ionicModelPtr_Type model):
M_surfaceVolumeRatio(2400.0),//cm^-1
M_membraneCapacitance(1.0),//uF
M_diffusionCoeff(1e-4),
M_ionicModelPtr(model),
M_initialTime(0.0),
M_timeStep(0.01),
M_endTime(100.0),
M_commPtr( new Epetra_MpiComm(MPI_COMM_WORLD) ),
M_meshPtr( new mesh_Type( M_commPtr ) ),
M_linearSolverPtr ( new LinearSolver() ),
M_globalSolution ( *(new vectorOfPtr_Type() ) ),
M_globalRhs ( *(new vectorOfPtr_Type() ) )
{
}




template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::HeartETAMonodomainSolver( 	list_Type list,
																		GetPot& dataFile,
																		ionicModelPtr_Type model):
M_surfaceVolumeRatio( list.get("surfaceVolumeRatio", 2400.0 ) ),
M_membraneCapacitance( list.get("membraneCapacitance", 1.0 ) ),
M_diffusionCoeff( list.get("diffusionCoeff", 0.001 ) ),
M_initialTime( list.get("initialTime", 0.0 ) ),
M_endTime( list.get("endTime", 100.0 ) ),
M_timeStep( list.get("diffusionCoeff", 0.01 ) ),
M_ionicModelPtr(model),
M_commPtr( new Epetra_MpiComm(MPI_COMM_WORLD) ),
M_meshPtr( new mesh_Type( M_commPtr ) ),
M_linearSolverPtr ( new LinearSolver() ),
M_globalSolution ( *(new vectorOfPtr_Type() ) ),
M_globalRhs ( *(new vectorOfPtr_Type() ) )
{
		setup(list.get("meshName", "lid16.mesh"), list.get("meshPath", "./"), dataFile, M_ionicModelPtr -> Size() );
}


template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::HeartETAMonodomainSolver( std::string meshName, std::string meshPath, GetPot& dataFile, ionicModelPtr_Type model):
	M_surfaceVolumeRatio(2400.0),//cm^-1
	M_membraneCapacitance(1.0),//uF
	M_diffusionCoeff(1e-4),
	M_ionicModelPtr(model),
	M_initialTime(0.0),
	M_timeStep(0.01),
	M_endTime(100.0),
	M_commPtr( new Epetra_MpiComm(MPI_COMM_WORLD) ),
	M_meshPtr( new mesh_Type( M_commPtr ) ),
	M_linearSolverPtr ( new LinearSolver() ),
	M_globalSolution ( *(new vectorOfPtr_Type() ) ),
	M_globalRhs ( *(new vectorOfPtr_Type() ) )
{
	setup(meshName, meshPath, dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::setup( std::string meshName, std::string meshPath, GetPot& dataFile, short int ionicSize)
{
	//partitioning the mesh
	MeshUtility::fillWithFullMesh( M_meshPtr, meshName, meshPath );


	M_ETFESpacePtr.reset( new ETFESpace_Type( M_meshPtr, &feTetraP1, M_commPtr) );
	M_feSpacePtr.reset(new FESpace< mesh_Type, MapEpetra >(M_meshPtr, "P1", 1, M_commPtr) );
	M_feSpacePtr.reset ( new feSpace_Type(M_meshPtr, "P1", 1, M_commPtr)  );
	M_massMatrixPtr.reset ( new matrix_Type( M_ETFESpacePtr->map() ) );
	M_stiffnessMatrixPtr.reset ( new matrix_Type( M_ETFESpacePtr->map() ) );
	M_globalMatrixPtr.reset ( new matrix_Type( M_ETFESpacePtr->map() ) );
	M_rhsPtr.reset ( new vector_Type( M_ETFESpacePtr->map(), Repeated ) );
	M_rhsPtrUnique.reset ( new vector_Type( *(M_rhsPtr), Unique ) );
	M_potentialPtr.reset ( new vector_Type( M_ETFESpacePtr->map() ) );





	//***********************//
	//  Setup Mass Matrix    //
	//***********************//
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_ETFESpacePtr -> mesh() ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   phi_i*phi_j
			   )
			   >> M_massMatrixPtr;

	}
	M_massMatrixPtr -> globalAssemble();

	//***********************//
	//Setup Stiffness Matrix //
	//***********************//
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_ETFESpacePtr -> mesh() ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   dot( grad(phi_i) , grad(phi_j) )
		   )
		   >> M_stiffnessMatrixPtr;

	}
	M_stiffnessMatrixPtr -> globalAssemble();

	//***********************//
	//Setup Global    Matrix //
	//***********************//
	(*M_globalMatrixPtr) = (*M_stiffnessMatrixPtr);
	(*M_globalMatrixPtr) *= M_diffusionCoeff;
	(*M_globalMatrixPtr) += ( (*M_massMatrixPtr) * ( 1.0 / M_timeStep ) );



	//***********************//
	//  Setup Linear Solver  //
	//***********************//



	prec_Type* precRawPtr;
	basePrecPtr_Type precPtr;
	precRawPtr = new prec_Type;
	precRawPtr->setDataFromGetPot( dataFile, "prec" );
	precPtr.reset( precRawPtr );
	Teuchos::RCP< Teuchos::ParameterList > solverParamList = Teuchos::rcp ( new Teuchos::ParameterList );
	solverParamList = Teuchos::getParametersFromXmlFile( "MonodomainSolverParamList.xml" );

	M_linearSolverPtr -> setCommunicator( M_commPtr );
	M_linearSolverPtr -> setParameters( *solverParamList );
	M_linearSolverPtr -> setPreconditioner( precPtr );
	M_linearSolverPtr -> setOperator( M_globalMatrixPtr );



	//***********************//
	//  Setup Exporter       //
	//***********************//
//	M_exporter -> setMeshProcId( M_meshPtr, M_comm -> MyPID() );
//	M_exporter -> setPrefix( "solution" );
//
//	M_exporter -> addVariable( ExporterData<mesh_Type>::ScalarField,  "potential", M_feSpacePtr, M_potential, UInt(0) );

	//**************************//
	//  Setup Initial condition //
	//**************************//
	(*M_potentialPtr) = 0.5;

	M_globalSolution.push_back( M_potentialPtr );
	for(int k = 1; k <  ionicSize; ++k )
		M_globalSolution.push_back( *(new vectorPtr_Type( new VectorEpetra( M_ETFESpacePtr -> map() ) ) ) );

	M_globalRhs.push_back( M_rhsPtr );
	for(int k = 1; k < ionicSize; ++k )
		M_globalRhs.push_back( *(new vectorPtr_Type( new VectorEpetra( M_ETFESpacePtr -> map() ) ) ) );
}
//	template<typename Mesh, typename SolverType>
//	HeartMonodomainSolver<Mesh, SolverType>::
//	HeartMonodomainSolver( const data_type&          dataType,
//					  FESpace<Mesh, MapEpetra>& uFESpace,
//					  BCHandler&                BCh_u,
//					  boost::shared_ptr<Epetra_Comm>&              comm ):
//		M_data                   ( dataType ),
//		M_uFESpace               ( uFESpace ),
//		M_comm                   ( comm ),
//		M_me                     ( M_comm->MyPID() ),
//		M_BChandlerElectric      ( &BCh_u ),
//		M_setBC                  ( true ),
//		M_localMap               ( M_uFESpace.map() ),
//		M_localMapVector         (M_localMap+M_localMap+M_localMap),
//		M_massMatrix             ( ),
//		M_stiffnessMatrix	     ( ),
//		M_matrNoBC               ( ),
//		M_rhsNoBC                ( M_localMap ),
//		M_solutionTransmembranePotential      ( M_localMap ),
//		M_fiberVector                  ( M_localMapVector, Repeated ),
//		M_residual               ( M_localMap ),
//		M_linearSolver           ( ),
//		M_preconditioner         ( ),
//		M_verbose                ( M_me == 0),
//		M_updated                ( false ),
//		M_reusePreconditioner    ( true ),
//		M_resetPreconditioner    ( true ),
//		M_maxIteration           ( -1 ),
//		M_recomputeMatrix        ( false ),
//		M_stiffnessElementaryMatrix ( M_uFESpace.fe().nbFEDof(), 1, 1 ),
//		M_massElementaryMatrix   ( M_uFESpace.fe().nbFEDof(), 1, 1 )
//	{
//
//		if (M_data.hasFibers() )
//			{
//				std::stringstream MyPID;
//				ifstream fibers(M_data.fibersFile().c_str());
//
//				std::cout << "fiber_file: " <<  M_data.fibersFile().c_str() << std::endl;
//				UInt NumGlobalElements= M_localMapVector.map(Repeated)->NumGlobalElements();
//				std::vector<Real> fiber_global_vector(NumGlobalElements);
//
//				for( UInt i=0; i< NumGlobalElements; ++i)
//					fibers>>fiber_global_vector[i];
//				 UInt NumMyElements = M_localMapVector.map(Repeated)->NumMyElements();
//				for(UInt j=0; j< NumMyElements; ++j)
//				{
//					UInt ig= M_localMapVector.map(Repeated)->MyGlobalElements()[j];
//					M_fiberVector[ig]= fiber_global_vector[ig];
//					}
//				std::cout << std::endl;
//				fiber_global_vector.clear();
//			}
//	};
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::setup( const GetPot& dataFile )
//	{
//
//		M_diagonalize = dataFile( "electric/space_discretization/diagonalize",  1. );
//
//		M_reusePreconditioner   = dataFile( "electric/prec/reuse", true);
//
//		M_linearSolver.setCommunicator(M_comm);
//
//		M_linearSolver.setDataFromGetPot( dataFile, "electric/solver" );
//
//		M_maxIteration = dataFile( "electric/solver/max_iter", -1);
//
//		std::string precType = dataFile( "electric/prec/prectype", "Ifpack");
//
//		M_preconditioner.reset( PRECFactory::instance().createObject( precType ) );
//		ASSERT(M_preconditioner.get() != 0, "monodomainSolver : Preconditioner not set");
//
//		M_preconditioner->setDataFromGetPot( dataFile, "electric/prec" );
//	}
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::buildSystem()
//	{
//		M_massMatrix.reset  ( new matrix_Type(M_localMap) );
//		M_stiffnessMatrix.reset( new matrix_Type(M_localMap) );
//
//		if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";
//
//		LifeChrono chrono;
//
//		LifeChrono chronoDer;
//		LifeChrono chronoStiff;
//		LifeChrono chronoMass;
//
//
//		LifeChrono chronoStiffAssemble;
//		LifeChrono chronoMassAssemble;
//		LifeChrono chronoZero;
//
//		M_comm->Barrier();
//
//		chrono.start();
//
//		//! Elementary computation and matrix assembling
//		//! Loop on elements
//
//		for ( UInt iVol = 0; iVol < M_uFESpace.mesh()->numVolumes(); iVol++ )
//		{
//			chronoZero.start();
//
//			M_stiffnessElementaryMatrix.zero();
//			M_massElementaryMatrix.zero();
//
//			chronoZero.stop();
//
//			chronoStiff.start();
//			switch(M_data.heartDiffusionFactor() )
//			{
//			case 0:
//				M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
//				if (M_data.hasFibers() )
//				{
//					stiff( M_data.longitudinalConductivity(),
//						   M_data.transversalConductivity(),
//						   M_fiberVector,
//						   M_stiffnessElementaryMatrix,
//						   M_uFESpace.fe(),
//						   M_uFESpace.dof(),
//						   0,
//						   0);
//				}
//				else
//				{
//					stiff( M_data.diffusivity(), M_stiffnessElementaryMatrix,  M_uFESpace.fe(), 0, 0 );
//				}
//				break;
//
//				case 1:
//					M_uFESpace.fe().updateFirstDerivQuadPt( M_uFESpace.mesh()->volumeList( iVol ) );
//					if (M_data.hasFibers() )
//					{
//						stiff( M_data.reducedConductivitySphere(),
//							   M_data.longitudinalConductivity(),
//							   M_data.transversalConductivity(),
//							   M_fiberVector,
//							   M_stiffnessElementaryMatrix,
//							   M_uFESpace.fe(),
//							   M_uFESpace.dof(), 0, 0);
//					}
//					else
//					{
//						stiff( M_data.reducedConductivitySphere(),
//							   M_data.diffusivity(), M_stiffnessElementaryMatrix,
//							   M_uFESpace.fe(),
//							   M_uFESpace.dof(),
//							   0,
//							   0);
//					}
//					break;
//			case 2:
//				M_uFESpace.fe().updateFirstDerivQuadPt( M_uFESpace.mesh()->volumeList( iVol ) );
//				if (M_data.hasFibers() )
//				{
//					stiff( M_data.reducedConductivityCylinder(),
//						   M_data.longitudinalConductivity(),
//						   M_data.transversalConductivity(),
//						   M_fiberVector,
//						   M_stiffnessElementaryMatrix,
//						   M_uFESpace.fe(),
//						   M_uFESpace.dof(),
//						   0,
//						   0);
//				}
//				else
//				{
//					stiff( M_data.reducedConductivityCylinder(),
//						   M_data.diffusivity(),
//						   M_stiffnessElementaryMatrix,
//						   M_uFESpace.fe(),
//						   M_uFESpace.dof(),
//						   0,
//						   0);
//				}
//				break;
//			case 3:
//				M_uFESpace.fe().updateFirstDerivQuadPt( M_uFESpace.mesh()->volumeList( iVol ) );
//				if (M_data.hasFibers() )
//				{
//					stiff( M_data.reducedConductivityBox(),
//						   M_data.longitudinalConductivity(),
//						   M_data.transversalConductivity(),
//						   M_fiberVector,
//						   M_stiffnessElementaryMatrix,
//						   M_uFESpace.fe(),
//						   M_uFESpace.dof(),
//						   0,
//						   0);
//				}
//				else
//				{
//					stiff( M_data.reducedConductivityBox(),
//						   M_data.diffusivity(),
//						   M_stiffnessElementaryMatrix,
//						   M_uFESpace.fe(),
//						   M_uFESpace.dof(),
//						   0,
//						   0);
//				}
//				break;
//			}
//			chronoStiff.stop();
//
//			chronoMass.start();
//			mass( 1., M_massElementaryMatrix, M_uFESpace.fe(), 0, 0 );
//			chronoMass.stop();
//
//
//			chronoStiffAssemble.start();
//			assembleMatrix( *M_stiffnessMatrix,
//							M_stiffnessElementaryMatrix,
//							M_uFESpace.fe(),
//							M_uFESpace.fe(),
//							M_uFESpace.dof(),
//							M_uFESpace.dof(),
//							0, 0,
//							0, 0);
//			chronoStiffAssemble.stop();
//
//
//			chronoMassAssemble.start();
//			assembleMatrix( *M_massMatrix,
//							M_massElementaryMatrix,
//							M_uFESpace.fe(),
//							M_uFESpace.fe(),
//							M_uFESpace.dof(),
//							M_uFESpace.dof(),
//							0, 0,
//							0, 0);
//			chronoMassAssemble.stop();
//
//		}
//
//		massCoefficient = M_data.volumeSurfaceRatio() * M_data.membraneCapacitance() / M_data.timeStep();
//
//		M_comm->Barrier();
//
//		chrono.stop();
//		if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
//
//		if (M_verbose) std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
//		chrono.start();
//
//		M_stiffnessMatrix->globalAssemble();
//		M_massMatrix->globalAssemble();
//
//		M_comm->Barrier();
//
//		M_matrNoBC.reset(new matrix_Type(M_localMap, M_stiffnessMatrix->meanNumEntries() ));
//
//		//! Computing 1/dt * M + K
//
//		*M_matrNoBC += *M_stiffnessMatrix;
//
//		*M_matrNoBC += *M_massMatrix*massCoefficient;
//
//		M_matrNoBC->globalAssemble();
//
//		chrono.stop();
//		if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;
//
//		if (false)
//			std::cout << "partial times:  \n"
//					  << " Der            " << chronoDer.diffCumul() << " s.\n"
//					  << " Zero           " << chronoZero.diffCumul() << " s.\n"
//					  << " Stiff          " << chronoStiff.diffCumul() << " s.\n"
//					  << " Stiff Assemble " << chronoStiffAssemble.diffCumul() << " s.\n"
//					  << " Mass           " << chronoMass.diffCumul() << " s.\n"
//					  << " Mass Assemble  " << chronoMassAssemble.diffCumul() << " s.\n"
//					  << std::endl;
//
//	}
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::
//	initialize( const source_Type& u0 )
//	{
//
//		vector_Type u(M_uFESpace.map());
//
//		M_uFESpace.interpolate(u0, u, 0.);
//
//		initialize(u);
//	}
//
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::
//	initialize( const Function& u0 )
//	{
//
//		 vector_Type u(M_uFESpace.map());
//		 M_uFESpace.interpolate(u0, u, 0.);
//
//		 initialize(u);
//	}
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::
//	initialize( const vector_Type& u0 )
//	{
//
//		M_solutionTransmembranePotential = u0;
//
//	}
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::
//	updatePDESystem(Real         alpha,
//					vector_Type& sourceVec )
//	{
//
//		LifeChrono chrono;
//
//		if (M_verbose)
//				std::cout << "  f-  Updating mass term and right hand side... "
//					  << std::flush;
//
//		chrono.start();
//
//		M_rhsNoBC = sourceVec;
//		M_rhsNoBC.globalAssemble();
//
//		chrono.stop();
//
//		if (M_verbose)
//			std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;
//
//		M_updated = false;
//
//		if (M_recomputeMatrix)
//			buildSystem();
//
//		if (M_verbose)
//			  std::cout << "  f-  Copying the matrices ...                 "
//						<< std::flush;
//
//		chrono.start();
//
//		M_matrNoBC.reset(new matrix_Type(M_localMap, M_stiffnessMatrix->meanNumEntries() ));
//
//		*M_matrNoBC += *M_stiffnessMatrix;
//
//		*M_matrNoBC += *M_massMatrix*alpha;
//
//
//		chrono.stop();
//		if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
//								 << std::flush;
//
//
//
//		M_updated = true;
//		M_matrNoBC->globalAssemble();
//
//	}
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::
//	updatePDESystem(vector_Type& sourceVec )
//	{
//
//		LifeChrono chrono;
//
//		if (M_verbose)
//				std::cout << "  f-  Updating right hand side... "
//					  << std::flush;
//
//		chrono.start();
//
//		M_rhsNoBC = sourceVec;
//		M_rhsNoBC.globalAssemble();
//
//		chrono.stop();
//
//		if (M_verbose)
//			std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;
//
//	}
//
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::PDEiterate( bcHandlerRaw_Type& bch )
//	{
//
//		LifeChrono chrono;
//		chrono.start();
//
//		matrixPtr_Type matrFull( new matrix_Type(*M_matrNoBC) );
//		vector_Type    rhsFull = M_rhsNoBC;
//
//		chrono.stop();
//
//		if (M_verbose)
//			std::cout << "done in " << chrono.diff() << " s.\n"
//					  << std::flush;
//
//		// boundary conditions update
//		M_comm->Barrier();
//		if (M_verbose) std::cout << "  f-  Applying boundary conditions ...         "
//				  << std::flush;
//
//		chrono.start();
//
//		applyBoundaryConditions( *matrFull, rhsFull, bch);
//
//		chrono.stop();
//
//		M_comm->Barrier();
//
//		if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
//
//		//! Solving the system
//		solveSystem( matrFull, rhsFull );
//
//		//  M_residual  = M_rhsNoBC;
//		//  M_residual -= *M_matrNoBC*M_solutionTransmembranePotential;
//
//	} // iterate()
//
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::solveSystem( matrixPtr_Type  matrFull,
//														  vector_Type&    rhsFull )
//	{
//		LifeChrono chrono;
//
//		if (M_verbose)
//			std::cout << "  f-  Setting up the solver ...                ";
//
//		chrono.start();
//		M_linearSolver.setMatrix(*matrFull);
//		chrono.stop();
//
//		if (M_verbose)
//			std::cout << "done in " << chrono.diff() << " s.\n"
//					  << std::flush;
//
//		if ( !M_reusePreconditioner || M_resetPreconditioner )
//		{
//			chrono.start();
//
//			if (M_verbose)
//				std::cout << "  f-  Computing the precond ...                ";
//
//			M_preconditioner->buildPreconditioner(matrFull);
//
//			Real condest = M_preconditioner->condest();
//
//			M_linearSolver.setPreconditioner(M_preconditioner);
//
//			chrono.stop();
//			if (M_verbose)
//			{
//				std::cout << "done in " << chrono.diff() << " s.\n";
//				std::cout << "Estimated condition number = " << condest << "\n" <<  std::flush;
//			}
//
//			M_resetPreconditioner = false;
//		}
//		else
//		{
//			if (M_verbose)
//				std::cout << "f-  Reusing  precond ...                \n" <<  std::flush;
//		}
//
//
//		chrono.start();
//
//		if (M_verbose)
//			std::cout << "f-  Solving system ...                                ";
//
//		Int numIter = M_linearSolver.solve(M_solutionTransmembranePotential, rhsFull);
//
//		chrono.stop();
//
//		if (M_verbose)
//		{
//			std::cout << "\ndone in " << chrono.diff()
//					  << " s. ( " << numIter << "  iterations. ) \n"
//					  << std::flush;
//		}
//
//
//		if (numIter > M_maxIteration)
//		{
//			M_resetPreconditioner = true;
//		}
//
//		M_comm->Barrier();
//
//	}
//
//
//
//	template<typename Mesh, typename SolverType>
//	void HeartMonodomainSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_Type&        matrix,
//																	  vector_Type&        rhs,
//																	  bcHandlerRaw_Type&  BCh )
//	{
//
//		// BC manage for the PDE
//		if ( !BCh.bcUpdateDone() )
//		{
//			BCh.bcUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
//		}
//
//		vector_Type rhsFull(M_rhsNoBC,Repeated, Zero);
//
//		bcManage( matrix, rhs, *M_uFESpace.mesh(), M_uFESpace.dof(),
//				  BCh, M_uFESpace.feBd(), 1.,M_data.time() );
//
//		rhs = rhsFull;
//		if ( BCh.hasOnlyEssential() && M_diagonalize )
//		{
//			matrix.diagonalize( 1*dim_u(),
//								M_diagonalize,
//								rhs,
//								0.);
//		}
//
//	} // applyBoundaryCondition


} // namespace LifeV


#endif //_MONODOMAINSOLVER_H_
