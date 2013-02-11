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

  @date 02-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @last update 02-2013

  This class provides interfaces to solve the monodomain equation
  ( reaction diffusion equation ) using the ETA framework.
  The solution can be performed using three different methods:
  -operator splitting method (at this point available only with forward Euler
  	  for the reaction step and backward Euler for the diffusion step. );
  -Ionic Currents Interpolation (at this point only forward Euler);
  -State Variable interpolation (at this point only forward Euler).
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


	//!Monodomain Solver
	/*!
	   The monodomain equation reads
	   \f \Chi

	*/


public:


    //! @name Type definitions
    //@{

	 typedef Mesh									mesh_Type;
	 typedef boost::shared_ptr< mesh_Type > 	 	meshPtr_Type;

	 typedef VectorEpetra							vector_Type;
	 typedef boost::shared_ptr<VectorEpetra> 	vectorPtr_Type;
	 typedef std::vector<vectorPtr_Type>			vectorOfPtr_Type;

	 typedef MatrixEpetra<Real> 					matrix_Type;
	 typedef boost::shared_ptr<matrix_Type>		matrixPtr_Type;

	 typedef boost::shared_ptr<Epetra_Comm>		commPtr_Type;

	 typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 > 						ETFESpace_Type;
	 typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > 	ETFESpacePtr_Type;
	 typedef FESpace< mesh_Type, MapEpetra >					 		  		feSpace_Type;
	 typedef boost::shared_ptr<feSpace_Type> 						 			feSpacePtr_Type;

	 typedef boost::shared_ptr<LinearSolver>						linearSolverPtr_Type;

	 typedef ExporterHDF5< mesh_Type >			exporter_Type;
	 typedef boost::shared_ptr<exporter_Type>						exporterPtr_Type;

	typedef LifeV::Preconditioner            	basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type> 	basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack    	prec_Type;
	typedef boost::shared_ptr<prec_Type>     	precPtr_Type;


	typedef HeartIonicModel 					superIonicModel;
	typedef boost::shared_ptr<IonicModel>	ionicModelPtr_Type;

	typedef Teuchos::ParameterList				list_Type;

    typedef boost::function< Real(const Real& t,
    								const Real& x,
    								const Real& y,
    								const Real& z,
    								const ID&   i ) > 	function_Type;

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
	HeartETAMonodomainSolver();

//	HeartETAMonodomainSolver( commPtr_Type 		comm	);
//
//	HeartETAMonodomainSolver( ionicModelPtr_Type model	);
//
//	HeartETAMonodomainSolver( ionicModelPtr_Type model,
//								 commPtr_Type 		comm	);

	HeartETAMonodomainSolver( list_Type 			list,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model	);

	HeartETAMonodomainSolver( list_Type 			list,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model,
								 commPtr_Type 		comm	);

	HeartETAMonodomainSolver( list_Type 			list,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model,
								 meshPtr_Type 		meshPtr	);

	HeartETAMonodomainSolver( std::string		meshName,
								 std::string 		meshPath,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model	);

	HeartETAMonodomainSolver( std::string 		meshName,
								 std::string 		meshPath,
								 GetPot& 			dataFile,
								 ionicModelPtr_Type model,
								 commPtr_Type 		comm	);

	HeartETAMonodomainSolver( GetPot& 			dataFile,
								 ionicModelPtr_Type model,
								 meshPtr_Type 		meshPtr	);


//    HeartETAMonodomainSolver(  const data_type& dataType,
//                      	  	  	  FESpace<Mesh, MapEpetra>& uFESpace,
//                      	  	  	  BCHandler& bcHandler,
//                      	  	  	  boost::shared_ptr<Epetra_Comm>& comm);

    //! Destructor
    virtual ~HeartETAMonodomainSolver() {}



    //@}

    //! @name Methods
    //@{

    inline const Real& surfaceVolumeRatio()	const { return M_surfaceVolumeRatio; 	}
    inline const Real& membraneCapacitance() 	const { return M_membraneCapacitance; 	}
    inline const Real& initialTime() 			const { return M_initialTime; 			}
    inline const Real& timeStep() 				const { return M_timeStep; 			}
    inline const Real& endTime() 				const { return M_endTime; 				}
    inline const VectorSmall<3>& diffusionTensor()		const { return M_diffusionTensor; 		}

    inline const std::string elementsOrder()	const { return M_elementsOrder; 		}


    inline const ionicModelPtr_Type 	ionicModelPtr		()	const { return M_ionicModelPtr; 		}
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
    inline const vectorPtr_Type		fiberPtr		    ()	const { return M_fiberPtr; 		}
    inline const vectorPtr_Type		appliedCurrentPtr	()	const { return M_appliedCurrentPtr;	}
    inline const linearSolverPtr_Type	linearSolverPtr	()	const { return M_linearSolverPtr; 		}
    inline const vectorOfPtr_Type&		globalSolution		()	const { return M_globalSolution; 		}
    inline const vectorOfPtr_Type&		globalRhs			()	const { return M_globalRhs; 			}


    inline void setSurfaceVolumeRatio	( const Real & p ) { this -> M_surfaceVolumeRatio = p; }
    inline void setMembraneCapacitance	( const Real & p ) { this -> M_membraneCapacitance = p;}
    inline void setInitialTime 			( const Real & p ) { this -> M_initialTime = p;  		}
    inline void setTimeStep 				( const Real & p ) { this -> M_timeStep = p;  			}
    inline void setEndTime 				( const Real & p ) { this -> M_endTime = p;  			}
    inline void setDiffusionTensor		( const VectorSmall<3> & p ) { this -> M_diffusionTensor = p;  	}


    inline void setIonicModelPtr 		( const ionicModelPtr_Type 	p 		)	 { this -> M_ionicModelPtr = p ; 		}
    inline void setCommPtr				( const commPtr_Type 		p		)	 { this -> M_commPtr = p ; 				}
    inline void setMeshPtr 				( const meshPtr_Type		p		)	 { this -> M_meshPtr = p ; 				}
    inline void setETFESpacePtr 			( const ETFESpacePtr_Type 	p		)	 { this -> M_ETFESpacePtr = p ; 		}
    inline void setFeSpacePtr			( const feSpacePtr_Type		p		)	 { this -> M_feSpacePtr = p ; 			}
    inline void setMassMatrixPtr			( const matrixPtr_Type		p		)	 { this -> M_massMatrixPtr = p ; 		}
    inline void setStiffnessMatrixPtr	( const matrixPtr_Type		p		)	 { this -> M_stiffnessMatrixPtr = p ;	}
    inline void setGlobalMatrixPtr		( const matrixPtr_Type		p		)	 { this -> M_globalMatrixPtr = p ;		}
    inline void setRhsPtr					( const vectorPtr_Type		p		)	 { this -> M_rhsPtr = p ; 				}
    inline void setRhsPtrUnique			( const vectorPtr_Type		p		)	 { this -> M_rhsPtrUnique = p ;
    																				   this -> M_globalRhs .at(0) = M_rhsPtrUnique; }
    inline void setPotentialPtr			( const vectorPtr_Type		p		)	 { this -> M_potentialPtr = p ;
    																				   this -> M_globalSolution.at(0) = M_potentialPtr; }
    inline void setAppliedCurrentPtr		( const vectorPtr_Type		p		)	 { this -> M_appliedCurrentPtr = p ;	}
    inline void setLinearSolverPtr		( const linearSolverPtr_Type	p	)	 { this -> M_linearSolverPtr = p ;		}
    inline void setGlobalSolution		( const vectorOfPtr_Type&	p		)	 { this -> M_globalSolution	= p; 		}
    inline void setGlobalRhs				( const vectorOfPtr_Type&	p		)	 { this -> M_globalRhs 		= p; 		}

    inline void setFiberPtr			( const vectorPtr_Type		p		)	 { this -> M_fiberPtr = p ; }

    inline void copyIonicModel 		( const ionicModelPtr_Type 	p	)	 { (*(M_ionicModelPtr) ) = *p ; 	}
    inline void copyComm				( const commPtr_Type 		p	)	 { (*(M_commPtr) ) = *p ; 			}
    inline void copyMesh 				( const meshPtr_Type		p	)	 { (*(M_meshPtr) ) = *p ; 			}
    inline void copyETFESpace 		( const ETFESpacePtr_Type 	p	)	 { (*(M_ETFESpacePtr) ) = *p ; 		}
    inline void copyFeSpace			( const feSpacePtr_Type		p	)	 { (*(M_feSpacePtr) ) = *p ; 		}
    inline void copyMassMatrix		( const matrixPtr_Type		p	)	 { (*(M_massMatrixPtr) ) = *p ; 	}
    inline void copyStiffnessMatrix	( const matrixPtr_Type		p	)	 { (*(M_stiffnessMatrixPtr) ) = *p ;}
    inline void copyGlobalMatrix		( const matrixPtr_Type		p	)	 { (*(M_globalMatrixPtr) ) = *p ;	}
    inline void copyRhs				( const vectorPtr_Type		p	)	 { (*(M_rhsPtr) ) = *p ; 			}
    inline void copyRhsUnique		( const vectorPtr_Type		p	)	 { (*(M_rhsPtrUnique) ) = *p ; 		}
    inline void copyPotential		( const vectorPtr_Type		p	)	 { (*(M_potentialPtr) ) = *p ; 		}
    inline void copyFiber				( const vectorPtr_Type		p	)	 { (*(M_fiberPtr) ) = *p ; 		}
    inline void copyAppliedCurrent	( const vectorPtr_Type		p	)	 { (*(M_appliedCurrentPtr) ) = *p ;	}
    inline void copyLinearSolver		( const linearSolverPtr_Type p	)	 { (*(M_linearSolverPtr) ) = *p ;	}
    inline void copyGlobalSolution	( const vectorOfPtr_Type&	p	)	 { for(int j = 0; j < M_ionicModelPtr -> Size(); j++ ) 	( *( M_globalSolution.at(j) ) )	= (*(p.at(j) ) ); }
    inline void copyGlobalRhs		( const vectorOfPtr_Type&	p	)	 { for(int j = 0; j < M_ionicModelPtr -> Size(); j++ ) 	( *( M_globalRhs.at(j) ) )	= (*(p.at(j) ) ); }

    //    inline void setGlobalRhs				( const vectorOfPtr_Type&	p		)	 { this -> M_globalRhs 		= p; 			}

    void setup( GetPot& dataFile, short int ionicSize);

    void setup(std::string meshName, std::string meshPath, GetPot& dataFile, short int ionicSize);

    void setupMassMatrix();
    void setupStiffnessMatrix();
    void setupStiffnessMatrix(VectorEpetra& fiber, VectorSmall<3> diffusion);
    void setupGlobalMatrix();


    void setupLinearSolver( GetPot dataFile );
    void setupLinearSolver( GetPot dataFile, list_Type list );

    void inline initializePotential() 				{ (*M_potentialPtr) 	*= 0; }
    void inline initializePotential(Real k)			{ (*M_potentialPtr)  	 = k; }
    void inline initializeAppliedCurrent()  		{ (*M_appliedCurrentPtr)*= 0; }
    void inline initializeAppliedCurrent(Real k)	{ (*M_appliedCurrentPtr) = k; }

    void setupGlobalSolution(short int ionicSize);
    void setupGlobalRhs(short int ionicSize);

    void setParameters(list_Type 	list);

    void inline partitionMesh( std::string	meshName, std::string 	meshPath)
    							{ 	MeshUtility::fillWithFullMesh( M_meshPtr, meshName, meshPath ); }

    void inline setPotentialFromFunction( function_Type f )
    							{ M_feSpacePtr -> interpolate( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type >( f ), *M_potentialPtr , 0); }
    void inline setAppliedCurrentFromFunction( function_Type f )
    							{ 	M_feSpacePtr -> interpolate( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type >( f ), *M_appliedCurrentPtr , 0); }

    void solveOneReactionStepFE();

    void inline updateRhs() {	 (*M_rhsPtrUnique) += (*M_massMatrixPtr) * (*M_potentialPtr) * ( 1.0 / M_timeStep ); }

    void solveOneDiffusionStepBE();

    void solveOneSplittingStep();
    void solveSplitting();

    void solveOneSplittingStep( exporter_Type& exporter, Real t );
    void solveSplitting( exporter_Type& exporter );

    void setupPotentialExporter(exporter_Type& exporter);
    void setupExporter(exporter_Type& exporter);

    void setupPotentialExporter(exporter_Type& exporter, std::string fileName);
    void setupExporter(exporter_Type& exporter, std::string fileName);

   // void setupPotentialRhsExporter(exporter_Type& exporter, std::string fileName);
   // void setupRhsExporter(exporter_Type& exporter, std::string fileName);


    void solveOneStepGatingVariablesFE();

    void computeRhsSVI();

    void computeRhsICI();

    void solveOneICIStep();

    void solveOneSVIStep();

    void solveICI();

    void solveSVI();

    void solveOneICIStep(exporter_Type& exporter, Real t);

    void solveOneSVIStep(exporter_Type& exporter, Real t);

    void solveICI(exporter_Type& exporter);

    void solveSVI(exporter_Type& exporter);

    void inline exportSolution(exporter_Type& exporter, Real t) { exporter.postProcess(t); }

    void inline exportRhs(exporter_Type& exporter, Real t) { exporter.postProcess(t); }
    //@}

private:

    void setParameters();
    void init();
	void init( commPtr_Type comm );
    void init(meshPtr_Type meshPtr);
	void init( ionicModelPtr_Type model);


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

	VectorSmall<3>		M_diffusionTensor;

	vectorPtr_Type			M_rhsPtr;
	vectorPtr_Type			M_rhsPtrUnique;
	vectorPtr_Type			M_potentialPtr;
	vectorPtr_Type			M_appliedCurrentPtr;

	linearSolverPtr_Type	M_linearSolverPtr;

	//exporterPtr_Type		M_exporter;

	vectorOfPtr_Type		M_globalSolution;
	vectorOfPtr_Type		M_globalRhs;

	std::string				M_elementsOrder;


	vectorPtr_Type			M_fiberPtr;



}; // class MonodomainSolver



//
// IMPLEMENTATION
//
// ===================================================
//! Constructors
// ===================================================
template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver()
{
	setParameters();
	init();
}


template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver( 	list_Type list,
								GetPot& dataFile,
								ionicModelPtr_Type model)
{
	init(model);
	setParameters(list);
	setup(list.get("meshName", "lid16.mesh"), list.get("meshPath", "./"), dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver( 	list_Type 			list,
								GetPot& 			dataFile,
								ionicModelPtr_Type 	model,
								commPtr_Type 		comm	)
{
	setParameters(list);
	init(comm);
	setup(list.get("meshName", "lid16.mesh"), list.get("meshPath", "./"), dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver( 	list_Type 			list,
								GetPot& 			dataFile,
								ionicModelPtr_Type 	model,
								meshPtr_Type 		meshPtr	)
{
	setParameters(list);
	init(meshPtr);
	setup(dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::
HeartETAMonodomainSolver( 	std::string 		meshName,
								std::string 		meshPath,
								GetPot& 			dataFile,
								ionicModelPtr_Type 	model	)
{
	setParameters();
	init(model);
	setup(meshName, meshPath, dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::
HeartETAMonodomainSolver( 	std::string 		meshName,
								std::string 		meshPath,
								GetPot& 			dataFile,
								ionicModelPtr_Type 	model,
								commPtr_Type 		comm	):
	M_ionicModelPtr(model)
{
	setParameters();
	init(comm);
	setup(meshName, meshPath, dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::
HeartETAMonodomainSolver( 	GetPot& 			dataFile,
								ionicModelPtr_Type 	model,
								meshPtr_Type 		meshPtr	):
	M_ionicModelPtr(model)
{
	setParameters();
	init(meshPtr);
	setup(dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setup(GetPot& dataFile, short int ionicSize)
{

	M_ETFESpacePtr.reset( new ETFESpace_Type( M_meshPtr, &feTetraP1, M_commPtr) );
	M_feSpacePtr.reset(new FESpace< mesh_Type, MapEpetra >(M_meshPtr, M_elementsOrder, 1, M_commPtr) );
	M_feSpacePtr.reset ( new feSpace_Type(M_meshPtr, M_elementsOrder, 1, M_commPtr)  );
	M_massMatrixPtr.reset ( new matrix_Type( M_ETFESpacePtr->map() ) );
	M_stiffnessMatrixPtr.reset ( new matrix_Type( M_ETFESpacePtr->map() ) );
	M_globalMatrixPtr.reset ( new matrix_Type( M_ETFESpacePtr->map() ) );
	M_rhsPtr.reset ( new vector_Type( M_ETFESpacePtr->map(), Repeated ) );
	M_rhsPtrUnique.reset ( new vector_Type( *(M_rhsPtr), Unique ) );
	M_potentialPtr.reset ( new vector_Type( M_ETFESpacePtr->map() ) );
	M_appliedCurrentPtr.reset ( new vector_Type( M_ETFESpacePtr->map() ) );


	/*************************/
	// INITIALIZE FIBERS   ***/
	/************************/
	boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
        ( new FESpace< mesh_Type, MapEpetra >( M_meshPtr, "P1", 3, M_commPtr) );

    M_fiberPtr.reset( new vector_Type( Space3D -> map() ) );

    int n1 = (*M_fiberPtr).epetraVector().MyLength();
    int d1 = n1 / 3;
    (*M_fiberPtr) *= 0;
	int j(0);
	for( int k(0); k<d1; k++){
	j = (*M_fiberPtr).blockMap().GID(k+d1);
	(*M_fiberPtr)[j]=1.0;
	}

	//***********************//
	//  Setup Mass Matrix    //
	//***********************//
	//setupMassMatrix();

	//***********************//
	//Setup Stiffness Matrix //
	//***********************//
	//setupStiffnessMatrix( *M_fiberPtr, M_diffusionTensor );

	//***********************//
	//Setup Global    Matrix //
	//***********************//
	//setupGlobalMatrix();



	//***********************//
	//  Setup Linear Solver  //
	//***********************//
	setupLinearSolver( dataFile );


	//**************************//
	//  Setup Initial condition //
	//**************************//
	initializePotential(0.0);
	initializeAppliedCurrent(0.0);
	setupGlobalSolution(ionicSize);
	setupGlobalRhs(ionicSize);
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setup( std::string meshName, std::string meshPath, GetPot& dataFile, short int ionicSize)
{
	//partitioning the mesh
	partitionMesh( meshName, meshPath );
	setup(dataFile, ionicSize);
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupMassMatrix()
{
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
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupStiffnessMatrix()
{
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_ETFESpacePtr -> mesh() ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   dot(  rotate( M_ETFESpacePtr, *M_fiberPtr, M_diffusionTensor ) * grad(phi_i) , grad(phi_j) )
		   )
		   >> M_stiffnessMatrixPtr;

	}
	M_stiffnessMatrixPtr -> globalAssemble();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupStiffnessMatrix(VectorEpetra& fiber, VectorSmall<3> diffusion)
{
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_ETFESpacePtr -> mesh() ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   dot( rotate( M_ETFESpacePtr, fiber, diffusion ) * grad(phi_i) , grad(phi_j) )
		   )
		   >> M_stiffnessMatrixPtr;

	}
	M_stiffnessMatrixPtr -> globalAssemble();
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupGlobalMatrix()
{
	(*M_globalMatrixPtr) = (*M_stiffnessMatrixPtr);
	//(*M_globalMatrixPtr) *= M_diffusionTensor;
	(*M_globalMatrixPtr) += ( (*M_massMatrixPtr) * ( 1.0 / M_timeStep ) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupLinearSolver( GetPot dataFile )
{
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
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupLinearSolver( GetPot dataFile, list_Type list )
{
	prec_Type* precRawPtr;
	basePrecPtr_Type precPtr;
	precRawPtr = new prec_Type;
	precRawPtr->setDataFromGetPot( dataFile, "prec" );
	precPtr.reset( precRawPtr );

	M_linearSolverPtr -> setCommunicator( M_commPtr );
	M_linearSolverPtr -> setParameters( list );
	M_linearSolverPtr -> setPreconditioner( precPtr );
	M_linearSolverPtr -> setOperator( M_globalMatrixPtr );
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupGlobalSolution(short int ionicSize)
{
	M_globalSolution.push_back( M_potentialPtr );
	for(int k = 1; k <  ionicSize; ++k )
		M_globalSolution.push_back( *(new vectorPtr_Type( new VectorEpetra( M_ETFESpacePtr -> map() ) ) ) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupGlobalRhs(short int ionicSize)
{
	M_globalRhs.push_back( M_rhsPtrUnique );
	for(int k = 1; k < ionicSize; ++k )
		M_globalRhs.push_back( *(new vectorPtr_Type( new VectorEpetra( M_ETFESpacePtr -> map() ) ) ) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init()
{
	M_linearSolverPtr.reset ( new LinearSolver() );
	M_globalSolution = *(new vectorOfPtr_Type() ) ;
	M_globalRhs =  *(new vectorOfPtr_Type() ) ;
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init( commPtr_Type comm  )
{
	init();
	M_commPtr = comm;
	M_meshPtr.reset( new mesh_Type( M_commPtr ) );
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init(meshPtr_Type meshPtr)
{
	init();
	M_meshPtr = meshPtr;
	M_commPtr = meshPtr -> comm();
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init(ionicModelPtr_Type model)
{
	init();
	M_commPtr.reset( new Epetra_MpiComm(MPI_COMM_WORLD)  );
	M_meshPtr.reset( new mesh_Type( M_commPtr ) );
	M_ionicModelPtr = model;

}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setParameters()
{
	M_surfaceVolumeRatio = 2400.0;
	M_membraneCapacitance= 1.0;
	M_diffusionTensor[0] = 0.001;
	M_diffusionTensor[1] = 0.001;
	M_diffusionTensor[2] = 0.001;
	M_initialTime		 = 0.0;
	M_endTime			 = 100.0;
	M_timeStep			 = 0.01;
	M_elementsOrder		 = "P1";

}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setParameters(list_Type list)
{
	M_surfaceVolumeRatio = list.get("surfaceVolumeRatio", 2400.0  );
	M_membraneCapacitance= list.get("membraneCapacitance", 1.0 );
	M_diffusionTensor[0] = list.get("longitudinalDiffusion", 0.001 );
	M_diffusionTensor[1] = list.get("transversalDiffusion", 0.001 );
	M_diffusionTensor[2] = M_diffusionTensor[1];
	M_initialTime		 = list.get("initialTime", 0.0 );
	M_endTime			 = list.get("endTime", 100.0 );
	M_timeStep			 = list.get("diffusionCoeff", 0.01 );
	M_elementsOrder		 = list.get("elementsOrder", "P1" );

}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneReactionStepFE()
{
	M_ionicModelPtr -> superIonicModel::computeRhs( M_globalSolution, *M_appliedCurrentPtr, M_globalRhs);

	for( UInt i = 0; i < M_ionicModelPtr -> Size() ; i++ )
		*( M_globalSolution.at(i) ) = *( M_globalSolution.at(i) ) + M_timeStep * ( *( M_globalRhs.at(i) ) );

}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneDiffusionStepBE()
{
	M_linearSolverPtr -> setRightHandSide( M_rhsPtrUnique );
	M_linearSolverPtr -> solve( M_potentialPtr );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSplittingStep()
{
		solveOneReactionStepFE();
		(*M_rhsPtrUnique)*=0;
		updateRhs();
		solveOneDiffusionStepBE();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSplittingStep(exporter_Type& exporter, Real t)
{
		solveOneSplittingStep();
		exportSolution(exporter, t);
}



template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSplitting()
{
	for( Real t = M_initialTime; t < M_endTime; )
	{
		solveOneSplittingStep();
		t = t + M_timeStep;
	}
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSplitting(exporter_Type& exporter)
{
	for( Real t = M_initialTime; t < M_endTime; )
	{
		t = t + M_timeStep;
		solveOneSplittingStep(exporter, t);

	}
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupPotentialExporter(exporter_Type& exporter)
{
	exporter.setMeshProcId( M_meshPtr, M_commPtr -> MyPID() );
	exporter.setPrefix("Potential");
	exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "Potential", M_feSpacePtr, M_potentialPtr, UInt(0) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupExporter(exporter_Type& exporter)
{
	exporter.setMeshProcId( M_meshPtr, M_commPtr -> MyPID() );
	exporter.setPrefix("Solution");
	std::string variableName;
	for( int i = 0; i < M_ionicModelPtr -> Size() ; i++ ){
		variableName = "Variable" + boost::lexical_cast<std::string>( i );
		exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  variableName, M_feSpacePtr, M_globalSolution.at(i), UInt(0) );
	}
}



template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupPotentialExporter(exporter_Type& exporter, std::string fileName)
{
	exporter.setMeshProcId( M_meshPtr, M_commPtr -> MyPID() );
	exporter.setPrefix(fileName);
	exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  "Potential", M_feSpacePtr, M_potentialPtr, UInt(0) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupExporter(exporter_Type& exporter, std::string fileName)
{
	exporter.setMeshProcId( M_meshPtr, M_commPtr -> MyPID() );
	exporter.setPrefix(fileName);
	std::string variableName;
	for( int i = 0; i < M_ionicModelPtr -> Size() ; i++ ){
		variableName = "Variable" + boost::lexical_cast<std::string>( i );
		exporter.addVariable( ExporterData<mesh_Type>::ScalarField,  variableName, M_feSpacePtr, M_globalSolution.at(i), UInt(0) );
	}
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneStepGatingVariablesFE()
{
	M_ionicModelPtr -> superIonicModel::computeRhs( M_globalSolution, M_globalRhs);

	for( UInt i = 1; i < M_ionicModelPtr -> Size() ; i++ )
		*( M_globalSolution.at(i) ) = *( M_globalSolution.at(i) ) + M_timeStep * ( *( M_globalRhs.at(i) ) );

}



template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
computeRhsICI()
{
	M_ionicModelPtr -> superIonicModel::computePotentialRhsICI( M_globalSolution, (*M_appliedCurrentPtr), M_globalRhs, (*M_massMatrixPtr) );
	updateRhs();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
computeRhsSVI()
{
	M_ionicModelPtr -> superIonicModel::computePotentialRhsSVI( M_globalSolution, (*M_appliedCurrentPtr), M_globalRhs, (*M_feSpacePtr) );
	updateRhs();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneICIStep()
{
	computeRhsICI();
	M_linearSolverPtr -> setRightHandSide( M_rhsPtrUnique );
	M_linearSolverPtr -> solve( M_potentialPtr );
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSVIStep()
{
	computeRhsSVI();
	M_linearSolverPtr -> setRightHandSide( M_rhsPtrUnique );
	M_linearSolverPtr -> solve( M_potentialPtr );
}




template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveICI()
{
	for( Real t = M_initialTime; t < M_endTime; )
	{
		t = t + M_timeStep;
		solveOneStepGatingVariablesFE();
		solveOneICIStep();
	}
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSVI()
{
	for( Real t = M_initialTime; t < M_endTime; )
	{
		t = t + M_timeStep;
		solveOneStepGatingVariablesFE();
		solveOneSVIStep();
	}
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneICIStep(exporter_Type& exporter, Real t)
{
	solveOneICIStep();
	exportSolution(exporter, t);
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSVIStep(exporter_Type& exporter, Real t)
{
	solveOneSVIStep();
	exportSolution(exporter, t);
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveICI(exporter_Type& exporter)
{
	for( Real t = M_initialTime; t < M_endTime; )
	{
		t = t + M_timeStep;
		solveOneStepGatingVariablesFE();
		solveOneICIStep( exporter, t);
	}
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSVI(exporter_Type& exporter)
{
	for( Real t = M_initialTime; t < M_endTime; )
	{
		t = t + M_timeStep;
		solveOneStepGatingVariablesFE();
		solveOneSVIStep( exporter, t);
	}
}




} // namespace LifeV


#endif //_MONODOMAINSOLVER_H_
