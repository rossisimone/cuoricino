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

#ifndef _EMSOLVER_H_
#define _EMSOLVER_H_

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

#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/GeneralizedActiveHolzapfelOgdenMaterial.hpp>
#include <lifev/em/solver/EMActiveStrainSolver.hpp>


namespace LifeV {

//! EMSolver - Class featuring the solution of the electromechanical problem with monodomain equation

template<typename Mesh , typename IonicModel>
class EMSolver {

	//!Monodomain Solver
	/*!
	 The monodomain equation reads
	 \f \Chi

	 */

public:

	//! @name Type definitions
	//@{

	typedef Mesh mesh_Type;
	typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

	typedef VectorEpetra vector_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

	typedef std::vector<vectorPtr_Type> vectorOfPtr_Type;

	typedef MatrixEpetra<Real> matrix_Type;
	typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<comm_Type> commPtr_Type;

	typedef ETFESpace<mesh_Type, MapEpetra, 3, 1> ETFESpace_Type;
	typedef boost::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 1> > ETFESpacePtr_Type;

	typedef ETFESpace<mesh_Type, MapEpetra, 3, 3> ETFESpaceVectorial_Type;
	typedef boost::shared_ptr< ETFESpaceVectorial_Type > ETFESpaceVectorialPtr_Type;

	typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
	typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

	typedef LinearSolver linearSolver_Type;
	typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;

//    typedef ExporterHDF5< mesh_Type >          exporter_Type;
//    typedef boost::shared_ptr<exporter_Type>                       exporterPtr_Type;
	typedef Exporter<mesh_Type> exporter_Type;    //                IOFile_Type;
	typedef boost::shared_ptr<exporter_Type> exporterPtr_Type; //                IOFilePtr_Type;

	typedef LifeV::Preconditioner basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack prec_Type;
	typedef boost::shared_ptr<prec_Type> precPtr_Type;

	typedef IonicModel ionicModel_Type;
	typedef ElectroIonicModel superIonicModel;
	typedef boost::shared_ptr<ionicModel_Type> ionicModelPtr_Type;

	typedef Teuchos::ParameterList list_Type;

	typedef boost::function<
			Real(const Real& t, const Real& x, const Real& y, const Real& z,
					const ID& i)> function_Type;

	typedef MatrixSmall<3, 3>                          matrixSmall_Type;

	typedef ElectroETAMonodomainSolver<mesh_Type, ionicModel_Type>					monodomainSolver_Type;
	typedef boost::shared_ptr<monodomainSolver_Type>					monodomainSolverPtr_Type;

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;
    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    typedef StructuralConstitutiveLawData			structureData_Type;
    typedef boost::shared_ptr<structureData_Type>			structureDataPtr_Type;

    typedef StructuralOperator< RegionMesh<LinearTetra> > structuralOperator_Type;
    typedef boost::shared_ptr< structuralOperator_Type > structuralOperatorPtr_Type;


    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >		physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

    typedef EMActiveStrainSolver<mesh_Type> activeStrain_Type;
    typedef boost::shared_ptr< activeStrain_Type >  activeStrainPtr_Type;




	//@}

	//! @name Constructors & Destructor
	//@{
	EMSolver();

	EMSolver( 	Teuchos::ParameterList parameterList,
			    const std::string data_file_name,
     			commPtr_Type comm );	//!Empty Constructor

	void setup();

	void update();
	/*!
	 */
//	ionicModelPtr_Type			M_ionicPtr;
	monodomainSolverPtr_Type	M_monodomainPtr;
	bool						M_usingDifferentMeshes;
	structureDataPtr_Type       M_solidDataPtr;
	structuralOperatorPtr_Type  M_solidPtr;
    bcInterfacePtr_Type         M_solidBC;
    activeStrainPtr_Type		M_activationPtr;
    Real						M_monodomainTimestep;
    Real						M_solidTimestep;

};

// ===================================================
//! Constructors
// ===================================================
template<typename Mesh, typename IonicModel>
EMSolver<Mesh, IonicModel>::EMSolver():
//	M_ionicPtr(),
	M_monodomainPtr(),
	M_usingDifferentMeshes(false),
	M_solidDataPtr(),
	M_solidPtr(),
	M_solidBC(),
	M_activationPtr(),
	M_monodomainTimestep(0.01),
	M_solidTimestep(1.0)
	{}

template<typename Mesh, typename IonicModel>
EMSolver<Mesh, IonicModel>::EMSolver( 	Teuchos::ParameterList parameterList,
		const std::string data_file_name, commPtr_Type comm )
{
	GetPot dataFile (data_file_name);
	//Initializing monodomain solver
	ionicModelPtr_Type ionicPtr( new IonicModel() );
    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");
	M_monodomainPtr.reset( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicPtr ) );
	M_monodomainPtr -> setInitialConditions();
	M_usingDifferentMeshes = false;

	//Initializing structural solver

	//Material data
	M_solidDataPtr.reset( new structureData_Type() );
	M_solidDataPtr -> setup(dataFile);

	//mesh
    std::string solidMeshName = parameterList.get ("solid_mesh_name", "no_solid_mesh");
    if(solidMeshName == "no_solid_mesh")
    	solidMeshName   = dataFile ( "solid/pace_discretization/mesh_file",   "no_solid_mesh" );

    if(solidMeshName == "no_solid_mesh" || solidMeshName == meshName )
    	{
    	M_usingDifferentMeshes = false;
    	}
    else M_usingDifferentMeshes = true;

    std::string solidMeshPath = parameterList.get ("solid_mesh_path", "");
    if(solidMeshPath == "")
    	solidMeshPath   = dataFile ( "solid/pace_discretization/mesh_dir",  "" );


    meshPtr_Type fullSolidMesh;
    meshPtr_Type localSolidMesh;
    if( M_usingDifferentMeshes  )
    {
    	fullSolidMesh.reset(new mesh_Type ( comm ) );
    	localSolidMesh.reset(new mesh_Type ( comm ) );
    	MeshUtility::fillWithFullMesh (localSolidMesh, fullSolidMesh,  solidMeshName,  solidMeshPath );
    }
    else
    {
    	fullSolidMesh = M_monodomainPtr -> fullMeshPtr();
    	localSolidMesh = M_monodomainPtr -> localMeshPtr();
    }

    //FESPACEs
    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh,
         														dOrder,
         														3,
         														localSolidMesh -> comm() ) );

    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, &feTetraP1, comm) );

    //boundary conditions
    M_solidBC.reset( new bcInterface_Type() );
    M_solidBC->createHandler();
    M_solidBC->fillHandler ( data_file_name, "solid" );

    //setup structural operator
    M_solidPtr -> setup (M_solidDataPtr,
            			dFESpace,
            			dETFESpace,
            			M_solidBC -> handler(),
            			comm);
    M_solidPtr -> setDataFromGetPot (dataFile);

    //activation
    M_activationPtr.reset( new activeStrain_Type(parameterList, comm) );

    //    solidFESpacePtr_Type aFESpace ( new solidFESpace_Type (M_monodomainPtr -> localMeshPtr(), dOrder, 1, comm) );
//    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (M_monodomainPtr -> localMeshPtr(), & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
//    scalarETFESpacePtr_Type aETFESpace ( new scalarETFESpace_Type (M_monodomainPtr -> localMeshPtr(), & (aFESpace->refFE() ), & (aFESpace->fe().geoMap() ), comm) );



}




} // namespace LifeV

#endif //_MONODOMAINSOLVER_H_
