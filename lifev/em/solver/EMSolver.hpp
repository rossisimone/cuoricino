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
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledScalar.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledScalar.hpp>
//#include <lifev/core/interpolation/RBFscalar.hpp>
#include <lifev/core/interpolation/RBFvectorial.hpp>


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

    typedef ExporterHDF5< mesh_Type >          exporter_Type;
    typedef boost::shared_ptr<exporter_Type>                       exporterPtr_Type;
//	typedef Exporter<mesh_Type> exporter_Type;    //                IOFile_Type;
//	typedef boost::shared_ptr<exporter_Type> exporterPtr_Type; //                IOFilePtr_Type;

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

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>                        FESpacePtr_Type;

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

    typedef RBFInterpolation<mesh_Type>           interpolation_Type;
    typedef boost::shared_ptr<interpolation_Type> interpolationPtr_Type;
    ///////////////////////////////////////////////////////////////////////////

	inline monodomainSolverPtr_Type    monodomainPtr()        { return 	M_monodomainPtr;}
	inline bool						usingDifferentMeshes(){ return M_usingDifferentMeshes; }
	inline structureDataPtr_Type       solidDataPtr()        	{ return M_solidDataPtr; }
	inline structuralOperatorPtr_Type  solidPtr()             	{ return M_solidPtr; }
	inline bcInterfacePtr_Type         solidBCPtr()              	{ return M_solidBCPtr; }
	inline activeStrainPtr_Type 	    activationPtr()     	{ return M_activationPtr;}
	inline Real		                    monodomainTimeStep() 	{ return M_monodomainTimeStep;}
	inline Real		                    solidTimeStep()       	{ return M_solidTimeStep; }

	inline void setMonodomainPtr(monodomainSolverPtr_Type p) { M_monodomainPtr = p;}
	inline void setMonodomainPtr(monodomainSolver_Type& p)   { *M_monodomainPtr = p;}
	inline void setUsingDifferentMeshes(bool p){  M_usingDifferentMeshes = p; }
	inline void setSolidDataPtr(structureDataPtr_Type p ) 	{ M_solidDataPtr = p; }
	inline void setSolidDataPtr(structureData_Type& p ) 	{ *M_solidDataPtr = p; }
	inline void setSolidPtr(structuralOperatorPtr_Type p)  	{ M_solidPtr = p; }
	inline void setSolidPtr(structuralOperator_Type& p)  	{ *M_solidPtr = p; }
	inline void setSolidBCPtr(bcInterfacePtr_Type p)        	{ M_solidBCPtr = p; }
	inline void setSolidBCPtr(bcInterface_Type& p)        	{ *M_solidBCPtr = p; }
	inline void setActivationPtr(activeStrainPtr_Type p)   	{ M_activationPtr = p;}
	inline void setActivationPtr(activeStrain_Type& p)   	{ *M_activationPtr = p;}
	inline void setMonodomainTimeStep(Real p) 	{ M_monodomainTimeStep = p;}
	inline void setSolidTimeStep(Real p)       	{ M_solidTimeStep = p; }


	//@}

	//! @name Constructors & Destructor
	//@{
	EMSolver();

	EMSolver(  Teuchos::ParameterList& parameterList,
				const std::string data_file_name,
				commPtr_Type comm );	//!Empty Constructor

	virtual ~EMSolver() {};

	void setup(Teuchos::ParameterList& parameterList,
				const std::string data_file_name,
				commPtr_Type comm,
				std::string parameterListName = "ParamList.xml");

	void setupMonodomainMatrix(Teuchos::ParameterList& parameterList);

	void updateMonodomainMatrix();

	void setupExporters(commPtr_Type comm, std::string dir = "./");

	void exportSolution(Real time = 0.0);

	void closeExporters();

	//void setupPreloadBC( Teuchos::ParameterList& parameterList );

	void exportSolidFibersDirection(commPtr_Type comm, std::string dir = "./" );
	void exportMonodomainFibersDirection(std::string dir = "./");
	void exportSolidSheetsDirection(commPtr_Type comm, std::string dir = "./" );
	void exportFibersAndSheetsFields(commPtr_Type comm, std::string dir = "./" );
	void exportActivationTime(commPtr_Type comm, std::string dir = "./" );

	inline void importSolidFibers(Teuchos::ParameterList& parameterList);
	inline void importSolidSheets(Teuchos::ParameterList& parameterList);
	inline void importMonodomainFibers(Teuchos::ParameterList& parameterList);

	void setFibersAndSheets(Teuchos::ParameterList& parameterList);

	void setupInterpolants(std::string parameterListName, Teuchos::ParameterList& parameterList, GetPot& dataFile);

	inline void registerActivationTime( Real time, Real threshold = 0.0)
	{
		M_monodomainPtr -> registerActivationTime(*M_activationTimePtr, time, threshold);
	}

	inline void setSolidFibers(vector_Type& fibers)
	{
		M_solidPtr -> material() -> setFiberVector(fibers);
	}
	inline void setSolidFibers(vectorPtr_Type fibers)
	{
		M_solidPtr -> material() -> setFiberVector(*fibers);
	}
	inline void setSolidSheets(vector_Type& sheets)
	{
		M_solidPtr -> material() -> setSheetVector(sheets);
	}
	inline void setSolidSheets(vectorPtr_Type sheets)
	{
		M_solidPtr -> material() -> setSheetVector(*sheets);
	}
	//monodomain -> setFiberPtr( electroFibers )
	inline void setMonodomainFibers(vectorPtr_Type fibers)
	{
		M_monodomainPtr -> setFiberPtr(fibers);
	}
	inline void setMonodomainFibers(vector_Type& fibers)
	{
		M_monodomainPtr -> setFiber(fibers);
	}


//	void update();

	inline void  setPotentialOnBoundary(Real value, UInt flag)
	{
		HeartUtility::setValueOnBoundary( *(M_monodomainPtr -> potentialPtr() ), M_monodomainPtr -> fullMeshPtr(), value, flag);
	}
	inline void  setPotentialFromFunction(function_Type f, Real time = 0.0)
	{
		M_monodomainPtr -> setPotentialFromFunction( f, time );
	}

		/*!
	 */
	monodomainSolverPtr_Type	M_monodomainPtr;
	bool						M_usingDifferentMeshes;
	structureDataPtr_Type       M_solidDataPtr;
	structuralOperatorPtr_Type  M_solidPtr;
    bcInterfacePtr_Type         M_solidBCPtr;
    activeStrainPtr_Type		M_activationPtr;
    Real						M_monodomainTimeStep;
    Real						M_solidTimeStep;
    exporterPtr_Type			M_monodomainExporterPtr;
    exporterPtr_Type			M_solidExporterPtr;
    exporterPtr_Type			M_activationExporterPtr;

    //Coarse To Fine ( C2F )
    vectorPtr_Type				M_monodomainDisplacementPtr;
    interpolationPtr_Type 		M_C2F;
    //Fine To Coarse ( F2C )
    vectorPtr_Type 				M_activationSolidPtr;
    interpolationPtr_Type 		M_F2C;
    meshPtr_Type				M_fullSolidMesh;

    vectorPtr_Type				M_activationTimePtr;


private:
//    void initSolid();
//    void initMonodomain();
//    void initActivation();

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
	M_solidBCPtr(),
	M_activationPtr(),
	M_monodomainTimeStep(0.01),
	M_solidTimeStep(1.0),
	M_monodomainExporterPtr(),
	M_solidExporterPtr(),
	M_activationExporterPtr(),
	M_monodomainDisplacementPtr(),
    M_C2F(),
    M_activationSolidPtr(),
    M_F2C(),
    M_fullSolidMesh(),
    M_activationTimePtr()
	{}

template<typename Mesh, typename IonicModel>
EMSolver<Mesh, IonicModel>::EMSolver( 	Teuchos::ParameterList& parameterList,
		const std::string data_file_name, commPtr_Type comm )
{
	if(comm->MyPID()==0)
	{
		std::cout << "\n==========================================";
		std::cout << "\n\t EM SOLVER: 'YOU ROCK!!!!' ";
		std::cout << "\n==========================================";
	}
	M_solidTimeStep = parameterList.get("emdt",1.0);

	GetPot dataFile (data_file_name);
	//Initializing monodomain solver
	if(comm->MyPID()==0)
	{
		std::cout << "\n==========================================";
		std::cout << "\n\t Initializing Monodomain Solver";
		std::cout << "\n==========================================";
	}

    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");
	ionicModelPtr_Type ionicPtr( new IonicModel() );
	M_monodomainPtr.reset( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicPtr ) );
	M_monodomainPtr -> setInitialConditions();
	M_monodomainPtr-> setParameters ( parameterList );

	M_usingDifferentMeshes = false;

	//Initializing structural solver

	if(comm->MyPID()==0)
	{
		std::cout << "\n==========================================";
		std::cout << "\n\t Initializing Solid Data and Solid mesh";
		std::cout << "\n==========================================";
	}

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


    M_fullSolidMesh.reset(new mesh_Type( comm ) );
    meshPtr_Type localSolidMesh(new mesh_Type( comm ) );
    if( M_usingDifferentMeshes  )
    {
    	localSolidMesh.reset(new mesh_Type ( comm ) );
    	MeshUtility::fillWithFullMesh (localSolidMesh, M_fullSolidMesh,  solidMeshName,  solidMeshPath );
    }
    else
    {
    	M_fullSolidMesh = M_monodomainPtr -> fullMeshPtr();
    	localSolidMesh = M_monodomainPtr -> localMeshPtr();
    }

    //FESPACEs
    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    FESpacePtr_Type dFESpace ( new FESpace_Type ( localSolidMesh,
         										   dOrder,
         										   3,
         										   localSolidMesh -> comm() ) );

    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, &feTetraP1, comm) );

    //boundary conditions
	if(comm->MyPID()==0)
	{
		std::cout << "\n==========================================";
		std::cout << "\n\t Creating BC handler";
		std::cout << "\n==========================================";
	}

    M_solidBCPtr.reset( new bcInterface_Type() );
    M_solidBCPtr->createHandler();
    M_solidBCPtr->fillHandler ( data_file_name, "solid" );

    //setup structural operator
	if(comm->MyPID()==0)
	{
		std::cout << "\n==========================================";
		std::cout << "\n\t Initializing Structure Solver";
		std::cout << "\n==========================================";
	}
    M_solidPtr.reset(new structuralOperator_Type() );
    M_solidPtr -> setup (M_solidDataPtr,
            			dFESpace,
            			dETFESpace,
            			M_solidBCPtr -> handler(),
            			comm);
    M_solidPtr -> setDataFromGetPot (dataFile);

    //activation
	if(comm->MyPID()==0)
	{
		std::cout << "\n==========================================";
		std::cout << "\n\t Initializing Activation Solver";
		std::cout << "\n==========================================";
	}
    M_activationPtr.reset( new activeStrain_Type(parameterList, dataFile, M_monodomainPtr -> localMeshPtr(), comm) );


    M_activationTimePtr.reset(new vector_Type( M_monodomainPtr -> feSpacePtr() -> map() ) );
    *M_activationTimePtr = -1.0;
}


template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::setup(Teuchos::ParameterList& parameterList,
											const std::string data_file_name,
											commPtr_Type comm,
											std::string parameterListName = "ParamList.xml")
{
	if(M_usingDifferentMeshes)
	{
		M_monodomainDisplacementPtr.reset( new vector_Type ( M_monodomainPtr -> displacementETFESpacePtr() -> map() ) );
		M_activationSolidPtr.reset( new vector_Type( M_solidPtr -> material() -> activationSpace() -> map() ) );
		GetPot dataFile(data_file_name);
    	setupInterpolants(parameterListName, parameterList, dataFile);
	}
	else
	{
		M_monodomainDisplacementPtr = M_solidPtr -> displacementPtr();
		M_activationSolidPtr = M_activationPtr -> gammafPtr();
	}

	setupMonodomainMatrix(parameterList);

}

template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::setupMonodomainMatrix(Teuchos::ParameterList& parameterList)
{
	bool lumpedMass = parameterList.get ("LumpedMass", true);
	if(lumpedMass) M_monodomainPtr -> setupLumpedMassMatrix();
	else M_monodomainPtr -> setupMassMatrix();

	M_monodomainPtr -> setDisplacementPtr( M_monodomainDisplacementPtr );
	M_monodomainPtr -> setupStiffnessMatrix();
	M_monodomainPtr -> setupGlobalMatrix();
}


template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::updateMonodomainMatrix()
{
	M_monodomainPtr -> setDisplacementPtr( M_monodomainDisplacementPtr );
	M_monodomainPtr -> setupStiffnessMatrix();
	M_monodomainPtr -> setupGlobalMatrix();
}

template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::setupExporters(commPtr_Type comm, std::string dir)
{
	if(comm->MyPID()==0)
	{
		std::cout << "\n==========================================";
		std::cout << "\n\t Setting up the exporters";
		std::cout << "\n==========================================";
	}
	M_monodomainExporterPtr.reset(new exporter_Type() );
	M_monodomainPtr -> setupExporter(*M_monodomainExporterPtr, "ElectroOutput", dir);
	M_activationExporterPtr.reset(new exporter_Type() );
	M_activationPtr -> setupExporter(*M_activationExporterPtr, comm, dir);

	M_solidExporterPtr.reset(new exporter_Type() );
	M_solidExporterPtr -> setMeshProcId( M_solidPtr -> mesh(), comm->MyPID());
	M_solidExporterPtr -> setPrefix( "StructureOutput" );
	M_solidExporterPtr -> setPostDir ( dir );
	M_solidExporterPtr -> addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", M_solidPtr -> dispFESpacePtr(), M_solidPtr -> displacementPtr(), UInt (0) );


	if(M_usingDifferentMeshes)
	{
		FESpacePtr_Type gfSolidFESpace ( new FESpace_Type ( M_solidPtr -> dispFESpace().mesh(),
															"P1", 	3,   comm ) );
		M_solidExporterPtr -> addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField,
											"gammaf",
											gfSolidFESpace,
											M_solidPtr -> material() -> gammaf(),
											UInt (0) );
		FESpacePtr_Type displacementMonodomainFESpace ( new FESpace_Type ( M_monodomainPtr -> localMeshPtr(),
																			"P1", 	3,   comm ) );
		M_monodomainExporterPtr -> addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField,
												"interpolated_displacement",
												displacementMonodomainFESpace,
												M_monodomainPtr -> displacementPtr(),
												UInt (0) );
	}
}


template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::exportSolution( Real time)
{
	M_monodomainExporterPtr -> postProcess(time);
	M_activationExporterPtr -> postProcess(time);
	M_solidExporterPtr -> postProcess(time);
}


template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::closeExporters()
{
	M_monodomainExporterPtr -> closeFile();
	M_activationExporterPtr -> closeFile();
	M_solidExporterPtr -> closeFile();
}

template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::exportSolidFibersDirection(commPtr_Type comm, std::string dir )
{
	exporter_Type exp;
	exp.setMeshProcId( M_solidPtr -> mesh(), comm->MyPID());
	exp.setPostDir ( dir );
	exp.setPrefix("SolidFibesrDirection");
	exp.addVariable(ExporterData<mesh_Type>::VectorField,
					"solid_fibers",
					M_solidPtr -> dispFESpacePtr(),
					M_solidPtr -> material() -> fiberVectorPtr(), UInt(0));
	exp.postProcess(0);
	exp.closeFile();
}


template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::exportMonodomainFibersDirection(std::string dir)
{
	M_monodomainPtr -> exportFiberDirection(dir);
}


template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::exportSolidSheetsDirection(commPtr_Type comm, std::string dir )
{
	exporter_Type exp;
	exp.setMeshProcId( M_solidPtr -> mesh(), comm->MyPID());
	exp.setPostDir ( dir );
	exp.setPrefix("SolidSheetsDirection");
	exp.addVariable(ExporterData<mesh_Type>::VectorField,
					"solid_sheets",
					M_solidPtr -> dispFESpacePtr(),
					M_solidPtr -> material() -> sheetVectorPtr(), UInt(0));
	exp.postProcess(0);
	exp.closeFile();
}

template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::exportFibersAndSheetsFields(commPtr_Type comm, std::string dir )
{
	exportSolidSheetsDirection( comm, dir);
	exportSolidFibersDirection( comm, dir);
	exportMonodomainFibersDirection(dir);
}


template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::exportActivationTime(commPtr_Type comm, std::string dir )
{
	exporter_Type exp;
	exp.setMeshProcId( M_monodomainPtr -> localMeshPtr(), comm->MyPID());
	exp.setPostDir ( dir );
	exp.setPrefix("ActivationTime");
	exp.addVariable(ExporterData<mesh_Type>::ScalarField,
					"activation_time",
					M_monodomainPtr -> feSpacePtr(),
					M_activationTimePtr, UInt(0));
	exp.postProcess(0);
	exp.closeFile();
}

template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::importSolidFibers(Teuchos::ParameterList& parameterList)
{
	std::string solidFibersFile = parameterList.get ("solid_fibers_file", "");
	std::string solidFibersField = parameterList.get ("solid_fibers_field", "");
    HeartUtility::importVectorField( M_solidPtr -> material() -> fiberVectorPtr(),
    		                         solidFibersFile,
    		                         solidFibersField,
    		                         M_solidPtr -> mesh() );
}
template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::importSolidSheets(Teuchos::ParameterList& parameterList)
{
	std::string solidSheetsFile = parameterList.get ("solid_sheets_file", "");
	std::string solidSheetsField = parameterList.get ("solid_sheets_field", "");
    HeartUtility::importVectorField( M_solidPtr -> material() -> sheetVectorPtr(),
    		                         solidSheetsFile,
    		                         solidSheetsField,
    		                         M_solidPtr -> mesh() );
}
template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::importMonodomainFibers(Teuchos::ParameterList& parameterList)
{
	std::string fibersFile = parameterList.get ("fibers_file", "");
	std::string fibersField = parameterList.get ("fibers_field", "");
	vectorPtr_Type monodomainFibers(new vector_Type( M_monodomainPtr -> displacementETFESpacePtr() -> map() ) );
	HeartUtility::importFibers(monodomainFibers, fibersFile, M_monodomainPtr -> localMeshPtr() );
//    HeartUtility::importVectorField( monodomainFibers,
//    		                         fibersFile,
//    		                         fibersField,
//    		                         M_monodomainPtr -> localMeshPtr() );
    setMonodomainFibers(monodomainFibers);
}

template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::setFibersAndSheets(Teuchos::ParameterList& parameterList)
{
	std::cout << "\nImporting Solid Fibers\n\n";
	this->importSolidFibers(parameterList);

	std::cout << "\nImporting Solid Sheets\n\n";
	importSolidSheets(parameterList);

	std::cout << "\nImporting Monodomain Fibers\n\n";
	importMonodomainFibers(parameterList);
}

template<typename Mesh, typename IonicModel>
void EMSolver<Mesh, IonicModel>::setupInterpolants(std::string parameterListName, Teuchos::ParameterList& parameterList, GetPot& dataFile)
{
	Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
	belosList = Teuchos::getParametersFromXmlFile ( parameterListName );

	int nFlags = 1;
	std::vector<int> flags (nFlags);
	flags[0] = -1;

	std::string c2f = parameterList.get ("c2f", "RBFrescaledVectorial");
	M_C2F.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( c2f ) );
	M_C2F->setup( M_fullSolidMesh,
				  M_solidPtr -> mesh(),
				  M_monodomainPtr -> fullMeshPtr(),
				  M_monodomainPtr -> localMeshPtr(),
				  flags);
	M_C2F -> setRadius( 2.0 * (double) MeshUtility::MeshStatistics::computeSize (* (M_fullSolidMesh) ).maxH );
	M_C2F -> setupRBFData ( M_solidPtr -> displacementPtr(), M_monodomainDisplacementPtr, dataFile, belosList);
	if(c2f == "RBFvectorial") M_C2F->setBasis("TPS");
	M_C2F->buildOperators();
	M_C2F->interpolate();
	M_C2F->solution (M_monodomainDisplacementPtr);



	std::string f2c = parameterList.get ("f2c", "RBFrescaledScalar");
	M_F2C.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( f2c ) );
	M_F2C->setup( M_monodomainPtr -> fullMeshPtr(),
				  M_monodomainPtr -> localMeshPtr(),
				  M_fullSolidMesh,
				  M_solidPtr -> mesh(),
				  flags);
	//WARNING
	std::cout<< "\nWARNING!!! Setting the Radius of interpolation using the full monodomain mesh.";
	std::cout<< "\nWARNING!!! You shoul use the full activation mesh, but it's not coded yet...";

	M_F2C -> setRadius( (double) MeshUtility::MeshStatistics::computeSize (* ( M_monodomainPtr -> fullMeshPtr()) ).maxH );
	M_F2C -> setupRBFData ( M_activationPtr -> gammafPtr(), M_activationSolidPtr , dataFile, belosList);
	M_F2C -> buildOperators();
	M_F2C->interpolate();
	M_F2C->solution (M_activationSolidPtr);

}




} // namespace LifeV

#endif //_MONODOMAINSOLVER_H_
