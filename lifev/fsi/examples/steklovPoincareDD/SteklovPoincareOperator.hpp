#ifndef STEKLOVPOINCAREOPERATOR_H
#define STEKLOVPOINCAREOPERATOR_H 1

// Includes for the interface
#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>

namespace LifeV
{

typedef RegionMesh<LinearTetra>      			mesh_Type;
typedef boost::shared_ptr<mesh_Type> 			meshPtr_Type;
typedef MeshPartitioner<mesh_Type >  			meshPartitioner_Type;


typedef FESpace<mesh_Type, MapEpetra >  		FESpace_Type;
typedef boost::shared_ptr<FESpace_Type> 		FESpacePtr_Type;

typedef MapEpetra 								map_Type;
typedef boost::shared_ptr<map_Type> 			mapPtr_Type;

typedef VectorEpetra 				   			vector_Type;
typedef boost::shared_ptr<vector_Type> 			vectorPtr_Type;

typedef MatrixEpetra<Real>             			matrix_Type;
typedef boost::shared_ptr<matrix_Type> 			matrixPtr_Type;

typedef RBFInterpolation<mesh_Type>             interpolation_Type;
typedef boost::shared_ptr<interpolation_Type>   interpolationPtr_Type;

typedef boost::shared_ptr<Epetra_Comm> 			communicatorPtr_Type;

class SteklovPoincareOperator
{
public:

	SteklovPoincareOperator ( );

    ~SteklovPoincareOperator ( );

    void buildTranferOperators(	const meshPtr_Type& fullMeshFluid, const meshPtr_Type& localMeshFluid,
								const meshPtr_Type& fullMeshSructure, const meshPtr_Type& localMeshSructure,
								const std::vector<int> flags, const GetPot data_file);

    void setUpData( const vectorPtr_Type& fluidVec, const vectorPtr_Type& structureVec);

    void transferFluidOnStructure (const vectorPtr_Type& vecFluid, vectorPtr_Type& vecStructure);

    void transferGammaFluidOnGammaStructure (const vectorPtr_Type& vecGammaFluid, vectorPtr_Type& vecGammaStructure);

    // Getters:

    mapPtr_Type fluidInterfaceMap()
    {
    	return M_interfaceMapFluid;
    }

    mapPtr_Type structureInterfaceMap()
    {
    	return M_interfaceMapStructure;
    }

private:

    // Private methods

    void createInterfaceMaps ( );


    // Private members

    mapPtr_Type M_interfaceMapStructure;
    mapPtr_Type M_interfaceMapFluid;

    mapPtr_Type M_interfaceStructureMapScalar;
    mapPtr_Type M_interfaceFluidMapScalar;

    // Transfer operators
    interpolationPtr_Type M_structureToFluid;
    interpolationPtr_Type M_fluidToStructure;
};

SteklovPoincareOperator::SteklovPoincareOperator()
{}

SteklovPoincareOperator::~SteklovPoincareOperator()
{}

void SteklovPoincareOperator::buildTranferOperators(const meshPtr_Type& fullMeshFluid, const meshPtr_Type& localMeshFluid,
													const meshPtr_Type& fullMeshSructure, const meshPtr_Type& localMeshSructure,
													const std::vector<int> flags, const GetPot data_file)
{
	M_structureToFluid.reset ( interpolation_Type::InterpolationFactory::instance().createObject (data_file("interpolation/interpolation_Type","none")));
	M_structureToFluid->setup(fullMeshSructure, localMeshSructure, fullMeshFluid, localMeshFluid, flags);

	M_fluidToStructure.reset ( interpolation_Type::InterpolationFactory::instance().createObject (data_file("interpolation/interpolation_Type","none")));
	M_fluidToStructure->setup(fullMeshFluid, localMeshFluid, fullMeshSructure, localMeshSructure, flags);
}

void SteklovPoincareOperator::setUpData( const vectorPtr_Type& fluidVec,
										 const vectorPtr_Type& structureVec)
{
	M_structureToFluid->setupRBFData (structureVec, fluidVec);
	M_structureToFluid->buildOperators ();

	M_fluidToStructure->setupRBFData (fluidVec, structureVec);
	M_fluidToStructure->buildOperators ();

	createInterfaceMaps();
}

void
SteklovPoincareOperator::createInterfaceMaps ( )
{
	M_structureToFluid->buildUnknownVectorialInterfaceMap();
	M_interfaceMapFluid.reset ( new map_Type(*M_structureToFluid->unknownVectorialInterfaceMap() ) );

	M_fluidToStructure->buildUnknownVectorialInterfaceMap();
	M_interfaceMapStructure.reset ( new map_Type( *M_fluidToStructure->unknownVectorialInterfaceMap() ) );
}

void
SteklovPoincareOperator::transferFluidOnStructure (const vectorPtr_Type& vecFluid, vectorPtr_Type& vecStructure)
{
	M_fluidToStructure->updateRhs (vecFluid);
	M_fluidToStructure->interpolate ( );
	M_fluidToStructure->solution ( vecStructure );
}

void
SteklovPoincareOperator::transferGammaFluidOnGammaStructure (const vectorPtr_Type& vecGammaFluid, vectorPtr_Type& vecGammaStructure)
{
	vectorPtr_Type fullFluid;
	vectorPtr_Type fullStructure;

	fullFluid.reset( new vector_Type ( *M_fluidToStructure->knownMap() ) );
	fullStructure.reset( new vector_Type ( *M_fluidToStructure->unknownMap() ) );

	fullFluid->subset ( *vecGammaFluid, *M_interfaceMapFluid, 0, 0 );
	transferFluidOnStructure ( fullFluid, fullStructure);

	*vecGammaStructure *= 0;
	vecGammaStructure->subset( *fullStructure, *M_interfaceMapStructure, 0, 0);
}

} // end namespace LifeV

#endif
