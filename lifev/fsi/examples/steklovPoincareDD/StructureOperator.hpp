#ifndef STRUCTUREOPERATOR_H
#define STRUCTUREOPERATOR_H 1

//datafile
#include <lifev/core/filter/GetPot.hpp>

//fespace
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

namespace LifeV
{

Real displ (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
	return 0.001;
}

Real fZer (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

typedef RegionMesh<LinearTetra>      mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
typedef MeshPartitioner<mesh_Type >  meshPartitioner_Type;

typedef FESpace<mesh_Type, MapEpetra >  FESpace_Type;
typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;

typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 > solidETFESpace_Type;
typedef boost::shared_ptr<solidETFESpace_Type>                solidETFESpacePtr_Type;

typedef VectorEpetra                    vector_Type;
typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
    
typedef TimeAdvance<vector_Type> 			timeAdvance_Type;
typedef boost::shared_ptr<timeAdvance_Type>	timeAdvancePtr_Type;

typedef boost::shared_ptr<Epetra_Comm> communicatorPtr_Type;

typedef BCHandler                                          bc_Type;
typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
typedef StructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

class StructureOperator
{
public:

    StructureOperator(boost::shared_ptr<Epetra_Comm>& comm);

    ~StructureOperator();
    
    void setup(const GetPot& dataFile);
    
    void setBC();
    
    void buildSystem(const GetPot& dataFile, const timeAdvancePtr_Type& timeAdvance);
    
    void iterate(const timeAdvancePtr_Type& timeAdvance);

    // Getters

    boost::shared_ptr<StructuralOperator<mesh_Type > > solver()
	{
    	return M_solid;
	}

    boost::shared_ptr<BCHandler> bc()
	{
    	return M_BCh;
	}

    FESpacePtr_Type feDisplacementSerial()
    {
    	return M_dFESpaceSerial;
    }

    FESpacePtr_Type feDisplacement()
    {
    	return M_dFESpace;
    }

    meshPtr_Type mesh()
    {
    	return M_fullMeshPtrStructure;
    }

private:
    
    // Methods
    void loadData(const GetPot& dataFile);
    
    void loadMesh();
    
    void partitionMesh();
    
    void createFESpaces();
    
    // Members
    communicatorPtr_Type M_comm;
    bool M_verbose;
    meshPtr_Type M_fullMeshPtrStructure;
    meshPtr_Type M_localMeshPtrStructure;
    boost::shared_ptr<StructuralConstitutiveLawData> M_dataStructure;
    boost::shared_ptr<meshPartitioner_Type> M_meshPartStructure;
    boost::shared_ptr<MeshData> M_meshDataStructure;
    FESpacePtr_Type M_dFESpace;
    FESpacePtr_Type M_dFESpaceSerial;
    std::string M_dOrder;
    solidETFESpacePtr_Type M_dETFESpace;
    boost::shared_ptr<BCHandler> M_BCh;
    boost::shared_ptr<StructuralOperator<mesh_Type > > M_solid;
    bcInterfacePtr_Type M_solidBCPtr;
};

StructureOperator::StructureOperator(boost::shared_ptr<Epetra_Comm>& comm):
M_comm (comm)
{
    M_verbose = M_comm->MyPID()==0;
}

StructureOperator::~StructureOperator()
{}

void StructureOperator::setup(const GetPot& dataFile)
{
    loadData(dataFile);
    loadMesh();
    partitionMesh();
    createFESpaces();
}

void StructureOperator::loadData(const GetPot& dataFile)
{
    M_dataStructure.reset(new StructuralConstitutiveLawData( ) );
    M_dataStructure->setup (dataFile);
    M_meshDataStructure.reset(new MeshData());
    M_meshDataStructure->setup (dataFile, "solid/space_discretization");
    
    M_dOrder = dataFile("solid/space_discretization/vel_order","P1");
    
    // Initialize the boost pointer of the structural operator
    M_solid.reset(new StructuralOperator<mesh_Type>());

    M_solidBCPtr.reset ( new bcInterface_Type() );
    M_solidBCPtr->createHandler();
    M_solidBCPtr->fillHandler ( "datanew", "solid" );
}

void StructureOperator::loadMesh()
{
    M_fullMeshPtrStructure.reset(new mesh_Type());
    readMesh(*M_fullMeshPtrStructure, *M_meshDataStructure);
    
    if (M_verbose)
        std::cout << "Mesh source: file(" << M_meshDataStructure->meshDir() << M_meshDataStructure->meshFile() << ")" << std::endl;
    
    if (M_verbose){
        std::cout << "Mesh statistics: " << std::endl;
        std::cout << "Mesh size (Hmax) : " << MeshUtility::MeshStatistics::computeSize(*M_fullMeshPtrStructure).maxH << std::endl;
        std::cout << "Mesh size (Hmin) : " << MeshUtility::MeshStatistics::computeSize(*M_fullMeshPtrStructure).minH << std::endl;
    }
    
}

void StructureOperator::partitionMesh()
{
    M_meshPartStructure.reset( new meshPartitioner_Type(M_fullMeshPtrStructure, M_comm));
    M_localMeshPtrStructure.reset(new mesh_Type(*M_meshPartStructure->meshPartition() ) );
}

void StructureOperator::createFESpaces()
{
    if (M_verbose)
        std::cout << std::endl << "[Creating the Structure FE space]" << std::endl;
    
    if (M_verbose)
        std::cout << "FE for the displacement: " << M_dOrder << "  ... " << std::flush;
            
    M_dFESpace.reset (new FESpace_Type (M_localMeshPtrStructure, M_dOrder, 3, M_comm) );
    
    M_dFESpaceSerial.reset (new FESpace_Type (M_fullMeshPtrStructure, M_dOrder, 3, M_comm) );

    if (M_verbose)
        std::cout << "ok.\n";
    
    M_dETFESpace.reset( new solidETFESpace_Type (*M_meshPartStructure, & (M_dFESpace->refFE() ), & (M_dFESpace->fe().geoMap() ), M_comm) );
}
    
void StructureOperator::setBC()
{
	vector <ID> compy (1);
	compy[0] = 1;

	BCFunctionBase bcf (fZer);
	BCFunctionBase disp (displ);

	M_BCh.reset(new BCHandler());
	M_BCh->addBC ("Top",   2,  Essential, Full, bcf, 3);
	M_BCh->addBC ("Base",  3, Essential, Full, bcf, 3);
	M_BCh->addBC ("Interface",  1, Essential, Component, disp, compy);
}
    
void StructureOperator::buildSystem(const GetPot& dataFile, const timeAdvancePtr_Type& timeAdvance)
{
    //setBC();
	BCFunctionBase disp (displ);
	M_solidBCPtr -> handler() -> addBC ( "Interface",  1, Essential, Full, disp, 3 );
	M_solid->setup (M_dataStructure, M_dFESpace, M_dETFESpace, M_solidBCPtr -> handler() /*M_BCh*/, M_comm);
    M_solid->setDataFromGetPot (dataFile);
    M_solid->buildSystem(timeAdvance->coefficientSecondDerivative(0)/(M_dataStructure->dataTime()->timeStep()*M_dataStructure->dataTime()->timeStep()));

    vectorPtr_Type a;
    a.reset( new vector_Type(M_solid->displacement(), Unique));
    M_solid->initialize ( a );
}    

void StructureOperator::iterate( const timeAdvancePtr_Type& timeAdvance)
{
    vector_Type rhs( M_dFESpace->map() );
    timeAdvance->updateRHSContribution ( 0.001 );
    Real timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative(0)/(M_dataStructure->dataTime()->timeStep()*M_dataStructure->dataTime()->timeStep());
    rhs += *M_solid->massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;
    M_solid->setRightHandSide ( rhs );

	M_solid->iterate ( M_solidBCPtr -> handler() /*M_BCh*/ );
}
    
} // end namespace LifeV

#endif
