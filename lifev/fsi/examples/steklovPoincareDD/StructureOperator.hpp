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

#include <lifev/core/fem/TimeAdvanceBDF.hpp>

namespace LifeV
{
    
typedef RegionMesh<LinearTetra>      mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
typedef MeshPartitioner<mesh_Type >  meshPartitioner_Type;

typedef FESpace<mesh_Type, MapEpetra >  FESpace_Type;
typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;

typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 > solidETFESpace_Type;
typedef boost::shared_ptr<solidETFESpace_Type>                solidETFESpacePtr_Type;

typedef VectorEpetra                    vector_Type;
typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
    
typedef boost::shared_ptr< TimeAdvanceBDF< vector_Type > > timeAdvancePtr_type;    
    
typedef boost::shared_ptr<Epetra_Comm> communicatorPtr_Type;

class StructureOperator
{
public:
    
    StructureOperator(boost::shared_ptr<Epetra_Comm>& comm);
    ~StructureOperator();
    
    void setup(const GetPot& dataFile);
    
    void setBC(const boost::shared_ptr<BCHandler> bc);
    
    void buildSystem(const GetPot& dataFile);
    
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
    timeAdvancePtr_type M_timeAdvance;
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
    
void StructureOperator::setBC(const boost::shared_ptr<BCHandler> bc)
{
    M_BCh.reset(new BCHandler(*bc));
}
    
void StructureOperator::buildSystem(const GetPot& dataFile)
{
    M_timeAdvance.reset(new TimeAdvanceBDF<vector_Type>() );
    M_timeAdvance->setup (M_dataStructure->dataTimeAdvance()->orderBDF() , 2);
    M_timeAdvance->setTimeStep (M_dataStructure->dataTime()->timeStep() );
    
    M_solid->setup (M_dataStructure, M_dFESpace, M_dETFESpace, M_BCh, M_comm);
    M_solid->setDataFromGetPot (dataFile);
    M_solid->buildSystem(M_timeAdvance->coefficientSecondDerivative(0)/(M_dataStructure->dataTime()->timeStep()*M_dataStructure->dataTime()->timeStep()));
}    
    
} // end namespace LifeV

#endif
