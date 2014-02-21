#ifndef FLUIDOPERATOR_H
#define FLUIDOPERATOR_H 1

#include <lifev/core/filter/GetPot.hpp>

//fespace
#include <lifev/core/fem/FESpace.hpp>

// fluid data
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/navier_stokes/solver/OseenSolver.hpp>

namespace LifeV
{
    
typedef RegionMesh<LinearTetra>      mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
typedef MeshPartitioner<mesh_Type >  meshPartitioner_Type;


typedef FESpace<mesh_Type, MapEpetra >  FESpace_Type;
typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;
    
typedef boost::shared_ptr<Epetra_Comm> communicatorPtr_Type;

class FluidOperator
{
public:
    
    FluidOperator(boost::shared_ptr<Epetra_Comm>& comm);
    
    ~FluidOperator();
    
    void setup(const GetPot& dataFile);
    
    void setBC(const boost::shared_ptr<BCHandler> bc);
    
    void buildSystem(const GetPot& dataFile);
    
    void buildInterfaceMap();

    // getters

    FESpacePtr_Type feVelocity()
    {
    	return M_uFESpace;
    }

    meshPtr_Type mesh()
    {
    	return M_fullMeshPtrFluid;
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
    meshPtr_Type M_fullMeshPtrFluid;
    meshPtr_Type M_localMeshPtrFluid;
    boost::shared_ptr<OseenData> M_oseenData;
    boost::shared_ptr<meshPartitioner_Type> M_meshPartFluid;
    boost::shared_ptr<MeshData> M_meshDataFluid;
    FESpacePtr_Type M_uFESpace;
    FESpacePtr_Type M_pFESpace;
    std::string M_uOrder;
    std::string M_pOrder;
    boost::shared_ptr<BCHandler> M_BCh;
    boost::shared_ptr<OseenSolver<mesh_Type> > M_fluid;
};

FluidOperator::FluidOperator(boost::shared_ptr<Epetra_Comm>& comm):
M_comm (comm)
{
    M_verbose = M_comm->MyPID()==0;
}

FluidOperator::~FluidOperator()
{}
  
void FluidOperator::setup(const GetPot& dataFile)
{
    loadData(dataFile);
    loadMesh();
    partitionMesh();
    // buildInterfaceMap();
    createFESpaces();
}
    
void FluidOperator::loadData(const GetPot& dataFile)
{
    M_oseenData.reset(new OseenData());
    M_oseenData->setup ( dataFile );
    
    M_meshDataFluid.reset(new MeshData());
    M_meshDataFluid->setup(dataFile, "fluid/space_discretization");
    
    M_uOrder = dataFile("fluid/space_discretization/vel_order","P1");
    M_pOrder = dataFile("fluid/space_discretization/pres_order","P1");
}

void FluidOperator::loadMesh()
{
    M_fullMeshPtrFluid.reset(new mesh_Type());
    readMesh(*M_fullMeshPtrFluid, *M_meshDataFluid);
    
    if (M_verbose)
        std::cout << "Mesh source: file(" << M_meshDataFluid->meshDir() << M_meshDataFluid->meshFile() << ")" << std::endl;
    
    if (M_verbose){
        std::cout << "Mesh statistics: " << std::endl;
        std::cout << "Mesh size (Hmax) : " << MeshUtility::MeshStatistics::computeSize(*M_fullMeshPtrFluid).maxH << std::endl;
        std::cout << "Mesh size (Hmin) : " << MeshUtility::MeshStatistics::computeSize(*M_fullMeshPtrFluid).minH << std::endl;
    }

}

void FluidOperator::partitionMesh()
{
    M_meshPartFluid.reset( new meshPartitioner_Type(M_fullMeshPtrFluid, M_comm));
    M_localMeshPtrFluid.reset(new mesh_Type(*M_meshPartFluid->meshPartition() ) );
}
    
void FluidOperator::buildInterfaceMap()
{
  /*
    for ( UInt i = 0; i < M_localMeshPtrFluid->numVertices(); ++i )
        if ( isInside (M_localMeshPtrFluid->point (i).markerID(), M_flags) )
            if (CheckVector->blockMap().LID (M_localMeshPtrFluid->point (i).id() ) != -1)
            {
                GID_nodes.insert (M_localMeshPtrFluid->point(i).id() );
            }
  */
}

/*    
bool FluidOperator::isInside (ID pointMarker, flagContainer_Type flags)
{
    int check = 0;
    if(flags[0]==-1)
        return true;
    else
    {
        for (UInt i = 0; i < flags.size(); ++i)
            if (pointMarker == flags[i])
            {
                ++check;
            }
        return (check > 0) ? true : false;
    }
}
*/

void FluidOperator::createFESpaces()
{
    if (M_verbose)
        std::cout << std::endl << "[Creating the fluid FE spaces]" << std::endl;
    
    if (M_verbose)
        std::cout << "FE for the velocity: " << M_uOrder << std::endl
        << "FE for the pressure: " << M_pOrder << std::endl;
    
    if (M_verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;
    
    M_uFESpace.reset (new FESpace_Type (M_localMeshPtrFluid, M_uOrder, 3, M_comm) );
    
    if (M_verbose)
        std::cout << "ok.\n" << "Building the pressure FE space ... " << std::flush;
    
    M_pFESpace.reset (new FESpace_Type (M_localMeshPtrFluid, M_pOrder, 1, M_comm) );
    
    if (M_verbose)
        std::cout << "ok.\n";
}
    
    
void FluidOperator::setBC(const boost::shared_ptr<BCHandler> bc)
{
    M_BCh.reset(new BCHandler(*bc));
}
    
void FluidOperator::buildSystem( const GetPot& dataFile )
{
    M_fluid.reset( new OseenSolver<mesh_Type>(M_oseenData, *M_uFESpace, *M_pFESpace, M_comm) );
    M_fluid->setUp (dataFile);
    M_fluid->buildSystem();
}
    
    
} // end namespace LifeV

#endif
