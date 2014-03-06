#ifndef FLUIDOPERATOR_H
#define FLUIDOPERATOR_H 1

#include <lifev/core/filter/GetPot.hpp>

//fespace
#include <lifev/core/fem/FESpace.hpp>

// fluid data
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/navier_stokes/solver/OseenSolverShapeDerivative.hpp>

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

	typedef BCHandler                                          		bc_Type;
	typedef boost::shared_ptr< bc_Type >                       		bcPtr_Type;
	typedef OseenSolverShapeDerivative< RegionMesh<LinearTetra> >   physicalSolverFluid_Type;
	typedef BCInterface3D< bc_Type, physicalSolverFluid_Type >		bcInterface_Type;
	typedef boost::shared_ptr< bcInterface_Type >           		bcInterfacePtr_Type;
    
    FluidOperator(boost::shared_ptr<Epetra_Comm>& comm);
    
    ~FluidOperator();
    
    void setup(const GetPot& dataFile);
    
    void buildSystem(const GetPot& dataFile);

    void iterate();

    // getters

    boost::shared_ptr<OseenData> data()
	{
    	return M_oseenData;
	}

    boost::shared_ptr<OseenSolverShapeDerivative<mesh_Type> > solver()
    {
    	return M_fluid;
    }

    FESpacePtr_Type feVelocitySerial()
    {
    	return M_uFESpaceSerial;
    }

    FESpacePtr_Type feVelocity()
    {
    	return M_uFESpace;
    }

    FESpacePtr_Type fePressure()
    {
    	return M_pFESpace;
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
    FESpacePtr_Type M_uFESpaceSerial;
    FESpacePtr_Type M_pFESpace;
    std::string M_uOrder;
    std::string M_pOrder;
    boost::shared_ptr<BCHandler> M_BCh;
    boost::shared_ptr<OseenSolverShapeDerivative<mesh_Type> > M_fluid;
    bcInterfacePtr_Type M_fluidBCPtr;
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
    createFESpaces();
    M_fluidBCPtr.reset ( new bcInterface_Type() );
    M_fluidBCPtr->createHandler();
    M_fluidBCPtr->fillHandler ( "dataFluid", "fluid" );
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
    M_uFESpaceSerial.reset (new FESpace_Type (M_fullMeshPtrFluid, M_uOrder, 3, M_comm) );
    
    if (M_verbose)
        std::cout << "ok.\n" << "Building the pressure FE space ... " << std::flush;
    
    M_pFESpace.reset (new FESpace_Type (M_localMeshPtrFluid, M_pOrder, 1, M_comm) );
    
    if (M_verbose)
        std::cout << "ok.\n";
}
    
void FluidOperator::buildSystem( const GetPot& dataFile )
{
    M_fluid.reset( new OseenSolverShapeDerivative<mesh_Type>(M_oseenData, *M_uFESpace, *M_pFESpace, M_comm) );
    M_fluid->setUp (dataFile);
    M_fluid->buildSystem();
}
    
void FluidOperator::iterate()
{
	M_fluid->iterate(*M_fluidBCPtr->handler());
}
    
} // end namespace LifeV

#endif
