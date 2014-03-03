// DAVIDE FORTI
// TO DO: - DOCUMENTATION

#undef HAVE_HDF5

#include <cassert>
#include <cstdlib>

#include <boost/timer.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


// LifeV includes
#include <lifev/core/LifeV.hpp>

// Include for the data
#include <lifev/fsi/solver/FSIData.hpp>

// Exporters
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

// Include Fluid operator (for the moment it can stay here)
#include "FluidOperator.hpp"

// Include Structure operator (for the moment it can stay here)
#include "StructureOperator.hpp"

// Includes for the interface
#include <lifev/core/fem/DOFInterface3Dto3D.hpp>

// Includes
#include "SteklovPoincareOperator.hpp"
#include <lifev/fsi/solver/HarmonicExtensionSolver.hpp>

// Boundary conditions
#include "boundaryConditions.hpp"

using namespace LifeV;

// Shortcuts
typedef VectorEpetra 							vector_Type;
typedef boost::shared_ptr<vector_Type> 			vectorPtr_Type;

typedef MapEpetra 					    		map_Type;
typedef boost::shared_ptr<map_Type> 			mapPtr_Type;

typedef boundaryConditions				     	bcFSIProblem_Type;
typedef boost::shared_ptr<bcFSIProblem_Type> 	bcFSIProblemPtr_Type;

typedef HarmonicExtensionSolver<mesh_Type>      meshMotion_Type;
typedef boost::shared_ptr<meshMotion_Type>      meshMotionPtr_Type;

typedef ExporterHDF5<mesh_Type>                 exporter_Type;
typedef boost::shared_ptr<exporter_Type>        exporterPtr_Type;

// Begin of the program
int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    bool verbose = Comm->MyPID() == 0;

    if(verbose)
        cout << "[Using parallel version ]" << endl;
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm() );
    std::cout << "[Using serial version ]" << std::endl;
#endif

    GetPot command_line (argc, argv);
    const bool check = command_line.search (2, "-c", "--check");

    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot data_file (data_file_name);

    // Data of the simulation
    boost::shared_ptr<FSIData> data( new FSIData);
    data->setup(data_file);

    // Instantiate the fluid and structure operators
    FluidOperator      fluid(Comm);
    StructureOperator  structure(Comm);

    // Call the setup of the operators
    fluid.setup(data_file);
    structure.setup(data_file);

    // ALE object (harmonic extension)
    meshMotionPtr_Type ale( new meshMotion_Type( *fluid.feVelocity(), Comm ) );
    ale->setUp(data_file);

    //////////////////////////
    // Build interface maps //
    //////////////////////////

    // Interface map of the fluid
    boost::shared_ptr<DOFInterface3Dto3D> DOFfluid;
    DOFfluid.reset(new DOFInterface3Dto3D(fluid.feVelocity()->refFE(), fluid.feVelocity()->dof()));

    DOFfluid->update( *fluid.feVelocity()->mesh(), data->fluidInterfaceFlag(),
    		          *fluid.feVelocity()->mesh(), data->fluidInterfaceFlag(),
        			   data->interfaceTolerance(), data->fluidInterfaceVertexFlag() );

    // Interface map of the structure
    boost::shared_ptr<DOFInterface3Dto3D> DOFstructure;
    DOFstructure.reset(new DOFInterface3Dto3D(structure.feDisplacement()->refFE(), structure.feDisplacement()->dof()));

    DOFstructure->update( *structure.feDisplacement()->mesh(), data->structureInterfaceFlag(),
    		              *structure.feDisplacement()->mesh(), data->structureInterfaceFlag(),
    					   data->interfaceTolerance() );

    // Important object that stores the dofwise connection at the interface
    boost::shared_ptr<DOFInterface3Dto3D> DOFstructureToFluid;
    DOFstructureToFluid.reset( new DOFInterface3Dto3D);

    DOFstructureToFluid->setup( fluid.feVelocitySerial()->refFE(), fluid.feVelocitySerial()->dof(),
    							structure.feDisplacementSerial()->refFE(), structure.feDisplacementSerial()->dof() );

    DOFstructureToFluid->update( *fluid.feVelocitySerial()->mesh(), data->fluidInterfaceFlag(),
    							 *structure.feDisplacementSerial()->mesh(), data->structureInterfaceFlag(),
    							 data->interfaceTolerance() );

    // Create Steklov Poincare operator
    SteklovPoincareOperator SPoperator(fluid.feVelocity(), structure.feDisplacement(),
    							       DOFstructureToFluid,	DOFstructure);

    // Interface map of the fluid
    mapPtr_Type fluidInterfaceMap;
    fluidInterfaceMap.reset ( new MapEpetra ( *SPoperator.createInterfaceMaps(DOFfluid->localDofMap(), fluid.feVelocity()->dof(), Comm, 0 ) ) );

    // Interface map of the structure
    mapPtr_Type structureInterfaceMap;
    structureInterfaceMap.reset ( new MapEpetra ( *SPoperator.createInterfaceMaps(DOFstructure->localDofMap(), structure.feDisplacement()->dof(), Comm, 1 ) ) );

    // Initialize the tranfer operator across the interface
    SPoperator.buildTranferOperators( fluidInterfaceMap, structureInterfaceMap,
    								  DOFstructureToFluid->localDofMap());

    ///////////////////////////
    //  Boundary Conditions  //
    ///////////////////////////

    bcFSIProblemPtr_Type bcFSI (new bcFSIProblem_Type);
    bcFSI->setup();
    fluid.setBC(bcFSI->fluidBC());
    structure.setBC(bcFSI->structureBC());

    /////////////////////////////
    //  Initialize exporters   //
    /////////////////////////////

    exporterPtr_Type exporterFluid ( new  exporter_Type ( data_file, "Fluid") );
    exporterPtr_Type exporterStructure ( new  exporter_Type ( data_file, "Structure") );

    UInt totalVelocityDofs (3*fluid.feVelocity()->dof().numTotalDof() );

    exporterFluid->setMeshProcId (fluid.feVelocity()->mesh(), Comm->MyPID() );
    exporterStructure->setMeshProcId (structure.feDisplacement()->mesh(), Comm->MyPID() );

    //////////////////////////
    //  Initialize systems  //
    //////////////////////////

    fluid.buildSystem(data_file);
    structure.buildSystem(data_file);
    ale->computeMatrix();

    //////////////////////////
    //  Initialize systems  //
    //////////////////////////

    vectorPtr_Type velAndPressure( new vector_Type( fluid.solver()->getMap(), exporterFluid->mapType() ) );
    vectorPtr_Type fluidDisp( new vector_Type( fluid.feVelocity()->map(), exporterFluid->mapType() ) );
    vectorPtr_Type structureDisp( new vector_Type( structure.feDisplacement()->map(), exporterFluid->mapType() ) );

    *velAndPressure *= 0;
    *fluidDisp      *= 0;
    *structureDisp  *= 0;

    exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f-velocity",fluid.feVelocity(), velAndPressure, UInt ( 0 ) );
    exporterFluid->addVariable ( ExporterData<mesh_Type>::ScalarField, "f-pressure",fluid.fePressure(), velAndPressure, UInt ( totalVelocityDofs ) );
    exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f-displacement", fluid.feVelocity(), fluidDisp, UInt (0) );

    exporterStructure->addVariable ( ExporterData<mesh_Type>::VectorField, "s-displacement", structure.feDisplacement(), structureDisp, UInt (0) );

    ///////////////////////////
    // Core of the algorithm //
    ///////////////////////////

    *structureDisp  += 100;
    //vectorPtr_Type structureDispGamma( new vector_Type( *structureInterfaceMap, exporterStructure->mapType() ) );
    //structureDispGamma->subset(*structureDisp, *structureInterfaceMap, 0, 0);

    //vectorPtr_Type fluidDispGamma( new vector_Type( *fluidInterfaceMap, exporterFluid->mapType() ) );
    SPoperator.transferStructureOnFluid(structureDisp, fluidDisp);

    exporterFluid->postProcess(0.0);
    exporterStructure->postProcess(0.0);

    //////////////////////////
    //  Closing simulation  //
    //////////////////////////


    exporterFluid->closeFile();
    exporterStructure->closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return 0;

}

