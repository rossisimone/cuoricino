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
#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>

// Includes
// #include "SteklovPoincareOperator.hpp"
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

typedef RBFInterpolation<mesh_Type>             interpolation_Type;
typedef boost::shared_ptr<interpolation_Type>   interpolationPtr_Type;

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

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

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

    int nFlags = 2;
    std::vector<int> flags (nFlags);
    flags[0] = 1;
    flags[1] = 20;

    interpolationPtr_Type RBFinterpolant;
    RBFinterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject (data_file("interpolation/interpolation_Type","none")));
    RBFinterpolant->setup(structure.mesh(), structure.feDisplacement()->mesh(),
    					  fluid.mesh(), fluid.feVelocity()->mesh(),
    					  flags);

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
    // InterMesh operators   //
    ///////////////////////////

    RBFinterpolant->setupRBFData (structureDisp, fluidDisp);
    RBFinterpolant->buildOperators ();

    // Interface map of the fluid

    // Interface map of the structure

    ///////////////////////////
    // Core of the algorithm //
    ///////////////////////////

    *structureDisp += 10;

    RBFinterpolant->updateRhs (structureDisp);
    RBFinterpolant->interpolate ();
    RBFinterpolant->solution (fluidDisp);

    exporterFluid->postProcess(0.0);
    exporterStructure->postProcess(0.0);


    *structureDisp = -1500;

    RBFinterpolant->updateRhs (structureDisp);
    RBFinterpolant->interpolate ();
    RBFinterpolant->solution (fluidDisp);

    exporterFluid->postProcess(1.0);
    exporterStructure->postProcess(1.0);

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

