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
#include <sys/stat.h>


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

// Includes
#include <lifev/fsi/solver/HarmonicExtensionSolver.hpp>

// Boundary conditions
#include <lifev/core/fem/BCVector.hpp>

// Steklov-Poincare operator
#include "SteklovPoincareOperator.hpp"

using namespace LifeV;

// Shortcuts
typedef boost::shared_ptr<map_Type> 				mapPtr_Type;
typedef HarmonicExtensionSolver<mesh_Type>      	meshMotion_Type;
typedef boost::shared_ptr<meshMotion_Type>      	meshMotionPtr_Type;
typedef ExporterHDF5<mesh_Type>                 	exporter_Type;
typedef boost::shared_ptr<exporter_Type>        	exporterPtr_Type;
typedef FESpace< RegionMesh<LinearTetra>, MapEpetra > FESpace_type;

using namespace LifeV;

Real initDisplacement (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
	switch (i)
	{
	case 0:
		return - 0.058549 * ( x - 0.5 );
		break;
	case 1:
		return  ( 0.256942 / 2 ) * y;
		break;
	case 2:
		return - 0.058549 * ( z + 0.5);
		break;
	default:
		ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
		return 0.;
		break;
	}
}

Real fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    bool verbose = Comm->MyPID() == 0;
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm() );
    bool verbose = true;
#endif

    if (verbose)
    {
    	try
    	{
    		if(argc!=9)
    			throw "Wrong code execution";
    	}
    	catch ( const char* Message )
    	{
    		std::cout   << "caught exception :  " << Message << "\n";
    		std::cout 	<< 	"Please use the code as it follows: " << std::endl;
    		std::cout	<< 	"mpirun -np X ./FSI_SteklovPoincare.exe  -df dataFluid -ds dataStructure -da dataAle -o Output" <<  std::endl;
    		std::cout	<< 	"./FSI_SteklovPoincare.exe  -df dataFluid -ds dataStructure -da dataAle -o Output" 				<<  std::endl;
    		return 0;
    	}
    }

    GetPot command_line (argc, argv);
    const bool check = command_line.search (2, "-c", "--check");

    const std::string data_file_name_fluid = command_line.follow ("dataFluid", 2, "-df", "--fileFluid");
    GetPot data_file_fluid (data_file_name_fluid);

    const std::string data_file_name_structure = command_line.follow ("dataStructure", 2, "-ds", "--fileStructure");
    GetPot data_file_structure (data_file_name_structure);

    const std::string data_file_name_ale = command_line.follow ("dataAle", 2, "-da", "--fileAle");
    GetPot data_file_ale (data_file_name_structure);

    std::string problem_folder = command_line.follow ("output", 2, "-o", "--output");

    // set the folder for the output
    if ( problem_folder.compare("./") )
    {
    	problem_folder += "/";
    	if ( Comm->MyPID() == 0 )
    		mkdir ( problem_folder.c_str(), 0777 );
    }

    // Instantiate the fluid and structure operators

    FluidOperator      fluid(Comm);
    StructureOperator  structure(Comm);

    // Call the setup of the operators
    fluid.setup(data_file_fluid);
    structure.setup(data_file_structure);
    meshMotionPtr_Type ale( new meshMotion_Type( *fluid.feVelocity(), Comm ) );
    ale->setUp(data_file_ale);

    //////////////////////////
    // Build interface maps //
    //////////////////////////

    int nFlags = 1;
    std::vector<int> flags (nFlags);
    flags[0] = 1;

    SteklovPoincareOperator FSI;
    FSI.buildTranferOperators ( fluid.mesh(), fluid.feVelocity()->mesh(),
    							structure.mesh(), structure.feDisplacement()->mesh(),
    							flags, data_file_fluid );

    ///////////////////////////
    //  Boundary Conditions  //
    ///////////////////////////

    BCHandler	BCh_ale;
    BCFunctionBase bcf (fZero);

    BCh_ale.addBC ("Edges", 20, EssentialVertices, Full, bcf, 3);
    BCh_ale.addBC ("Edges", 2,  EssentialVertices, Full, bcf, 3);
    BCh_ale.addBC ("Base",  3,  EssentialVertices, Full, bcf, 3);
    boost::shared_ptr<BCVector> BCVectorALE;
    boost::shared_ptr<BCVector> BCVectorStructureDisp;

    fluid.buildSystem(data_file_fluid);

    /////////////////////////////
    //  Initialize exporters   //
    /////////////////////////////

    exporterPtr_Type exporterFluid ( new  exporter_Type ( data_file_fluid, "Fluid") );
    exporterPtr_Type exporterStructure ( new  exporter_Type ( data_file_structure, "Structure") );
    UInt totalVelocityDofs (3*fluid.feVelocity()->dof().numTotalDof() );

    exporterFluid->setMeshProcId (fluid.feVelocity()->mesh(), Comm->MyPID() );
    exporterStructure->setMeshProcId (structure.feDisplacement()->mesh(), Comm->MyPID() );

    vectorPtr_Type velAndPressure( new vector_Type( fluid.solver()->getMap(), exporterFluid->mapType() ) );
    vectorPtr_Type fluidDisp( new vector_Type( fluid.feVelocity()->map(), exporterFluid->mapType() ) );
    vectorPtr_Type weakStressFluid( new vector_Type( fluid.feVelocity()->map(), exporterFluid->mapType() ) );
    vectorPtr_Type structureDisp( new vector_Type( structure.feDisplacement()->map(), exporterStructure->mapType() ) );

    *velAndPressure  *= 0;
    *fluidDisp       *= 0;
    *weakStressFluid *= 0;
    *structureDisp   *= 0;

    exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f-velocity",fluid.feVelocity(), velAndPressure, UInt ( 0 ) );
    exporterFluid->addVariable ( ExporterData<mesh_Type>::ScalarField, "f-pressure",fluid.fePressure(), velAndPressure, UInt ( totalVelocityDofs ) );
    exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f-weakStress",fluid.feVelocity(), weakStressFluid, UInt ( 0 ) );
    exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f-displacement", fluid.feVelocity(), fluidDisp, UInt (0) );
    exporterStructure->addVariable ( ExporterData<mesh_Type>::VectorField, "s-displacement", structure.feDisplacement(), structureDisp, UInt (0) );


    ///////////////////////////
    // Time Advance Objects  //
    ///////////////////////////

    // Create the objects
    timeAdvancePtr_Type 	fluidTA;
    timeAdvancePtr_Type 	structureTA;
    timeAdvancePtr_Type 	aleTA;

    fluidTA.reset( TimeAdvanceFactory::instance().createObject( data_file_fluid("fluid/time_discretization/method", "BDF") ) );
    structureTA.reset(TimeAdvanceFactory::instance().createObject( data_file_structure("solid/time_discretization/method", "BDF") ) );
    aleTA.reset(TimeAdvanceFactory::instance().createObject( data_file_ale("mesh_motion/time_discretization/method", "BDF") ) );

    fluidTA->setup ( fluid.data()->dataTimeAdvance()->orderBDF(), 1);
    fluidTA->setTimeStep ( fluid.data()->dataTime()->timeStep() );

    std::vector<Real> parameters (2);
    parameters[0]  = data_file_structure ("solid/time_discretization/theta", 0.25);
    parameters[1]  = data_file_structure ("solid/time_discretization/zeta", 0.5);
    UInt order = data_file_structure ("solid/time_discretization/BDF_order", 1);

    structureTA->setup ( order , 2);
    structureTA->setTimeStep ( structure.data()->dataTime()->timeStep() );
    structure.buildSystem(data_file_structure, structureTA);

    aleTA->setup ( fluid.data()->dataTimeAdvance()->orderBDF(), 1 );
    aleTA->setTimeStep ( fluid.data()->dataTime()->timeStep() );

    // Initialize the time advance objects

    std::vector<vectorPtr_Type> velInit;
    std::vector<vectorPtr_Type> dispFluidInit;
    std::vector<vectorPtr_Type> dispStructureInit;

    for(UInt i = 0; i < fluidTA->size(); ++i)
    	velInit.push_back(fluid.solver()->solution());

    fluidTA->setInitialCondition ( velInit);

    for(UInt i = 0; i < aleTA->size(); ++i)
    	dispFluidInit.push_back ( fluidDisp );

    aleTA->setInitialCondition ( dispFluidInit ) ;

    vectorPtr_Type initialDisplacement (new vector_Type (structure.solver()->displacement(), Unique) );
    structure.feDisplacement()->interpolate ( static_cast<FESpace_type::function_Type> ( initDisplacement ), *initialDisplacement, 0.0 );

    for(UInt i = 0; i < structureTA->size(); ++i)
    	dispStructureInit.push_back ( initialDisplacement );

    structureTA->setInitialCondition ( dispStructureInit ) ;

    structureTA->updateRHSContribution ( structure.data()->dataTime()->timeStep() );

    //////////////////////////
    //  Initialize systems  //
    //////////////////////////

    vectorPtr_Type a;
    a.reset( new vector_Type(structure.solver()->displacement(), Unique));
    structure.solver()->initialize ( initialDisplacement );

    ///////////////////////////////////////////
    // Initialize the interface displacement //
    ///////////////////////////////////////////

    // Set the vectors that will be used by intermesh operators
    FSI.setUpData(fluidDisp, structureDisp);

    // Diplacement of the interface defined on the fluid domain
    vectorPtr_Type lambdaF;
    lambdaF.reset ( new vector_Type ( *FSI.fluidInterfaceMap() ) );
    *lambdaF += 0.001;

    vector_Type fluid_vector ( fluid.feVelocity()->map() );
    fluid_vector.subset(*lambdaF, *FSI.fluidInterfaceMap(), 0, 0);

    ///////////////////////////
    // Core of the algorithm //
    ///////////////////////////

    // Initial guess, lambda equal to zero

    /////////////////////////////////////////////
    // EvalResidual: S_f(lambda) + S_s(lambda) //
    /////////////////////////////////////////////

    // 1) ALE problem

    if(Comm->MyPID()==0)
    {
    	std::cout << "\n\n";
    	std::cout << "[***************************************]" << std::endl;
    	std::cout << "       [Solving the ALE problem] " << std::endl;
    	std::cout << "[***************************************]" << std::endl;
    }

    BCVectorALE.reset(new BCVector(fluid_vector, fluid.feVelocity()->dof().numTotalDof(),(UInt)0));
    BCh_ale.addBC ("Interface", 1,  Essential, Full, *BCVectorALE, 3);

    ale->iterate (BCh_ale);
    *fluidDisp = ale->disp();

    aleTA->updateRHSFirstDerivative ( fluid.data()->dataTime()->timeStep() );
    aleTA->shiftRight( *fluidDisp );
    aleTA->extrapolation( *fluidDisp );

    // 2) Move the fluid mesh

    if(Comm->MyPID()==0)
    {
    	std::cout << "\n\n";
    	std::cout << "[***************************************]" << std::endl;
    	std::cout << "       [Moving the fluid mesh] " << std::endl;
    	std::cout << "[***************************************]" << std::endl;
    }

    vector_Type meshDisp (*fluidDisp, Repeated);
    fluid.feVelocity()->mesh()->meshTransformer().moveMesh(meshDisp, fluid.feVelocity()->dof().numTotalDof());
    fluid.solver()->setRecomputeMatrix(true);

    // WARNING: HERE WE COMPUTE THE FLUID MESH VELOCITY
    // vel is a vector of velocity of the fluid domain with the map of fluid displacement that in principle can be different from
    // the one of the fluid velocity. For example if the ale order is 1 and the order of the fluid velocity is 2 we need to interpolate
    // Hence pay attention if vel and beta have different maps.
    vector_Type vel ( aleTA->firstDerivative() );
    vector_Type beta ( fluid.feVelocity()->map() );
    fluidTA->extrapolation( beta );
    beta -= vel;

    // 3) Solve the fluid problem

    if(Comm->MyPID()==0)
    {
    	std::cout << "\n\n";
    	std::cout << "[***************************************]" << std::endl;
    	std::cout << "       [Solve the fluid problem] " << std::endl;
    	std::cout << "[***************************************]" << std::endl;
    }

    // Add the possibility of conservative formulation
    vector_Type rhs ( velAndPressure->map() );
    rhs *= 0;

    if (fluid.solver()->recomputeMatrix())
    {
    	double alpha = fluidTA->coefficientFirstDerivative (0) / fluid.data()->dataTime()->timeStep();
    	fluid.solver()->updateSystem ( alpha, beta, rhs );
    }
    else
    {
    	fluid.solver()->updateRightHandSide ( rhs );
    }

    rhs = fluid.solver()->matrixMass() * fluidTA->rhsContributionFirstDerivative();
    fluid.solver()->updateRightHandSide ( rhs );

    fluid.iterate();

    *weakStressFluid = fluid.solver()->residual();

    *velAndPressure = *fluid.solver()->solution();

    // 4) Solve the structure problem

    if(Comm->MyPID()==0)
    {
    	std::cout << "\n\n";
    	std::cout << "[***************************************]" << std::endl;
    	std::cout << "       [Solve the structure problem] " << std::endl;
    	std::cout << "[***************************************]" << std::endl;
    }

    vectorPtr_Type lambdaS;
    lambdaS.reset ( new vector_Type ( *FSI.structureInterfaceMap() ) );
    FSI.transferGammaFluidOnGammaStructure(lambdaF, lambdaS);

    vector_Type structure_vector ( structureDisp->map() );
    structure_vector.subset(*lambdaS, *FSI.structureInterfaceMap(), 0, 0);

    lambdaS->spy("lambdaS");

    BCVectorStructureDisp.reset(new BCVector(structure_vector, structure.feDisplacement()->dof().numTotalDof(),(UInt)0));
    structure.bc()->addBC ("Interface", 1, EssentialVertices, Full, *BCVectorStructureDisp, 3);

    structure.iterate(structureTA);
    *structureDisp = structure.solver()->displacement();

    ///////////////////////////////
    // PostProcessing the result //
    ///////////////////////////////

    exporterFluid->postProcess ( structure.data()->dataTime()->timeStep() );
    exporterStructure->postProcess ( fluid.data()->dataTime()->timeStep() );

    //////////////////////////
    //  Closing simulation  //
    //////////////////////////

    exporterFluid->closeFile();
    exporterStructure->closeFile();

    if(Comm->MyPID()==0)
    {
    	std::cout << "\n\n";
    	std::cout << "[***************************************]" << std::endl;
    	std::cout << "       [End of the simulation] " << std::endl;
    	std::cout << "[***************************************]" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return 0;

}

