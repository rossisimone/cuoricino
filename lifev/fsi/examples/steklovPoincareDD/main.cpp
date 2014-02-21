// DAVIDE FORTI - TO DO DOCUMENTATION

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

using namespace LifeV;

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
    FluidOperator     fluid(Comm);
    StructureOperator structure(Comm);

    // Call the setup of the operators
    fluid.setup(data_file);
    structure.setup(data_file);

    // Build interface map
    boost::shared_ptr<DOFInterface3Dto3D> DOFfluidToSolid;
    DOFfluidToSolid.reset(new DOFInterface3Dto3D());

    DOFfluidToSolid->setup(	structure.feDisplacement()->refFE(), structure.feDisplacement()->dof(),
    						fluid.feVelocity()->refFE(), fluid.feVelocity()->dof());

    DOFfluidToSolid->update( *structure.mesh(), data->structureInterfaceFlag(), *fluid.mesh(), data->fluidInterfaceFlag(),
        					 data->interfaceTolerance(), data->fluidInterfaceVertexFlag() );

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return 0;

}

