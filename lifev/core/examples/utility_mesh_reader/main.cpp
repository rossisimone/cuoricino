/*
 * test_genAlpha.cpp
 *
 *  Created on: Jul 27, 2010
 *      Author: uvilla
 */

#include <Epetra_ConfigDefs.h>
#include <Epetra_Comm.h>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>

using namespace LifeV;

// Do not edit
int main(int argc, char **argv)
{
	using namespace LifeV;
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	std::cout<< "MPI Initialization\n";
#endif

#ifdef EPETRA_MPI
	boost::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	boost::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm());
#endif
	bool verbose = comm->MyPID()==0;

	std::string dataFileName;
	GetPot command_line(argc, argv);
	dataFileName = command_line.follow("data", 2, "-f", "--file");
	GetPot dataFile(dataFileName);

    //Create the mesh data and read and partitioned the mesh
    boost::shared_ptr< RegionMesh <LinearTetra> > meshFullPtr ( new RegionMesh <LinearTetra> ( comm ) );
    std::string meshName = dataFile ("mesh/mesh_file", "cube4x4.mesh");
    std::string meshPath =  dataFile ("mesh/mesh_dir", "./");
    std::string meshOrder=  "P1";
    bool isPartitioned = false;


    //MeshUtility::fillWithMesh( meshFullPtr, isPartitioned, meshName, meshPath, "P1" );
    MeshUtility::fillWithMesh( meshFullPtr,isPartitioned, meshName, meshPath, meshOrder);

    boost::shared_ptr< RegionMesh <LinearTetra> > meshStructPtr ( new RegionMesh <LinearTetra> ( comm ) );
    MeshUtility::fillWithStructuredMesh( meshStructPtr,
    							 1,
    							 5,
    							 5,
    							 5,
    							 true,
    							 1.0,
    							 1.0,
    							 1.0,
    							 0.0,
    							 0.0,
    							 0.0 );
	comm.reset();

#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout<< "MPI Finalization \n";
#endif

	return 0;
}
