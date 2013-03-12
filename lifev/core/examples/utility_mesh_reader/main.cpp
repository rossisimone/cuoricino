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

	std::string dataFileName;
	GetPot command_line(argc, argv);
	dataFileName = command_line.follow("data", 2, "-f", "--file");
	GetPot dataFile(dataFileName);

	if( comm -> MyPID() == 0 ) cout << "\n\nReading and partitioning the cube mesh without saving the global mesh: ... ";
    //Create the mesh data and read and partitioned the mesh
    boost::shared_ptr< RegionMesh <LinearTetra> > meshPtr ( new RegionMesh <LinearTetra> ( comm ) );
    std::string meshName = dataFile ("mesh/mesh_file", "cube4x4.mesh");
    std::string meshPath =  dataFile ("mesh/mesh_dir", "./");
    bool isPartitioned = false;
    MeshUtility::fillWithMesh( meshPtr, isPartitioned, meshName, meshPath );
	if( comm -> MyPID() == 0 ) cout << "... DONE! ";


    //create the mesh data and read and save both the global mesh and the partitioned mesh
	if( comm -> MyPID() == 0 ) cout << "\n\nReading and partitioning the cube mesh saving the global mesh: ... ";
	boost::shared_ptr< RegionMesh <LinearTetra> > meshFullPtr ( new RegionMesh <LinearTetra> ( comm ) );
    boost::shared_ptr< RegionMesh <LinearTetra> > meshLocalPtr ( new RegionMesh <LinearTetra> ( comm ) );
    MeshUtility::fillWithFullMesh( meshFullPtr, meshLocalPtr, meshName, meshPath );
	if( comm -> MyPID() == 0 ) cout << "... DONE! ";

    //create a 3D structured mesh
	if( comm -> MyPID() == 0 ) cout << "\n\nCreating a structured mesh: ... ";
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
    							 0.0 );\


	if( comm -> MyPID() == 0 ) cout << "... DONE!\n\n ";

	comm.reset();

#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout<< "MPI Finalization \n";
#endif

	return 0;
}
