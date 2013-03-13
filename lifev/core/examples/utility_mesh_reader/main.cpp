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
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>

using namespace LifeV;

// Do not edit
int main (int argc, char** argv)
{
    using namespace LifeV;
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::cout << "MPI Initialization\n";
#endif

#ifdef EPETRA_MPI
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_SerialComm() );
#endif

    std::string dataFileName;
    GetPot command_line (argc, argv);
    dataFileName = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

    if ( comm -> MyPID() == 0 )
    {
        cout << "\n\nReading and partitioning the cube mesh without saving the global mesh: ... \n";
    }
    //Create the mesh data and read and partitioned the mesh
    boost::shared_ptr< RegionMesh <LinearTetra> > meshPtr ( new RegionMesh <LinearTetra> ( comm ) );
    std::string meshName = dataFile ("mesh/mesh_file", "cube4x4.mesh");
    std::string meshPath =  dataFile ("mesh/mesh_dir", "./");
    std::string meshOrder =  dataFile ("mesh/mesh_order", "P1");
    bool isPartitioned = false;
    MeshUtility::fillWithMesh ( meshPtr, isPartitioned, meshName, meshPath, meshOrder );
    if ( comm -> MyPID() == 0 )
    {
        cout << "... DONE! ";
    }


    //create the mesh data and read and save both the global mesh and the partitioned mesh
    if ( comm -> MyPID() == 0 )
    {
        cout << "\n\nReading and partitioning the cube mesh saving the global mesh: ... \n";
    }
    boost::shared_ptr< RegionMesh <LinearTetra> > meshFullPtr ( new RegionMesh <LinearTetra> ( comm ) );
    boost::shared_ptr< RegionMesh <LinearTetra> > meshLocalPtr ( new RegionMesh <LinearTetra> ( comm ) );
    MeshUtility::fillWithMesh ( meshLocalPtr, meshFullPtr, isPartitioned, meshName, meshPath, "P1" );
    if ( comm -> MyPID() == 0 )
    {
        cout << "... DONE! ";
    }

    //create a 3D structured mesh
    if ( comm -> MyPID() == 0 )
    {
        cout << "\n\nCreating a structured mesh without saving the full mesh: ... \n";
    }
    boost::shared_ptr< RegionMesh <LinearTetra> > meshStructPtr ( new RegionMesh <LinearTetra> ( comm ) );
    MeshUtility::fillWithStructuredMesh ( meshStructPtr,
                                          1,
                                          std::vector<UInt> (3, 5),
                                          true,
                                          std::vector<UInt> (3, 1),
                                          std::vector<UInt> (3, 0) );


    if ( comm -> MyPID() == 0 )
    {
        cout << "... DONE!\n\n ";
    }


    //create a 3D structured mesh saving the full mesh
    if ( comm -> MyPID() == 0 )
    {
        cout << "\n\nCreating a structured mesh saving the full mesh: ... \n";
    }
    boost::shared_ptr< RegionMesh <LinearTetra> > meshLocalStructPtr ( new RegionMesh <LinearTetra> ( comm ) );
    boost::shared_ptr< RegionMesh <LinearTetra> > meshFullStructPtr ( new RegionMesh <LinearTetra> ( comm ) );
    MeshUtility::fillWithStructuredMesh (meshStructPtr,
                                         meshLocalStructPtr,
                                         1,
                                         std::vector<UInt> (3, 5),
                                         true,
                                         std::vector<UInt> (3, 1),
                                         std::vector<UInt> (3, 0) );


    if ( comm -> MyPID() == 0 )
    {
        cout << "... DONE!\n\n ";
    }


    comm.reset();

#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout << "MPI Finalization \n";
#endif

    return 0;
}
