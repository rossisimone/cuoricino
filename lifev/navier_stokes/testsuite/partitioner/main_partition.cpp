//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Test of the serial mesh partitioning - part 1

    @author Radu Popescu <radu.popescu@epfl.ch>
    @date 02-07-2010

    Partition a mesh using a single (MPI) process and save mesh partitions
    to a HDF5 file.
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/filter/ExporterHDF5Mesh3D.hpp>

#include <iostream>
#include <string>


using namespace LifeV;

int main( int argc, char** argv )
{

#ifdef HAVE_HDF5

    boost::shared_ptr<Epetra_Comm> comm;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    comm.reset(new Epetra_MpiComm(MPI_COMM_WORLD) );
    std::cout << "% using MPI version" << std::endl;
#else
    comm.reset( new Epetra_SerialComm() );
    std::cout << "% using serial version" << std::end;
#endif

    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    OseenData oseenData;
    oseenData.setup(dataFile);

    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr(new RegionMesh<LinearTetra>);
    readMesh(*fullMeshPtr, meshData);

    MeshPartitioner<RegionMesh<LinearTetra> > meshPart;
    meshPart.setup(4, comm);

    meshPart.attachUnpartitionedMesh(fullMeshPtr);
    meshPart.doPartitionGraph();
    meshPart.doPartitionMesh();

    // Release the original mesh from the MeshPartitioner object and delete the RegionMesh object
    meshPart.releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    ExporterHDF5Mesh3D<RegionMesh<LinearTetra> > HDF5Output(dataFile, meshPart.meshPartition(), "cylinderPart",
                                                              comm->MyPID());
    HDF5Output.addPartitionGraph(meshPart.elementDomains(), comm);
    HDF5Output.addMeshPartitionAll(meshPart.meshPartitions(), comm);
    HDF5Output.postProcess(0);
    HDF5Output.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

#endif // HAVE_HDF5

    return(EXIT_SUCCESS);
}

