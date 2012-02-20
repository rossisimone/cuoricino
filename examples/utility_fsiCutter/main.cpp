/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author: Radu Popescu <radu.popescu@epfl.ch>
  Copyright (C) 2010 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file main.cpp
   \author Radu Popescu <radu.popescu@epfl.ch>
   \date 2010-08-24
 */

/*
  Test of the offline partitioning for FSI simulations
*/

#include <iostream>
#include <string>

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

#include <boost/shared_ptr.hpp>

#include <life/lifecore/LifeV.hpp>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitionerOfflineFSI.hpp>
#include <life/lifefilters/ExporterHDF5Mesh3D.hpp>

using namespace LifeV;

typedef RegionMesh<LinearTetra> Mesh;

int main( int argc, char** argv )
{
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

    int fluidInterfaceFlag = dataFile("interface/fluid_flag", 1);
    int solidInterfaceFlag = dataFile("interface/structure_flag", fluidInterfaceFlag);
    Real interfaceTolerance = dataFile("interface/tolerance", 0.0);

    int fluidInterfaceVertexFlag;
    int solidInterfaceVertexFlag;
    int vertexFlag;
    vertexFlag = dataFile("interface/edgeFlag",      -1);
    vertexFlag = dataFile("interface/fluid_vertex_flag", vertexFlag);

    fluidInterfaceVertexFlag = vertexFlag;

    solidInterfaceVertexFlag = dataFile("interface/structure_vertex_flag", -1);

    std::string fluidOrder(dataFile("fluid/space_discretization/vel_order", "P1"));
    std::string solidOrder(dataFile("solid/space_discretization/order", "P1"));

    boost::shared_ptr<MeshData> fluidMeshData(new MeshData);
    fluidMeshData->setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<Mesh> uncutFluidMesh(new Mesh);
    readMesh(*uncutFluidMesh, *fluidMeshData);

    boost::shared_ptr<MeshData> solidMeshData(new MeshData);
    solidMeshData->setup(dataFile, "solid/space_discretization");

    boost::shared_ptr<Mesh> uncutSolidMesh(new Mesh);
    readMesh(*uncutSolidMesh, *solidMeshData);

    boost::shared_ptr<MeshPartitionerOfflineFSI<Mesh> >
    cutter(new MeshPartitionerOfflineFSI<Mesh>);
    cutter->setup(uncutFluidMesh, uncutSolidMesh, 4, 4, fluidOrder, solidOrder,
                  fluidInterfaceFlag, solidInterfaceFlag, interfaceTolerance,
                  fluidInterfaceVertexFlag, solidInterfaceVertexFlag, comm);
    cutter->showMe();

    cutter->execute();

    ExporterHDF5Mesh3D<Mesh> fluidOutput(dataFile, uncutFluidMesh, "FSIFluidPartitions",
                                         comm->MyPID());
    fluidOutput.addPartitionGraph(cutter->fluidGraph(), comm);
    fluidOutput.addMeshPartitionAll(cutter->fluidPartitions(), comm);
    fluidOutput.addDOFInterface(cutter->dofStructureToHarmonicExtension(),
                                "dofStructureToHarmonicExtension",
                                cutter->fluidInterfaceFlag(),
                                cutter->solidInterfaceFlag(),
                                comm);
    fluidOutput.postProcess(0);
    fluidOutput.closeFile();

    ExporterHDF5Mesh3D<Mesh> solidOutput(dataFile, uncutSolidMesh, "FSISolidPartitions",
                                         comm->MyPID());
    solidOutput.addPartitionGraph(cutter->solidGraph(), comm);
    solidOutput.addMeshPartitionAll(cutter->solidPartitions(), comm);
    solidOutput.addDOFInterface(cutter->dofStructureToHarmonicExtension(),
                                "dofStructureToHarmonicExtension",
                                cutter->fluidInterfaceFlag(),
                                cutter->solidInterfaceFlag(),
                                comm);
    solidOutput.postProcess(0);
    solidOutput.closeFile();


#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return(EXIT_SUCCESS);
}

