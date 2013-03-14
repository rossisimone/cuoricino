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
    @brief 0D test with the Negroni Lascano model of 1996.

    @date 01−2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
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



#include <fstream>
#include <string>

#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/heart/solver/HeartETAMonodomainSolver.hpp>
#include <lifev/heart/solver/HeartIonicSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace LifeV;


Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    //********************************************//
    // Starts the chronometer.                    //
    //********************************************//
    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList list1 = * ( Teuchos::getParametersFromXmlFile ( "ParamList1.xml" ) );
    Teuchos::ParameterList list2 = * ( Teuchos::getParametersFromXmlFile ( "ParamList2.xml" ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    std::string meshName1 = list1.get ("mesh_name", "lid16.mesh");
    std::string meshPath1 = list1.get ("mesh_path", "./");
    std::string meshName2 = list2.get ("mesh_name", "lid64.mesh");
    std::string meshPath2 = list2.get ("mesh_path", "./");

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // Load the mesh                              //
    //********************************************//
    meshPtr_Type mesh1 ( new mesh_Type ( comm ) );
    meshPtr_Type fullMesh1 ( new mesh_Type ( comm ) );
    MeshUtility::fillWithFullMesh (mesh1, fullMesh1, meshName1, meshPath1);

    meshPtr_Type mesh2 ( new mesh_Type ( comm ) );
    meshPtr_Type fullMesh2 ( new mesh_Type ( comm ) );
    MeshUtility::fillWithFullMesh (mesh2, fullMesh2, meshName2, meshPath2);


    //********************************************//
    // Create fe spaces                           //
    //********************************************//
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > coarse
    ( new FESpace< mesh_Type, MapEpetra > ( mesh1, "P1", 3, comm ) );

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > fine
    ( new FESpace< mesh_Type, MapEpetra > ( mesh2, "P1", 3, comm ) );


    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    VectorSmall<3> fibers;
    fibers[0] =  list1.get ("fiber_X", std::sqrt (2) / 2.0 );
    fibers[1] =  list1.get ("fiber_Y", std::sqrt (2) / 2.0 );
    fibers[2] =  list1.get ("fiber_Z", 0.0 );




    //********************************************//
    // Fill the fiber direction in an EpetraVector//
    //********************************************//
    vectorPtr_Type fiber1 ( new vector_Type ( coarse -> map() ) );

    int n1 = (*fiber1).epetraVector().MyLength();
    int d1 = n1 / 3;
    (*fiber1) *= 0;
    int i (0);
    int j (0);
    int k (0);
    for ( int l (0); l < d1; l++)
    {

        i = (*fiber1).blockMap().GID (l);
        j = (*fiber1).blockMap().GID (l + d1);
        k = (*fiber1).blockMap().GID (l + 2 * d1);
        (*fiber1) [i] = fibers[0];
        (*fiber1) [j] = fibers[1];
        (*fiber1) [k] = fibers[2];
    }

    //********************************************//
    // Create the new fiber direction in the finer//
    //              mesh                          //
    //********************************************//
    vectorPtr_Type fiber2 ( new vector_Type ( fine -> map() ) );

    //********************************************//
    // Interpolate the fiber direction in the fine//
    // mesh                                       //
    //********************************************//

    int n2 = (*fiber2).epetraVector().MyLength();
    int d2 = n2 / 3;
    (*fiber2) *= 0;
    int i2 (0);
    int j2 (0);
    int k2 (0);
    for ( int l (0); l < d2; l++)
    {

        i2 = (*fiber2).blockMap().GID (l);
        j2 = (*fiber2).blockMap().GID (l + d2);
        k2 = (*fiber2).blockMap().GID (l + 2 * d2);
        (*fiber2) [i2] = fibers[0];
        (*fiber2) [j2] = fibers[1];
        (*fiber2) [k2] = fibers[2];
    }


    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< mesh_Type > exporter1;
    exporter1.setMeshProcId ( mesh1, comm -> MyPID() );
    exporter1.setPrefix ("CoarseFiberDirection");
    exporter1.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", coarse, fiber1, UInt (0) );
    exporter1.postProcess (0);
    exporter1.closeFile();


    ExporterHDF5< mesh_Type > exporter2;
    exporter2.setMeshProcId ( mesh2, comm -> MyPID() );
    exporter2.setPrefix ("FineFiberDirection");
    exporter2.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", fine, fiber2, UInt (0) );
    exporter2.postProcess (0);
    exporter2.closeFile();


    if ( comm->MyPID() == 0 )
    {
        cout << "\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
