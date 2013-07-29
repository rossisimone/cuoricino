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
    @brief Test of the MeshTrasformer class

    @author Antonio Cervone <ant.cervone@gmail.com>
    @contributor
    @maintainer Antonio Cervone <ant.cervone@gmail.com>

    @date 2013-03-19

 */

// ===================================================
//! Includes
// ===================================================
// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/RegionMeshStructured.hpp>
#include <lifev/core/mesh/MeshTransformer.hpp>

using namespace LifeV;

typedef LinearTriangle               elem_Type;
typedef RegionMesh<elem_Type>        mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

int main ( int argc, char** argv )
{
    // verbosity
    bool verbose = true;

    // communicator
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    verbose = ( comm->MyPID() == 0 );
#else
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm );
#endif

    // number of elements for each direction in the mesh
    const UInt numMeshElem = 1;

    // init mesh
    meshPtr_Type mesh ( new mesh_Type ( comm ) );

    // build mesh on the unit square
    regularMesh ( *mesh, 1, std::vector<UInt> ( 2, numMeshElem ), false,
                  Vector3D ( 1.0, 1.0, 1.0 ),
                  Vector3D ( 0.0, 0.0, 0.0 ) );

    // get the trasformer handler
    MeshUtility::MeshTransformer<mesh_Type>& transformer = mesh->meshTransformer();

    // save the points before moving them
    transformer.savePoints();

    // set first transformation
    Vector3D scale ( 2.0, 3.0, 1.0 );
    Vector3D rotate ( 0.0, 0.0, M_PI );
    Vector3D translate ( 2.0, 3.0, 0.0 );

    // apply transformation
    transformer.transformMesh ( scale, rotate, translate );

    // get back to original shape
    transformer.transformMesh ( Vector3D ( 0.5, 1. / 3., 1.0 ), Vector3D ( 0.0, 0.0, M_PI ), Vector3D ( 1.0, 1.0, 0.0 ) );

    // check that the points correspond
    for ( UInt i = 0; i < mesh->numPoints(); i++ )
    {
        Vector3D origPoint = castToVector3D ( transformer.pointInitial ( i ).coordinates() );
        Vector3D  newPoint = castToVector3D ( mesh->point ( i ).coordinates() );
        if ( verbose )
        {
            std::cout << "origPoint " << i << ": " << origPoint << std::endl;
            std::cout << " newPoint " << i << ": " <<  newPoint << std::endl;
        }
        if ( ( origPoint - newPoint).norm() > 1.e-6 )
        {
            return EXIT_FAILURE;
        }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}
