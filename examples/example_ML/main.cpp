/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-16

  Copyright (C) 2005 EPFL

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
   \author Gilles Fourestey ( gilles.fourestey@epfl.ch )
   \date 2008-10-16
 */


/*
  Lid-driven cavity problem
*/

// includes and whatnot

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

#include <life/lifecore/LifeV.hpp>

#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>

#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>
#include "MLTester.hpp"

#include <iostream>


const int UPWALL   = 2;
const int WALL     = 1;
const int SLIPWALL = 20;


using namespace LifeV;

typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_Type;

typedef OseenSolver< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
typedef boost::shared_ptr<vector_Type>                   vector_ptrtype;

Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

Real uLid(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
    case 0:
        return 1.0;
        break;
    case 1:
        return 0.0;
        break;
    case 2:
        return 0.0;
        break;
    }
    return 0;
}


int
main( int argc, char** argv )
{

    //
    // Mpi Communicator definition ( see http://www.mpi-forum.org/docs/docs.html for documentation ).
    // The communicator (paralell or sequential) is then given to Epetra.
    // This is standard and can be "copy/pasted"
    //

    boost::shared_ptr<Epetra_Comm> comm;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    comm.reset(new Epetra_MpiComm(MPI_COMM_WORLD));
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
        int ntasks = 0;
//        int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
        std::cout << "My PID = " << comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
#else
    comm.reset(new Epetra_SerialComm);
    cout << "% using serial Version" << endl;
#endif

    // a flag to see who's the leader for output purposes
    bool verbose = comm->MyPID() == 0;

    // We now proceed to the data file. Its name can be given using the
    // -f or --file argument after the name of launch program.
    // By default, it's data.

    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    boost::shared_ptr<OseenData> oseenData;
    oseenData->setup( dataFile );

    // Now for the boundary conditions :
    // BCHandler is the class that stores the boundary conditions. Here we will
    // set 3 boundary conditions :
    // top               : (ux, uy, uz) = (1., 0., 0.) essential BC
    // left, right, down : (ux, uy, uz) = (0., 0., 0.) essential BC
    // front and rear    : uz = 0 essential BC

    BCHandler bcH;

    std::vector<ID> zComp(1);
    zComp[0] = 3;

    BCFunctionBase uIn  ( boost::bind(&uLid, _1, _2, _3, _4, _5) );
    BCFunctionBase uZero( zero_scalar );

    // boundary conditions definition.
    // the first two are classical essential or dirichlet conditions
    bcH.addBC( "Upwall",   UPWALL,   Essential, Full,      uIn,   3 );
    bcH.addBC( "Wall",     WALL,     Essential, Full,      uZero, 3 );
    // this bc is imposed only on some compenants, that is the ones given in zComp
    // Here it's the thirs, ie z, in order to have u.n = 0
    bcH.addBC( "Slipwall", SLIPWALL, Essential, Component, uZero, zComp );

    // Read the mesh
    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra>);
    readMesh(*fullMeshPtr, meshData);

    // Partition the mesh
    MeshPartitioner< RegionMesh<LinearTetra> >   meshPart(fullMeshPtr, comm);

    // Now we proceed with the FESpace definition
    // here we decided to use P2/P1 elements

    // first the velocity FE space

    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;

    boost::shared_ptr<FESpace< RegionMesh<LinearTetra>, MapEpetra > > uFESpacePtr(
                    new FESpace< RegionMesh<LinearTetra>, MapEpetra >(meshPart,"P2",3,comm) );

    // then the pressure FE space

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

    boost::shared_ptr<FESpace< RegionMesh<LinearTetra>, MapEpetra > > pFESpacePtr(
                    new FESpace< RegionMesh<LinearTetra>, MapEpetra >(meshPart,"P1",1,comm) );


    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalVelDof   = uFESpacePtr->map().map(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpacePtr->map().map(Unique)->NumGlobalElements();


    if (verbose) std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;


    // now that the FE spaces are built, we proceed to the NS solver constrution
    // we will use oseen here

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    MLTester< RegionMesh<LinearTetra> > fluid (oseenData,
                                                 *uFESpacePtr,
                                                 *pFESpacePtr,
                                                 comm);


    // this is the total map ( velocity + pressure ). it will be used to create
    // vectors to strore the solutions

    MapEpetra fullMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;

    // Now, the fluid solver is set up using the data file
    fluid.setUp(dataFile);
    // the we build the constant matrices
    fluid.buildSystem();
    // creating the default list
    fluid.createDefaultList(dataFile, "fluid/prec");


    // finally, let's create an exporter in order to view the results
    // here, we use the ensight exporter

    ExporterEnsight<RegionMesh<LinearTetra> > ensight( dataFile, meshPart.meshPartition(), "cavity", comm->MyPID());

    // we have to define a variable that will store the solution
    vector_ptrtype velAndPressure ( new vector_Type(*fluid.solution(), Repeated ) );

    // and we add the variables to be saved
    // the velocity
    ensight.addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",
                         uFESpacePtr, velAndPressure, UInt(0) );

    // and the pressure
    ensight.addVariable( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "pressure",
                         pFESpacePtr, velAndPressure, UInt(3*uFESpacePtr->dof().numTotalDof() ) );

    // everything is ready now
    // a little barrier to synchronize the processes
    MPI_Barrier(MPI_COMM_WORLD);


    // Initialization


//    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
//    Real tFinal = oseenData->dataTime()->endTime();

    // bdf object to store the previous solutions

    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup(oseenData->dataTime()->orderBDF());


    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

    oseenData->dataTime()->setTime(t0);

    // advection speed (beta) and rhs definition using the full map
    // (velocity + pressure)
    vector_Type beta( fullMap );
    vector_Type rhs ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    beta *= 0.;
    rhs  *= 0.;

    // updating the system with no mass matrix, advection and rhs set to zero,
    // that is the stokes problem
    fluid.updateSystem(0, beta, rhs );


    // iterating the solver in order to produce the solution
    fluid.testML( bcH );
    //fluid.iterate( bcH );

// Finalizing the MPI session

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return( EXIT_SUCCESS );
}


