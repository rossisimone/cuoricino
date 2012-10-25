/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
             Christian Vergara <>
       Date: 2008-10-30

  Copyright (C) 2008 EPFL

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
   \file cylinder.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2008-06-13
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

#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>

#include <boost/shared_ptr.hpp>
#include <vector>

#include <life/lifesolver/OseenSolver.hpp>

#include "cylinder.hpp"
#include <iostream>
#include <math.h>

#define LM 1

using namespace LifeV;


//cylinder

#define TUBE20_MESH_SETTINGS

#ifdef TUBE20_MESH_SETTINGS
const int INLET    = 2;
const int WALL     = 1;
const int SLIPWALL = 20;
const int OUTLET   = 3;
#endif



Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

Real flux(const Real& t, const ID& i)
{

    // we have to impose (q*n)
    // here : n = ( 0, 0, 1)
    Real const pi(3.14159265);

    switch (i)
    {
    case 0:
        return 0.0;
        break;
    case 2:
        if ( t <= 1 )
            return sin(pi*t);
        return 0.0;
        break;
    case 1:
        return 0.0;
        break;
    }
    return 0;
}



struct Cylinder::Private
{
    Private() :
            nu(1), D(1), period(1),
            lambda()
    {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_Type;

    std::string data_file_name;

    Real Re;
    Real nu;  /**< viscosity (in m^2/s) */
    Real D;
    Real period;
    boost::shared_ptr< std::vector<Real> > lambda;

    void setLambda(boost::shared_ptr< std::vector<Real> >& lam) { lambda = lam;}

    boost::shared_ptr<Epetra_Comm>   comm;
    /**
     * get the characteristic velocity
     *
     * @return the characteristic velocity
     */
    Real Ubar() const { return nu*Re/D; }

    double Um_2d() const { return 3*Ubar()/2; }
    /**
    * Poiseuille inflow
    *
    * Define the velocity profile at the inlet for the 3D cylinder
    */
    Real poiseuille( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& /*z*/,
                     const ID&   id ) const
    {
        double r = std::sqrt(x*x + y*y);

        if (id == 2)
            return Um_2d()*2*((D/2.)*(D/2.) - r*r);

        return 0.;
    }

    fct_Type getU_pois()
    {
        fct_Type f;
        f = boost::bind(&Cylinder::Private::poiseuille, this, _1, _2, _3, _4, _5);
        return f;
    }


    /**
     * u3d 3D velocity profile.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real flux3d( const ID&   id, const Real t  ) const
    {
        return (Ubar() * flux( period*t, id ));
    }

    fct_Type get_flux3d()
    {
        fct_Type f;
        f = boost::bind(&Cylinder::Private::flux3d, this, _1, _2);
        return f;
    }

    /**
     * u3d 3D lagrabge multiplier.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real lambda3d( const Real& /*t*/,
                   const Real& /* x */,
                   const Real& /* y */,
                   const Real& /* z */,
                   const ID&   /*id*/ ) const
    {
        return ( 1.33333 );
        //            return ( (*lambda)[0] );
    }

    fct_Type get_lambda3d()
    {
        fct_Type f;
        f = boost::bind(&Cylinder::Private::lambda3d, this, _1, _2, _3, _4, _5);
        return f;
    }


};

Cylinder::Cylinder( int argc,
                    char** argv )
        :
        d( new Private )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    d->data_file_name = data_file_name;

    d->Re = dataFile( "fluid/problem/Re", 1. );
    d->nu = dataFile( "fluid/physics/viscosity", 1. ) /
            dataFile( "fluid/physics/density", 1. );
    d->period  = dataFile( "fluid/problem/period", 20. );
    d->D  = dataFile( "fluid/problem/D", 1. );

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm.reset( new Epetra_SerialComm() );
#endif

    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "Re = " << d->Re << std::endl
                  << "nu = " << d->nu << std::endl
                  << "period  = " << d->period  << std::endl
                  << "D  = " << d->D  << std::endl;
    }
}

void
Cylinder::run()

{

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef OseenSolver< mesh_Type >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

    // Reading from data file
    //
    GetPot dataFile( d->data_file_name.c_str() );

    int me       = d->comm->MyPID();
    bool verbose = (me == 0);

    // Boundary conditions
    BCHandler bcH;
    BCFunctionBase uZero( zero_scalar );
    std::vector<ID> zComp(1);
    zComp[0] = 3;

    BCFunctionBase uPois(  d->getU_pois() );


//     boost::shared_ptr< std::vector<Real> > lambda;
//     lambda.reset( new std::vector<Real>(1) );
//     (*lambda)[1] = INLET;
//     d-> setLambda(lambda);

    BCFunctionBase lambdaIn  (  d->get_lambda3d() );


#ifdef TUBE20_MESH_SETTINGS
    //cylinder
#if LM
    //bcH.addBC( "InletLM",    INLET,    Flux,        Full,   lambdaIn, 3 );
    bcH.addBC( "InletLM",    INLET,    Flux,   Full,      lambdaIn, 3 );
#else
    bcH.addBC( "Inlet",    INLET,    Essential,   Full,      uPois, 3 );
#endif
    bcH.addBC( "Outlet",   OUTLET,   Natural,     Full,    uZero, 3 );
    //bcH.addBC( "Outlet",   OUTLET,   Flux,     Full,    lambdaIn, 3 );
    bcH.addBC( "Wall",     WALL,     Essential,   Full,    uZero, 3 );
    bcH.addBC( "Slipwall", SLIPWALL, Essential,   Full,    uZero, 3 );

#endif


    // fluid solver
    boost::shared_ptr<OseenData> oseenData(new OseenData());
    oseenData->setup( dataFile );

    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr(new mesh_Type);
    readMesh(*fullMeshPtr, meshData);

    boost::shared_pt<mesh_Type> localMeshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart(fullMeshPtr, d->comm);
        localMeshPtr = meshPart.meshPartition();
    }

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "    Time discretization order " << oseenData->dataTime()->orderBDF() << std::endl;

    //oseenData->meshData()->setMesh(localMeshPtr);

    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;
    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    feSpacePtr_Type uFESpacePtr( new feSpace_Type(localMeshPtr,uOrder,3,d->comm) );

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;
    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");
    feSpacePtr_Type pFESpacePtr( new feSpace_Type(localMeshPtr,pOrder,1,d->comm) );

    if (verbose)
        std::cout << "ok." << std::endl;



    UInt totalVelDofs   = uFESpacePtr->map().map(Unique)->NumGlobalElements();
    UInt totalPressDofs = pFESpacePtr->map().map(Unique)->NumGlobalElements();
    UInt totalDofs      = totalVelDofs + totalPressDofs;


#if LM
    bcH.setOffset("InletLM", totalDofs);
    //bcH.setOffset("Outlet", totalDofs + 1);
#endif

    // Lagrange multipliers for flux imposition
    std::vector<int> lagrangeMultipliers(0);


    if  (d->comm->MyPID() == 0) // Adding lagrange multipliers in the first processor
    {

        lagrangeMultipliers.resize(1); // just one flux to impose (n);
        lagrangeMultipliers[0] = 1;    // just take a numbering of the lagrange multipliers [1:n];
        //            lagrangeMultipliers[1] = 2;    // just take a numbering of the lagrange multipliers [1:n];

    }


    if (verbose) std::cout << "Total Velocity Dofs = " << totalVelDofs << std::endl;
    if (verbose) std::cout << "Total Pressure Dofs = " << totalPressDofs << std::endl;
    if (verbose) std::cout << "Total Dofs          = " << totalDofs << std::endl;

    UInt myTotalVelDofs   = uFESpacePtr->map().map(Unique)->NumMyElements();
    UInt myTotalPressDofs = pFESpacePtr->map().map(Unique)->NumMyElements();
    UInt myTotalDofs      = myTotalVelDofs + myTotalPressDofs;


    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "        " << me << " has " << myTotalDofs << " dofs " << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose) std::cout << "Calling the fluid constructor  ... " << std::flush;

#if LM
    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpacePtr,
                                    *pFESpacePtr,
                                    d->comm,
                                    lagrangeMultipliers.size());
#else
    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpacePtr,
                                    *pFESpacePtr,
                                    d->comm);
#endif

    if (verbose) std::cout << "ok." << std::endl;

    MapEpetra fullMap(fluid.getMap());


    // intialization

    LifeChrono chrono;
    chrono.start();

    fluid.setUp(dataFile);
    fluid.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);

    // Initialization

    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();


    // bdf object to store the previous solutions

    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup(oseenData->dataTime()->orderBDF());

    // initialization with stokes solution

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

    oseenData->dataTime()->setTime(t0);

    vector_Type beta( fullMap );
    vector_Type rhs ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    beta *= 0.;
    rhs  *= 0.;

    fluid.updateSystem(0, beta, rhs );
    fluid.iterate( bcH );

    chrono.stop();
    if (verbose) std::cout << "\n \n -- Total time = " << chrono.diff() << std::endl << std::endl;


    double fluxin  = fluid.flux(2);
    double fluxout = fluid.flux(3);

    if (verbose)
    {
        std::cout << " Inlet Flux  = " << fluxin  << std::endl;
        std::cout << " Outlet Flux = " << fluxout << std::endl;
    }

//    fluid.postProcess();


    bdf.bdfVelocity().setInitialCondition( *fluid.solution() );

    fluid.resetPreconditioner();

    ExporterEnsight<mesh_Type > ensight( dataFile, localMeshPtr, "cylinder", d->comm->MyPID());

    vectorPtr_Type velAndPressure ( new vector_Type(*fluid.solution(), Repeated ) );

    ensight.addVariable( ExporterData<mesh_Type>::VectorField, "velocity", uFESpacePtr,
                         velAndPressure, UInt(0) );

    ensight.addVariable( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                         velAndPressure, UInt(3*uFESpacePtr->dof().numTotalDof() ) );
    ensight.postProcess( 0 );

    // Temporal loop

    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        oseenData->dataTime()->setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time "<< oseenData->dataTime()->time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        Real alpha = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / oseenData->dataTime()->timeStep();
	//Matteo
	//beta = bdf.bdfVelocity().extrapolation();
        bdf.bdfVelocity().extrapolation(beta);
        bdf.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
        rhs  = fluid.matrixMass()*bdf.bdfVelocity().rhsContributionFirstDerivative();


        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdfVelocity().shiftRight( *fluid.solution() );

//         if (((iter % save == 0) || (iter == 1 )))
//         {
        *velAndPressure = *fluid.solution();
        ensight.postProcess( time );
//        fluid.postProcess();
//         }


        MPI_Barrier(MPI_COMM_WORLD);


        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//////////////////////


