/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-19

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
   \file ethiersteiman.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-19
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

//#include "life/lifesolver/NavierStokesSolver.hpp"
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifesolver/ADRData.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>
//#include <life/lifefilters/ExporterHDF5.hpp>

#include <life/lifesolver/AdvectionDiffusionReactionSolver.hpp>
#include <life/lifesolver/OseenSolver.hpp>

#include <iostream>

#include "masstransport.hpp"
#include "ud_functions.hpp"


using namespace LifeV;


// const int INLET    = 102;
// const int WALL     = 101;
// const int SLIPWALL = 120;
// const int OUTLET   = 103;

const int INLET    = 2;
const int WALL     = 1;
const int SLIPWALL = 20;
const int OUTLET   = 3;


Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}



std::set<UInt> parseList( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    std::string::size_type commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );
        setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
        stringList = stringList.substr( commaPos+1 );
    }
    setList.insert( atoi( stringList.c_str() ) );
    return setList;
}



struct MassTransport::Private
{
    Private() :
            nu       ( 1. ),
            D        ( 1. ),
            uBar     ( 1. ),
            centered ( true ),
            steady   ( false )
    {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    double         Re;

    std::string    data_file_name;

    double         nu;       /**< viscosity (in m^2/s) */
    double         uBar;     /**< height and width of the domain (in m) */
    double         D;        /**< diameter of the cylinder (in m) */
    bool           centered; /**< true if the cylinder is at the origin */


    bool           steady;
    Epetra_Comm*   comm;

    /**
     * get the characteristic velocity
     *
     * @return the characteristic velocity
     */
    double Ubar() const { return uBar; }

    double Um_3d() const { return 1.*Ubar()/1.; }

    /**
     * u3d 3D velocity profile.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real u3d( const Real& /* t */,
              const Real& x,
              const Real& y,
              const Real& z,
              const ID&   id ) const
    {
        if ( id == 3 )
        {
            if ( centered )
            {
                return 1.;
                return Um_3d() * (D + x)*(D - x) * (D + y)*(D - y) / pow(D,4);
            }
            else
            {
                return 16 * Um_3d() * y * x * (D-y) * (D-x) / pow(D,4);
            }
        }
        else
        {
            return 0;
        }
    }

    fct_type getU_3d()
    {
        fct_type f;
        f = boost::bind(&MassTransport::Private::u3d, this, _1, _2, _3, _4, _5);
        return f;
    }

    Real cc3d( const Real& /* t */,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   id ) const
    {
        return 100.;
    }


    fct_type getCC_3d()
    {
        fct_type f;
        f = boost::bind(&MassTransport::Private::cc3d, this, _1, _2, _3, _4, _5);
        return f;
    }



};








MassTransport::MassTransport( int argc,
                              char** argv ):
        d( new Private )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    d->data_file_name = data_file_name;

    d->Re       = dataFile( "fluid/problem/Re", 1. );

    d->uBar     = dataFile( "fluid/problem/uBar", .1 );
    d->D        = dataFile( "fluid/problem/D", .1 );
    d->centered = (bool)dataFile( "fluid/problem/centered", 1 );

    d->nu       = d->uBar*d->D/d->Re;


    std::cout << "Re    = " << d->Re << std::endl;
    std::cout << "nu    = " << d->nu << std::endl;
    std::cout << "uBar  = " << d->uBar << std::endl;
    std::cout << "D     = " << d->D  << std::endl;


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    d->comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm = new Epetra_SerialComm();
#endif

    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running.\n" << std::endl;
//         std::cout << "Re = " << d->Re << std::endl
//                   << "nu = " << d->nu << std::endl;
    }
}

void
MassTransport::run()
{
    typedef OseenSolver< RegionMesh<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type>                   vector_ptrtype;

    // Reading from data file
    //

    GetPot dataFile( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);


    // Problem definition
    typedef MassTransport Problem;


    //
    // the fluid solver
    //

    if (verbose) std::cout << "\n*********************** fluid solver setup\n\n" << std::flush;


    // Boundary conditions

    BCHandler      bcH;
    BCFunctionBase uIn  (  d->getU_3d() );
    BCFunctionBase uZero(  fZero );

    //cylinder
    bcH.addBC( "Inlet",    INLET,    Essential,   Full,      uIn, 3 );
    bcH.addBC( "Outlet",   OUTLET,   Natural  ,   Full,      uZero, 3 );
    bcH.addBC( "wall",     WALL,     Essential,   Full,      uZero, 3 );
    bcH.addBC( "spliwall", SLIPWALL, Essential,   Full,      uZero, 3 );

    OseenData<RegionMesh<LinearTetra> > oseenData;
    oseenData.setup( dataFile );

    oseenData.viscosity( d->nu/dataFile( "fluid/physics/density", 1. ));


    partitionMesh< RegionMesh<LinearTetra> >   meshPart(*oseenData.meshData()->mesh(), *d->comm);

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Time discretization order " << oseenData.dataTime()->orderBDF() << std::endl;

    oseenData.meshData()->setMesh(meshPart.meshPartition());

    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;

    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    FESpace< RegionMesh<LinearTetra>, MapEpetra > uFESpace(meshPart,uOrder,3,*d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");
    FESpace< RegionMesh<LinearTetra>, MapEpetra > pFESpace(meshPart,pOrder,1,*d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalVelDof   = uFESpace.map().getMap(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpace.map().getMap(Unique)->NumGlobalElements();


    if (verbose) std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;

    if (verbose) std::cout << "Calling the fluid solver constructor ... ";

    OseenSolver< RegionMesh<LinearTetra> > fluid (oseenData,
                                                    uFESpace,
                                                    pFESpace,
                                                    *d->comm);
    MapEpetra fullFluidMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;

    fluid.setUp(dataFile);
    fluid.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);


    //
    // The ADR solver
    //

    if (verbose) std::cout << "\n*********************** adr solver setup \n\n" << std::flush;

    // the boundary conditions

    BCHandler bcADR;
    BCFunctionBase uADR  (  d->getCC_3d() );

    //cylinder
    bcADR.addBC( "Inlet",    INLET,    Essential,   Full,      uADR,  1 );
    bcADR.addBC( "Outlet",   OUTLET,   Natural  ,   Full,      uZero, 1 );
    bcADR.addBC( "wall",     WALL,     Essential,   Full,      uZero, 1 );
    bcADR.addBC( "Splipwall",SLIPWALL, Essential,   Full,      uZero, 1 );

    DataADR<RegionMesh<LinearTetra> > dataADR( dataFile );

    std::string adrOrder =  dataFile( "adr/space_discretization/order", "P1");

    const ReferenceFE*    refFE_adr(0);
    const QuadratureRule* qR_adr(0);
    const QuadratureRule* bdQr_adr(0);

    if ( adrOrder.compare("P1") == 0 )
    {
        if (verbose) std::cout << "P1 velocity " << std::flush;
        refFE_adr = &feTetraP1;
        qR_adr    = &quadRuleTetra4pt; // DoE 5
        bdQr_adr  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( uOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 velocity ";
        refFE_adr = &feTetraP2;
        qR_adr    = &quadRuleTetra15pt;  // DoE 2
        bdQr_adr  = &quadRuleTria3pt;   // DoE 2
    }

    dataADR.setMesh(meshPart.meshPartition());

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Time discretization order " << dataADR.orderBDF() << std::endl;

    FESpace< RegionMesh<LinearTetra>, MapEpetra > adrFESpace(meshPart,
                                                               *refFE_adr,
                                                               *qR_adr,
                                                               *bdQr_adr,
                                                               1,
                                                               *d->comm);

    UInt totalADRDof   = adrFESpace.map().getMap(Unique)->NumGlobalElements();

    if (verbose) std::cout << "Calling the ADR solver constructor ... \n";

    ADRSolver< RegionMesh<LinearTetra> > adr (dataADR,
                                                adrFESpace,
                                                uFESpace,
                                                *d->comm);

    adr.setUp(dataFile);
    adr.buildSystem();

    // Initialization

    if (verbose) std::cout << "\n*********************** solutions initialization \n\n" << std::flush;


    Real dt     = oseenData.dataTime()->timeStep();
    Real t0     = oseenData.dataTime()->initialTime();
    Real tFinal = oseenData.dataTime()->endTime();


    // bdf object to store the previous solutions

    TimeAdvanceBDFNavierStokes<vector_type> bdf(oseenData.dataTime()->orderBDF());

    // initialization with stokes solution

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Computing the fluid initial solution ... " << std::endl << std::endl;

    oseenData.dataTime()->setTime(t0);

    vector_type betaFluid( uFESpace.map() );
    vector_type rhsFluid ( fullFluidMap );

    MPI_Barrier(MPI_COMM_WORLD);

    betaFluid *= 0.;
    rhsFluid  *= 0.;

    std::string const proj =  dataFile( "fluid/space_discretization/initialization", "proj");
    bool const L2proj( !proj.compare("proj") );

//     fluid.initialize( Problem::uexact, Problem::pexact );
    //fluid.initialize( zero_scalar, zero_scalar);

    betaFluid.subset(fluid.solution());


    if (L2proj)
    {
        fluid.updateSystem( 0., betaFluid, rhsFluid );
        fluid.iterate(bcH);
    }


    vector_type velpressure ( fluid.solution(), Repeated );

//     vel.subset(velpressure);
//     press.subset(velpressure, uFESpace.dim()*uFESpace.fieldDim());


//     fluid.updateSystem(0, beta, rhs );
//     fluid.iterate();


    bdf.bdfVelocity().setInitialCondition( fluid.solution() );

    fluid.resetPreconditioner();

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    exporter.reset( new Ensight<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "rclux", d->comm->MyPID()) );
    // hdf5 exporter, still under development
    //    exporter.reset( new Hdf5exporter<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "ethiersteinman", d->comm->MyPID()) );

    //    Ensight<RegionMesh<LinearTetra> > exporter( dataFile, meshPart.meshPartition(), "ethiersteinman", d->comm->MyPID());


    // adr solver initial solution


    dataADR.setTime(t0);

    betaFluid *= 0.;

    MapEpetra fullAdrMap(adr.getMap());
    vector_type rhsADR ( fullAdrMap );

    rhsADR  *= 0.;
    adrFESpace.interpolateBC(bcADR, rhsADR, 0.);

    rhsADR = adr.matrMass()*rhsADR;
    adr.updateSystem(1., betaFluid, rhsADR);
    adr.iterate(bcADR);


    adr.resetPrec();

    // post processing setup


    vector_ptrtype velAndPressure ( new vector_type(fluid.solution(), Repeated ) );
    vector_ptrtype concentration  ( new vector_type(adr.solution(), Repeated ) );

    exporter->addVariable( ExporterData::VectorField, "velocity", velAndPressure,
                           UInt(0), uFESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::ScalarField, "pressure", velAndPressure,
                           UInt(3*uFESpace.dof().numTotalDof() ),
                           UInt(  pFESpace.dof().numTotalDof() ) );

    exporter->addVariable( ExporterData::ScalarField, "concentration", concentration,
                           UInt(0), UInt(adrFESpace.dof().numTotalDof()));




//    exporter->postProcess( 0 );

    fluid.postProcess();
    adr.postProcess();


    // Temporal loop

    Chrono chrono;
    int iter = 1;


    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        oseenData.dataTime()->setTime(time);

        if (verbose)
        {
            std::cout << "\n\n*********************** Temporal loop : " << std::flush;
            std::cout << "we are now at time "<< oseenData.dataTime()->time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        double alphaNS = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / oseenData.dataTime()->timeStep();

//         betaFluid = bdf.bdfVelocity().extrapolation();

        betaFluid.subset(bdf.bdfVelocity().extrapolation());
        rhsFluid  = fluid.matrixMass()*bdf.bdfVelocity().time_der( oseenData.dataTime()->timeStep() );


//         std::cout << "alphaNS " << alphaNS << std::endl;
//         std::cout << "norm beta " << betaFluid.Norm2() << std::endl;
//         std::cout << "norm rhs  " << rhsFluid.Norm2() << std::endl;

        fluid.updateSystem( alphaNS, betaFluid, rhsFluid );
        fluid.iterate( bcH );

        bdf.bdfVelocity().shiftRight( fluid.solution() );

        velpressure   = fluid.solution();

//         vector_type vel  (uFESpace.map(), Repeated);
//        vector_type press(pFESpace.map(), Repeated);

//         vel.subset  (velpressure);
//        press.subset(velpressure, uFESpace.dim()*uFESpace.fieldDim());


        double alphaADR = 1./oseenData.dataTime()->timeStep();

//        adrFESpace.interpolateBC(bcADR, rhsADR, 0.);
        rhsADR  = adr.matrMass()*adr.solution();
        rhsADR *= 1./oseenData.dataTime()->timeStep();

//        betaFluid *= 0.;
        adr.updateSystem( alphaADR, betaFluid, rhsADR );
        adr.iterate(bcADR);


//         {
        *velAndPressure = fluid.solution();
        *concentration  = adr.solution();
//        exporter->postProcess( time );


        fluid.postProcess();
        adr.postProcess();
//         }


        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//////////////////////


