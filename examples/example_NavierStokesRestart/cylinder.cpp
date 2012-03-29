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
   \file cylinder.cpp
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
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

#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/ExporterHDF5.hpp>
#endif
#include <life/lifefilters/ExporterEnsight.hpp>

#include <life/lifesolver/OseenSolver.hpp>

#include "cylinder.hpp"
#include "exactSolution.hpp"

#include <iostream>



using namespace LifeV;



const int INLET       = 2;
const int WALL        = 1;
const int OUTLET      = 3;
const int RINGIN      = 20;
const int RINGOUT     = 30;




Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

Real u2(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

  //Parameters
  Real R = 0.5; //radius
  Real mu = 0.035; //Dynamic viscosity (The density = 1.0 )
  Real Re = 300; //Reynols
  Real L = 10; //Length
  Real Vavg = ( Re * mu ) / ( 2 * R ); //Characteristic velocity
  Real Vmax = 2 * Vavg; //Maximum velocity
  Real DeltaP = - ( 8 * Vavg * mu * L ) / ( R * R ); //Pressure drop

  /*
  std::cout << "The average velocity: "<< Vavg << std::endl;
  std::cout << "The maximum velocity: "<< Vmax << std::endl;
  std::cout << "The length: "<< L << std::endl;
  std::cout << "The viscosity is: "<< mu << std::endl;
  std::cout << "The Reynolds is: "<< Re << std::endl;
  std::cout << "The Applied pressure drop is: "<< DeltaP << std::endl;
  */

  switch (i)
    {
    case 0:
      return 0.0;
      break;
    case 1:
      return 0.0;
      break;
    case 2:
      //return -(1/(4*mu))*(DeltaP/L)*(R*R- (x*x + y*y));
      return -DeltaP;
    }
    return 0;
}

Real zeroVectorial(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
    case 0:
        return 0.0;
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


void
postProcessFluxesPressures( OseenSolver< RegionMesh<LinearTetra> >& nssolver,
                            BCHandler& bcHandler,
                            const LifeV::Real& t, bool _verbose )
{
    LifeV::Real Q, P;
    UInt flag;

    for ( BCHandler::bcBaseIterator_Type it = bcHandler.begin();
            it != bcHandler.end(); ++it )
    {
        flag = it->flag();

        Q = nssolver.flux(flag);
        P = nssolver.pressure(flag);

        if ( _verbose )
        {
            std::ofstream outfile;
            std::stringstream filenamess;
            std::string filename;

            // file name contains the label
            filenamess << flag;
            // writing down fluxes
            filename = "flux_label" + filenamess.str() + ".m";
            outfile.open(filename.c_str(),std::ios::app);
            outfile << Q << " " << t << "\n";
            outfile.close();
            // writing down pressures
            filename = "pressure_label" + filenamess.str() + ".m";
            outfile.open(filename.c_str(),std::ios::app);
            outfile << P << " " << t << "\n";
            outfile.close();
            // reset ostringstream
            filenamess.str("");
        }
    }

}


struct Cylinder::Private
{
    Private() :
            //check(false),
            nu(1),
            //rho(1),
            H(1), D(1)
            //H(20), D(1)
            //H(0.41), D(0.1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_Type;

    double Re;

    std::string data_file_name;

    double      nu;  /**< viscosity (in m^2/s) */
    //const double rho; /**< density is constant (in kg/m^3) */
    double      H;   /**< height and width of the domain (in m) */
    double      D;   /**< diameter of the cylinder (in m) */
    bool        centered; /**< true if the cylinder is at the origin */

    std::string initial_sol;

    boost::shared_ptr<Epetra_Comm>   comm;
    /**
     * get the characteristic velocity
     *
     * @return the characteristic velocity
     */
    double Ubar() const { return nu*Re/D; }

    /**
     * get the magnitude of the profile velocity
     *
     *
     * @return the magnitude of the profile velocity
     */
    double Um_3d() const { return 9*Ubar()/4; }

    double Um_2d() const { return 3*Ubar()/2; }


    /**
     * u3d 3D velocity profile.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real u3d( const Real& /* t */,
              const Real& /* x */,
              const Real& y,
              const Real& z,
              const ID&   id ) const
    {
        if ( id == 0 )
        {
            if ( centered )
            {
                return Um_3d() * (H+y)*(H-y) * (H+z)*(H-z) / pow(H,4);
            }
            else
            {
                return 16 * Um_3d() * y * z * (H-y) * (H-z) / pow(H,4);
            }
        }
        else
        {
            return 0;
        }
    }

    fct_Type getU_3d()
    {
        fct_Type f;
        f = boost::bind(&Cylinder::Private::u3d, this, _1, _2, _3, _4, _5);
        return f;
    }

    /**
     * u2d flat 2D velocity profile.
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real u2d( const Real& t,
              const Real& /*x*/,
              const Real& /*y*/,
              const Real& /*z*/,
              const ID&   id ) const
    {

        switch (id)
        {
        case 0: // x component
            return 0.0;
            break;
        case 2: // z component
            if ( t <= 0.003 )
                return 1.3332e4;
            //      return 0.01;
            return 0.0;
            break;
        case 1: // y component
            return 0.0;
            //      return 1.3332e4;
            //    else
            //      return 0.0;
            break;
        }
        return 0;
    }

    fct_Type getU_2d()
    {
        fct_Type f;
        f = boost::bind(&Cylinder::Private::u2d, this, _1, _2, _3, _4, _5);
        return f;
    }

    /**
     * one flat (1,1,1)
     *
     * Define the velocity profile at the inlet for the 2D cylinder
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


    Real oneU( const Real& /*t*/,
               const Real& /*x*/,
               const Real& /*y*/,
               const Real& /*z*/,
               const ID&   /*id*/ ) const
    {
        //            if (id == 3)
        return 10.;

        return 0.;
    }

    fct_Type getU_one()
    {
        fct_Type f;
        f = boost::bind(&Cylinder::Private::oneU, this, _1, _2, _3, _4, _5);
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

    d->Re          = dataFile( "fluid/problem/Re", 1. );
    d->nu          = dataFile( "fluid/physics/viscosity", 1. ) /
                     dataFile( "fluid/physics/density", 1. );
    d->H           = 20.;//dataFile( "fluid/problem/H", 20. );
    d->D           =               dataFile( "fluid/problem/D", 1. );
    d->centered    = (bool)        dataFile( "fluid/problem/centered", 0 );
    d->initial_sol = (std::string) dataFile( "fluid/problem/initial_sol", "stokes");
    std::cout << d->initial_sol << std::endl;


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    int ntasks = 0;
    d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "Re = " << d->Re << std::endl
                  << "nu = " << d->nu << std::endl
                  << "H  = " << d->H  << std::endl
                  << "D  = " << d->D  << std::endl;
    }
//    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm.reset( new Epetra_SerialComm() );
#endif

}

void
Cylinder::run()
{
    typedef RegionMesh<LinearTetra>                        Mesh;

    typedef OseenSolver< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                   vectorPtr_Type;
    typedef FESpace< Mesh, MapEpetra >                       feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>                  feSpacePtr_Type;

    // Reading from data file
    //
    GetPot dataFile( d->data_file_name );

//    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (d->comm->MyPID() == 0);

    // Boundary conditions
    BCHandler bcH;

    std::vector<ID> compxy(2);
    compxy[0]=0;compxy[1]=1;

    BCFunctionBase uIn  (  d->getU_2d() );
    BCFunctionBase uOne (  d->getU_one() );
    BCFunctionBase uPois(  d->getU_pois() );

    BCFunctionBase uZero( zero_scalar );
    BCFunctionBase zeroVect(  zeroVectorial );
    BCFunctionBase Stress(  u2 );

    bcH.addBC( "Inlet",   INLET,   Natural,   Full,   Stress, 3 );
    bcH.addBC( "DirichletInlet",    INLET,   Essential, Component, uZero, compxy );
    bcH.addBC( "Ringin",   RINGIN,   Essential,     Full,     uZero  , 3 );
    bcH.addBC( "Ringout",  RINGOUT,  Essential,     Full,     uZero  , 3 );
    bcH.addBC( "Outlet",   OUTLET,   Natural,   Full,  zeroVect, 3 );
    bcH.addBC( "DirichletOutlet",   OUTLET,  Essential, Component, uZero, compxy );
    bcH.addBC( "Wall",     WALL,     Essential,   Full,     uZero, 3 );

    
    // File where the errors are registered
    std::ofstream out_norm;
    std::ofstream out_normP;
    if (verbose)
    {
        out_norm.open("normVel.txt");
        out_norm << "  time   "
        <<"  L2_errorVel    "
        <<"  H1_errorVel    "
        <<"  L2_relErrorVel "
        <<"  H1_relErrorVel \n";
        out_norm.close();

        out_normP.open("normP.txt");
        out_normP << "  time   "
        <<"  L2_errorP    "
        <<"  H1_errorP    "
        <<"  L2_relErrorP "
        <<"  H1_relErrorP \n";
        out_normP.close();


    }

    int numLM = 0; //Because there are not fluxes BCs

    boost::shared_ptr<OseenData> oseenData(new OseenData() );
    oseenData->setup( dataFile );

    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<Mesh> fullMeshPtr(new Mesh);
    readMesh(*fullMeshPtr, meshData);

    MeshPartitioner< Mesh >   meshPart(fullMeshPtr, d->comm);

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Time discretization order " << oseenData->dataTime()->orderBDF() << std::endl;

    //oseenData.meshData()->setMesh(meshPart.meshPartition());

    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;

    feSpacePtr_Type uFESpacePtr( new feSpace_Type( meshPart,uOrder,3,d->comm) );

    if (verbose)
        std::cout << "ok." << std::endl;


    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

    feSpacePtr_Type pFESpacePtr( new feSpace_Type(meshPart, pOrder, 1, d->comm) );

    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalVelDof   = uFESpacePtr->map().map(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpacePtr->map().map(Unique)->NumGlobalElements();


    if (verbose) std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    OseenSolver< RegionMesh<LinearTetra> > fluid (oseenData,
                                                    *uFESpacePtr,
                                                    *pFESpacePtr,
                                                    d->comm, numLM);
    MapEpetra fullMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;

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

    // exactSolution
    AnalyticalSolVelocity uExact;
    AnalyticalSolPressure pExact;
    
    vector_Type u (*fluid.solution(),Repeated);
    vector_Type velocityComp(uFESpacePtr->map(),Repeated);
    vector_Type pressureComp(pFESpacePtr->map(),Repeated);

    Real H1_errorVel(0);
    Real H1_relErrorVel(0);
    Real L2_errorVel(0);
    Real L2_relErrorVel(0);

    Real H1_errorP(0);
    Real H1_relErrorP(0);
    Real L2_errorP(0);
    Real L2_relErrorP(0);


    vector_Type beta( fullMap );
    vector_Type rhs ( fullMap );


    std::string expFileName = dataFile( "exporter/filename", "fluid");
    LifeV::ExporterHDF5<Mesh> exporter( dataFile, meshPart.meshPartition(), expFileName, d->comm->MyPID());

    vectorPtr_Type velAndPressure ( new vector_Type(*fluid.solution(), exporter.mapType() ) );

    exporter.addVariable( ExporterData<Mesh>::VectorField, "velocity", uFESpacePtr,
                          velAndPressure, UInt(0) );

    exporter.addVariable( ExporterData<Mesh>::ScalarField, "pressure", pFESpacePtr,
                          velAndPressure, UInt(3*uFESpacePtr->dof().numTotalDof()) );



    // initialization with stokes solution

    if (d->initial_sol == "stokes")
    {
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

        oseenData->dataTime()->setTime(t0);

        MPI_Barrier(MPI_COMM_WORLD);

        beta *= 0.;
        rhs  *= 0.;

        fluid.updateSystem(0, beta, rhs );
        fluid.iterate( bcH );

//    fluid.postProcess();

        *velAndPressure = *fluid.solution();
        exporter.postProcess( 0 );
        fluid.resetPreconditioner();
    }
    else    if (d->initial_sol == "restart")
    {
        // if (verbose)
        //     std::cout << "  f- Restarting the solver at time " << t0 << " ... " << std::flush;

        std::string start    = dataFile( "fluid/importer/start", "00000");
        std::string filename = dataFile("fluid/importer/filename", "cylinder");

        LifeV::ExporterHDF5<Mesh> importer( dataFile, filename);
        importer.setMeshProcId(uFESpacePtr->mesh(), d->comm->MyPID());

        importer.addVariable( ExporterData<Mesh>::VectorField,
                              "velocity",
                              uFESpacePtr,
                              velAndPressure,
                              UInt ( 0 ) );

        importer.addVariable( ExporterData<Mesh>::ScalarField,
                              "pressure",
                              pFESpacePtr,
                              velAndPressure,
                              3*uFESpacePtr->dof().numTotalDof() );

        exporter.setTimeIndex(importer.importFromTime(0.0));


        //
#if 0
        vectorPtr_Type vel      (new LifeV::VectorEpetra(uFESpacePtr->map(), importer.mapType()));
        vectorPtr_Type pressure (new LifeV::VectorEpetra(pFESpacePtr->map(), importer.mapType()));

        LifeV::ExporterData initSolVel(LifeV::ExporterData<Mesh>::VectorField,
                                       std::string("velocity." + start),
                                       vel,
                                       UInt(0),
                                       uFESpacePtr->dof().numTotalDof(),
                                       UInt(0),
                                       LifeV::ExporterData<Mesh>::Node );

        LifeV::ExporterData initSolPress(LifeV::ExporterData<Mesh>::ScalarField,
                                         std::string("pressure." + start),
                                         pressure,
                                         UInt(0),
                                         pFESpacePtr->dof().numTotalDof(),
                                         UInt(0),
                                         LifeV::ExporterData<Mesh>::Node  );

        importer.rd_var(initSolVel  );
        importer.rd_var(initSolPress);

        velAndPressure->subset(*vel,           vel->getMap(), (UInt) 0,                            0);
        velAndPressure->subset(*pressure, pressure->getMap(), (UInt) 0, uFESpacePtr->dof().numTotalDof());

        std::cout << vel->Norm2() << std::endl;
        std::cout << pressure->Norm2() << std::endl;

        //
#endif

        std::cout << "ok." << std::endl;
        exporter.postProcess( 0. );

        double norm = velAndPressure->norm2();
        if (verbose)
            std::cout << "   f- restart solution norm = " << norm << std::endl;
        fluid.initialize(*velAndPressure);

    }


    bdf.bdfVelocity().setInitialCondition( *fluid.solution() );
    MPI_Barrier(MPI_COMM_WORLD);

    // Temporal loop

    LifeChrono chrono;
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

        double alpha = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / oseenData->dataTime()->timeStep();

	bdf.bdfVelocity().extrapolation(beta);
        bdf.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
        rhs  = fluid.matrixMass()*bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdfVelocity().shiftRight( *fluid.solution() );

        *velAndPressure = *fluid.solution();
	u = *fluid.solution();

        exporter.postProcess( time );

	//computation of the errors
	velocityComp.subset(u,0); //Extracting the velocity
	pressureComp.subset(u, uFESpacePtr->dim() );  //Extracting the pressure
	
	/*
	std::cout << "Norm of the velocity field "<< velocityComp.norm2() << std::endl;
	std::cout << "Norm of the pressure field "<< pressureComp.norm2() << std::endl;
	*/

        L2_errorVel = uFESpacePtr->l2Error(uexact, velocityComp, time ,&L2_relErrorVel);
	H1_errorVel = uFESpacePtr->h1Error(uExact, velocityComp, time ,&H1_relErrorVel);

        L2_errorP = pFESpacePtr->l2Error(pexact, pressureComp, time ,&L2_relErrorP);
	H1_errorP = pFESpacePtr->h1Error(pExact, pressureComp, time ,&H1_relErrorP);


        //save the norm
        if (verbose)
        {
            out_norm.open("normVel.txt", std::ios::app);
            out_norm << time             << "   "
            << L2_errorVel       << "   "
            << H1_errorVel      << "   "
            << L2_relErrorVel << "   "
            << H1_relErrorVel << "\n";
            out_norm.close();


            out_normP.open("normP.txt", std::ios::app);
            out_normP << time             << "   "
            << L2_errorP       << "   "
            << H1_errorP      << "   "
            << L2_relErrorP << "   "
            << H1_relErrorP << "\n";
            out_normP.close();

        }
	

        MPI_Barrier(MPI_COMM_WORLD);

	

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//////////////////////


