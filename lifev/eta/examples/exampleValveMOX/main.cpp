/* -*- mode: c++ -*-

   This file is part of the LifeV Applications.

   Author(s): Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   Date: 2011-01-27
   Aymen Laadhari <aymen.laadhari@epfl.ch>
   Date: 2012-12-11
   Marco Fedele <fedele.marco@gmail.com>
   Date: 2013-10-25

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
   \author Aymen Laadhari <aymen.laadhari@epfl.ch>
   \author Marco Fedele <fedele.marco@gmail.com>
   \date 2013-10-25
*/



// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <unistd.h> //to make pause

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
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/GradientRecovery.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
//#include <lifev/core/mesh/RegionMesh3D.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshData.hpp>
//#include <lifev/level_set/fem/LevelSetQRAdapter.hpp>
#include <lifev/core/fem/PostProcessingBoundary.hpp>
#include <lifev/navier_stokes/solver/StabilizationIP.hpp>

#include <iostream>
#include <fstream> // fstream object is used to save data into file

#include "levelset_valve.hpp"
#include "physiological_bc.hpp"


using namespace LifeV;



/*******************************************
***** typedef, macro, global variables *****
********************************************/

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetraStructured<Real> matrix_block_type;
typedef VectorEpetraStructured vector_block_type;
typedef MatrixEpetra<Real> matrix_type;
typedef VectorEpetra vector_type;


#define HEAVISIDE_PSI_1     ( value(1.0) - eval( smooth_heaviside, value(ETlsFESpace, psiNonSolution) ) )
#define DIRAC_PHI_1         ( eval( smooth_dirac, value(ETlsFESpace, phiNonSolution)) )
#define HEAVISIDE_PSI_2     ( value(1.0) - eval( smooth_heaviside, value(ETlsFESpace, psiLeftSolution) ) )
#define DIRAC_PHI_2         ( eval( smooth_dirac, value(ETlsFESpace, phiLeftSolution)) )
#define HEAVISIDE_PSI_3     ( value(1.0) - eval( smooth_heaviside, value(ETlsFESpace, psiRightSolution) ) )
#define DIRAC_PHI_3         ( eval( smooth_dirac, value(ETlsFESpace, phiRightSolution)) )

//#define BDF2_TIME           1       // to put in the dataFile


// Physics
bool leafletsPosition = false;
Real angleMin = 5.*M_PI/180.; //5 degrees
Real angleMax = 75.*M_PI/180.; //75 degrees
Real angle(angleMin); //Initially AV has closed position
Real dotAngle = 0.0; //mixte variable needed for RK2
Real thresholdAngleCoeff(0.001);

// levelset
//Real epsilon(0.15);
std::string statusValve("moving");
//bool useSmoothHeaviside(true);




void computeAngle ( const Real & Pin, const Real & Pout, const Real & Qin, const Real & currentTime, 
                    const Real & dt, const std::string inletBC, Real & angle, Real &  angleCoeff )
{
    // Compute leaflets angle
    //if ( leafletsPosition == false )
    //{
        // solve angle: Hein method => 2nd order Range-Kutta
        Real const KpI( 40.126032/(unitsFactor*pRescaleFactor) );
        Real const KfI( 50.0 );
        Real const KqI( 2.e-3*unitsFactor*unitsFactor*unitsFactor/uRescaleFactor );
        Real const KvI( 7.e-3*unitsFactor*unitsFactor*unitsFactor/uRescaleFactor );
        Real psiStar, angleStar;

        Real M1 = KpI*(Pin-Pout)+KqI*Qin;//Qout
        Real M2 = KvI*Qin;
        if (Pin<Pout)
            M2 = 0.0;

        psiStar = dotAngle + dt * (-KfI*dotAngle+M1*std::cos(angle)+M2*std::sin(2.0*angle));
        angleStar = angle+dt*dotAngle;

        dotAngle = dotAngle + 0.5*dt*( (-KfI*dotAngle+M1*std::cos(angle)+M2*std::sin(2.0*angle)) + psiStar );
        angle = angle + 0.5*dt*(dotAngle+angleStar);

        if (angle<= angleMin)
        {
            angle = angleMin;
        }
        else if (angle>=angleMax)
        {
            angle = angleMax;
        }
    //}

    if (statusValve == "closed")
        angleCoeff = 0;
    else if (statusValve == "open")
        angleCoeff = 1;
    else if (statusValve == "onoff")
    {
        if (inletBC == "switch")
        {
            if (angleCoeff == 0 && linearFluxIn(currentTime) > 0)
                angleCoeff = 1;
            else if (angleCoeff == 1 && linearFluxIn(currentTime) == 0)
                angleCoeff = 0;
        }
        else
        {
            if (angleCoeff == 0 && (Pin - Pout) > 0)
                angleCoeff = 1;
            else if ( angleCoeff == 1 && Qin < 0 )
                angleCoeff = 0;
        }
    }
    //statusValve = "moving" is default!
    else
        angleCoeff = ( angle - angleMin ) / ( angleMax - angleMin );
};





int main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_SerialComm);
#endif

    // a flag to see who's the leader for output purposes
    bool verbose = (Comm->MyPID() == 0);


    /********************************
    ******* Reading data file *******
    *********************************/
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    unitsFactor = dataFile("mesh/units_factor", 1.);
    Real resistance = dataFile("levelset/resistance", 1.e6);
    epsilon = dataFile("levelset/epsilon", 0.35115);
    statusValve = dataFile("levelset/status_valve", "moving");
    useSmoothHeaviside = dataFile("levelset/smooth_heaviside", true);
    bool useGlobalPsi = dataFile("levelset/global_psi", false);
    bool useGlobalLS = dataFile("levelset/global_ls", false);
    Real viscosity = dataFile("physics/viscosity", 0.035);
    Real density = dataFile("physics/density", 1.0);
    angleMin = dataFile("physics/angle_min", 0.08726);//5.*M_PI/180.);
    angleMax = dataFile("physics/angle_max", 1.308);//75.*M_PI/180.);
    angle = dataFile("physics/angle", angleMin);
    std::string uOrder = dataFile("problem/space_discretization/vel_order", "P1");
    std::string pOrder = dataFile("problem/space_discretization/press_order", "P1");
    std::string lsOrder = dataFile("problem/space_discretization/ls_order", "P1");
    //UInt bdfOrder = dataFile("problem/time/bdf_order", 2);
    Real currentTime = dataFile("problem/time/initial", 0.0);
    Real finalTime = dataFile("problem/time/final", 1.0);
    Real dtDiastole = dataFile("problem/time/diastole_dt", 0.1);
    Real dtSistole = dataFile("problem/time/sistole_dt", 0.1);
    Real dtMovingValve = dataFile("problem/time/movingvalve_dt", 0.01);
    UInt const wallLabel = dataFile("problem/bc/wall_label", 1);
    UInt const inletLabel = dataFile("problem/bc/inlet_label", 3);
    UInt const outletLabel = dataFile("problem/bc/outlet_label", 2);
    std::string inletBC = dataFile("problem/bc/inlet", "p");
    std::string outletBC = dataFile("problem/bc/outlet", "p");
    pRescaleFactor = dataFile("problem/bc/p_rescale_factor", 1.0);
    uRescaleFactor = dataFile("problem/bc/u_rescale_factor", 1.0);
    flowrateCorrection = dataFile("problem/bc/flowrate_correction", 1.0);
    Real gammaBeta = dataFile ( "problem/ipstab/gammaBeta", 1.0);
    Real gammaDiv = dataFile ( "problem/ipstab/gammaDiv", 0.2);
    Real gammaPress = dataFile ( "problem/ipstab/gammaPress", 0.5);
    UInt exportEach = dataFile("exporter/each", 1);

    if (statusValve == "open")
        angleCoeff = 1;

    if (verbose)
    {
        std::cout << std::endl;
        std::cout << " ### Simulation data ### " << std::endl;
        std::cout << " mesh:\t\t\t " << dataFile("mesh/mesh_file", "undefined") << std::endl;
        std::cout << " levelset:\t\t R = " << resistance << ", eps = " << epsilon << \
            ", statusValve = " << statusValve << ", units factor = " << unitsFactor << std::endl;
        std::cout << "\t\t\t use smooth heaviside = " << useSmoothHeaviside << \
            ", use global psi = " << useGlobalPsi << ", use global levelset = " << useGlobalLS << std::endl;
        std::cout << " physics:\t\t viscosity = " << viscosity << ", density = " << density << std::endl;
        std::cout << " FE spaces:\t\t velocity = " << uOrder << ", pressure = " << pOrder << \
            ", levelset = " << lsOrder << std::endl;
        std::cout << " Labels:\t\t inlet = " << inletLabel << ", outlet = " << outletLabel << \
            ", wall = " << wallLabel << std::endl;
        std::cout << " BC:\t\t\t inlet = " << inletBC << ", outlet = " << outletBC << std::endl;
        std::cout << " rescale factor:\t pressure = " << pRescaleFactor << ", velocity = " << uRescaleFactor << std::endl;
        std::cout << " IP stabilization:\t gammaBeta = " << gammaBeta << ", gammaDiv = " << gammaDiv \
            << ", gammaPress = " << gammaPress << std::endl;
        std::cout << " Time:\t\t\t initial = " << currentTime << ", final = " << finalTime << std::endl;
        std::cout << " dt:\t\t\t diastole = " << dtDiastole << ", sistole = " << dtSistole \
            << ", moving valve = " << dtMovingValve << std::endl;
        std::cout << std::endl;
    }


    /*************************************************************************
    ******* Initializing (mesh, FE spaces, solution, solver, exporter) *******
    **************************************************************************/
    if (verbose) std::cout << " Reading and partitioning the mesh " << std::endl;

    // Load the mesh
    MeshData dataMesh;
    dataMesh.setup(dataFile, "mesh");
    boost::shared_ptr < mesh_Type > fullMeshPtr(new mesh_Type);
    readMesh(*fullMeshPtr, dataMesh);
    // Partition the mesh
    MeshPartitioner< mesh_Type >   meshPart(fullMeshPtr, Comm);
    // Free the global mesh
    fullMeshPtr.reset();

    if (verbose) std::cout << " Building FESpaces  " << std::endl;

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace( new FESpace<mesh_Type, MapEpetra> (meshPart, uOrder, 3, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > pFESpace( new FESpace<mesh_Type, MapEpetra> (meshPart, pOrder, 1, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > lsFESpace( new FESpace<mesh_Type, MapEpetra> (meshPart, lsOrder, 1, Comm) );

    if (verbose) std::cout << std::endl << " ### Dof Summary ###: " <<  std::endl;
    if (verbose) std::cout << " Velocity  : " << uFESpace->map().map(Unique)->NumGlobalElements() << std::endl;
    if (verbose) std::cout << " Pressure  : " << pFESpace->map().map(Unique)->NumGlobalElements() << std::endl;
    if (verbose) std::cout << " Level set : " << lsFESpace->map().map(Unique)->NumGlobalElements() << std::endl;
    if (verbose) std::cout << " Building EA FESpaces  " << std::endl;

    boost::shared_ptr<ETFESpace< mesh_Type,MapEpetra,3,3> > ETuFESpace(new ETFESpace<mesh_Type,MapEpetra,3,3>(meshPart,&(uFESpace->refFE()),Comm));
    boost::shared_ptr<ETFESpace< mesh_Type,MapEpetra,3,1> > ETpFESpace(new ETFESpace<mesh_Type,MapEpetra,3,1>(meshPart,&(pFESpace->refFE()),Comm));
    boost::shared_ptr<ETFESpace< mesh_Type,MapEpetra,3,1> > ETlsFESpace(new ETFESpace<mesh_Type,MapEpetra,3,1>(meshPart,&(lsFESpace->refFE()),Comm));

    if (verbose) std::cout << " Initial conditions " << std::endl;
    vector_type phiNonSolution(ETlsFESpace->map(),Unique);
    vector_type phiLeftSolution(ETlsFESpace->map(), Unique);
    vector_type phiRightSolution(ETlsFESpace->map(), Unique);
    vector_type psiNonSolution(ETlsFESpace->map(), Unique);
    vector_type psiLeftSolution(ETlsFESpace->map(), Unique);
    vector_type psiRightSolution(ETlsFESpace->map(), Unique);
    vector_type phiGlobalSolution(ETlsFESpace->map(), Unique);
    vector_type psiGlobalSolution(ETlsFESpace->map(), Unique);

    if (useGlobalLS)
    {
        lsFESpace->interpolate(phiGlobalFct, phiGlobalSolution, 0.0); //phi_global
        lsFESpace->interpolate(psiGlobalFct, psiGlobalSolution, 0.0); //psi_global
    }
    else
    {
        lsFESpace->interpolate(phiNonFct, phiNonSolution, 0.0); //phi_0
        lsFESpace->interpolate(phiLeftFct, phiLeftSolution, 0.0); //phiII_0
        lsFESpace->interpolate(phiRightFct, phiRightSolution, 0.0); //phiIII_0
        if (useGlobalPsi)
        {
            lsFESpace->interpolate(psiGlobalFct, psiNonSolution, 0.0); //psi
            lsFESpace->interpolate(psiGlobalFct, psiLeftSolution, 0.0); //psiII
            lsFESpace->interpolate(psiGlobalFct, psiRightSolution, 0.0); //psiIII
        }
        else
        {
            lsFESpace->interpolate(psiNonFct, psiNonSolution, 0.0); //psi
            lsFESpace->interpolate(psiLeftFct, psiLeftSolution, 0.0); //psiII
            lsFESpace->interpolate(psiRightFct, psiRightSolution, 0.0); //psiIII
        }
    }

    //vector_type LSSolutionOld(phiNonSolution, Repeated);

    vector_type pressureSolution(ETpFESpace->map(), Unique);
    pFESpace->interpolate(initPressureFct, pressureSolution, currentTime);
    vector_type velocitySolution(ETuFESpace->map(), Unique);
    uFESpace->interpolate(zeroFct, velocitySolution, 0.0);
    vector_type velocitySolutionOld(velocitySolution, Repeated);
#ifdef BDF2_TIME
    vector_type velocitySolutionOldOld(velocitySolution, Repeated);
#endif

    vector_block_type NSSolution(ETuFESpace->map() | ETpFESpace->map(), Unique);
    NSSolution *= 0.0;
    NSSolution.subset( pressureSolution, ETpFESpace->map(), UInt(0), 3*uFESpace->dof().numTotalDof() );
    vector_block_type NSSolutionOld(NSSolution, Unique);
    //NSSolutionOld *= 0.0;

    if (verbose) std::cout << " Building the solvers " << std::endl;
    SolverAztecOO NSSolver;
    NSSolver.setCommunicator(Comm);
    NSSolver.setDataFromGetPot(dataFile,"solver");
    NSSolver.setupPreconditioner(dataFile,"prec");

    if (verbose) std::cout << " Building the exporter " << std::endl;
    ExporterHDF5<mesh_Type> exporter ( dataFile, meshPart.meshPartition(), "solution", Comm->MyPID());
    exporter.setMultimesh(false);
    boost::shared_ptr<vector_type> phiNonExported( new vector_type (phiNonSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> psiNonExported( new vector_type (psiNonSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> phiLeftExported( new vector_type (phiLeftSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> psiLeftExported( new vector_type (psiLeftSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> phiRightExported( new vector_type (phiRightSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> psiRightExported( new vector_type (psiRightSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> phiGlobalExported( new vector_type (phiGlobalSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> psiGlobalExported( new vector_type (psiGlobalSolution, exporter.mapType()) );
    boost::shared_ptr<vector_type> NSExported( new vector_type (NSSolution, exporter.mapType()) );

    const UInt PressureOffset( 3 * uFESpace->dof().numTotalDof());

    if (useGlobalLS)
    {
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "phi-global", lsFESpace, phiGlobalExported , UInt(0));
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "psi-global", lsFESpace, psiGlobalExported, UInt(0));
    }
    else
    {
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "phi-non", lsFESpace, phiNonExported , UInt(0));
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "psi-non", lsFESpace, psiNonExported, UInt(0));
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "phi-left", lsFESpace, phiLeftExported , UInt(0));
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "psi-left", lsFESpace, psiLeftExported, UInt(0));
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "phi-right", lsFESpace, phiRightExported , UInt(0));
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "psi-right", lsFESpace, psiRightExported, UInt(0));
    }
    exporter.addVariable( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace, NSExported , UInt(0));
    exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace, NSExported , PressureOffset);

    if (verbose) std::cout << " Exporting the initial condition " << std::endl;
    exporter.postProcess(currentTime);


    /***************************
    ********* Time loop ********
    ****************************/
    Real dt(dtDiastole), dtOld(dtDiastole);
    bool dtDecreased(false);
    UInt niter(0), dtDecreased_iter(0);
    Real max_dtChangeFactor(dtDiastole/dtSistole);
    if (statusValve == "moving")
        max_dtChangeFactor = dtDiastole/dtMovingValve;
    Real PinBC = -1.*inletPressureFct(currentTime, 0, 0, 0, 0);
    Real PoutBC = -1.*outletPressureFct(currentTime, 0, 0, 0, 0);
    Real QinBC = linearFluxIn(currentTime) / (unitsFactor * unitsFactor * unitsFactor) * uRescaleFactor;
    Real Pin(PinBC), Pout(PoutBC), Qin(0), Qout(0), Plv(PinBC), Pao(PoutBC), Qlv(0);

    std::ofstream file;
    file.open("solution.dat");
    file << "t,PinBC,PoutBC,QinBC,Pin,Pout,Qin,Qout,Plv,Pao,Qlv" << std::endl;

    while  ( currentTime < finalTime )
    {

        LifeChrono ChronoIteration, ChronoItem;
        ChronoIteration.start();
        niter += 1;
        dtOld = dt;

        if (verbose)
        {
            std::cout << std::endl <<std::endl << "-----------------------------------------------------------------" << std::endl;
            std::cout << "\tIter " << niter << std::endl;
            std::cout << "-----------------------------------------------------------------" << std::endl << std::endl;
        }

        /*********************************
        ***** Status Valve Algorithm *****
        **********************************/

        if (verbose) std::cout << "[Status valve] computations ... " << std::flush;
        ChronoItem.start();

        if ( (angleCoeff<=thresholdAngleCoeff) ) //ie. angle<=angleMin+eps
        {
            if (PinBC - PoutBC > 0.)
                leafletsPosition = false;
            else
                leafletsPosition = true;
        }
        else if (angleCoeff>=1.-thresholdAngleCoeff) //ie. angle>=angleMax-eps
        {
            if (PinBC - PoutBC < 0.)
                leafletsPosition = false;
            else
                leafletsPosition = true;
        }
        else if ( (angleCoeff>thresholdAngleCoeff) && (angleCoeff<1.-thresholdAngleCoeff) )
            leafletsPosition = false;
        else
            leafletsPosition = true;
        if ( leafletsPosition == false )
            dt = dtMovingValve;
        currentTime += dt;
        PinBC = -1.*inletPressureFct(currentTime, 0, 0, 0, 0);
        PoutBC = -1.*outletPressureFct(currentTime, 0, 0, 0, 0);
        QinBC = linearFluxIn(currentTime) / (unitsFactor * unitsFactor * unitsFactor) * uRescaleFactor;
        computeAngle(Plv, Pao, Qlv, currentTime, dt, inletBC, angle, angleCoeff);

        // update dt and using angleCoeff value
        currentTime -= dt;
        if (angleCoeff == 0)
            dt = dtDiastole;
        else if (angleCoeff == 1)
            dt = dtSistole;
        else
            dt = dtMovingValve;
        if (dtDecreased && dt > dtOld)
            dt = dtOld;
        if ( dt < dtOld ) // to avoid up and down in few iterations
        {
            dtDecreased = true;
            dtDecreased_iter = niter;
        }
        if( dtDecreased && niter > (dtDecreased_iter + max_dtChangeFactor) )
        {
            dtDecreased = false;
        }
        // update PinBC, PoutBC and angle
        currentTime += dt;
        PinBC = -1.*inletPressureFct(currentTime, 0, 0, 0, 0);
        PoutBC = -1.*outletPressureFct(currentTime, 0, 0, 0, 0);
        QinBC = linearFluxIn(currentTime) / (unitsFactor * unitsFactor * unitsFactor) * uRescaleFactor;
        computeAngle(Plv, Pao, Qlv, currentTime, dt, inletBC, angle, angleCoeff);

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;

        if (verbose)
        {
            std::cout << "\tangle coefficient = " << angleCoeff << ",\tleaflets position = " << leafletsPosition <<
                ",\tangle [" << angleMin*180./M_PI << "," << angleMax*180./M_PI << "] = " << angle*180./M_PI << std::endl;
            std::cout << "\tTime = " << currentTime << ",\tdt = " << dt << std::endl;
            std::cout << "\tBC: Pin = " << PinBC << ", Pout = " << PoutBC << ", Pin - Pout = " << PinBC - PoutBC \
                << ", Qin = " << QinBC << std::endl << std::endl;
        }

        if (verbose) std::cout << "[Status Valve] update levelset with new angle ... " << std::flush;
        ChronoItem.start();

        if (useGlobalLS)
        {
            lsFESpace->interpolate(phiGlobalFct, phiGlobalSolution, 0.0); //interpolate phi_ global with the new angle
            lsFESpace->interpolate(psiGlobalFct, psiGlobalSolution, 0.0);
        }
        else
        {
            lsFESpace->interpolate(phiNonFct, phiNonSolution,0.0); //interpolate phi_0 with the new angle
            lsFESpace->interpolate(phiLeftFct, phiLeftSolution,0.0); //phiII_0
            lsFESpace->interpolate(phiRightFct, phiRightSolution,0.0); //phiIII_0
            if (useGlobalPsi)
            {
                lsFESpace->interpolate(psiGlobalFct, psiNonSolution, 0.0); //psi
                lsFESpace->interpolate(psiGlobalFct, psiLeftSolution, 0.0); //psiII
                lsFESpace->interpolate(psiGlobalFct, psiRightSolution, 0.0); //psiIII
            }
            else
            {
                lsFESpace->interpolate(psiNonFct, psiNonSolution, 0.0); //psi
                lsFESpace->interpolate(psiLeftFct, psiLeftSolution, 0.0); //psiII
                lsFESpace->interpolate(psiRightFct, psiRightSolution, 0.0); //psiIII
            }
        }

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


        /***********************************
        ******* Navier-Stokes Matrix *******
        ************************************/
        if (verbose) std::cout << "[Navier-Stokes] Assembling the matrix ... " << std::flush;
        ChronoItem.start();

#ifdef BDF2_TIME
        vector_type velocityExtrapolated(velocitySolutionOld*2-velocitySolutionOldOld, Repeated);
        Real alpha(1.5);
        vector_type velocityBdfRHS(velocitySolutionOld*2-velocitySolutionOldOld*0.5, Repeated);
#else
        vector_type velocityExtrapolated(velocitySolutionOld, Repeated);
        Real alpha(1.0);
        vector_type velocityBdfRHS(velocitySolutionOld, Repeated);
#endif

        boost::shared_ptr<matrix_block_type> NSMatrix(new matrix_block_type( ETuFESpace->map() | ETpFESpace->map() ));
        *NSMatrix *= 0.0;

        {
            boost::shared_ptr<SmoothDeltaFct> smooth_dirac(new SmoothDeltaFct);
            boost::shared_ptr<SmoothHeavisideFct> smooth_heaviside(new SmoothHeavisideFct);
            boost::shared_ptr<ValveInterfaceFunctor> valve_interface(new ValveInterfaceFunctor);

            using namespace ExpressionAssembly;

            if (useGlobalLS)
            {
                integrate
                (
                    elements(ETuFESpace->mesh()), // Mesh
                    //adapt(ETlsFESpace,LSSolutionOld, uFESpace->qr()), // QR //SAMUEL
                    uFESpace->qr(), //AYMEN
                    ETuFESpace,
                    ETuFESpace,
                    // NS
                    value(density*unitsFactor*unitsFactor*unitsFactor) *
                    (
                        value(alpha/dt) * dot(phi_i,phi_j) +
                        dot(grad(phi_j) * value(ETuFESpace,velocityExtrapolated),phi_i)
                    )
                    + value(viscosity*unitsFactor) * dot(grad(phi_i),grad(phi_j))
                    // Penalization force on leaflets
                    + value(resistance) * eval( valve_interface, X ) * dot(phi_j,phi_i)
                )
                >> NSMatrix->block(0,0);
            }
            else
            {
                integrate
                (
                    elements(ETuFESpace->mesh()), // Mesh
                    //adapt(ETlsFESpace,LSSolutionOld, uFESpace->qr()), // QR //SAMUEL
                    uFESpace->qr(), //AYMEN
                    ETuFESpace,
                    ETuFESpace,
                    // NS
                    value(density*unitsFactor*unitsFactor*unitsFactor)
                    * ( value(alpha/dt) * dot(phi_i,phi_j)
                        + dot(grad(phi_j)*value(ETuFESpace,velocityExtrapolated),phi_i)
                        )
                    + value(viscosity*unitsFactor) * dot(grad(phi_i),grad(phi_j))
                    // Penalization force on leaflets
                    + value(resistance) * (
                        HEAVISIDE_PSI_1 * DIRAC_PHI_1 + // Leaflet I
                        HEAVISIDE_PSI_2 * DIRAC_PHI_2 + // Leaflet II
                        HEAVISIDE_PSI_3 * DIRAC_PHI_3 // Leaflet III
                        )
                    * dot(phi_j,phi_i)
                )
                >> NSMatrix->block(0,0);
            }

            integrate
            (
                elements(ETuFESpace->mesh()), // Mesh
                uFESpace->qr(),
                ETuFESpace,
                ETpFESpace,
                value(-1.0)*phi_j*div(phi_i)
            )
            >> NSMatrix->block(0,1);

            integrate
            (
                elements(ETuFESpace->mesh()), // Mesh
                uFESpace->qr(),
                ETpFESpace,
                ETuFESpace,
                phi_i*div(phi_j)
            )
            >> NSMatrix->block(1,0);

        }

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


        /*****************************
        ******* IP-Stab matrix *******
        ******************************/
        if ( gammaBeta != 0 || gammaDiv != 0 || gammaPress != 0 )
        {
            if (verbose) std::cout << "[Navier-Stokes] Adding IP stabilization to the matrix ... " << std::flush;
            ChronoItem.start();

            MapEpetra fullMap ( ETuFESpace->map() + ETpFESpace->map() );
            details::StabilizationIP<mesh_Type, DOF> M_ipStabilization;
            M_ipStabilization.setFeSpaceVelocity ( *uFESpace );
            M_ipStabilization.setViscosity ( viscosity*unitsFactor );
            M_ipStabilization.setGammaBeta ( gammaBeta );
            M_ipStabilization.setGammaDiv ( gammaDiv );
            M_ipStabilization.setGammaPress ( gammaPress );
            boost::shared_ptr<matrix_type> stabMatrix (new matrix_type ( fullMap ));
            M_ipStabilization.apply ( *stabMatrix, velocityExtrapolated, false );
            stabMatrix->globalAssemble();
            //if(currentTime == 0.26)
            //    stabMatrix->spy("IPMatrix");
            *NSMatrix += *stabMatrix;
            ChronoItem.stop();
            if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        if (verbose) std::cout << "[Navier-Stokes] Closing the matrix ... " << std::flush;
        ChronoItem.start();
        NSMatrix->globalAssemble();
        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


        /********************************
        ******* Navier-Stokes RHS *******
        *********************************/
        if (verbose) std::cout << "[Navier-Stokes] Assembling the rhs ... " << std::flush;
        ChronoItem.start();
        vector_block_type NSRhs( ETuFESpace->map() | ETpFESpace->map(), Unique );
        NSRhs *= 0.0;

        {
            using namespace ExpressionAssembly;

            integrate
            (
                elements(ETuFESpace->mesh()), // Mesh
                uFESpace->qr(),
                ETuFESpace,
                value(density*unitsFactor*unitsFactor*unitsFactor)
                * dot( value(1/dt)*value(ETuFESpace,velocityBdfRHS), phi_i)
            )
            >> NSRhs.block(0);

        }// end namespace ExpressionAssembly

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;
        if (verbose) std::cout << "[Navier-Stokes] Closing the rhs ... " << std::flush;
        ChronoItem.start();
        vector_block_type NSRhsUnique( NSRhs, Unique );
        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


        /*******************************
        ******* Navier-Stokes BC *******
        ********************************/
        if (verbose) std::cout << "[Navier-Stokes] Applying boundary conditions: " << std::flush;
        ChronoItem.start();

        BCHandler NSBCHandler;
        BCFunctionBase ZeroBCFct( zeroFct );
        BCFunctionBase inletPressureBCFct( inletPressureFct );
        BCFunctionBase inletFixedPressureBCFct( inletFixedPressureFct );
        BCFunctionBase outletPressureBCFct( outletPressureFct );
        BCFunctionBase outletFixedPressureBCFct( outletFixedPressureFct );
        BCFunctionBase uInlet( aortaVelIn );

        if (inletBC == "p")
        {
            if (verbose) std::cout << "inlet pressure, " << std::flush;
            NSBCHandler.addBC("Inlet", inletLabel, Natural, Normal, inletPressureBCFct);
        }
        else if (inletBC == "pconst")
        {
            if (verbose) std::cout << "inlet fixed pressure, " << std::flush;
            NSBCHandler.addBC("Inlet", inletLabel, Natural, Normal, inletFixedPressureBCFct);
        }
        if (inletBC == "v")
        {
            if (verbose) std::cout << "inlet velocity, " << std::flush;
            NSBCHandler.addBC("Inlet", inletLabel, Essential, Full, uInlet, 3);
        }
        if (inletBC == "switch")
        {
            if ( linearFluxIn(currentTime)>0 )
            {
                if (verbose) std::cout << "inlet velocity, " << std::flush;
                NSBCHandler.addBC("Inlet", inletLabel, Essential, Full, uInlet, 3);
            }
            else
            {
                if (verbose) std::cout << "inlet pressure, " << std::flush;
                NSBCHandler.addBC("Inlet", inletLabel, Natural, Normal, inletPressureBCFct);
            }
        }

        if (outletBC == "p")
        {
            if (verbose) std::cout << "outlet pressure ... " << std::flush;
            NSBCHandler.addBC("Outlet", outletLabel, Natural, Normal, outletPressureBCFct);
        }
        if (outletBC == "pconst")
        {
            if (verbose) std::cout << "outlet fixed pressure ... " << std::flush;
            NSBCHandler.addBC("Outlet", outletLabel, Natural, Normal, outletFixedPressureBCFct);
        }
        if (outletBC == "nostress")
        {
            if (verbose) std::cout << "outlet nostress ... " << std::flush;
            NSBCHandler.addBC("Outlet", outletLabel, Natural, Full, ZeroBCFct, 3);
        }

        NSBCHandler.addBC("Wall", wallLabel, Essential, Full, ZeroBCFct, 3);

        //Update the FE BCs
        NSBCHandler.bcUpdate( *meshPart.meshPartition(), uFESpace->feBd(), uFESpace->dof() );
        boost::shared_ptr<matrix_type> NSMatrixNoBlock(new matrix_type( uFESpace -> map() + pFESpace -> map(), NSMatrix->matrixPtr() ));
        bcManage(*NSMatrixNoBlock, NSRhsUnique,
                 *uFESpace->mesh(), uFESpace->dof(),
                 NSBCHandler, uFESpace->feBd(), 1.0, currentTime);
        //NSMatrix->diagonalize(3*ETuFESpace->dof().numTotalDof(),1.0,NSRhsUnique,0.0);

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


        /*******************************************
        ******* Navier-Stokes Solving System *******
        ********************************************/
        if (verbose) std::cout << "[Navier-Stokes] Solving the system " << std::endl;
        NSSolver.setMatrix(*NSMatrixNoBlock);
        NSSolver.solveSystem(NSRhsUnique, NSSolution, NSMatrixNoBlock);
        if (verbose) std::cout << "[Navier-Stokes] Time advancing ... " << std::flush;
        ChronoItem.start();

        NSSolutionOld = NSSolution;
#ifdef BDF2_TIME
        velocitySolutionOldOld = velocitySolutionOld;
#endif
        velocitySolutionOld.subset(NSSolutionOld);
        pressureSolution.subset( NSSolution, 3*uFESpace->dof().numTotalDof() );

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl << std::endl;


        /******************************************************************************
        ******* Navier-Stokes Computing Value to next iteration valve algorithm *******
        *******************************************************************************/

        if (verbose) std::cout << "[Navier-Stokes] computing some variables at the end of the iteration  ... " << std::flush;
        ChronoItem.start();

        {
            boost::shared_ptr<PostProcessingBoundary<mesh_Type> > M_postProcessing;
            MapEpetra M_localMap ( uFESpace->map() + pFESpace->map() );//Constructor
            M_postProcessing.reset(
                new PostProcessingBoundary<mesh_Type>( uFESpace->mesh(),
                                                       &uFESpace->feBd(),
                                                       &uFESpace->dof(),
                                                       &pFESpace->feBd(),
                                                       &pFESpace->dof(),
                                                       M_localMap )
                ); //create the object
            Qin = -M_postProcessing->flux( velocitySolutionOld, inletLabel );
            Qout = M_postProcessing->flux( velocitySolutionOld, outletLabel );
            Pin = M_postProcessing->average( pressureSolution, inletLabel ) [0];
            Pout = M_postProcessing->average( pressureSolution, outletLabel ) [0];
            Plv = Pin;//(2.517647, 7.147916, -19.128582)
            Pao = Pout;
            Qlv = Qin;
        }

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;

        if (verbose)
        {
            std::cout << "\tQin = " << Qin << ",\tQout = " << Qout << "\tPin = " << Pin << ",\tPout = " << Pout << std::endl;
            std::cout << "\tQlv = " << Qlv << ",\tPlv = " << Plv << "\tPao = " << Pao << std::endl;
        }

        /***********************************************
        ******* Navier-Stokes Exporting Solution *******
        ************************************************/
        if (verbose) std::cout << "[Navier-Stokes] Exporting: " << std::flush;
        ChronoItem.start();

        file << currentTime << "," << PinBC << "," << PoutBC << "," << QinBC << "," <<
            Pin << "," << Pout << "," << Qin << "," << Qout << "," << Plv << "," << Pao << "," << Qlv << std::endl;

        if (useGlobalLS)
        {
            *phiGlobalExported = phiGlobalSolution;
            *psiGlobalExported = psiGlobalSolution;
        }
        else
        {
            *phiNonExported  = phiNonSolution;
            *psiNonExported = psiNonSolution;
            *phiLeftExported = phiLeftSolution;
            *psiLeftExported = psiLeftSolution;
            *phiRightExported = phiRightSolution;
            *psiRightExported = psiRightSolution;
        }
        *NSExported = NSSolutionOld;

        if (niter%exportEach == 0)
            exporter.postProcess(currentTime);

        ChronoItem.stop();
        if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;

        ChronoIteration.stop();
        if (verbose) std::cout << std::endl << " Total iteration time: " << ChronoIteration.diff() << " s" << std::endl;

    } /****** end time loop ******/

    exporter.closeFile();
    file.close();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return( EXIT_SUCCESS );
}
