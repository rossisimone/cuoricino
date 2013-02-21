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
    @brief main for the test_heart

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
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

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/heart/solver/HeartMonodomainSolver.hpp>
#include <lifev/heart/solver/HeartIonicSolver.hpp>
#include <lifev/heart/solver/HeartAlievPanfilov.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
//#include "heart.hpp"

using namespace LifeV;

Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

    typedef HeartMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type   vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;

    LifeChrono chronoinitialsettings;
    LifeChrono chronototaliterations;
    chronoinitialsettings.start();

    Real normu;
    Real meanu;
    Real minu;

    //! Construction of data classes
    UInt ion_model;

    boost::shared_ptr<HeartFunctors> M_heart_fct;

    /*********************/
    //old constructor
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //! Pointer to access functors
    M_heart_fct.reset (new HeartFunctors ( dataFile) );
    ion_model = dataFile ("electric/physics/ion_model", 6);
    M_heart_fct->M_comm.reset (new Epetra_MpiComm ( MPI_COMM_WORLD ) );

    if (!M_heart_fct->M_comm->MyPID() )
    {
        std::cout << "My PID = " << M_heart_fct->M_comm->MyPID() << std::endl;
    }
    /****************************/


    HeartIonicData dataIonic (M_heart_fct->M_dataFile);

    MeshData meshData;
    meshData.setup (M_heart_fct->M_dataFile, "electric/space_discretization");
    boost::shared_ptr<mesh_Type > fullMeshPtr (new mesh_Type ( M_heart_fct->M_comm ) );
    readMesh (*fullMeshPtr, meshData);
    bool verbose = (M_heart_fct->M_comm->MyPID() == 0);


    const ReferenceFE*    refFE_u;
    const QuadratureRule* qR_u;
    const QuadratureRule* bdQr_u;


    //! Construction of the partitioned mesh
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_heart_fct->M_comm);

    //! Initialization of the FE type and quadrature rules for both the variables
    refFE_u = &feTetraP1;
    qR_u    = &quadRuleTetra15pt;
    bdQr_u  = &quadRuleTria3pt;


    //! Construction of the FE spaces
    if (verbose)
    {
        std::cout << "Building the potential FE space ... " << std::flush;
    }

    feSpacePtr_Type uFESpacePtr ( new feSpace_Type ( meshPart,
                                                     *refFE_u,
                                                     *qR_u,
                                                     *bdQr_u,
                                                     1,
                                                     M_heart_fct->M_comm) );


    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    UInt totalUDof  = uFESpacePtr->map().map (Unique)->NumGlobalElements();
    if (verbose)
    {
        std::cout << "Total Potential DOF = " << totalUDof << std::endl;
    }

    if (verbose)
    {
        std::cout << "Calling the ionic model constructor ... ";
    }
    boost::shared_ptr< HeartAlievPanfilov < mesh_Type > > ionicModel;

    if (verbose)
    {
        std::cout << "Ion Model = Aliev-Panfilov" << std::endl << std::flush;
    }
    ionicModel.reset (new HeartAlievPanfilov< mesh_Type > (dataIonic,
                                                           *meshPart.meshPartition(),
                                                           *uFESpacePtr,
                                                           *M_heart_fct->M_comm) );



    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    ionicModel->initialize( );


    std::cout << "buildsystem ok" << std::endl;
    //! Initialization
    Real dt     = dataIonic.timeStep();
    Real t0     = 0;
    Real tFinal = dataIonic.endTime();
    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose)
    {
        std::cout << "Setting the initial solution ... " << std::endl << std::endl;
    }
    dataIonic.setTime (t0);

    if (verbose)
    {
        std::cout << " ok " << std::endl;
    }

    //! Setting generic Exporter postprocessing
    boost::shared_ptr< Exporter<mesh_Type > > exporter;
    std::string const exporterType =  M_heart_fct->M_dataFile ( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<mesh_Type > ( M_heart_fct->M_dataFile,
                                                        "heart_ionicModel" ) );
        exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
        exporter->setMeshProcId ( meshPart.meshPartition(), M_heart_fct->M_comm->MyPID() );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<mesh_Type > ( M_heart_fct->M_dataFile,
                                                             meshPart.meshPartition(),
                                                             "heart",
                                                             M_heart_fct->M_comm->MyPID() ) );
        }
        else
        {
            exporter.reset ( new ExporterEnsight<mesh_Type > ( M_heart_fct->M_dataFile,
                                                               meshPart.meshPartition(),
                                                               "heart",
                                                               M_heart_fct->M_comm->MyPID() ) );
        }
    }


    vectorPtr_Type Uptr ( new vector_Type (ionicModel->getPotential(), Repeated ) );
    vectorPtr_Type rptr ( new vector_Type (ionicModel->getRecoveryVariable(), Repeated ) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential", uFESpacePtr,
                            Uptr, UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,  "recovery variable", uFESpacePtr,
                            rptr, UInt (0) );


    exporter->postProcess ( 0 );

    MPI_Barrier (MPI_COMM_WORLD);
    chronoinitialsettings.stop();

    //! Temporal loop
    LifeChrono chrono;
    Int iter = 1;
    chronototaliterations.start();
    for ( Real time = t0 + dt ; time <= tFinal + dt / 2.; time += dt, iter++)
    {
        dataIonic.setTime (time);
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time " << dataIonic.time() << " s. " << std::endl;
            std::cout << std::endl;
        }
        chrono.start();
        MPI_Barrier (MPI_COMM_WORLD);
        ionicModel->solveIonicModel ( *M_heart_fct, dataIonic.timeStep(), time );

        normu = ionicModel->getPotential().norm2();
        ionicModel->getPotential().epetraVector().MeanValue (&meanu);
        ionicModel->getPotential().epetraVector().MaxValue (&minu);

        if (verbose)
        {
            std::cout << "norm u " << normu << std::endl;
            std::cout << "mean u " << meanu << std::endl;
            std::cout << "max u " << minu << std::endl << std::flush;
        }

        *Uptr = ionicModel->getPotential();
        *rptr = ionicModel->getRecoveryVariable();


        exporter->postProcess ( time );
        MPI_Barrier (MPI_COMM_WORLD);
        chrono.stop();
        if (verbose)
        {
            std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
        }
        chronototaliterations.stop();
    }

    if (verbose)
    {
        std::cout << "Total iterations time " << chronototaliterations.diff() << " s." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total initial settings time " << chronoinitialsettings.diff() << " s." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total execution time " << chronoinitialsettings.diff() + chronototaliterations.diff() << " s." << std::endl;
    }


    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
