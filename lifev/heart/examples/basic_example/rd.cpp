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

#include "heart.hpp"

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;

Real init_u1 ( const Real& /* t */,
               const Real& /* x */,
               const Real& /* y */,
               const Real& /* z */,
               const ID& /* i */ )
{
    return 5.85;
}
Real init_u2 ( const Real& /* t */,
               const Real& /* x */,
               const Real& /* y */,
               const Real& /* z */,
               const ID& /* i */ )
{
    return 0.12;
}

// ===================================================
//! Constructors
// ===================================================

ReactionDiffusion::ReactionDiffusion ( Int argc,
                                       char** argv )
{
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //! Pointer to access functors
    M_heart_fct.reset (new ReactionDiffusionFunctors ( dataFile) );
    M_heart_fct->M_comm.reset (new Epetra_MpiComm ( MPI_COMM_WORLD ) );

    if (!M_heart_fct->M_comm->MyPID() )
    {
        std::cout << "My PID = " << M_heart_fct->M_comm->MyPID() << std::endl;
    }
}


// ===================================================
//! Methods
// ===================================================

void
ReactionDiffusion::run()
{
    typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

    LifeChrono chronoinitialsettings;
    LifeChrono chronototaliterations;
    chronoinitialsettings.start();
    Real normu;
    Real meanu;
    Real minu;

    //! Construction of data classes

    ReactionDiffusionData _data (M_rd_fct);


    MeshData meshData;
    meshData.setup (M_rd_fct->M_dataFile, "space_discretization");
    boost::shared_ptr<mesh_Type > fullMeshPtr (new mesh_Type);
    readMesh (*fullMeshPtr, meshData);
    bool verbose = (M_heart_fct->M_comm->MyPID() == 0);

    //! Boundary conditions handler and function
    BCFunctionBase uZero ( zero_scalar );
    BCHandler bcH;
    bcH.addBC ("All", ALL, Natural,  Full,   uZero,  1 );

    const ReferenceFE*    refFE_u;
    const QuadratureRule* qR_u;
    const QuadratureRule* bdQr_u;


    //! Construction of the partitioned mesh
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_rd_fct->M_comm);
    std::string u1order =  M_rd_fct->M_dataFile ( "space_discretization/u1_order", "P1");

    //! Initialization of the FE type and quadrature rules for both the variables
    if ( u1order.compare ("P1") == 0 )
    {
        if (verbose)
        {
            std::cout << "u1 and u1 : P1 " << std::flush;
        }
        refFE_u = &feTetraP1;
        qR_u    = &quadRuleTetra15pt;
        bdQr_u  = &quadRuleTria3pt;
    }
    else
    {
        cout << "\n " << uOrder << " finite element not implemented yet \n";
        exit (1);
    }
    //! Construction of the FE spaces
    if (verbose)
    {
        std::cout << "Building the species FE space ... " << std::flush;
    }

    feSpacePtr_Type uFESpacePtr ( new feSpace_Type ( meshPart,
                                                     *refFE_u,
                                                     *qR_u,
                                                     *bdQr_u,
                                                     1,
                                                     M_rd_fct->M_comm) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }
    UInt totalUDof  = uFESpacePtr->map().map (Unique)->NumGlobalElements();

    if (verbose)
    {
        std::cout << "DOF per species = " << totalUDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Calling the constructor ... ";
    }

    ReactionDiffusionSolver< mesh_Type > RDModel (_data, *uFESpacePtr, bcH, M_rd_fct->M_comm);


    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }
    MapEpetra fullMap (RDModel.getMap() );
    vector_Type rhs ( fullMap);
    RDModel.setup ( M_rd_fct->M_dataFile );
    std::cout << "setup ok" << std::endl;

    if (verbose)
    {
        std::cout << "Calling the model constructor ... ";
    }
    boost::shared_ptr< ReactionDiffusionSolver< mesh_Type > > RDModel;

    RDModel.initialize ( M_rd_fct->init_u1(), M_rd_fct->init_u2() );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    //! Building time-independent part of the system
    RDModel.buildSystem( );
    std::cout << "buildsystem ok" << std::endl;
    //! Initialization
    Real dt     = _data.timeStep();
    Real t0     = 0;
    Real tFinal = _data.endTime();
    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose)
    {
        std::cout << "Setting the initial solution ... " << std::endl << std::endl;
    }
    _data.setTime (t0);
    RDModel.resetPreconditioner();
    if (verbose)
    {
        std::cout << " ok " << std::endl;
    }

    //! Setting generic Exporter postprocessing
    boost::shared_ptr< Exporter<mesh_Type > > exporter;
    std::string const exporterType =  M_rd_fct->M_dataFile ( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<mesh_Type > ( M_rd_fct->M_dataFile,
                                                        "rd" ) );
        exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
        exporter->setMeshProcId ( meshPart.meshPartition(), M_rd_fct->M_comm->MyPID() );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<mesh_Type > ( M_rd_fct->M_dataFile,
                                                             meshPart.meshPartition(),
                                                             "rd",
                                                             M_rd_fct->M_comm->MyPID() ) );
        }
        else
        {
            exporter.reset ( new ExporterEnsight<mesh_Type > ( M_rd_fct->M_dataFile,
                                                               meshPart.meshPartition(),
                                                               "rd",
                                                               M_rd_fct->M_comm->MyPID() ) );
        }
    }

    vectorPtr_Type U1ptr ( new vector_Type (RDModel.solutionU1(), Repeated ) );
    vectorPtr_Type U2ptr ( new vector_Type (RDModel.solutionU2(), Repeated ) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,  "u1", uFESpacePtr,
                            U1ptr, UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,  "u2", uFESpacePtr,
                            U2ptr, UInt (0) );
    exporter->postProcess ( 0 );
    MPI_Barrier (MPI_COMM_WORLD);
    chronoinitialsettings.stop();

    //! Temporal loop
    LifeChrono chrono;
    Int iter = 1;
    chronototaliterations.start();
    for ( Real time = t0 + dt ; time <= tFinal + dt / 2.; time += dt, iter++)
    {
        _data.setTime (time);
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time " << _data.time() << " s. " << std::endl;
            std::cout << std::endl;
        }
        chrono.start();
        MPI_Barrier (MPI_COMM_WORLD);
        computeRhs ( rhs, RDModel,  _data );
        RDModel.updatePDESystem ( rhs );
        RDModel.PDEiterate ( bcH );
        normu1 = RDModel.solutionU1().norm2();
        normu2 = RDModel.solutionU2().norm2();
        RDModel.solutionU1().epetraVector().MeanValue (&meanu1);
        RDModel.solutionU2().epetraVector().MeanValue (&meanu2);
        if (verbose)
        {
            std::cout << "norm u1 " << normu1 << std::endl;
            std::cout << "norm u2 " << normu2 << std::endl;
            std::cout << "mean u1 " << meanu1 << std::endl;
            std::cout << "mean u2 " << meanu2 << std::endl << std::flush;
        }

        *U1ptr = RDModel.solutionU1();
        *U2ptr = RDModel.solutionU2();

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
}


void ReactionDiffusion::computeRhs ( vector_Type& rhs,
                                     ReactionDiffusionSolver< mesh_Type >& RDModel,
                                     ReactionDiffusionData& data )
{
    bool verbose = (M_rd_fct->M_comm->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type u1VecRep (RDModel.solutionU1(), Repeated);
    vector_Type u2VecRep (RDModel.solutionU2(), Repeated);

    VectorElemental elvec_u1 ( RDModel.U1FESpace().fe().nbFEDof(), 1 ),
                    elvec_u2 ( RDModel.U2FESpace().fe().nbFEDof(), 1 ),
                    elvec_Reaction1 ( RDModel.U1FESpace().fe().nbFEDof(), 1 ),
                    elvec_Reaction2 ( RDModel.U1FESpace().fe().nbFEDof(), 1 );
    for (UInt iVol = 0; iVol < RDModel.U1FESpace().mesh()->numVolumes(); ++iVol)
    {
        RDModel.U1FESpace().fe().updateJacQuadPt ( RDModel.U1FESpace().mesh()->volumeList ( iVol ) );
        elvec_u1.zero();
        elvec_u2.zero();
        elvec_Reaction1.zero();
        elvec_Reaction2.zero();
        UInt eleIDu = RDModel.U1FESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) RDModel.U1FESpace().fe().nbFEDof();
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = RDModel.U1FESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u1.vec() [ iNode ] = u1VecRep[ig];
            elvec_u2.vec() [ iNode ] = u2VecRep[ig];
        }

        UInt eleID = RDModel.U1FESpace().fe().currentLocalId();

        // Computation of the Reaction terms
        Real u1_ig, u2_ig;
        for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
        {
            u1_ig = u2_ig = 0.;
            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {
                u1_ig += elvec_u1 ( i ) * uFESpace.fe().phi ( i, ig );
            }
            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {
                u2_ig += elvec_u2 ( i ) * uFESpace.fe().phi ( i, ig );
            }

            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {
                elvec_Reaction1 ( i ) += (- (beta + 1.0) * u1_ig + pow (u1_ig, 2) * u2_ig + alpha) *
                                         uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
                elvec_Reaction2 ( i ) += (beta * u1_ig - pow (u1_ig, 2) * u2_ig) *
                                         uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
            }
        }

        UInt totalUDof  = RDModel.U1FESpace().map().map (Unique)->NumGlobalElements();

        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = RDModel.U1FESpace().dof().localToGlobalMap ( eleIDu, iNode );
            rhs.sumIntoGlobalValues (ig, elvec_Reaction1.vec() [iNode] );
            rhs.sumIntoGlobalValues (ig + totalUDof, elvec_Reaction2.vec() [iNode] );
        }
    }
    rhs.globalAssemble();

    rhs += RDModel.matrMass() * RDModel.BDFU1().time_der (data.timeStep() );

    MPI_Barrier (MPI_COMM_WORLD);

    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}

