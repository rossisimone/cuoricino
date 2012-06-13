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
/**
    @file
    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author Anna Scotti <anna.scotti@mail.polimi.it>

    @date 2012-03-30
*/

// ===================================================
//! Includes
// ===================================================

#include "darcy.hpp"
#include "user_fun.hpp"

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

// ===================================================
//!              Standard functions
// ===================================================

Real UOne( const Real& /* t */,
           const Real& /* x */,
           const Real& /* y */,
           const Real& /* z */,
           const ID&   /* icomp */)
{
    return 1.;
}

Real UZero( const Real& /* t */,
            const Real& /* x */,
            const Real& /* y */,
            const Real& /* z */,
            const ID&   /* icomp */)
{
    return 0.;
}


// ===================================================
//!                  Private Members
// ===================================================

struct darcy::Private
{
    Private() {}

    // Policy for scalar functions
    typedef boost::function<Real ( const Real&, const Real&,
                                   const Real&, const Real&, const ID& )>
    fct_type;

    // Policy for vector functions
    typedef boost::function<Vector ( const Real&, const Real&,
                                     const Real&, const Real&, const ID& )>
    Vfct_type;

    // Policy for matrix functions
    typedef boost::function<Matrix ( const Real&, const Real&,
                                     const Real&, const Real&,
                                     const std::vector<Real>& )>
    Mfct_type;

    std::string    data_file_name;
    std::string    discretization_section;

    boost::shared_ptr<Epetra_Comm>   comm;

    // Function Types

    fct_type getUOne ( )
    {
        fct_type f;
        f = boost::bind( &UOne, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getUZero ( )
    {
        fct_type f;
        f = boost::bind( &UZero, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getAnalyticalSolution ( )
    {
        fct_type f;
        f = boost::bind( &dataProblem::analyticalSolution, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getAnalyticalFlux ( )
    {
        fct_type f;
        f = boost::bind( &dataProblem::analyticalFlux, _1, _2, _3, _4, _5 );
        return f;
    }

    Mfct_type getInversePermeability ( )
    {
        Mfct_type m;
        m = boost::bind( &dataProblem::inversePermeability, _1, _2, _3, _4 , _5 );
        return m;
    }

    fct_type getSource ( )
    {
        fct_type f;
        f = boost::bind( &dataProblem::source, _1, _2, _3, _4, _5 );
        return f;
    }

    Vfct_type getVectorSource ( )
    {
        Vfct_type f;
        f = boost::bind( &dataProblem::vectorSource, _1, _2, _3, _4, _5 );
        return f;
    }

};

// ===================================================
//!                  Constructors
// ===================================================

darcy::darcy( int argc,
              char** argv )
        : Members( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "darcy";

#ifdef EPETRA_MPI
    Members->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    Members->comm.reset( new Epetra_SerialComm() );
#endif

}

// ===================================================
//!                      Methods
// ===================================================

Real
darcy::run()
{
    using boost::dynamic_pointer_cast;
    using namespace Structured2DLabel;

    typedef LinearTriangle geoElement_Type;

    typedef RegionMesh< geoElement_Type > regionMesh_Type;
    typedef boost::shared_ptr< regionMesh_Type > regionMeshPtr_Type;

    typedef SolverAztecOO solver_Type;

    typedef DarcySolver< regionMesh_Type, solver_Type > darcyLinearSolver_Type;
    typedef boost::shared_ptr< darcyLinearSolver_Type > darcyLinearSolverPtr_Type;

    typedef darcyLinearSolver_Type::vector_Type vector_Type;
    typedef boost::shared_ptr< vector_Type > vectorPtr_Type;

    typedef FESpace< regionMesh_Type, MapEpetra > feSpace_Type;
    typedef boost::shared_ptr< feSpace_Type > feSpacePtr_Type;

    LifeChrono chronoTotal;
    LifeChrono chronoReadAndPartitionMesh;
    LifeChrono chronoBoundaryCondition;
    LifeChrono chronoFiniteElementSpace;
    LifeChrono chronoProblem;
    LifeChrono chronoProcess;
    LifeChrono chronoError;

    // Start chronoTotal for measure the total time for the computation
    chronoTotal.start();

    // Reading from data file
    GetPot dataFile( Members->data_file_name.c_str() );

    // Create the leader process, i.e. the process with MyPID equal to zero
    bool isLeader = ( Members->comm->MyPID() == 0 );


    //
    // The Darcy Solver
    //

    if ( isLeader )
        std::cout << "The Darcy solver" << std::endl << std::flush;

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes
    chronoReadAndPartitionMesh.start();

    // Create the data file
    DarcyData<regionMesh_Type> darcyData;

    // Set up the data
    darcyData.setup( dataFile );

    // Create the mesh file handler
    MeshData meshData;

    // Set up the mesh file
    meshData.setup( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the the mesh
    boost::shared_ptr<regionMesh_Type> fullMeshPtr( new regionMesh_Type );

    // Select if the mesh is structured or not
    if ( meshData.meshType() != "structured" )
    {
        // Set up the mesh
        readMesh( *fullMeshPtr, meshData );
    }
    else
    {
        // Set up the structured mesh
        regularMesh2D( *fullMeshPtr, 0,
                       dataFile( ( Members->discretization_section + "/space_discretization/nx" ).data(), 4 ),
                       dataFile( ( Members->discretization_section + "/space_discretization/ny" ).data(), 4 ),
                       dataFile( ( Members->discretization_section + "/space_discretization/verbose" ).data(), false ),
                       dataFile( ( Members->discretization_section + "/space_discretization/lx" ).data(), 1. ),
                       dataFile( ( Members->discretization_section + "/space_discretization/ly" ).data(), 1. ) );
    }

    // Create the partitioner
    MeshPartitioner< regionMesh_Type >  meshPart;

    // Partition the mesh using ParMetis
    meshPart.doPartition ( fullMeshPtr, Members->comm );

    // Stop chronoReadAndPartitionMesh
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    if ( isLeader )
        std::cout << "Time for read and partition the mesh " <<
                  chronoReadAndPartitionMesh.diff() << std::endl << std::flush;

    // Create the boundary conditions

    // Start chronoBoundaryCondition for measure the total time for create the boundary conditions
    chronoBoundaryCondition.start();

    BCFunctionBase dirichletBDfun, neumannBDfun1, neumannBDfun2;
    BCFunctionRobin robinBDfun;

    dirichletBDfun.setFunction ( dataProblem::dirichlet );
    neumannBDfun1.setFunction  ( dataProblem::neumann1 );
    neumannBDfun2.setFunction  ( dataProblem::neumann2 );
    // dp/dn = first_parameter + second_parameter * p
    robinBDfun.setFunctions_Robin( dataProblem::robin,
                                   Members->getUOne() );

    BCHandler bcDarcy;

    bcDarcy.addBC(    "Top",    TOP,    Natural,    Full,    neumannBDfun2, 1 );
    bcDarcy.addBC( "Bottom", BOTTOM,    Robin,      Scalar,  robinBDfun       );
    bcDarcy.addBC(   "Left",   LEFT,    Essential,  Scalar,  dirichletBDfun   );
    bcDarcy.addBC(  "Right",  RIGHT,    Natural,    Full,    neumannBDfun1, 1 );

    // Stop chronoBoundaryCondition
    chronoBoundaryCondition.stop();

    // The leader process print chronoBoundaryCondition
    if ( isLeader )
    {
        std::cout << "Time for create the boundary conditions handler " <<
                  chronoBoundaryCondition.diff() << std::endl << std::flush;

    }

    // Create the solution spaces

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces
    chronoFiniteElementSpace.start();

    // Primal solution parameters
    const ReferenceFE*    refFE_primal ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_primal    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_primal  ( static_cast<QuadratureRule*>(NULL) );

    refFE_primal = &feTriaP0;
    qR_primal    = &quadRuleTria4pt;
    bdQr_primal  = &quadRuleSeg1pt;

    // Dual solution parameters
    const ReferenceFE*    refFE_dual ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_dual    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_dual  ( static_cast<QuadratureRule*>(NULL) );

    refFE_dual = &feTriaRT0;
    qR_dual    = &quadRuleTria4pt;
    bdQr_dual  = &quadRuleSeg1pt;

    // Interpolate of dual solution parameters
    const ReferenceFE*    refFE_dualInterpolate ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_dualInterpolate    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_dualInterpolate  ( static_cast<QuadratureRule*>(NULL) );

    refFE_dualInterpolate = &feTriaP0;
    qR_dualInterpolate    = &quadRuleTria4pt;
    bdQr_dualInterpolate  = &quadRuleSeg1pt;

    // Hybrid solution parameters
    // hybrid
    const ReferenceFE*    refFE_hybrid ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_hybrid    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_hybrid  ( static_cast<QuadratureRule*>(NULL) );

    refFE_hybrid = &feTriaRT0Hyb;
    qR_hybrid    = &quadRuleTria4pt;
    bdQr_hybrid  = &quadRuleSeg1pt;

    // dual dot outward unit normal
    const ReferenceFE*    refFE_VdotN ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_VdotN    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_VdotN  ( static_cast<QuadratureRule*>(NULL) );

    refFE_VdotN = &feTriaRT0VdotNHyb;
    qR_VdotN    = &quadRuleTria4pt;
    bdQr_VdotN  = &quadRuleSeg1pt;

    // Finite element space of the primal variable
    feSpacePtr_Type p_FESpacePtr( new feSpace_Type( meshPart,
                                                    *refFE_primal,
                                                    *qR_primal,
                                                    *bdQr_primal,
                                                    1,
                                                    Members->comm ) );

    // Finite element space of the dual variable
    feSpacePtr_Type u_FESpacePtr( new feSpace_Type( meshPart,
                                                    *refFE_dual,
                                                    *qR_dual,
                                                    *bdQr_dual,
                                                    1,
                                                    Members->comm ) );

    // Finite element space of the interpolation of dual variable
    feSpacePtr_Type uInterpolate_FESpacePtr( new feSpace_Type( meshPart,
                                                               *refFE_dualInterpolate,
                                                               *qR_dualInterpolate,
                                                               *bdQr_dualInterpolate,
                                                               2,
                                                               Members->comm ) );

    // Vector for the interpolated dual solution
    vectorPtr_Type dualInterpolated( new vector_Type ( uInterpolate_FESpacePtr->map(), Repeated ) );

    // Finite element space of the hybrid variable
    feSpacePtr_Type hybrid_FESpace( new feSpace_Type( meshPart,
                                                      *refFE_hybrid,
                                                      *qR_hybrid,
                                                      *bdQr_hybrid,
                                                      1,
                                                      Members->comm ) );

    // Finite element space of the  outward unit normal variable
    feSpacePtr_Type VdotN_FESpace( new feSpace_Type( meshPart,
                                                     *refFE_VdotN,
                                                     *qR_VdotN,
                                                     *bdQr_VdotN,
                                                     1,
                                                     Members->comm ) );

    // Stop chronoFiniteElementSpace
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace
    if ( isLeader )
        std::cout << "Time for create the finite element spaces " <<
                  chronoFiniteElementSpace.diff() << std::endl << std::flush;

    // Start chronoProblem for measure the total time for create the problem
    chronoProblem.start();

    // Instantiation of the DarcySolver class
    darcyLinearSolverPtr_Type darcySolver;

    darcySolver.reset( new darcyLinearSolver_Type( darcyData, *p_FESpacePtr,
                                                   *u_FESpacePtr, *hybrid_FESpace,
                                                   *VdotN_FESpace,
                                                   Members->comm ) );

    // Stop chronoProblem
    chronoProblem.stop();

    // The leader process print chronoProblem
    darcySolver->getDisplayer().leaderPrint( "Time for create the problem ",
                                             chronoProblem.diff(), "\n" );

    // Process the problem

    // Start chronoProcess for measure the total time for the simulation
    chronoProcess.start ();

    // Setup phase for the linear solver
    darcySolver->setup ();

    // Set the source term
    darcySolver->setSource ( Members->getSource() );

    // Set the vector source term
    darcySolver->setVectorSource ( Members->getVectorSource() );

    // Create the inverse permeability
    inversePermeability < regionMesh_Type > invPerm ( Members->getInversePermeability(),
                                                 *p_FESpacePtr );

    // Set the inverse of the permeability
    darcySolver->setInversePermeability ( invPerm );

    // Set the boudary conditions
    darcySolver->setBC ( bcDarcy );

    // Set the exporter for the solution
    boost::shared_ptr< Exporter< regionMesh_Type > > exporter;

    // Shared pointer used in the exporter for the primal solution
    vectorPtr_Type primalExporter;

    // Shared pointer used in the exporter for the dual solution
    vectorPtr_Type dualExporter;

    // Type of the exporter
    std::string const exporterType =  dataFile( "exporter/type", "ensight");

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare("hdf5") == 0 )
    {
        exporter.reset( new ExporterHDF5< regionMesh_Type > ( dataFile, dataFile( "exporter/file_name", "PressureVelocity" ) ) );
    }
    else
#endif
    {
        if ( exporterType.compare("none") == 0 )
        {
            exporter.reset( new ExporterEmpty< regionMesh_Type > ( dataFile, dataFile( "exporter/file_name", "PressureVelocity" ) ) );
        }
        else
        {
            exporter.reset( new ExporterEnsight< regionMesh_Type > ( dataFile, dataFile( "exporter/file_name", "PressureVelocity" ) ) );
        }
    }

    // Set directory where to save the solution
    exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

    exporter->setMeshProcId( meshPart.meshPartition(), Members->comm->MyPID() );

    // Set the exporter primal pointer
    primalExporter.reset( new vector_Type ( *( darcySolver->primalSolution() ),
                                            exporter->mapType() ) );

    // Add the primal variable to the exporter
    exporter->addVariable( ExporterData< regionMesh_Type >::ScalarField,
                           dataFile( "exporter/name_primal", "Pressure" ),
                           p_FESpacePtr,
                           primalExporter,
                           static_cast<UInt>( 0 ),
                           ExporterData< regionMesh_Type >::SteadyRegime,
                           ExporterData< regionMesh_Type >::Cell );

    // Set the exporter dual pointer
    dualExporter.reset( new vector_Type ( *dualInterpolated, exporter->mapType() ) );

    // Add the variable to the exporter
    exporter->addVariable( ExporterData< regionMesh_Type >::VectorField,
                           dataFile( "exporter/name_dual", "Velocity" ),
                           uInterpolate_FESpacePtr,
                           dualExporter,
                           static_cast<UInt>( 0 ),
                           ExporterData< regionMesh_Type >::SteadyRegime,
                           ExporterData< regionMesh_Type >::Cell );

    // Display the total number of unknowns
    darcySolver->getDisplayer().leaderPrint( "Number of unknowns : ",
                                             hybrid_FESpace->map().map(Unique)->NumGlobalElements(), "\n" );

    // Export the partitioning
    exporter->exportPID( meshPart );

    // Solve the problem

    // Build the linear system and the right hand side
    darcySolver->buildSystem();

    // Solve the linear system
    darcySolver->solve();

    // Post process of the primal and dual variables
    darcySolver->computePrimalAndDual();

    // Save the solution

    // Copy the primal solution to the exporter
    *primalExporter = *( darcySolver->primalSolution() );

    // Interpolate the dual vector field spammed as Raviart-Thomas into a P0 vector field
    *dualInterpolated = uInterpolate_FESpacePtr->feToFEInterpolate( *u_FESpacePtr,
                                                                    *( darcySolver->dualSolution() ) );

    // Copy the dual interpolated solution to the exporter
    *dualExporter = *dualInterpolated;

    // Save the solution into the exporter
    exporter->postProcess( static_cast<Real>(0) );

    // Stop chronoProcess
    chronoProcess.stop();

    // The leader process print chronoProcess
    darcySolver->getDisplayer().leaderPrint( "Time for process ",
                                             chronoProcess.diff(), "\n" );

    // Compute the errors

    // Start chronoError for measure the total time for computing the errors.
    chronoError.start();

    // Compute the error L2 norms
    Real primalL2Norm(0), exactPrimalL2Norm(0), primalL2Error(0), primalL2RelativeError(0);
    Real dualL2Norm(0), exactDualL2Norm(0), dualL2Error(0), dualL2RelativeError(0);

    // Norms and errors for the pressure
    darcySolver->getDisplayer().leaderPrint( "\nPRESSURE ERROR\n" );

    // Compute the L2 norm for the primal solution
    primalL2Norm = p_FESpacePtr->l2Norm( *( darcySolver->primalSolution() ) );

    // Display the L2 norm for the primal solution
    darcySolver->getDisplayer().leaderPrint( " L2 norm of primal unknown:            ",
                                             primalL2Norm, "\n" );

    // Compute the L2 norm for the analytical primal
    exactPrimalL2Norm = p_FESpacePtr->l2NormFunction( Members->getAnalyticalSolution(),
                                                      darcyData.dataTime()->endTime() );

    // Display the L2 norm for the analytical primal
    darcySolver->getDisplayer().leaderPrint( " L2 norm of primal exact:              ",
                                             exactPrimalL2Norm, "\n" );

    // Compute the L2 error for the primal solution
    primalL2Error = p_FESpacePtr->l2ErrorWeighted( Members->getAnalyticalSolution(),
                                                   *( darcySolver->primalSolution() ),
                                                   Members->getUOne(),
                                                   darcyData.dataTime()->endTime() );

    // Display the L2 error for the primal solution
    darcySolver->getDisplayer().leaderPrint( " L2 error of primal unknown:           ",
                                             primalL2Error, "\n" );

    // Compute the L2 realative error for the primal solution
    primalL2RelativeError = primalL2Error / exactPrimalL2Norm;

    // Display the L2 relative error for the primal solution
    darcySolver->getDisplayer().leaderPrint( " L2 relative error of primal unknown:  ",
                                             primalL2RelativeError, "\n" );

    // Norms and errors for the interpolated dual
    darcySolver->getDisplayer().leaderPrint( "\n\nINTERPOLATED DARCY VELOCITY ERROR\n" );

    // Compute the L2 norm for the interpolated dual solution
    dualL2Norm = uInterpolate_FESpacePtr->l2Norm( *dualInterpolated );

    // Display the L2 norm for the interpolated dual solution
    darcySolver->getDisplayer().leaderPrint( " L2 norm of dual unknown:              ",
                                             dualL2Norm, "\n" );

    // Compute the L2 norm for the analytical dual
    exactDualL2Norm = uInterpolate_FESpacePtr->l2NormFunction( Members->getAnalyticalFlux(),
                                                               darcyData.dataTime()->endTime() );

    // Display the L2 norm for the analytical dual
    darcySolver->getDisplayer().leaderPrint( " L2 norm of dual exact:                ",
                                             exactDualL2Norm, "\n" );

    // Compute the L2 error for the dual solution
    dualL2Error = uInterpolate_FESpacePtr->l2Error( Members->getAnalyticalFlux(),
                                                    *dualInterpolated,
                                                    darcyData.dataTime()->endTime(),
                                                    NULL );

    // Display the L2 error for the dual solution
    darcySolver->getDisplayer().leaderPrint( " L2 error of dual unknown:             ",
                                             dualL2Error, "\n" );

    // Compute the L2 relative error for the dual solution
    dualL2RelativeError = dualL2Error / exactDualL2Norm;

    // Display the L2 relative error for the dual solution
    darcySolver->getDisplayer().leaderPrint( " L2 relative error of Dual unknown:    ",
                                             dualL2RelativeError, "\n" );

    // Stop chronoError
    chronoError.stop();

    // The leader process print chronoError
    darcySolver->getDisplayer().leaderPrint( "Time for compute errors ",
                                             chronoError.diff(), "\n" );

    // Stop chronoTotal
    chronoTotal.stop();

    // The leader process print chronoTotal
    darcySolver->getDisplayer().leaderPrint( "Total time for the computation ",
                                             chronoTotal.diff(), "\n" );

    // Return the error, needed for the succes/failure of the test
    return primalL2Error;

} // run
