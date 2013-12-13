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
    @brief 1D Test RogersMcCulloch, with a pacing protocol; Map of the APD and deltaAPD; pseudo-ECG

    @date 22.04.2013
    @author Simone Rossi <simone.rossi@epfl.ch> & Marie Dupraz <dupraz.marie@gmail.com>

    @contributor
    @mantainer Marie Dupraz <dupraz.marie@gmail.com>
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
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroIonicSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/electrophysiology/solver/IonicModels/IonicFitzHughNagumo.hpp>
#include <lifev/electrophysiology/solver/ElectroFunctors.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
// ---------------------------------------------------------------
// In order to use the ETA framework, a special version of the
// FESpace structure must be used. It is called ETFESpace and
// has basically the same role as the FESpace.
// ---------------------------------------------------------------

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

//--------------------------------------------------------
// For the pseudo- ECG
//--------------------------------------------------------
#include <lifev/core/function/Norm.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

// ---------------------------------------------------------------
// The most important file to include is the Integrate.hpp file
// which contains all the definitions required to perform the
// different integrations.
// ---------------------------------------------------------------

//#include <lifev/eta/expression/Integrate.hpp>
//
//#include <lifev/eta/expression/ExpressionDot.hpp>

using std::cout;
using std::endl;
using namespace LifeV;

// Choice of the fibers direction : ||.||=1
Real Fibers (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& i)
{
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    Real zmin = 0.0;
    Real zmax = monodomainList.get ("domain_Z", 1. ); // 1.0;
    Real L = zmax - zmin;
    Real thetamin = M_PI / 3.;
    Real thetamax = -M_PI / 3.;

    Real ztheta = (L - z) / L;
    Real theta = (thetamax - thetamin) * ztheta + thetamin;

    switch (i)
    {
        case 0:
            return  std::cos (theta); // x_fib; //y/N; //std::sqrt (2.0) / 2.0; //
            break;
        case 1:
            return  std::sin (theta); // y_fib; //-x/N; //std::sqrt (2.0) / 2.0; //
            break;
        case 2:
            return 0.0;
            break;
        default:
            return 0.0;
            break;
    }
}

Real Stimulus2 (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if ( x <= 0.1 )
    {
        return 80.0;
    }
    else if ( x <= 0.2)
    {
        return 80.0 * ( 0.2 - x ) / (0.1);
    }
    else
    {
        return 0.0;
    }
}
Real Cut (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    if ( y <= 0.05 )
    {
        return 1.0;
    }
    else if ( y <= 0.1)
    {
        return ( 0.1 - y ) / (0.05);
    }
    else
    {
        return 0.0;
    }
}

Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    {
        boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
        if ( Comm->MyPID() == 0 )
        {
            cout << "% using MPI" << endl;
        }

        //********************************************//
        // Starts the chronometer.                    //
        //********************************************//
        LifeChrono chronoinitialsettings;

        chronoinitialsettings.start();

        typedef RegionMesh<LinearTetra>            mesh_Type;
        typedef boost::shared_ptr<mesh_Type>       meshPtr_Type;
        typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;
        typedef FESpace< mesh_Type, MapEpetra >    feSpace_Type;
        typedef boost::shared_ptr<feSpace_Type>    feSpacePtr_Type;
        typedef boost::function < Real (const Real& /*t*/,
                                        const Real &   x,
                                        const Real &   y,
                                        const Real& /*z*/,
                                        const ID&   /*i*/ ) > function_Type;

        typedef ElectroETAMonodomainSolver< mesh_Type, IonicFitzHughNagumo > monodomainSolver_Type;
        typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
        typedef VectorEpetra                        vector_Type;
        typedef MatrixEpetra<Real>                  matrix_Type;
        typedef LifeV::Preconditioner               basePrec_Type;
        typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
        typedef LifeV::PreconditionerML             prec_Type;
        typedef boost::shared_ptr<prec_Type>        precPtr_Type;

        boost::shared_ptr< Exporter<mesh_Type > > exporter;

        //********************************************//
        // Import parameters from an xml list. Use    //
        // Teuchos to create a list from a given file //
        // in the execution directory.                //
        //********************************************//

        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Importing parameters list...";
        }
        Teuchos::ParameterList FHNParameterList = * ( Teuchos::getParametersFromXmlFile ( "FitzHughNagumoParameters.xml" ) );
        Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );

        if ( Comm->MyPID() == 0 )
        {
            std::cout << " Done!" << endl;
        }

        //********************************************//
        // Creates a new model object representing the//
        // model from Fitz-Hugh Nagumo.  The          //
        // model input are the parameters. Pass  the  //
        // parameter list in the constructor          //
        //********************************************//
        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Building Constructor for Fitz-Hugh Nagumo Model with parameters ... ";
        }
        boost::shared_ptr<IonicFitzHughNagumo>  model ( new IonicFitzHughNagumo (FHNParameterList) );
        if ( Comm->MyPID() == 0 )
        {
            std::cout << " Done!" << endl;
        }

        //********************************************//
        // We need the GetPot datafile for to setup   //
        // the preconditioner.                        //
        //********************************************//
        GetPot command_line (argc, argv);
        const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
        GetPot dataFile (data_file_name);

        // +-----------------------------------------------+
        // |               Loading the mesh                |
        // +-----------------------------------------------+
        if ( Comm->MyPID() == 0 )
        {
            std::cout << std::endl << "[Loading the mesh]" << std::endl;
        }

        meshPtr_Type fullMeshPtr ( new mesh_Type ( Comm ) );

        VectorSmall<3> meshDim;
        meshDim[0] =  monodomainList.get ("meshDim_X", 10 );
        meshDim[1] =  monodomainList.get ("meshDim_Y", 10 );
        meshDim[2] =  monodomainList.get ("meshDim_Z", 10 );
        VectorSmall<3> domain;
        domain[0] =  monodomainList.get ("domain_X", 1. );
        domain[1] =  monodomainList.get ("domain_Y", 1. );
        domain[2] =  monodomainList.get ("domain_Z", 1. );

        // Building the mesh from the source
        regularMesh3D ( *fullMeshPtr,
                        1,
                        meshDim[0], meshDim[1], meshDim[2],
                        false,
                        domain[0], domain[1], domain[2],
                        0.0, 0.0, 0.0 );

        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Mesh size  : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).maxH << std::endl;
        }
        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Partitioning the mesh ... " << std::endl;
        }
        meshPtr_Type meshPtr;
        {
            MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, Comm );
            meshPtr = meshPart.meshPartition();
        }
        fullMeshPtr.reset(); //Freeing the global mesh to save memory

        //********************************************//
        // We create three solvers to solve with:     //
        // 1) Operator Splitting method               //
        // 2) Ionic Current Interpolation             //
        // 3) State Variable Interpolation            //
        //********************************************//
        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Building Monodomain Solvers... ";
        }

        monodomainSolverPtr_Type splitting ( new monodomainSolver_Type ( dataFile, model, meshPtr ) );
        const feSpacePtr_Type FESpacePtr =  splitting->feSpacePtr(); //FE Space
        if ( Comm->MyPID() == 0 )
        {
            std::cout << " Splitting solver done... ";
        }

        //********************************************//
        // Setting up the initial condition form      //
        // a given function.                          //
        //********************************************//
        if ( Comm->MyPID() == 0 )
        {
            cout << "\nInitializing potential:  " ;
        }

        function_Type f = &Stimulus2;
        splitting -> setPotentialFromFunction ( f ); //initialize potential

        //setting up initial conditions
        * ( splitting -> globalSolution().at (1) ) = FHNParameterList.get ("W0", 0.011);

        if ( Comm->MyPID() == 0 )
        {
            cout << "Done! \n" ;
        }

        //*******************************************//
        // Setting up the pseudo-ECG                 //
        //*******************************************//
        Real ecg_position_X = monodomainList.get ("ecg_position_X", 1.);
        Real ecg_position_Y = monodomainList.get ("ecg_position_Y", 1.);
        Real ecg_position_Z = monodomainList.get ("ecg_position_Z", 1.);
        vector_Type ecgDistance = (* (splitting->globalSolution().at (0) ) ) ;
        ecgDistance *= 0.;
        FESpacePtr->interpolate ( static_cast<function_Type> ( Norm::f ), ecgDistance, 0.0 );
        Norm::setPosition ( ecg_position_X , ecg_position_Y , ecg_position_Z ); // Set electrode position

        Real pseudoEcgReal (0.);
        Real Global_pseudoEcgReal (0.);
        vector_Type solutionLaplacian = ( (* (splitting->globalSolution().at (0) ) ) );
        vector_Type pseudoEcgVec = ( (* (splitting->globalSolution().at (0) ) ) );

        // Discrete Laplacian matrix
        // setting up the assembler
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        adrAssembler.setup ( FESpacePtr, FESpacePtr );
        // define the matrices
        boost::shared_ptr<matrix_Type> systemMatrixL ( new matrix_Type ( FESpacePtr->map() ) );
        boost::shared_ptr<matrix_Type> systemMatrixM ( new matrix_Type ( FESpacePtr->map() ) );
        boost::shared_ptr<vector_Type> rhs_Laplacian ( new vector_Type ( FESpacePtr->map() ) );
        boost::shared_ptr<vector_Type> pseudoEcgVec_ptr ( new vector_Type ( FESpacePtr->map() ) );

        // fill the matrix
        adrAssembler.addDiffusion ( systemMatrixL, -1.0 );
        adrAssembler.addMass ( systemMatrixM, 1.0 );
        // closed
        systemMatrixL->globalAssemble();
        systemMatrixM->globalAssemble();

        // uncomment to check the matrices with MATLAB
        //    matrix_Type LaplacianMatrix ( *systemMatrixL );
        //    matrix_Type MassMatrix ( *systemMatrixM );
        //    LaplacianMatrix.spy("matriceL_check");
        //    MassMatrix.spy("matriceM_check");

        //********************************************//
        // Setting up the time data                   //
        //********************************************//
        splitting -> setParameters ( monodomainList );

        //********************************************//
        // Create a fiber direction                   //
        //********************************************//
        VectorSmall<3> fibers;
        fibers[0] =  monodomainList.get ("fiber_X", std::sqrt (2.0) / 2.0 );
        fibers[1] =  monodomainList.get ("fiber_Y", std::sqrt (2.0) / 2.0 );
        fibers[2] =  monodomainList.get ("fiber_Z", 0.0 );

        splitting ->setupFibers (fibers);

        //    function_Type Fiber_fct;
        //    Fiber_fct = &Fibers;
        //    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
        //    ( new FESpace< mesh_Type, MapEpetra > ( splitting->localMeshPtr(), "P1", 3, Comm) );
        //    vectorPtr_Type fibers_vect(new vector_Type( Space3D->map() ));
        //    Space3D -> interpolate(static_cast< FESpace < RegionMesh<LinearTetra >, MapEpetra >::function_Type>(Fiber_fct), *fibers_vect,0. );
        //
        //    splitting ->setFiberPtr(fibers_vect);

        //********************************************//
        // Create the global matrix: mass + stiffness //
        //********************************************//
        splitting -> setupLumpedMassMatrix();
        splitting -> setupStiffnessMatrix();
        splitting -> setupGlobalMatrix();

        //********************************************//
        // Creating exporters to save the solution,
        // APD, deltaAPD, pseudo-ECG
        //********************************************//
        ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;

        string filenameSplitting =  monodomainList.get ("OutputFile", "" );
        filenameSplitting += "RM_spiral";

        splitting -> setupPotentialExporter ( exporterSplitting, filenameSplitting );

        splitting -> exportSolution ( exporterSplitting, 0);


        std::ofstream output  ("ecg_output_fibre.txt");

        //********************************************//
        // Solver initialization for the discrete Laplacian //
        //********************************************//
        if (  Comm->MyPID() == 0 )
        {
            std::cout << std::endl << "[Solvers initialization]" << std::endl;
        }
        prec_Type* precRawPtr;
        basePrecPtr_Type precPtr;
        precRawPtr = new prec_Type;
        precRawPtr->setDataFromGetPot ( dataFile, "prec" );
        precPtr.reset ( precRawPtr );
        if (  Comm->MyPID() == 0 )
        {
            std::cout << "Setting up LinearSolver (Belos)... " << std::flush;
        }
        Teuchos::RCP< Teuchos::ParameterList > belosList2 = Teuchos::rcp ( new Teuchos::ParameterList );
        belosList2 = Teuchos::getParametersFromXmlFile ( "SolverParamList2.xml" );
        LinearSolver linearSolver2;
        linearSolver2.setCommunicator ( Comm );
        linearSolver2.setParameters ( *belosList2 );
        linearSolver2.setPreconditioner ( precPtr );
        if (  Comm->MyPID() == 0 )
        {
            std::cout << "done" << std::endl;
        }
        linearSolver2.showMe();

        //********************************************//
        // Solving the monodomain problem             //
        //********************************************//
        if ( Comm->MyPID() == 0 )
        {
            cout << "\nstart solving:  " ;
        }

        Real dt = monodomainList.get ("timeStep", 0.1);
        Real deltat = monodomainList.get ("SaveStep", 1.0);
        Real TF = monodomainList.get ("endTime", 150.0);
        Real TCut1 = monodomainList.get ("TCut", 35.0) - dt / 2.0;
        Real TCut2 = monodomainList.get ("TCut", 35.0) + dt / 2.0;

        int iter ( (deltat / dt) );
        int k (0);
        Int sz = 0;
        sz = (* (splitting->globalSolution().at (0) ) ).size();

        for ( Real t = 0.0; t < TF; )
        {

            //--------------------------------------
            // ECG initialization
            pseudoEcgReal = 0.;
            pseudoEcgVec_ptr.reset ( new vector_Type ( FESpacePtr->map(), Unique ) );
            rhs_Laplacian.reset ( new vector_Type ( FESpacePtr->map(), Unique ) );
            //-----------------------------------------------------------------

            t = t + dt;

            k++;
            if (k % iter == 0)
            {
                splitting -> solveOneSplittingStep (exporterSplitting, t);
            }
            else
            {
                splitting -> solveOneSplittingStep();
            }

            // ECG : discrete laplacian of the solution
            (*rhs_Laplacian) = (*systemMatrixL) * (* (splitting->globalSolution().at (0) ) );

            linearSolver2.setOperator ( systemMatrixM );
            linearSolver2.setRightHandSide ( rhs_Laplacian );
            linearSolver2.solve ( pseudoEcgVec_ptr );

            pseudoEcgVec = (*pseudoEcgVec_ptr) / ecgDistance;

            for (int i = 0; i <= sz - 1; i++)
            {
                if ( (* (splitting->globalSolution().at (0) ) ).isGlobalIDPresent (i) )
                {
                    pseudoEcgReal += pseudoEcgVec[i];
                }
            }

            MPI_Allreduce (&pseudoEcgReal, &Global_pseudoEcgReal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // rapporte a une variable connue de tous les procs

            if ( t >= TCut1 && t <= TCut2)
            {
                function_Type g = &Cut;
                vectorPtr_Type M_Cut (new VectorEpetra ( splitting->feSpacePtr()->map() ) );
                const feSpacePtr_Type feSpace =  splitting->feSpacePtr();
                feSpacePtr_Type* feSpace_noconst = const_cast< feSpacePtr_Type* > (&feSpace);
                (*feSpace_noconst)->interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( g ), *M_Cut , 0);
                * (splitting->globalSolution().at (0) ) = * (splitting->globalSolution().at (0) ) * (*M_Cut);
            }

            if (  Comm->MyPID() == 0 )
            {
                output << Global_pseudoEcgReal << "\n";
            }
        }
        exporterSplitting.closeFile();
        output.close();

        //********************************************//
        // Saving Fiber direction to file             //
        //********************************************//
        splitting -> exportFiberDirection();

        chronoinitialsettings.stop();
        std::cout << "\n\n\nElapsed time : " << chronoinitialsettings.diff() << std::endl;

        if ( Comm->MyPID() == 0 )
        {
            cout << "\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
        }

        MPI_Barrier (MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
