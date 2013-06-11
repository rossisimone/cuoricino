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
    @brief 3D test with the Rogers-McCulloch model.

    @date 01−2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor Marie Dupraz <dupraz.marie@gmail.com>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>

    Use to run some test : choice of the fibers direction, right mesh size, ...
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
#include <lifev/core/mesh/MeshUtility.hpp>
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


Real Stimulus2 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    Real pacingSite_X = monodomainList.get ("pacingSite_X", 0.);
    Real pacingSite_Y = monodomainList.get ("pacingSite_Y", 0.);
    Real pacingSite_Z = monodomainList.get ("pacingSite_Z", 0. );
    Real stimulusRadius = 0.1; // monodomainList.get ("stimulusRadius", 0.1);

    if (  ( ( x - pacingSite_X ) * ( x - pacingSite_X ) +  ( y - pacingSite_Y ) * ( y - pacingSite_Y ) +  ( z - pacingSite_Z ) * ( z - pacingSite_Z )  )
                    <= ( stimulusRadius * stimulusRadius ) )
        return 80.0;
    else
        return 0.0;
}

//// Choice of the fibers direction : ||.||=1
//Real Fibers (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
//{
//    Real zmin = 0.0;
//    Real zmax = 0.5;
//    Real L = zmax - zmin;
//    Real thetamin = M_PI /3.;
//    Real thetamax = -M_PI /3.;
//
//    Real ztheta = (L-z)/L;
//    Real theta = (thetamax - thetamin) * ztheta + thetamin;
//
//        switch(i){
//            case 0:
//                return  std::cos(theta); // x_fib; //y/N; //std::sqrt (2.0) / 2.0; //
//                break;
//            case 1:
//                return  std::sin(theta); // y_fib; //-x/N; //std::sqrt (2.0) / 2.0; //
//                break;
//            case 2:
//                return 0.0;
//                break;
//            default:
//                return 0.0;
//                break;
//    }
//}


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
    LifeChrono chrono;

    chronoinitialsettings.start();

    typedef RegionMesh<LinearTetra>            mesh_Type;
    typedef boost::shared_ptr<mesh_Type>       meshPtr_Type;
    typedef VectorEpetra vector_Type;
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

    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > >    filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                                  hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                                              hdf5FilterPtr_Type;
    typedef ExporterData<mesh_Type>                           exporterData_Type;
    typedef Exporter< mesh_Type >                             IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >                  IOFilePtr_Type;
    typedef ExporterHDF5< mesh_Type >                         hdf5IOFile_Type;

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


    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if ( Comm->MyPID() == 0 )
    {
        std::cout << std::endl << "[Loading the mesh]" << std::endl;
    }

    meshPtr_Type fullMeshPtr ( new mesh_Type ( Comm ) );

    VectorSmall<3> meshDim;
    meshDim[0] =  monodomainList.get ("meshDim_X", 32 );
    meshDim[1] =  monodomainList.get ("meshDim_Y", 32 );
    meshDim[2] =  monodomainList.get ("meshDim_Z", 2 );
    VectorSmall<3> domain;
    domain[0] =  monodomainList.get ("domain_X", 1. );
    domain[1] =  monodomainList.get ("domain_Y", 1. );
    domain[2] =  monodomainList.get ("domain_Z", 0.0625 );

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
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

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

    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }

//    //***********************************************//
//    //              RESTART Protocol                 //
//    //***********************************************//
//
//    string filenameStart =  monodomainList.get ("StartFile", "FHNStartSplitting" );
//    std::string sol ( filenameStart ); // name of the file from which we want to restart
//
//    // Setting up the importer
//    filterPtr_Type importer ( new hdf5Filter_Type( ) );
//    importer->setMeshProcId ( splitting -> localMeshPtr(), Comm -> MyPID() );
//    importer->setPrefix ( sol );
//
//    // Import the value of the potential and gating variable
//    monodomainSolver_Type::vectorPtr_Type newSol;
//    newSol.reset(new VectorEpetra ( splitting -> feSpacePtr() -> map(), LifeV::Unique ));
//
//    exporterData_Type ImportDataPotential (exporterData_Type::ScalarField, "Variable0.00030", splitting -> feSpacePtr(),
//                                              newSol, UInt (0), exporterData_Type::UnsteadyRegime);
//
//    monodomainSolver_Type::vectorPtr_Type newGating;
//    newGating.reset(new VectorEpetra ( splitting -> feSpacePtr() -> map(), LifeV::Unique ));
//
//    exporterData_Type ImportDataGating (exporterData_Type::ScalarField, "Variable1.00030", splitting -> feSpacePtr(),
//                                                  newGating, UInt (0), exporterData_Type::UnsteadyRegime);
//
//    importer->readVariable ( ImportDataPotential );
//    importer->readVariable( ImportDataGating );
//
//    //********************************************//
//    // Setting up the initial condition form      //
//    // a given function.                          //
//    //********************************************//
//    if ( Comm->MyPID() == 0 )
//    {
//        cout << "\nInitializing potential:  " ;
//    }
//
//    splitting-> setPotentialPtr(newSol);
//    * ( splitting -> globalSolution().at (1) ) = *newGating;

    //****************************************************************************************//

    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nInitializing potential:  " ;
    }

    //Compute the potential at t0
    function_Type f = &Stimulus2;
    splitting -> setPotentialFromFunction ( f ); //initialize potential
//    ICI -> setPotentialFromFunction ( f ); //initialize potential

    //setting up initial conditions
    * ( splitting -> globalSolution().at (1) ) = FHNParameterList.get ("W0", 0.011);
//    * ( ICI -> globalSolution().at (1) ) = FHNParameterList.get ("W0", 0.011);

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the time data                   //
    //********************************************//
    splitting -> setParameters ( monodomainList );
//    ICI -> setParameters ( monodomainList );

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    VectorSmall<3> fibers;
    fibers[0] =  monodomainList.get ("fiber_X", std::sqrt (2.0) / 2.0 );
    fibers[1] =  monodomainList.get ("fiber_Y", std::sqrt (2.0) / 2.0 );
    fibers[2] =  monodomainList.get ("fiber_Z", 0.0 );

    splitting->setupFibers(fibers);
//    ICI->setupFibers(fibers);
//
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

//    chrono.start();
//    ICI -> setupLumpedMassMatrix();
//    ICI -> setupStiffnessMatrix();
//    ICI -> setupGlobalMatrix();
//    chrono.stop();
//    timematrixICI = chrono.globalDiff(*Comm);


    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;
//    ExporterHDF5< RegionMesh <LinearTetra> > exporterICI;
//    ExporterHDF5< RegionMesh <LinearTetra> > exporterSVI;
//
    string filenameSplitting =  monodomainList.get ("OutputFile", "FHN" );
    filenameSplitting += "Splitting";
//    string filenameICI =  monodomainList.get ("OutputFile", "FHN");
//    filenameICI += "ICI";
//    string filenameSVI =  monodomainList.get ("OutputFile", "FHN" );
//    filenameSVI += "SVI";
//
    splitting -> setupExporter ( exporterSplitting, filenameSplitting );
//    ICI -> setupExporter ( exporterICI, filenameICI );
//    SVI -> setupExporter ( exporterSVI, filenameSVI );
//
    splitting -> exportSolution ( exporterSplitting, 0);
//    ICI -> exportSolution ( exporterICI, 0);
//    SVI -> exportSolution ( exporterSVI, 0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

    Real Savedt = monodomainList.get ("saveStep", 0.1);

    splitting   -> solveSplitting (exporterSplitting, Savedt);
    exporterSplitting.closeFile();


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
//    ofstream output("norm.txt");
//
//    Real diff_norm (0.);
//    Real normS (0.);
//    Real normSNL (0.);
//
//    for ( Real t = 0.0; t < TF; ){
//        t = t + dt;
//        splitting   -> solveOneSplittingStep();
//        normS = (*(splitting->globalSolution().at(0))).norm2();
//        SplittingNoLumping -> solveOneSplittingStep();
//        normSNL = (*(SplittingNoLumping->globalSolution().at(0))).norm2();
//
//
//        diff_norm = (*(splitting->globalSolution().at(0))-*(SplittingNoLumping->globalSolution().at(0))).norm2();
//        output << normS << " " << normSNL << " " << diff_norm << " " << (diff_norm)/normS << std::endl;
//
//    }
//
//    if ( Comm->MyPID() == 0 )
//    {
//        ofstream output_time(monodomainList.get ("OutputFile", "Times").c_str());
//        output_time << "Time matrix splitting : " << timematrixsplitting << std::endl;
//        output_time << "Time matrix ICI : " << timematrixICI << std::endl;
//        output_time << "Time solving splitting lumped : " << timesolvesplitting << std::endl;
//        output_time << "Time solving splitting NL: " << timesolvesplittingNOL << std::endl;
//    }
//
//    vector_Type finalSolution (*(splitting->globalSolution().at(0)));
//    vector_Type GlobalfinalSolution (finalSolution);
//
//    std::cout << finalSolution.size() << std::endl;
//    std::cout << GlobalfinalSolution.size() << std::endl;
//
//    MPI_Allreduce(&finalSolution,&GlobalfinalSolution,3267,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); // rapporte a une variable connue de tous les procs
////
////    if ( Comm->MyPID() == 0 ){
////        for (int j=0; j<(*(splitting->globalSolution().at(0))).size(); j++){
////            output_solution << GlobalfinalSolution[j] << " " << GlobalfinalSolution[j] <<"\n";
////        }
////    }
//
//    ofstream output_solution("output_solution.txt");
//    for (int j=0; j<(*(splitting->globalSolution().at(0))).size(); j++){
//        if( (*(splitting->globalSolution().at(0))).isGlobalIDPresent(j) ){
//            output_solution << (*(splitting->globalSolution().at(0)))[j] << " "
//                        << (*(splitting->globalSolution().at(1)))[j] <<"\n";
//        }
//    }
