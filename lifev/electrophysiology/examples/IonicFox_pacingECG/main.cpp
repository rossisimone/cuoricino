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
    @brief 3D test with the Fox model 2002.

    @date 04−2013
    @author Marie Dupraz <dupraz.marie@gmail.com>

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

#include <lifev/electrophysiology/solver/IonicModels/IonicFox.hpp>
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

Real InitialPotential (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -94.7;
}

Real PacingProtocol ( const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID&   /*id*/)
{
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    Real pacingSite_X = monodomainList.get ("pacingSite_X", 0.);
    Real pacingSite_Y = monodomainList.get ("pacingSite_Y", 0.);
    Real stimulusRadius = monodomainList.get ("stimulusRadius", 0.1);
    Real stimulusValue = monodomainList.get ("stimulusValue", 80.);

    Real returnValue1;

    // --- Pacing protocol parameters -----------------------------------

    std::vector<double> returnPeriods;
    std::vector<double> returnStimulusTime;

    if ( ( ( ( x - pacingSite_X ) * ( x - pacingSite_X ) +  ( y - pacingSite_Y ) * ( y - pacingSite_Y )  )
        <= ( stimulusRadius * stimulusRadius ) ) ){
        returnValue1 = stimulusValue;
    }else{
        returnValue1 = 0.;
    }

    Real returnValue = returnValue1;
    return returnValue;
}

//Real Stimulus2 (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
//{
//    if ( sqrt( x*x + y*y ) <= 0.1 )
//    	return -14.7;
//    else
//    	return -94.7;
//}

Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
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

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >    feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>    feSpacePtr_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) > function_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, IonicFox > monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
    typedef VectorEpetra                        vector_Type;
    typedef MatrixEpetra<Real>                  matrix_Type;
    typedef LifeV::Preconditioner               basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
    typedef LifeV::PreconditionerML             prec_Type;
    typedef boost::shared_ptr<prec_Type>        precPtr_Type;

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList FoxParameterList = * ( Teuchos::getParametersFromXmlFile ( "FoxParameters.xml" ) );
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
        std::cout << "Building Constructor for Fox Model with default parameters ... ";
    }
    boost::shared_ptr<IonicFox>  model ( new IonicFox ( ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = monodomainList.get ("mesh_path", "./");

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

    monodomainSolverPtr_Type splitting ( new monodomainSolver_Type ( meshName, meshPath, dataFile, model ) );
    const feSpacePtr_Type FESpacePtr =  splitting->feSpacePtr(); //FE Space
//    monodomainSolverPtr_Type ICI ( new monodomainSolver_Type ( dataFile, model, splitting ->  localMeshPtr() ) );

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

    //Compute the potential at t0
    function_Type pacing = &PacingProtocol;
    function_Type f = &InitialPotential;
//    function_Type f = &Stimulus2;

    splitting -> setPotentialFromFunction ( f ); //initialize potential

    // APD calculation variables
    Int sz = 0;
    sz = (*(splitting->globalSolution().at(0))).size();
    Real threshold = monodomainList.get ("threshold", 10.);
    Real trep = 0.;
    vector_Type tact = (*(splitting->globalSolution().at(0)));
    tact*=0;
    vector_Type apd = tact;
    vector_Type delta_apd = tact;
    vector_Type previouspotential = (*(splitting->globalSolution().at(0)));

    //Compute the potential at t0
//    function_Type f = &Stimulus2;
//    splitting -> setPotentialFromFunction ( f ); //initialize potential
//    ICI -> setPotentialFromFunction ( f ); //initialize potential


    //setting up initial conditions
//    * ( splitting -> globalSolution().at (0) ) = FoxParameterList.get ("InitialPotential", -94.7);
    * ( splitting -> globalSolution().at (1) ) = FoxParameterList.get ("gatingM",2.4676e-4);
    * ( splitting -> globalSolution().at (2) ) = FoxParameterList.get ("gatingH", 0.99869);
    * ( splitting -> globalSolution().at (3) ) = FoxParameterList.get ("gatingJ", 0.99887);
    * ( splitting -> globalSolution().at (4) ) = FoxParameterList.get ("gatingXKR", 0.229);
    * ( splitting -> globalSolution().at (5) ) = FoxParameterList.get ("gatingXKS", 0.0001);
    * ( splitting -> globalSolution().at (6) ) = FoxParameterList.get ("gatingXT0", 3.742e-5);
    * ( splitting -> globalSolution().at (7) ) = FoxParameterList.get ("gatingYT0", 1.0);
    * ( splitting -> globalSolution().at (8) ) = FoxParameterList.get ("gatingF", 0.983);
    * ( splitting -> globalSolution().at (9) ) = FoxParameterList.get ("gatingD", 0.0001);
    * ( splitting -> globalSolution().at (10) ) = FoxParameterList.get ("gatingFCA", 0.942);
    * ( splitting -> globalSolution().at (11) ) = FoxParameterList.get ("concentrationCAIN", 0.0472);
    * ( splitting -> globalSolution().at (12) ) = FoxParameterList.get ("concentrationCASR", 320.0);

//    * ( ICI -> globalSolution().at (1) ) = FoxParameterList.get ("gatingM",2.4676e-4);
//    * ( ICI -> globalSolution().at (2) ) = FoxParameterList.get ("gatingH", 0.99869);
//    * ( ICI -> globalSolution().at (3) ) = FoxParameterList.get ("gatingJ", 0.99887);
//    * ( ICI -> globalSolution().at (4) ) = FoxParameterList.get ("gatingXKR", 0.229);
//    * ( ICI -> globalSolution().at (5) ) = FoxParameterList.get ("gatingXKS", 0.0001);
//    * ( ICI -> globalSolution().at (6) ) = FoxParameterList.get ("gatingXT0", 3.742e-5);
//    * ( ICI -> globalSolution().at (7) ) = FoxParameterList.get ("gatingYT0", 1.0);
//    * ( ICI -> globalSolution().at (8) ) = FoxParameterList.get ("gatingF", 0.983);
//    * ( ICI -> globalSolution().at (9) ) = FoxParameterList.get ("gatingD", 0.0001);
//    * ( ICI -> globalSolution().at (10) ) = FoxParameterList.get ("gatingFCA", 0.942);
//    * ( ICI -> globalSolution().at (11) ) = FoxParameterList.get ("concentrationCAIN", 0.0472);
//    * ( ICI -> globalSolution().at (12) ) = FoxParameterList.get ("concentrationCASR", 320.0);

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the pacing protocol             //
    //********************************************//
    Real pacingPeriod = monodomainList.get ("pacingPeriod", 500.);
    Real pacingPeriodMin = monodomainList.get ("pacingPeriodMin", 400.);
    Real pacingDelta = monodomainList.get ("pacingDelta", 0.);
    Real stimulusStart = monodomainList.get ("stimulusStart", 0.);
    Real stimulusStop = monodomainList.get ("stimulusStop", 0.05);
    Real stimulusNumber = monodomainList.get ("stimulusNumber", 1);

    std::vector<double> returnPeriods;
    std::vector<double> returnStimulusTime;
    int i(0);
    int control(0);

    Real NumberPacingPeriods ( (pacingPeriod-pacingPeriodMin)/pacingDelta );
    //--- Pacing method
    if ( pacingDelta >0 ){    // IF pacing
        for(int k=0; k <= NumberPacingPeriods-1; k++ ){
            for(i = stimulusNumber * k; i <= stimulusNumber * (k+1)-1; i++){
                returnPeriods.push_back ( pacingPeriod - k * pacingDelta );
                if ( i==0 ){
                    returnStimulusTime.push_back( returnPeriods[0] );
                }else{
                    returnStimulusTime.push_back ( returnStimulusTime[i-1]+returnPeriods[i] );
                }
            }
        }
    }
    //*******************************************//
    // Setting up the pseudo-ECG                 //
    //*******************************************//
    Real ecg_position_X = monodomainList.get ("ecg_position_X", 1.);
    Real ecg_position_Y = monodomainList.get ("ecg_position_Y", 1.);
    Real ecg_position_Z = monodomainList.get ("ecg_position_Z", 0.5);
    vector_Type ecgDistance = (*(splitting->globalSolution().at(0))) ;
    ecgDistance *= 0.;
    FESpacePtr->interpolate ( static_cast<function_Type> ( Norm::f ), ecgDistance, 0.0 );
    Norm::setPosition ( ecg_position_X , ecg_position_Y , ecg_position_Z ); // Set electrode position

    Real pseudoEcgReal(0.);
    Real Global_pseudoEcgReal(0.);
    vector_Type solutionLaplacian = ( (*(splitting->globalSolution().at(0))) );
    vector_Type pseudoEcgVec = ( (*(splitting->globalSolution().at(0))) );

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
//    ICI -> setParameters ( monodomainList );

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    VectorSmall<3> fibers;
    fibers[0] =  monodomainList.get ("fiber_X", std::sqrt (2.0) / 2.0 );
    fibers[1] =  monodomainList.get ("fiber_Y", std::sqrt (2.0) / 2.0 );
    fibers[2] =  monodomainList.get ("fiber_Z", 0.0 );

    splitting ->setupFibers (fibers);
//    ICI->setupFibers(fibers);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    splitting -> setupLumpedMassMatrix();
    splitting -> setupStiffnessMatrix();
    splitting -> setupGlobalMatrix();
//
//    ICI -> setupLumpedMassMatrix();
//    ICI -> setupStiffnessMatrix();
//    ICI -> setupGlobalMatrix();
//
//    monodomainSolverPtr_Type SVI ( new monodomainSolver_Type ( *ICI ) );

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;
    string filenameSplitting =  monodomainList.get ("OutputFile", "Fox" );
    filenameSplitting += "Splitting";
    splitting -> setupExporter ( exporterSplitting, filenameSplitting );
//    splitting -> exportSolution ( exporterSplitting, 0);

    vectorPtr_Type APDptr ( new vector_Type (apd, Repeated ) );
    exporterSplitting.addVariable ( ExporterData<mesh_Type>::ScalarField,  "apd", FESpacePtr,
                            APDptr, UInt (0) );

    vectorPtr_Type DELTA_APDptr ( new vector_Type (delta_apd, Repeated ) );
    exporterSplitting.addVariable ( ExporterData<mesh_Type>::ScalarField,  "delta_apd", FESpacePtr,
                            DELTA_APDptr, UInt (0) );

    *APDptr = apd;
    *DELTA_APDptr = delta_apd;

    std::ofstream output  ("ecg_output.txt");

//    ExporterHDF5< RegionMesh <LinearTetra> > exporterICI;
//    ExporterHDF5< RegionMesh <LinearTetra> > exporterSVI;

//    string filenameICI =  monodomainList.get ("OutputFile", "Fox");
//    filenameICI += "ICI";
//    string filenameSVI =  monodomainList.get ("OutputFile", "Fox" );
//    filenameSVI += "SVI";

//    ICI -> setupExporter ( exporterICI, filenameICI );
//    SVI -> setupExporter ( exporterSVI, filenameSVI );

//    splitting -> exportSolution ( exporterSplitting, 0);
//    ICI -> exportSolution ( exporterICI, 0);
//    SVI -> exportSolution ( exporterSVI, 0);

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
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

//    Real dt = monodomainList.get ("timeStep", 0.1);
//    Real TF = monodomainList.get ("endTime", 150.0);
//    Real Savedt = monodomainList.get ("saveStep", 1.0);

//    splitting   -> solveSplitting ( exporterSplitting, Savedt );
//    exporterSplitting.closeFile();

//    ICI         -> solveICI ( exporterICI, Savedt );
//    exporterICI.closeFile();
//
//    SVI         -> solveSVI ( exporterSVI, Savedt );
//    exporterSVI.closeFile();
//
    Real dt = monodomainList.get ("timeStep", 0.1);
   Real deltat = monodomainList.get ("SaveStep", 1.0);
   Real TF = monodomainList.get ("endTime", 150.0);

   int iter((deltat / dt));
   int k(0);

   for ( Real t = 0.0; t < TF; ){

       // APD calculation
       previouspotential = (*(splitting->globalSolution().at(0)));
       //--------------------------------------
       // ECG initialization
       pseudoEcgReal = 0.;
       pseudoEcgVec_ptr.reset ( new vector_Type ( FESpacePtr->map(), Unique ) );
       rhs_Laplacian.reset( new vector_Type ( FESpacePtr->map(), Unique ) );
       //-----------------------------------------------------------------

       control = 0;
       t = t + dt;

       if( pacingDelta > 0 ){
           for(i = 0; i<= NumberPacingPeriods * stimulusNumber - 1; i++){
               if ( control < 1 ){
                   if ( ( t >=stimulusStart &&   t <=  stimulusStop + dt) ||
                       (t >= (returnStimulusTime[i] + stimulusStart) && t <= (returnStimulusTime[i] + stimulusStop + dt)) ){
                       splitting -> setAppliedCurrentFromFunction ( pacing );
                       control = control + 1;
                   }else{
                       splitting -> initializeAppliedCurrent();
                   }
               }
           }
       }

       k++;
       if (k % iter == 0)
           splitting -> solveOneSplittingStep(exporterSplitting, t);
       else
           splitting -> solveOneSplittingStep();
//      splitting -> solveOneSplittingStep (exporterSplitting, t);

       // ECG : discrete laplacian of the solution
       (*rhs_Laplacian) = (*systemMatrixL) * (*(splitting->globalSolution().at(0)));

       linearSolver2.setOperator ( systemMatrixM );
       linearSolver2.setRightHandSide ( rhs_Laplacian );
       linearSolver2.solve ( pseudoEcgVec_ptr );

       pseudoEcgVec = (*pseudoEcgVec_ptr)/ecgDistance;

//      // APD calculation
       for (int i = 0; i <= sz-1; i++){
           if( (*(splitting->globalSolution().at(0))).isGlobalIDPresent(i) ){

               if ( ( previouspotential[i] < threshold ) && ( (*(splitting->globalSolution().at(0)))[i] >= threshold ) ){
                   tact[i] = t -((-threshold + previouspotential[i])/((*(splitting->globalSolution().at(0)))[i] - previouspotential[i]) )*dt;
               }else if ( ( previouspotential[i] >= threshold ) && ( (*(splitting->globalSolution().at(0)))[i] < threshold ) ){
                   trep = t -((-threshold + previouspotential[i])/( (*(splitting->globalSolution().at(0)))[i] - previouspotential[i]) )*dt;
                   delta_apd[i] = (trep - tact[i]) - apd[i];
                   apd[i] = trep - tact[i];
               }

   //          // Pseudo-ECG summation
               pseudoEcgReal += pseudoEcgVec[i];
           }
       }
       *APDptr = apd;
       *DELTA_APDptr = delta_apd;

       MPI_Allreduce(&pseudoEcgReal,&Global_pseudoEcgReal,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); // rapporte a une variable connue de tous les procs

       if (  Comm->MyPID() == 0 ){
           output << Global_pseudoEcgReal << "\n";
       }
   }

//
    exporterSplitting.closeFile();
//
//
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
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
