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
    @3D monodomain model with TenTusscher dynamics.

    @date 05−2013
    @author Ricardo Ruiz <ricardo.ruiz@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscherFE.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

using namespace LifeV;

Real Stimulus2 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    if ( x<= 0.15 && y<=0.15 && z<=0.15 && t <= 2.0)
        return 30.0;
    else
        return -85.23;
}

Real PacingProtocol ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 2;

    Real returnValue1;

    if ( std::abs( x - pacingSite_X ) <= stimulusRadius && std::abs( z - pacingSite_Z ) <= stimulusRadius && std::abs( y - pacingSite_Y ) <= stimulusRadius){
        returnValue1 = stimulusValue;
    }else{
        returnValue1 = 0.;
    }

    return returnValue1;
}


Int main ( Int argc, char** argv )
{
  typedef RegionMesh<LinearTetra>                         mesh_Type;
  typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
  typedef boost::function < Real (const Real &   t,
                                  const Real &   x,
                                  const Real &   y,
                                  const Real &   z,
                                  const ID&   /*i*/ ) >   function_Type;
  typedef IonicTenTusscherFE					ionicModel_Type;
  typedef boost::shared_ptr< ionicModel_Type >  ionicModelPtr_Type;
  typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type >        monodomainSolver_Type;
  typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
  typedef VectorEpetra				vector_Type;
  typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
  typedef MatrixEpetra<Real> matrix_Type;
  typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;


#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }


    //********************************************//
    // Import parameters from an xml list. 
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    if ( comm->MyPID() == 0 )
      {
        std::cout << "Mesh Loading...";
      }

    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = monodomainList.get ("mesh_path", "./");

    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    if ( comm->MyPID() == 0 )
      {
	std::cout << " Done!" << endl;
      }



    if ( comm->MyPID() == 0 )
      {
        std::cout << "Building Constructor for the tenTusscher model with parameters ... ";
      }
    Teuchos::ParameterList parameterList  = * ( Teuchos::getParametersFromXmlFile ( "TenTusscherParameters.xml" ) );
    ionicModelPtr_Type  ionicModel ( new ionicModel_Type(  parameterList ) );
    if ( comm->MyPID() == 0 )
      {
        std::cout << " Done!" << endl;
      }
    // Initial pacing
    function_Type pacing = &PacingProtocol;

    //********************************************//
    // We solve with the Operator Splitting method               //
     //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Operator Splitting Solver for the Monodomain Problem... ";
    }

    monodomainSolverPtr_Type monodomain ( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicModel ) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }




    //********************************************//
    // initial conditions
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        cout << "\nInitializing potential:  " ;
    }

    //monodomain -> setInitialConditions();
    function_Type Stim = &Stimulus2;
    monodomain -> setPotentialFromFunction( Stim );
    //HeartUtility::setValueOnBoundary( *(monodomain -> potentialPtr() ), monodomain -> fullMeshPtr(), 1.0, 30 );

    //setting up initial conditions for TenTusscher
      * ( monodomain -> globalSolution().at (1) ) = 0.00172;
      * ( monodomain -> globalSolution().at (2) ) = 0.7045;
      * ( monodomain -> globalSolution().at (3) ) = 0.7444;
      * ( monodomain -> globalSolution().at (4) ) = 0.00621;
      * ( monodomain -> globalSolution().at (5) ) = 0.4712;
      * ( monodomain -> globalSolution().at (6) ) = 0.0095;
      * ( monodomain -> globalSolution().at (7) ) = 0.00003373;
      * ( monodomain -> globalSolution().at (8) ) = 0.7888;
      * ( monodomain -> globalSolution().at (9) ) = 0.9953; //?
      * ( monodomain -> globalSolution().at (10) ) = 2.42e-8;
      * ( monodomain -> globalSolution().at (11) ) = 0.999998;
      * ( monodomain -> globalSolution().at (12) ) = 0.9755; //??
      * ( monodomain -> globalSolution().at (13) ) = 0.000126; // Intracellular Ca
      * ( monodomain -> globalSolution().at (14) ) = 3.64;     // Sarcoplasmic Ca
      * ( monodomain -> globalSolution().at (15) ) = 8.604;    // Intracellular Na
      * ( monodomain -> globalSolution().at (16) ) = 136.89;   // Intracellular Ka

    if ( comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the data                   //
    //********************************************//
    monodomain -> setParameters ( monodomainList );
 
    //********************************************//
    // Create a fiber direction                   //
    //********************************************//

    VectorSmall<3> fibers;
    fibers[0] =  monodomainList.get ("fiber_X", 0.0 );
    fibers[1] =  monodomainList.get ("fiber_Y", 0.0 );
    fibers[2] =  monodomainList.get ("fiber_Z", 1.0 );

    monodomain ->setupFibers (fibers);

    //boost::shared_ptr<VectorEpetra> fiber ( new VectorEpetra ( Space3D -> map() ) );
    //HeartUtility::setupFibers ( *fiber, 0.0, 0.0, 1.0 );
    //splitting -> setFiberPtr(fiber);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    monodomain -> setupMassMatrix();
    monodomain -> setupStiffnessMatrix();
    monodomain -> setupGlobalMatrix();


    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterMonodomain;

    monodomain -> setupExporter ( exporterMonodomain, monodomainList.get ("OutputFile", "TenTussch" ));
    monodomain -> exportSolution ( exporterMonodomain, 0.0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

    Real Savedt = monodomainList.get ("saveStep", 0.1);
    monodomain   -> solveSplitting ( exporterMonodomain, Savedt );

    exporterMonodomain.closeFile();

    if ( comm->MyPID() == 0 )
    {
        cout << "\nT ======Finished====== \n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
