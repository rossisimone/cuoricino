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

Real Stimulus2 (const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    if ( x<= 0.05 )
        return 1.0;
    else if( x<= 0.1)
        return 1.0*( 0.1 - x )/(0.05);
    else
        return 0.0;
}


Int main ( Int argc, char** argv )
{
  typedef RegionMesh<LinearTetra>                         mesh_Type;
  typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
  typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
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

     
    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = monodomainList.get ("mesh_path", "./");
    meshPtr_Type mesh ( new mesh_Type ( comm ) );

    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);



    if ( comm->MyPID() == 0 )
      {
        std::cout << "Building Constructor for the tenTusscher model with parameters ... ";
      }
    ionicModelPtr_Type  ionicModel ( new ionicModel_Type() );
    if ( comm->MyPID() == 0 )
      {
        std::cout << " Done!" << endl;
      }

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

    function_Type f = &Stimulus2;
    monodomain -> setPotentialFromFunction ( f );
  
    //setting up initial conditions
    * ( monodomain -> globalSolution().at (1) ) = 1.0;
    * ( monodomain -> globalSolution().at (2) ) = 1.0;
    * ( monodomain -> globalSolution().at (3) ) = 0.021553043080281;
    * ( monodomain -> globalSolution().at (4) ) = 0.0;

    //! Or simply 
    //monodomain -> setInitialConditions();
    //monodomain -> setPotentialFromFunction( Vlid );
    //HeartUtility::setValueOnBoundary( *(monodomain -> potentialPtr() ), monodomain -> fullMeshPtr(), 1.0, 30 );


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
    fibers[0] =  monodomainList.get ("fiber_X", std::sqrt (2) / 2.0 );
    fibers[1] =  monodomainList.get ("fiber_Y", std::sqrt (2) / 2.0 );
    fibers[2] =  monodomainList.get ("fiber_Z", 0.0 );

    monodomain ->setupFibers (fibers);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    monodomain -> setupLumpedMassMatrix();
    monodomain -> setupStiffnessMatrix();
    monodomain -> setupGlobalMatrix();


    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;

    string filenameSplitting =  monodomainList.get ("OutputFile", "TenTusscher" );
    filenameSplitting += "Monodomain";
    monodomain -> setupExporter ( exporterSplitting, filenameSplitting );
    monodomain -> exportSolution ( exporterSplitting, 0.0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

    Real Savedt = monodomainList.get ("saveStep", 0.1);
    monodomain   -> solveSplitting ( exporterSplitting, Savedt );
    exporterSplitting.closeFile();

    if ( comm->MyPID() == 0 )
    {
        cout << "\nT ======Finished====== \n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
