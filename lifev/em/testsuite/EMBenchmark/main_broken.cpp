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
    @brief Electrophysiology benchmark of Niederer et al. 2011

    @date 01−2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
//#ifdef EPETRA_MPI
//
//#else
//#include <Epetra_Serialcomm.h>
//#endif

#include <mpi.h>
#include <Epetra_Mpicomm.h>
//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"



#include <fstream>
#include <string>

#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/em/solver/EMSolver.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <sys/stat.h>

using namespace LifeV;

Real PacingProtocolMM ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 10;

    Real returnValue;

    if ( std::abs( x - pacingSite_X ) <= stimulusRadius
    		 &&
    	 std::abs( z - pacingSite_Z ) <= stimulusRadius
    	 	 &&
    	 std::abs( y - pacingSite_Y ) <= stimulusRadius
    	 	 &&
    	 t <= 2)
    {
    	returnValue = stimulusValue;
    }
    else{
        returnValue = 0.;
    }

    return returnValue;
}

Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_comm>  comm ( new Epetra_Mpicomm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    //********************************************//
    // Starts the chronometer.                    //
    //********************************************//

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

   	typedef ElectroIonicModel	                                     ionicModel_Type;
    typedef boost::shared_ptr<ionicModel_Type>                       ionicModelPtr_Type;
    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type > monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >               monodomainSolverPtr_Type;

    typedef VectorEpetra				                             vector_Type;
    typedef boost::shared_ptr<vector_Type>                           vectorPtr_Type;

    typedef EMSolver<mesh_Type, ionicModel_Type>	emSolver_Type;
    typedef boost::shared_ptr<emSolver_Type> emSolverPtr_Type;

    LifeChrono chronoinitialsettings;

    if ( comm->MyPID() == 0 )
      	chronoinitialsettings.start();

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

	//********************************************//
	// We need the GetPot datafile to setup       //
	//                                            //
	//********************************************//
	GetPot command_line(argc, argv);
	const string data_file_name = command_line.follow("data", 2, "-f",
			"--file");
	GetPot dataFile(data_file_name);


    //********************************************//
    // We create three solvers to solve with:     //
    // 1) Operator Splitting method               //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

	emSolverPtr_Type emSolverPtr( new emSolver_Type(monodomainList, data_file_name, comm));


    if ( comm->MyPID() == 0 )
    {
        std::cout << " solver done... ";
    }
    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nInitializing potential and gating variables:  " ;
    }


    function_Type stimulus;
   	stimulus = &PacingProtocolMM;
    emSolverPtr -> monodomainPtr() -> setAppliedCurrentFromFunction(stimulus, 0.0);

    if ( comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nSetting fibers:  " ;
    }

    VectorSmall<3> fibers;
    fibers[0]=0.0;
    fibers[1]=0.0;
    fibers[2]=1.0;
    emSolverPtr -> monodomainPtr() -> setupFibers( fibers );

    if ( comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    emSolverPtr -> setup(monodomainList, data_file_name, comm);

	emSolverPtr -> setupExporters(comm, problemFolder);


    //********************************************//
    // Activation time						      //
    //********************************************//
	emSolverPtr -> registerActivationTime(0.0, 0.8);
	emSolverPtr -> exportSolution(comm, 0.0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nstart solving:  " ;
    }

    Real dt = monodomainList.get ("timeStep", 0.1);
    Real TF = monodomainList.get ("endTime", 150.0);
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;
    Int k(0);

    Real timeReac = 0.0;
    Real timeDiff = 0.0;
    Real timeReacDiff = 0.0;
    LifeChrono chrono;

    std::string solutionMethod = monodomainList.get ("solutionMethod", "splitting");


    for ( Real t = 0.0; t < TF; )
    {
    	emSolverPtr -> monodomainPtr() -> setAppliedCurrentFromFunction ( stimulus, t );

		for(int j(0); j<reactionSubiter; j++)
			emSolverPtr -> monodomainPtr() -> solveOneReactionStepFE(reactionSubiter);
		//solve diffusion step
		emSolverPtr -> solveOneDiffusionStep();

        //register activation time
        k++;
        t = t + dt;



        if( k % iter == 0 )
        {
        	emSolverPtr -> registerActivationTime(t, 0.8);
        	emSolverPtr -> exportSolution(comm, t);
        }
        if ( comm->MyPID() == 0 )
        	std::cout<<"\n\n\nActual time : "<<t<<std::endl<<std::endl<<std::endl;
    }

	emSolverPtr -> exportActivationTime(comm, problemFolder);

	emSolverPtr -> closeExporters(comm);
    if ( comm->MyPID() == 0 )
          std::cout << "Exporting fibers: " << std::endl;

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    emSolverPtr -> monodomainPtr() -> exportFiberDirection(problemFolder);


//    if ( comm->MyPID() == 0 )
//    {
//    	chronoinitialsettings.stop();
//    	std::cout << "\n\n\nTotal lapsed time : " << chronoinitialsettings.diff() << std::endl;
//        if( solutionMethod == "splitting" )
//        {
//			std::cout<<"Diffusion time : "<<timeDiff<<std::endl;
//			std::cout<<"Reaction time : "<<timeReac<<std::endl;
//        }
//        else if( solutionMethod == "ICI" )
//        {
//        	std::cout<<"Solution time : "<<timeReacDiff<<std::endl;
//        }
//        else if( solutionMethod == "SVI" )
//        {
//        	std::cout<<"Solution time : "<<timeReacDiff<<std::endl;
//        }

        std::cout << "\n\nThank you for using EMSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
 //   }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();

    return ( EXIT_SUCCESS );
}
