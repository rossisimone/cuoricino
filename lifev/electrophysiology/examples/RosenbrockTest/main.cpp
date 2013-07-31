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
    @brief test of Rosenbrock performance.

    @date 27-04-2013
    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @contributor
    @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
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
#include <lifev/electrophysiology/solver/IonicModels/IonicJafriRiceWinslow.hpp>
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


Real Stimulus1 (const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    return 80.0 * ( 0.5 + 0.5 * ( std::tanh ( - 300 * ( ( x - 0.4 ) * ( x - 0.6 ) + ( y - 0.4 ) * ( y - 0.6 ) ) ) ) );
}

Real Stimulus2 (const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    if ( x<= 0.1 )
    	return 80.0;
    else if( x<= 0.2)
    	return 80.0*( 0.2 - x )/(0.1);
    else
    	return 0.0;
}

Real Stimulus3 (const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    if ( y<= 0.1 )
    	return 80.0;
    else if( y<= 0.2)
    	return 80.0*( 0.2 - y )/(0.1);
    else
    	return 0.0;
}

Real Cut (const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    if ( y<= 4.5 )
    	return 1.0;
    else if( y<= 5.5)
    	return ( 0.55 - y/10.0 )/(0.1);
    else
    	return 0.0;
}

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

    if ( Comm->MyPID() == 0 )
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

    typedef ElectroETAMonodomainSolver< mesh_Type, IonicJafriRiceWinslow > monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList JRWParameterList = * ( Teuchos::getParametersFromXmlFile ( "JafriRiceWinslowParameters.xml" ) );
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
        std::cout << "Building Constructor for JafriRiceWinslow model with parameters ... ";
    }
    boost::shared_ptr<IonicJafriRiceWinslow>  model ( new IonicJafriRiceWinslow (JRWParameterList) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    model->showMe();


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
    function_Type f = &Stimulus2;
    splitting -> setPotentialFromFunction ( f ); //initialize potential

    //setting up initial conditions
    //* ( splitting -> globalSolution().at (0) ) = - 84.1638;
    * ( splitting -> globalSolution().at (1) ) = 0.0328302;
    * ( splitting -> globalSolution().at (2) ) = 0.988354;
    * ( splitting -> globalSolution().at (3) ) = 0.992540;
    * ( splitting -> globalSolution().at (4) ) = 0.000928836;
    * ( splitting -> globalSolution().at (5) ) = 10.2042;
    * ( splitting -> globalSolution().at (6) ) = 143.727;
    * ( splitting -> globalSolution().at (7) ) = 5.4;
    * ( splitting -> globalSolution().at (8) ) = 9.94893e-11;
    * ( splitting -> globalSolution().at (9) ) = 1.243891;
    * ( splitting -> globalSolution().at (10) ) = 1.36058e-4;
    * ( splitting -> globalSolution().at (11) ) = 1.17504;
    * ( splitting -> globalSolution().at (12) ) = 0.762527;
    * ( splitting -> globalSolution().at (13) ) = 1.19168e-3;
    * ( splitting -> globalSolution().at (14) ) = 6.30613e-9;
    * ( splitting -> globalSolution().at (15) ) = 0.236283;
    * ( splitting -> globalSolution().at (16) ) = 0.997208;
    * ( splitting -> globalSolution().at (17) ) = 6.38897e-5;
    * ( splitting -> globalSolution().at (18) ) = 1.535e-9;
    * ( splitting -> globalSolution().at (19) ) = 1.63909e-14;
    * ( splitting -> globalSolution().at (20) ) = 6.56337e-20;
    * ( splitting -> globalSolution().at (21) ) = 9.84546e-21;
    * ( splitting -> globalSolution().at (22) ) = 2.72826e-3;
    * ( splitting -> globalSolution().at (23) ) = 6.99215e-7;
    * ( splitting -> globalSolution().at (24) ) = 6.71989e-11;
    * ( splitting -> globalSolution().at (25) ) = 2.87031e-15;
    * ( splitting -> globalSolution().at (26) ) = 4.59752e-20;
    * ( splitting -> globalSolution().at (27) ) = 0.0;
    * ( splitting -> globalSolution().at (28) ) = 0.998983;
    * ( splitting -> globalSolution().at (29) ) = 0.00635;
    * ( splitting -> globalSolution().at (30) ) = 0.13598;

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

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

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    splitting -> setupLumpedMassMatrix();
    splitting -> setupStiffnessMatrix();
    splitting -> setupGlobalMatrix();

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;

    splitting -> setupExporter ( exporterSplitting, monodomainList.get ("OutputFile", "Splitting") );

    splitting -> exportSolution ( exporterSplitting, 0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nstart solving: \n\n " ;
    }

    Real TF     = JRWParameterList.get( "endTime", 5.0 );
    Real dt     = JRWParameterList.get( "timeStep", 5.77e-5 );
    Real TCut1 = JRWParameterList.get ("TCut", 35.0) - dt/2.0;
    Real TCut2 = JRWParameterList.get ("TCut", 35.0) + dt/2.0;
    Int iter = JRWParameterList.get( "savedt", 1.0) / dt;
    Int meth = JRWParameterList.get ("meth", 1.0);
    Real dt_min = 1e-10;
    Int k(0),j(0);
    Int nodes;

    Real timeReac = 0.0;
    Real timeDiff = 0.0;
    LifeChrono chrono;

    if (meth <= 1.0)
    {

        monodomainSolver_Type::vectorPtr_Type dtVec ( new VectorEpetra ( splitting->feSpacePtr() -> map(), LifeV::Unique ) );
        ExporterHDF5<mesh_Type> Exp;
        Exp.setMeshProcId ( splitting -> localMeshPtr(), splitting -> commPtr() -> MyPID() );
        Exp.setPrefix (monodomainList.get ("OutputTimeSteps", "TimeSteps"));
        Exp.addVariable ( ExporterData<mesh_Type>::ScalarField,  "dt", splitting->feSpacePtr(), dtVec, UInt (0) );

		//splitting   -> solveSplitting ( exporterSplitting );

        cout<<"Starting for...\n";
		for ( Real t = 0.0; t < TF-1e-8; )
		{
			t = t + dt;

			cout<<"Done!\nStarting reaction...";
			chrono.start();
			if(meth==1)
				splitting->solveOneReactionStepROS3P(dtVec, dt_min);
				//splitting->solveOneReactionStepROS3P();
			else
				splitting->solveOneReactionStepFE();
			chrono.stop();
			timeReac += chrono.diff();

			(*splitting->rhsPtrUnique()) *= 0.0;
			splitting->updateRhs();

			cout<<"Done!\nStarting diffusion...";
			chrono.start();
			splitting->solveOneDiffusionStepBE();
			chrono.stop();
			timeDiff += chrono.diff();

			if( k % iter == 0 )
			{
				splitting -> exportSolution (exporterSplitting, t);
				Exp.postProcess (t);
			}

			nodes = dtVec->epetraVector().MyLength();
			j = dtVec->blockMap().GID (0);
			dt_min = (*dtVec)[j];
			for(int i=1; i<nodes; i++)
			{
				j = dtVec->blockMap().GID (i);
				if(dt_min>(*dtVec)[j])
					dt_min = (*dtVec)[j];
			}

			k++;

			if( t >= TCut1 && t<=TCut2)
			{
				function_Type g = &Cut;
				vectorPtr_Type M_Cut(new VectorEpetra( splitting->feSpacePtr()->map() ));
				const feSpacePtr_Type feSpace =  splitting->feSpacePtr();
				feSpacePtr_Type* feSpace_noconst = const_cast< feSpacePtr_Type* >(&feSpace);
				(*feSpace_noconst)->interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( g ), *M_Cut , 0);
				*(splitting->globalSolution().at(0)) = *(splitting->globalSolution().at(0))*(*M_Cut);
				//*(splitting->globalSolution().at(1)) = *(splitting->globalSolution().at(1))*(*M_Cut);
			}

			if ( Comm->MyPID() == 0 )
				std::cout<<"\n\n\nActual time : "<<t<<std::endl<<std::endl<<std::endl;

		}

		Exp.closeFile();

    }
    else
    {
    	for ( Real t = 0.0; t < TF; )
    	{
			t = t + dt;

			chrono.start();
			splitting -> solveOneStepGatingVariablesFE();
			chrono.stop();
			timeReac += chrono.diff();

			chrono.start();
			if(meth==2)
				splitting -> solveOneICIStep ();
			else
				splitting -> solveOneSVIStep ();
			chrono.stop();
			timeDiff += chrono.diff();

			if( k % iter == 0 )
				splitting -> exportSolution (exporterSplitting, t);

			k++;

			if( t >= TCut1 && t<=TCut2)
			{
				function_Type g = &Cut;
				vectorPtr_Type M_Cut(new VectorEpetra( splitting->feSpacePtr()->map() ));
				const feSpacePtr_Type feSpace =  splitting->feSpacePtr();
				feSpacePtr_Type* feSpace_noconst = const_cast< feSpacePtr_Type* >(&feSpace);
				(*feSpace_noconst)->interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( g ), *M_Cut , 0);
				*(splitting->globalSolution().at(0)) = *(splitting->globalSolution().at(0))*(*M_Cut);
				//*(splitting->globalSolution().at(1)) = *(splitting->globalSolution().at(1))*(*M_Cut);
			}

			if ( Comm->MyPID() == 0 )
				std::cout<<"\n\n\nActual time : "<<t<<std::endl<<std::endl<<std::endl;
    	}
    }

    exporterSplitting.closeFile();



    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    splitting -> exportFiberDirection();

    if ( Comm->MyPID() == 0 )
    {
    	chronoinitialsettings.stop();
    	std::cout << "\n\n\nTotal elapsed time : " << chronoinitialsettings.diff() << std::endl;
    	std::cout<<"Diffusion/ICI/SVI time : "<<timeDiff<<std::endl;
    	std::cout<<"Reaction/GatingVar time : "<<timeReac<<std::endl;

        cout << "\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
