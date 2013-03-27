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
 *  @file
 *  @brief File containing the Integrated Heart example
 *
 *  @date 2012-09-25
 *  @author Toni Lassila <toni.lassila@epfl.ch>
 *  @maintainer Toni Lassila <toni.lassila@epfl.ch>
 *
*/

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <cassert>
#include <cstdlib>

#include <boost/timer.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV includes
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/fsi/solver/FSISolver.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/fsi/solver/FSIMonolithicGI.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

#include "ud_functions.hpp"
#include "boundaryConditions.hpp"
#include "flowConditions.hpp"
//#include "lumpedHeart.hpp"

class Problem
{
public:

    typedef boost::shared_ptr<LifeV::FSISolver> fsi_solver_ptr;

    typedef LifeV::FSIOperator::data_Type                          data_Type;
    typedef LifeV::FSIOperator::dataPtr_Type                       dataPtr_Type;

    typedef LifeV::FSIOperator::vector_Type        vector_Type;
    typedef LifeV::FSIOperator::vectorPtr_Type     vectorPtr_Type;

    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;

    typedef LifeV::ExporterEnsight<LifeV::FSIOperator::mesh_Type>  ensightFilter_Type;
    typedef boost::shared_ptr<ensightFilter_Type>                 ensightFilterPtr_Type;
#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<LifeV::FSIOperator::mesh_Type>  hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                  hdf5FilterPtr_Type;
#endif
    typedef LifeV::FactorySingleton<LifeV::Factory<LifeV::FSIOperator,  std::string> > FSIFactory_Type;

    Problem ( GetPot const& data_file ) :
        M_Tstart (0.),
        M_saveEvery (5)
        //M_returnValue(EXIT_FAILURE)

    {
        using namespace LifeV;

        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "linearVenantKirchhoff", &FSIOperator::createVenantKirchhoffLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "exponential", &FSIOperator::createExponentialMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "neoHookean", &FSIOperator::createNeoHookeanMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "nonLinearVenantKirchhoff", &FSIOperator::createVenantKirchhoffNonLinear );

        std::cout << "register MonolithicGE : " << FSIMonolithicGE::S_register << std::endl;
        std::cout << "register MonolithicGI : " << FSIMonolithicGI::S_register << std::endl;

        M_data = dataPtr_Type ( new data_Type() );
        M_data->setup ( data_file );

        /* Set physical parameters from getpot, see FlowConditions.cpp for details */
        FlowConditions::setParamsFromGetPot (data_file);

        M_fsi = fsi_solver_ptr ( new FSISolver( ) );
        MPI_Barrier ( MPI_COMM_WORLD );

        M_fsi->setData ( M_data );
        M_fsi->FSIOper()->setDataFile ( data_file );

        // Setting FESpace and DOF

        std::string  fluidMeshPartitioned    =  data_file ( "problem/fluidMeshPartitioned", "none" );
        std::string  solidMeshPartitioned    =  data_file ( "problem/solidMeshPartitioned", "none" );

#ifdef HAVE_HDF5
        if ( fluidMeshPartitioned.compare ( "none" ) )
        {
            FSIOperator::meshFilter_Type fluidMeshFilter ( data_file, fluidMeshPartitioned );
            fluidMeshFilter.setComm ( M_fsi->FSIOper()->worldComm() );
            FSIOperator::meshFilter_Type solidMeshFilter ( data_file, solidMeshPartitioned );
            solidMeshFilter.setComm ( M_fsi->FSIOper( )->worldComm( ) );
            M_fsi->FSIOper( )->partitionMeshes ( fluidMeshFilter, solidMeshFilter );
            M_fsi->FSIOper( )->setupFEspace( );
            M_fsi->FSIOper( )->setupDOF ( fluidMeshFilter );
            fluidMeshFilter.closeFile( );
            solidMeshFilter.closeFile( );
        }
        else
#endif
        {
            M_fsi->FSIOper( )->partitionMeshes( );
            M_fsi->FSIOper( )->setupFEspace( );
            M_fsi->FSIOper( )->setupDOF( );
        }

        Debug ( 10000 ) << "Setting up the FESpace and DOF \n";

        MPI_Barrier ( MPI_COMM_WORLD );

#ifdef DEBUG
        Debug ( 10000 ) << "Setting up the BC \n";
#endif

        /* Start at end-diastolic configuration with both valves closed */

        M_aorticValveIsOpen = false;
        M_mitralValveIsOpen = false;

        M_fsi->setFluidBC ( BCh_monolithicFlux ( M_aorticValveIsOpen, M_mitralValveIsOpen ) );
        M_fsi->setSolidBC ( BCh_monolithicSolid ( *M_fsi->FSIOper( ) ) );

        M_fsi->setup();

        M_fsi->setFluidBC ( BCh_monolithicFluid ( *M_fsi->FSIOper( ) ) );
        M_fsi->setHarmonicExtensionBC ( BCh_harmonicExtension ( *M_fsi->FSIOper( ) ) );

        dynamic_cast<LifeV::FSIMonolithic*> (M_fsi->FSIOper().get() )->mergeBCHandlers();

#ifdef DEBUG
        Debug ( 10000 ) << "BC set\n";
#endif

        std::string const exporterType =  data_file ( "exporter/type", "ensight" );
        std::string const fluidName    =  data_file ( "exporter/fluid/filename", "fluid" );
        std::string const solidName    =  data_file ( "exporter/solid/filename", "solid" );

#ifdef HAVE_HDF5
        if (exporterType.compare ("hdf5") == 0)
        {
            M_exporterFluid.reset ( new  hdf5Filter_Type ( data_file, fluidName) );
            M_exporterSolid.reset ( new  hdf5Filter_Type ( data_file, solidName) );
        }
        else
#endif
        {
            if (exporterType.compare ("none") == 0)
            {
                M_exporterFluid.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), fluidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID() ) );
                M_exporterSolid.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), solidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID() ) );
            }
            else
            {
                M_exporterFluid.reset ( new  ensightFilter_Type ( data_file, fluidName) );
                M_exporterSolid.reset ( new  ensightFilter_Type ( data_file, solidName) );
            }
        }


        // Load using Ensight/HDF5
        M_saveEvery = data_file ("exporter/saveEvery", 5);

        //        std::string restartType(data_file("importer/restartType","newSimulation"));
        M_fsi->initialize();

        M_velAndPressure.reset ( new vector_Type ( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluid->mapType() ) );
        M_fluidDisp.reset     ( new vector_Type ( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ) );
        M_solidDisp.reset ( new vector_Type ( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ) );


        M_exporterFluid->setMeshProcId (M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID() );
        M_exporterSolid->setMeshProcId (M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID() );

        M_exporterFluid->addVariable ( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-velocity",
                                       M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressure, UInt (0) );
        M_exporterFluid->addVariable ( ExporterData<FSIOperator::mesh_Type>::ScalarField, "f-pressure",
                                       M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressure,
                                       UInt (3 * M_fsi->FSIOper()->uFESpace().dof().numTotalDof() ) );
        M_exporterFluid->addVariable ( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-displacement",
                                       M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDisp, UInt (0) );

        M_exporterSolid->addVariable ( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-displacement",
                                       M_fsi->FSIOper()->dFESpacePtr(), M_solidDisp, UInt (0) );

        //M_fsi->FSIOper()->fluid().setupPostProc(); //this has to be called if we want to initialize the postProcess

        // Initialize the lumped parameter valves
        FC0.initParameters ( *M_fsi->FSIOper(), 3);

        M_data->dataFluid()->dataTime()->setInitialTime ( M_Tstart );
        M_data->dataFluid()->dataTime()->setTime ( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setInitialTime ( M_Tstart );
        M_data->dataSolid()->dataTime()->setTime ( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->timeDataALE()->setInitialTime ( M_Tstart );
        M_data->timeDataALE()->setTime ( M_data->dataFluid()->dataTime()->initialTime() );
    }

    /*!
      This routine runs the temporal loop
     */
    void
    run()
    {
        boost::timer _overall_timer;

        LifeV::UInt iter = 1;
        //        LifeV::UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

        dynamic_cast<LifeV::FSIMonolithic*> (M_fsi->FSIOper().get() )->enableStressComputation (1);

        vectorPtr_Type solution ( new vector_Type ( (*M_fsi->FSIOper()->couplingVariableMap() ) ) );

        M_fsi->FSIOper()->extrapolation ( *solution );

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(), M_data->dataSolid()->dataTime()->updateTime(), ++iter)
        {
            //    M_fsi->FSIOper()->displayer().leaderPrintMax( "FSI-  Mitral flow rate:  ", M_fsi->FSIOper()->fluid().flux(INLET, M_fsi->displacement()) );
            //    M_fsi->FSIOper()->displayer().leaderPrintMax( "FSI-  Aortic flow rate:  ", M_fsi->FSIOper()->fluid().flux(OUTLET, M_fsi->displacement()) );
            //    M_fsi->FSIOper()->displayer().leaderPrintMax( "FSI-  Aortic valve area: ", M_fsi->FSIOper()->fluid().area(OUTLET) );

            /* Renew the parameters in the lumped parameter valve model.The average pressure on the aortic valve
             * needs to be passed to the valve to set the correct outgoing flux. */

            FC0.renewLumpedParameters (OUTLET, M_fsi->FSIOper()->fluid().pressure (OUTLET, M_fsi->displacement() ) );

            /*
                    if ( M_aorticValveIsOpen && M_fsi->FSIOper()->fluid().flux(OUTLET, M_fsi->displacement()) < 0)
                    {
                        // Close the aortic valve
                        M_aorticValveIsOpen = false;
                        std::cout<<"Closing aortic valve"<<endl;
                        //M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->FSIOper(), M_aorticValveIsOpen, M_mitralValveIsOpen));
                        M_fsi->setFluidBC( LifeV::BCh_monolithicFlux( M_aorticValveIsOpen, M_mitralValveIsOpen ) );
                    }
                    else{
                        if ( (!M_aorticValveIsOpen) && ( M_fsi->FSIOper()->fluid().pressure(OUTLET, M_fsi->displacement())
                                        > LifeV::FlowConditions::outPressure(M_data->dataFluid()->dataTime()->time(),-1,-1,-1,3) ) )
                        {
                            // Open the aortic valve
                            M_aorticValveIsOpen=true;
                            std::cout << "Opening aortic valve" << endl;
                            //M_fsi->setFluidBC( BCh_monolithicFluid( *M_fsi->FSIOper() ) );
                            M_fsi->setFluidBC( LifeV::BCh_monolithicFlux( M_aorticValveIsOpen, M_mitralValveIsOpen ) );
                        }
                    }
            */
            FC0.renewParameters ( *M_fsi, 3 );

            boost::timer _timer;

            if (iter % M_saveEvery == 0)
            {
                M_fsi->FSIOper()->exportSolidDisplacement (*M_solidDisp);
                M_fsi->FSIOper()->exportFluidVelocityAndPressure (*M_velAndPressure);
                M_exporterSolid->postProcess ( M_data->dataFluid()->dataTime()->time() );

                *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();
                M_exporterFluid->postProcess ( M_data->dataFluid()->dataTime()->time() );
            }

            // This is just the previous solution. Should use the extrapolation from time advance
            M_fsi->FSIOper()->extrapolation ( *solution );

            M_fsi->iterate ( solution );

            if (M_data->method().compare ("monolithicGI") == 0)
            {
                // Saving the solution
                if ( M_data->dataFluid()->domainVelImplicit() )
                {
                    // shift_right of the solution of all the time advance classes in the FSIOperator
                    M_fsi->FSIOper()->updateSolution ( *solution );
                    // This resets M_uk to *solution, but indeed is superseeded in the call of evalResidual
                }
                else
                {
                    // This resets M_uk to *solution, but indeed is superseeded in the call of evalResidual
                    updateSolutionDomainVelocityFalse ( solution );
                }
            }
            else //Case when monolithicGE is used
            {
                M_fsi->FSIOper()->updateSolution ( *solution );
            }
            M_fsi->FSIOper()->displayer().leaderPrintMax ("[fsi_run] Iteration ", iter);
            M_fsi->FSIOper()->displayer().leaderPrintMax (" was done in : ", _timer.elapsed() );

        }

        std::cout << "Total computation time = "
                  << _overall_timer.elapsed() << "s" << "\n";



    }

private:

    //    void initializeStokes(  GetPot const& data_file);
    //    void restartFSI(  GetPot const& data_file);

    //    void checkCEResult(const LifeV::Real& time);
    //    void checkGCEResult(const LifeV::Real& time);

    void updateSolutionDomainVelocityFalse ( const vectorPtr_Type solution );

    fsi_solver_ptr M_fsi;
    dataPtr_Type   M_data;

    filterPtr_Type M_exporterSolid;
    filterPtr_Type M_exporterFluid;
    filterPtr_Type M_importerSolid;
    filterPtr_Type M_importerFluid;
    vectorPtr_Type M_velAndPressure;
    vectorPtr_Type M_fluidDisp;
    vectorPtr_Type M_solidDisp;

    std::vector<vectorPtr_Type> M_solidStencil;
    std::vector<vectorPtr_Type> M_fluidStencil;
    std::vector<vectorPtr_Type> M_ALEStencil;

    LifeV::FlowConditions FC0;
    //LifeV::LumpedHeart LH;

    LifeV::Real    M_Tstart;
    LifeV::UInt           M_saveEvery;

    bool M_aorticValveIsOpen;
    bool M_mitralValveIsOpen;

    /*struct RESULT_CHANGED_EXCEPTION
    {
    public:
        RESULT_CHANGED_EXCEPTION(LifeV::Real time)
        {
            std::cout<<"Some modifications led to changes in the l2 norm of the solution at time"<<time<<std::endl;
        }
    };*/
};

struct FSIChecker
{
    FSIChecker ( GetPot const& _data_file ) :
        data_file ( _data_file )
    {}

    int operator() ()
    {
        boost::shared_ptr<Problem> fsip;

        try
        {
            fsip = boost::shared_ptr<Problem> ( new Problem ( data_file ) );
            fsip->run();
        }
        catch ( std::exception const& _ex )
        {
            std::cout << "caught exception :  " << _ex.what() << "\n";
        }
        return EXIT_FAILURE;

    }

    GetPot                data_file;
    LifeV::Vector         disp;
};


namespace LifeV
{

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}
}

int main (int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#else
    std::cout << "% using serial Version" << std::endl;
#endif

    GetPot command_line (argc, argv);

    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot data_file (data_file_name);
    FSIChecker _sp_check ( data_file );
    int returnValue = _sp_check();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return returnValue;

}

void Problem::updateSolutionDomainVelocityFalse ( const vectorPtr_Type solution )
{

    //The solution for the fluid and solid is at the current time
    //The solution for the ALE is at the previous time

    //I build the vector I want to push in the exporters
    //In this case, the vector is ( f^n+1, s^n+1, df^n )
    LifeV::UInt nDofsALE = M_fsi->FSIOper()->mmFESpace().fieldDim() * M_fsi->FSIOper()->mmFESpace().dof().numTotalDof();

    //Extract the previous solution
    vector_Type previousSolution ( M_fsi->FSIOper()->solution() );
    vector_Type previousDisplacement ( M_fsi->FSIOper()->mmFESpace().map() );
    previousDisplacement *= 0.0;

    LifeV::UInt sizeOfSolutionVector = previousSolution.size();
    LifeV::UInt offsetStartCopying = sizeOfSolutionVector - nDofsALE;

    previousDisplacement.subset (previousSolution,  offsetStartCopying );

    //After having saved the previous displacement we can push the current solution
    M_fsi->FSIOper()->updateSolution ( *solution );

    M_fsi->FSIOper()->ALETimeAdvance()->shiftRight ( previousDisplacement );
    M_fsi->FSIOper()->ALETimeAdvance()->updateRHSFirstDerivative ( M_data->dataFluid()->dataTime()->timeStep() );

}
