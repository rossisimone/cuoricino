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
 *  @brief File containing the Multiscale Model FSI3D
 *
 *  @date 19-04-2010
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/models/MultiscaleModelFSI3DActivated.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"


namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelFSI3DActivated::MultiscaleModelFSI3DActivated() :
    super				           (),
    M_monodomain                   (),
    M_importerElectro              (),
    M_exporterElectro              ()
{
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelFSI3DActivated::setupData ( const std::string& fileName )
{
	super::setupData(fileName);

	//    std::string xmlpath = dataFile("electrophysiology/monodomain_xml_path","./");
	//    std::string xmlfile = dataFile("electrophysiology/monodomain_xml_file","MonodomainSolverParamList.xml");

		std::string xmlpath = "./"
		std::string xmlfile = 		"MonodomainSolverParamList.xml";
	    Teuchos::RCP< Teuchos::ParameterList > monodomainList = Teuchos::rcp ( new Teuchos::ParameterList );
	    monodomainList = Teuchos::getParametersFromXmlFile ( xmlpath + xmlfile );

	    ionicModelPtr_Type  ionicModel ( new ionicModel_Type() );

	    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
	    std::string meshPath = monodomainList.get ("mesh_path", "./");

	    monodomainSolverPtr_Type splitting ( new monodomainSolver_Type ( meshName, meshPath, dataFile, model ) );
}

void
MultiscaleModelFSI3DActivated::setupModel()
{

	super::setupModel();







}

void
MultiscaleModelFSI3DActivated::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3DActivated::buildModel() \n";
#endif

    // Display data
    //    if ( M_comm->MyPID() == 0 )
    //        M_data->showMe();

    M_FSIoperator->buildSystem();
    M_FSIoperator->updateSystem();

    // Update BCInterface solver variables
    updateBC();
}

void
MultiscaleModelFSI3DActivated::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3DActivated::updateModel() \n";
#endif

    // Update System
    M_FSIoperator->updateSystem();

    // Update BCInterface solver variables
    updateBC();

    // TODO This is a workaround of Paolo Crosetto to make it possible to perform subiterations
    // in the multiscale when using 3D FSI models. In the future this should be replaced with
    // a more proper implementation.
    M_FSIoperator->precPtrView()->setRecompute ( 1, true );
    M_nonLinearRichardsonIteration = 0;
}

void
MultiscaleModelFSI3DActivated::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3DActivated::solveModel() \n";
#endif

    displayModelStatus ( "Solve" );

    // TODO This is a workaround of Paolo Crosetto to make it possible to perform subiterations
    // in the multiscale when using 3D FSI models. In the future this should be replaced with
    // a more proper implementation.
    if ( M_nonLinearRichardsonIteration != 0 )
    {
        M_FSIoperator->resetRHS();
        M_FSIoperator->updateRHS();
        M_FSIoperator->applyBoundaryConditions();
    }

    // Non-linear Richardson solver
    UInt maxSubIterationNumber = M_data->maxSubIterationNumber();
    M_FSIoperator->extrapolation ( *M_stateVariable );

    NonLinearRichardson ( *M_stateVariable, *M_FSIoperator,
                          M_data->absoluteTolerance(), M_data->relativeTolerance(),
                          maxSubIterationNumber, M_data->errorTolerance(),
                          M_data->NonLinearLineSearch(),
                          M_nonLinearRichardsonIteration,
                          1
                        );

    // TODO This is a workaround of Paolo Crosetto to make it possible to perform subiterations
    // in the multiscale when using 3D FSI models. In the future this should be replaced with
    // a more proper implementation.
    M_FSIoperator->precPtrView()->setRecompute ( 1, false );
    M_nonLinearRichardsonIteration = 1;

#ifdef FSI_WITH_BOUNDARYAREA
    for ( UInt j ( 0 ); j < M_boundaryFlagsAreaPerturbed.size(); ++j )
    {
        M_boundaryFlagsAreaPerturbed[j] = false;
    }
#endif
}

void
MultiscaleModelFSI3DActivated::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3DActivated::updateSolution() \n";
#endif

    M_FSIoperator->updateSolution ( *M_stateVariable );
}

void
MultiscaleModelFSI3DActivated::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3DActivated::saveSolution() \n";
#endif

    exportFluidSolution();
    if ( M_FSIoperator->isFluid() )
    {
        M_exporterFluid->postProcess ( M_data->dataFluid()->dataTime()->time() );
    }

    exportSolidSolution();
    if ( M_FSIoperator->isSolid() )
    {
        M_exporterSolid->postProcess ( M_data->dataSolid()->dataTime()->time() );
    }

#ifdef HAVE_HDF5
    if ( M_data->dataFluid()->dataTime()->isLastTimeStep() )
    {
        if ( M_FSIoperator->isFluid() )
        {
            M_exporterFluid->closeFile();
        }
        if ( M_FSIoperator->isSolid() )
        {
            M_exporterSolid->closeFile();
        }
    }
#endif

}

void
MultiscaleModelFSI3DActivated::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();

        std::cout << "FSI method          = " << M_data->method() << std::endl << std::endl;

        std::cout << "Velocity FE order   = " << M_FSIoperator->dataFluid()->uOrder() << std::endl
                  << "Pressure FE order   = " << M_FSIoperator->dataFluid()->pOrder() << std::endl
                  << "Structure FE order  = " << M_FSIoperator->dataSolid()->order() << std::endl << std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_FSIoperator->uFESpace().dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_FSIoperator->pFESpace().dof().numTotalDof() << std::endl
                  << "Harmonic ext. DOF   = " << M_FSIoperator->mmFESpace().dof().numTotalDof() << std::endl
                  << "Structure DOF       = " << M_FSIoperator->dFESpace().dof().numTotalDof() << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->uFESpace().mesh() ) ).maxH << std::endl
                  << "Fluid mesh meanH    = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->uFESpace().mesh() ) ).meanH << std::endl
                  << "Solid mesh maxH     = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->dFESpace().mesh() ) ).maxH << std::endl
                  << "Solid mesh meanH    = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->dFESpace().mesh() ) ).meanH << std::endl << std::endl;
    }
}

Real
MultiscaleModelFSI3DActivated::checkSolution() const
{
    return M_stateVariable->norm2();
}



} // Namespace multiscale
} // Namespace LifeV
