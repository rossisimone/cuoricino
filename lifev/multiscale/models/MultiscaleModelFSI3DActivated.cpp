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
    M_fiber						   (),
    M_fullSolidMesh				   (),
    M_monodomain                   (),
    M_importerElectro              (),
    M_exporterElectro              (),
    M_gammaf					   (),
    M_usingDifferentMeshes		   (false),
    M_oneWayCoupling			   (true),
    M_activationCenter             (3),
    M_activationRadius			   (),
    M_activationMarker			   (),
    M_dataFileName				   (),
    M_gammafSolid                  (),
    M_displacementMonodomain       (),
    M_activationSpacePtr			(),
    M_minCalciumLikeVariable		(0.21),
    M_maxCalciumLikeVariable		(0.85)
{
	M_coarseToFineInterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( "RBFlocallyRescaledVectorial" ) );
	M_fineToCoarseInterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( "RBFlocallyRescaledScalar" ) );
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelFSI3DActivated::setupData ( const std::string& fileName )
{

    // FSI setup
	M_dataFileName = fileName;
    super::setupData(M_dataFileName);
    M_fullSolidMesh.reset( new mesh_Type( M_FSIoperator -> solidMesh() ) );
    // Meshes and monodomain data

    GetPot dataFile(M_dataFileName);

    std::string xmlpath = dataFile("electrophysiology/monodomain_xml_path","./");
    std::string xmlfile = dataFile("electrophysiology/monodomain_xml_file","MonodomainSolverParamList.xml");
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( xmlpath + xmlfile ) );

    minimalModelPtr_Type  ionicModel ( new minimalModel_Type() );

    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = monodomainList.get ("mesh_path", "./");

    // Exporters

    const std::string exporterType = dataFile ( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( exporterType.compare ( "hdf5" ) == 0 )
    {
        M_exporterElectro.reset ( new hdf5IOFile_Type() );
    }
    else
#endif
        M_exporterElectro.reset ( new ensightIOFile_Type() );

    M_monodomain.reset ( new monoSolver_Type ( meshName, meshPath, dataFile, ionicModel ) );
    M_monodomain -> setParameters ( monodomainList );

    // Activation parameters for the initial condition

    M_activationCenter[0] = monodomainList.get ("activation_center_X", 0.);
    M_activationCenter[1] = monodomainList.get ("activation_center_Y", 0.);
    M_activationCenter[2] = monodomainList.get ("activation_center_Z", 0.);
    M_activationRadius = monodomainList.get ("activation_radius", 1.);
    M_activationMarker = monodomainList.get ("activation_marker", 0);

    M_usingDifferentMeshes = monodomainList.get("using_different_meshes", false);
    M_oneWayCoupling	   = monodomainList.get("one_way_coupling", true);
    // Electrics solution exporter

    M_exporterElectro -> setDataFromGetPot ( dataFile );
    std::string prefix = multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) +  "_Electro_" + number2string ( multiscaleProblemStep );
    M_exporterElectro->setPostDir ( multiscaleProblemFolder );
    M_monodomain -> setupExporter ( *M_exporterElectro, prefix);

    // Fiber directions

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( M_monodomain -> localMeshPtr(), "P1", 3, M_monodomain -> commPtr() ) );
    M_fiber.reset ( new vector_Type ( Space3D -> map() ) );
    std::string nm = monodomainList.get("fiber_file","FiberDirection") ;
    HeartUtility::importFibers( M_fiber, nm, M_monodomain -> localMeshPtr() );
    M_monodomain-> setFiberPtr(M_fiber);

    // Activation function
    M_displacementMonodomain.reset( new vector_Type( Space3D -> map() ) );

    Space3D.reset();
    M_gammaf.reset( new vector_Type( M_monodomain -> potentialPtr() -> map() ) );
}


LifeV::Real
MultiscaleModelFSI3DActivated::activationFunction(const Real& t, const Real& x, const Real& y, const Real& z, const LifeV::ID& i)
{
    return std::exp( -( std::pow(x-M_activationCenter[0],2) + std::pow(y-M_activationCenter[1],2) + std::pow(z-M_activationCenter[2],2) ) / std::pow(M_activationRadius,2) );
}

void
MultiscaleModelFSI3DActivated::setupModel()
{

	super::setupModel();
	boost::shared_ptr<mesh_Type> solidLocalMeshPtr( new super::mesh_Type( super::solver() -> solidLocalMesh() ) );
	M_activationSpacePtr.reset( new FESpace<mesh_Type, MapEpetra>( solidLocalMeshPtr, "P1", 1,  M_comm ) );
	M_gammafSolid.reset( new vector_Type( M_activationSpacePtr -> map() ) );
    M_exporterSolid->addVariable ( IOData_Type::ScalarField, "activation", M_activationSpacePtr , M_gammafSolid, static_cast<UInt> ( 0 ) );


    M_monodomain ->	setInitialConditions();


	HeartUtility::importFibers( super::solver() -> solid().material() -> fiberVector(),
    							super::solver() -> solid().material() -> materialData() -> fileFiberDirections(),
    							solidLocalMeshPtr );

    HeartUtility::setValueOnBoundary( *(M_monodomain -> potentialPtr() ), M_monodomain -> fullMeshPtr(), 1.0, M_activationMarker );

    function_Type f( boost::bind ( &MultiscaleModelFSI3DActivated::activationFunction, this, _1, _2, _3, _4, _5 ) );
    vectorPtr_Type smoother( new vector_Type( M_monodomain -> potentialPtr() -> map() ) );
    M_monodomain -> feSpacePtr() -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *smoother , 0);
    (*smoother) *= *(M_monodomain -> potentialPtr() );
    //M_monodomain -> setPotentialPtr(smoother);
    M_monodomain -> setPotential(*smoother);

    //setting up initial conditions

    super::solver() -> solid().material() -> setGammaf( *M_gammaf );
    if(!M_oneWayCoupling) M_monodomain -> setDisplacementPtr( M_displacementMonodomain );


    setupInterpolant();



}

void
MultiscaleModelFSI3DActivated::buildModel()
{
	super::buildModel();
	if( M_usingDifferentMeshes )
	{
	M_fineToCoarseInterpolant -> buildOperators();

	assert( M_fullSolidMesh != 0 && "Burn in hell!!!");
		//if( !M_oneWayCoupling ) M_coarseToFineInterpolant -> buildOperators();
	}

	M_monodomain -> setupLumpedMassMatrix();
    M_monodomain -> setupStiffnessMatrix();
    M_monodomain -> setupGlobalMatrix();

}

void
MultiscaleModelFSI3DActivated::updateModel()
{
	super::updateModel();
}

void
MultiscaleModelFSI3DActivated::solveModel()
{

    if ( M_nonLinearRichardsonIteration == 0 )
    {
    	if(!M_oneWayCoupling)
    	{
			if( M_usingDifferentMeshes )
			{
				M_coarseToFineInterpolant -> updateRhs( M_solidDisplacement );
				M_coarseToFineInterpolant -> interpolate();
				M_coarseToFineInterpolant -> solution ( M_displacementMonodomain );
			}
			else M_displacementMonodomain = M_solidDisplacement;

			M_monodomain -> setDisplacementPtr( M_displacementMonodomain );
			M_monodomain -> setupStiffnessMatrix();
			M_monodomain -> setupGlobalMatrix();
    	}
        // TODO: Better handling of different time steps (HeartETAMonodomainSolver uses ms for time)
        Real timeStep = base::globalData() -> dataTime() -> timeStep();
        Real tn       = base::globalData() -> dataTime() -> time() - timeStep;

        M_monodomain -> setInitialTime( 1000.0 * tn );
        M_monodomain -> setEndTime( 1000.0 * (tn + timeStep) );
        M_monodomain ->solveSplitting();

        *M_gammaf = *( M_monodomain -> globalSolution().at(3) );
        //M_gammaf.reset( new  vector_Type( *( M_monodomain -> globalSolution().at(3) ) ) );


        if( M_maxCalciumLikeVariable < M_gammaf -> maxValue() )
        		M_maxCalciumLikeVariable = M_gammaf -> maxValue();
        if( M_minCalciumLikeVariable > M_gammaf -> minValue() )
        		M_minCalciumLikeVariable = M_gammaf -> minValue();
        Real beta = -0.3;

        HeartUtility::rescaleVector( *M_gammaf, M_minCalciumLikeVariable, M_maxCalciumLikeVariable, beta);
        //rescaling parameters for gammaf with minimal model
        //M_gammafSolid.reset( new vector_Type( M_activationSpacePtr -> map() ) );
        if( M_usingDifferentMeshes )
        {
        M_fineToCoarseInterpolant -> updateRhs( M_gammaf );
        M_fineToCoarseInterpolant -> interpolate();
        M_fineToCoarseInterpolant -> solution ( M_gammafSolid );
        }
        else M_gammafSolid = M_gammaf;




        *M_gammafSolid *= ( M_gammafSolid -> operator <=(0.0) );

        super::solver() -> solid().material() -> setGammaf( *M_gammafSolid );

    }

    super::solveModel();



}

void
MultiscaleModelFSI3DActivated::updateSolution()
{
	super::updateSolution();
}

void
MultiscaleModelFSI3DActivated::saveSolution()
{
	super::saveSolution();
	M_monodomain -> exportSolution( *M_exporterElectro, base::globalData() -> dataTime() -> time() );
    if ( super::data() ->dataFluid()->dataTime()->isLastTimeStep() )
    {
        	M_exporterElectro->closeFile();
    }
}

void
MultiscaleModelFSI3DActivated::showMe()
{
	super::showMe();
	if ( M_comm->MyPID() == 0 )
    {
        std::cout << "Ionic model: Minimal Model" << std::endl << std::endl;
        std::cout << "Electrophysiology model: Monodomain" << std::endl << std::endl;
        std::cout << "Electrophysiology DOF = " << ( M_monodomain -> ionicModelPtr() -> Size() ) * M_monodomain -> feSpacePtr()  -> dof().numTotalDof() << std::endl << std::endl;
    }


}

Real
MultiscaleModelFSI3DActivated::checkSolution() const
{
    super::checkSolution();
}


void MultiscaleModelFSI3DActivated::setupInterpolant()
{
	    GetPot dataFile( M_dataFileName );
	    std::string xmlpath = dataFile("electrophysiology/monodomain_xml_path","./");
	    std::string xmlfile = dataFile("electrophysiology/monodomain_xml_file","MonodomainSolverParamList.xml");

	    Teuchos::RCP< Teuchos::ParameterList > monodomainList = Teuchos::rcp ( new Teuchos::ParameterList );
	    monodomainList = Teuchos::getParametersFromXmlFile (  xmlpath + xmlfile );

	    Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( xmlpath + xmlfile ) );
	    std::string meshName = list.get ("solid_mesh_name", "lid16.mesh");
	    std::string meshPath = list.get ("solid_mesh_path", "./");
		boost::shared_ptr<mesh_Type> solidLocalMeshPtr( new super::mesh_Type( super::solver() -> solidLocalMesh() ) );
		MeshUtility::fillWithFullMesh(solidLocalMeshPtr, M_fullSolidMesh, meshName, meshPath);



	    if(M_usingDifferentMeshes)
	    {
	    M_fineToCoarseInterpolant -> setup( M_monodomain -> fullMeshPtr(), M_monodomain -> localMeshPtr(),
	    		M_fullSolidMesh, M_FSIoperator -> solidLocalMeshPtr(),  std::vector<int>(1,-1) );
	    M_fineToCoarseInterpolant -> setRadius( (double) MeshUtility::MeshStatistics::computeSize ( * (M_monodomain -> fullMeshPtr() ) ).maxH );
	    M_fineToCoarseInterpolant -> setupRBFData ( M_gammaf, M_gammafSolid, dataFile, monodomainList);

			if( !M_oneWayCoupling )
			{
			M_coarseToFineInterpolant -> setup( M_fullSolidMesh, M_FSIoperator -> solidLocalMeshPtr(),
												M_monodomain -> fullMeshPtr(), M_monodomain -> localMeshPtr(),
												std::vector<int>(1,-1) );
//			M_coarseToFineInterpolant -> setup( M_fullSolidMesh, M_FSIoperator -> solidLocalMeshPtr(),
//												M_monodomain -> fullMeshPtr(), M_monodomain -> localMeshPtr(),
//												std::vector<int>(1,-1) );
			//M_coarseToFineInterpolant -> setRadius( (double) MeshUtility::MeshStatistics::computeSize ( * (M_fullSolidMesh) ).maxH );
			//if( M_monodomain -> displacementPtr() ) return 0.0;
			assert( M_monodomain -> displacementPtr() != 0 && "Burn in hell!!!");
			M_coarseToFineInterpolant -> setupRBFData ( M_solidDisplacement , M_monodomain -> displacementPtr(), dataFile, monodomainList);
     	 M_coarseToFineInterpolant -> buildOperators();
			}
	    }
}

} // Namespace multiscale
} // Namespace LifeV
