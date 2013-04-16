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
    M_exporterElectro              (),
    M_fiber						   (),
    M_gammaf						()
{
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelFSI3DActivated::setupData ( const std::string& fileName )
{
	super::setupData(fileName);

	GetPot dataFile(fileName);

    	std::string xmlpath = dataFile("electrophysiology/monodomain_xml_path","./");
    	std::string xmlfile = dataFile("electrophysiology/monodomain_xml_file","MonodomainSolverParamList.xml");
    	Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
//	    monodomainList = Teuchos::getParametersFromXmlFile ( xmlpath + xmlfile );
	  //  monodomainList = Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" );

    	minimalModelPtr_Type  ionicModel ( new minimalModel_Type() );

	    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
	    std::string meshPath = monodomainList.get ("mesh_path", "./");

	    M_monodomain.reset ( new monoSolver_Type ( meshName, meshPath, dataFile, ionicModel ) );
        M_monodomain -> setParameters ( monodomainList );


        M_exporterElectro -> setDataFromGetPot ( dataFile );
        std::string prefix = multiscaleProblemPrefix  + number2string ( M_ID ) +  "electrophysiology" + "_" + number2string ( multiscaleProblemStep );
        M_exporterElectro->setPostDir ( multiscaleProblemFolder );
  //      M_monodomain -> setupExporter ( *M_exporterElectro, prefix);

        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
        ( new FESpace< mesh_Type, MapEpetra > ( M_monodomain -> localMeshPtr(), "P1", 3, M_monodomain -> commPtr() ) );
        M_fiber.reset ( new vector_Type ( Space3D -> map() ) );
    	std::string nm = monodomainList.get("fiber_file","FiberDirection") ;
        HeartUtility::importFibers( M_fiber, nm, M_monodomain -> localMeshPtr() );
        M_monodomain -> copyFiber(M_fiber);


        Space3D.reset();
        M_gammaf.reset( new vector_Type( M_monodomain -> potentialPtr() -> map() ) );
}

void
MultiscaleModelFSI3DActivated::setupModel()
{

	super::setupModel();


    HeartUtility::setValueOnBoundary( *(M_monodomain -> potentialPtr() ), M_monodomain -> fullMeshPtr(), 1.0, 43 );
    HeartUtility::setValueOnBoundary( *(M_monodomain -> potentialPtr() ), M_monodomain -> fullMeshPtr(), 1.0, 45 );

//    function_Type f = &smoothing;
//    vectorPtr_Type smoother( new vector_Type( electro -> potentialPtr() -> map() ) );
//    electro -> feSpacePtr() -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *smoother , 0);
//    (*smoother) *= *(electro -> potentialPtr() );
//    electro -> setPotentialPtr(smoother);

    //setting up initial conditions
    * ( M_monodomain -> globalSolution().at (1) ) = 1.0;
    * ( M_monodomain -> globalSolution().at (2) ) = 1.0;
    * ( M_monodomain -> globalSolution().at (3) ) = 0.021553043080281;


    super::solver() -> solid().material() -> setFiberVector( *M_fiber );
    super::solver() -> solid().material() -> setGammaf( *M_gammaf );



}

void
MultiscaleModelFSI3DActivated::buildModel()
{
	super::buildModel();

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
	M_monodomain -> setInitialTime( base::globalData() -> dataTime() -> time() );
	M_monodomain -> setEndTime( ( M_monodomain -> initialTime() ) + base::globalData() -> dataTime() -> timeStep() );
	M_monodomain ->solveSplitting();


    M_gammaf.reset( new  vector_Type( *( M_monodomain -> globalSolution().at(3) ) ) );

    //rescaling parameters for gammaf with minimal model
	Real maxCalciumLikeVariable = 0.838443;
    Real minCalciumLikeVariable = 0.021553;
    Real beta = -0.3;

    HeartUtility::rescaleVector( *M_gammaf, minCalciumLikeVariable, maxCalciumLikeVariable, beta);
    super::solver() -> solid().material() -> setGammaf( *M_gammaf );


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
//	M_monodomain -> exportSolution( *M_exporterElectro, base::globalData() -> dataTime() -> time() );
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



} // Namespace multiscale
} // Namespace LifeV
