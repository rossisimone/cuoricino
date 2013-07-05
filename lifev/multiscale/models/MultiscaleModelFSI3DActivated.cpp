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

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <boost/typeof/typeof.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

namespace LifeV
{

namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelFSI3DActivated::MultiscaleModelFSI3DActivated() :
    super                          (),
    M_fiber                        (),
    M_fullSolidMesh                (),
    M_monodomain                   (),
    M_importerElectro              (),
    M_exporterElectro              (),
    M_gammaf                       (),
    M_gammas                       (),
    M_gamman                       (),
    M_gammasSolid                  (),
    M_gammanSolid                  (),
    M_usingDifferentMeshes         (false),
    M_oneWayCoupling               (true),
    M_activationCenter             (3),
    M_activationRadius             (),
    M_activationMarker             (),
//    M_endocardiumMarker				(),
//    M_epicardiumMarker				(),
    M_dataFileName                 (),
    M_gammafSolid                  (),
    M_displacementMonodomain       (),
    M_monodomainDisplacementETFESpace(),
    M_activationSpacePtr           (),
    M_minCalciumLikeVariable       (0.21),
    M_maxCalciumLikeVariable       (0.85),
    M_activationSolver             (),
    M_coarseToFineInterpolant		(),
    M_fineToCoarseInterpolant		(),
	M_rescalingVector				(),
	M_activationETFESpace			(),
	M_orthotropicActivationFactor	(3),
	M_activationOperator			(),
	M_preloadVector					(),
	M_preloadInTime					(false)
{
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelFSI3DActivated::setupData ( const std::string& fileName )
{

	// FSI setup
    M_dataFileName = fileName;
    super::setupData (M_dataFileName);
    M_fullSolidMesh.reset ( new mesh_Type ( M_FSIoperator -> solidMesh() ) );

    // Meshes and monodomain data

    GetPot dataFile (M_dataFileName);

    std::string xmlpath = dataFile ("electrophysiology/monodomain_xml_path", "./");
    std::string xmlfile = dataFile ("electrophysiology/monodomain_xml_file", "MonodomainSolverParamList.xml");
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( xmlpath + xmlfile ) );

    minimalModelPtr_Type  ionicModel ( new minimalModel_Type() );

    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = monodomainList.get ("mesh_path", "./");





    std::map< std::string, FSI3DActivated_ActivationModelType > activationModelTypeMap;
    activationModelTypeMap["Algebraic"]           = Algebraic;
    activationModelTypeMap["SimpleODE"]           = SimpleODE;
    activationModelTypeMap["StretchDependentODE"] = StretchDependentODE;

    std::map< std::string, FSI3DActivated_ActivationType > activationTypeMap;
    activationTypeMap["TransverselyIsotropic"] = TransverselyIsotropic;
    activationTypeMap["Orthotropic"]           = Orthotropic;

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
//    M_endocardiumMarker = monodomainList.get ("endocardium_marker", 0);
//    M_epicardiumMarker = monodomainList.get ("epicardium_marker", 0);


    M_usingDifferentMeshes = monodomainList.get ("using_different_meshes", false);
    M_oneWayCoupling       = monodomainList.get ("one_way_coupling", true);

    // Electrics solution exporter

    M_exporterElectro -> setDataFromGetPot ( dataFile );
    std::string prefix = multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) +  "_Electro_" + number2string ( multiscaleProblemStep );
    M_exporterElectro->setPostDir ( multiscaleProblemFolder );
    M_monodomain -> setupExporter ( *M_exporterElectro, prefix);

    // Fiber directions
	boost::shared_ptr<Epetra_Comm> comm( M_monodomain -> commPtr() );
    M_monodomainDisplacementETFESpace.reset
    ( new vectorialETFESpace_Type ( M_monodomain -> localMeshPtr(),  &feTetraP1,  comm ) );
    M_fiber.reset ( new vector_Type ( M_monodomainDisplacementETFESpace -> map() ) );
    std::string nm = monodomainList.get ("fiber_file", "FiberDirection") ;
    HeartUtility::importFibers ( M_fiber, nm, M_monodomain -> localMeshPtr() );
    M_monodomain-> setFiberPtr (M_fiber);



    // Activation function

    M_displacementMonodomain.reset ( new vector_Type ( M_monodomainDisplacementETFESpace -> map() ) );

    M_gammaf.reset( new vector_Type ( M_monodomain -> potentialPtr() -> map() ) );
    M_gammas.reset( new vector_Type( M_gammaf -> map() ) );
    M_gamman.reset( new vector_Type( M_gammaf -> map() ) );

    M_activationModelType 	= activationModelTypeMap[dataFile ( "electrophysiology/activation_model", "SimpleODE" )];
    M_activationType 		= activationTypeMap[dataFile ( "electrophysiology/activation_type", "TransverselyIsotropic" )];

    M_preloadInTime			= dataFile("solid/physics/preload",0);
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >          FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>                        FESpacePtr_Type;

    FESpacePtr_Type aFESpace ( new FESpace_Type (M_monodomain -> localMeshPtr(), "P1", 1, M_monodomain -> commPtr()) );
//    scalarETFESpace_Type cippalippa(M_monodomain -> localMeshPtr(), &feTetraP1, M_monodomain -> localMeshPtr() -> comm() );

//    boost::shared_ptr<mesh_Type> mesh(M_monodomain -> localMeshPtr());
        M_activationETFESpace.reset ( new scalarETFESpace_Type ( M_monodomain -> localMeshPtr(), &feTetraP1,  comm ) );


    //=====================
    //*********************
    // Activation solver  *
    //*********************
    //=====================
    if(M_activationModelType==StretchDependentODE)
    {
	    prec_Type* precRawPtr;
		basePrecPtr_Type precPtr;
		precRawPtr = new prec_Type;
		precRawPtr->setDataFromGetPot (dataFile, "electrophysiology/prec");
		precPtr.reset (precRawPtr);

		Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp (
																	  new Teuchos::ParameterList);

		xmlpath = dataFile ("electrophysiology/monodomain_xml_path", "./");
		xmlfile = dataFile ("electrophysiology/monodomain_xml_file", "MonodomainSolverParamList.xml");

		solverParamList = Teuchos::getParametersFromXmlFile (xmlpath + xmlfile);

		M_activationSolver.reset ( new LinearSolver() );
		M_activationSolver->setCommunicator ( M_monodomain -> commPtr() );
		M_activationSolver->setParameters ( *solverParamList );
		M_activationSolver->setPreconditioner ( precPtr );

	    M_activationOperator.reset(new matrix_Type( M_monodomain -> massMatrixPtr() -> map() ) ) ;

	  	{
	  		using namespace ExpressionAssembly;

	  		integrate(elements(M_monodomain -> localMeshPtr() ), M_monodomain -> feSpacePtr() -> qr(), M_monodomain -> ETFESpacePtr(),
	  				M_monodomain -> ETFESpacePtr(), phi_i * phi_j) >> M_activationOperator;
	  	}

	  	M_activationOperator -> globalAssemble();

	  	M_activationSolver->setOperator( M_activationOperator );

	  	M_preloadVector.reset( new vector_Type( M_fiber -> map() ) );
    }
//=======================================
// setup INTERPOLATION
//=======================================
    M_coarseToFineInterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( monodomainList.get ("disp_interpolation_type","RBFrescaledVectorial" ) ) );
    M_fineToCoarseInterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( monodomainList.get ("gammaf_interpolation_type","RBFrescaledScalar" ) ) );

}


LifeV::Real
MultiscaleModelFSI3DActivated::activationFunction (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const LifeV::ID& /*i*/)
{
    return std::exp ( - ( std::pow (x - M_activationCenter[0], 2) + std::pow (y - M_activationCenter[1], 2) + std::pow (z - M_activationCenter[2], 2) ) / std::pow (M_activationRadius, 2) );
}

void
MultiscaleModelFSI3DActivated::setupModel()
{

    super::setupModel();

    boost::shared_ptr<mesh_Type> solidLocalMeshPtr ( new super::mesh_Type ( super::solver() -> solidLocalMesh() ) );

    M_activationSpacePtr.reset ( new FESpace<mesh_Type, MapEpetra> ( solidLocalMeshPtr, "P1", 1,  M_comm ) );

    M_gammafSolid.reset ( new vector_Type ( M_activationSpacePtr -> map() ) );
    M_gammasSolid.reset( new vector_Type( M_gammafSolid -> map() ) );
    M_gammanSolid.reset( new vector_Type( M_gammafSolid -> map() ) );

    M_exporterSolid->addVariable ( IOData_Type::ScalarField, "activation", M_activationSpacePtr , M_gammafSolid, static_cast<UInt> ( 0 ) );

    M_monodomain -> setInitialConditions();


    HeartUtility::setValueOnBoundary ( * (M_monodomain -> potentialPtr() ), M_monodomain -> fullMeshPtr(), 1.0, M_activationMarker );

    function_Type f ( boost::bind ( &MultiscaleModelFSI3DActivated::activationFunction, this, _1, _2, _3, _4, _5 ) );
    vectorPtr_Type smoother ( new vector_Type ( M_monodomain -> potentialPtr() -> map() ) );
    M_monodomain -> feSpacePtr() -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *smoother , 0);
    (*smoother) *= * (M_monodomain -> potentialPtr() );
    M_monodomain -> setPotential (*smoother);

    //setting up initial conditions

    super::solver() -> solid().material() -> setGammaf ( *M_gammaf );
    if (!M_oneWayCoupling)
    {
        M_monodomain -> setDisplacementPtr ( M_displacementMonodomain );
    }
    setupInterpolant();

    //SETUP RESCALING VECTOR

    GetPot dataFile ( M_dataFileName );
    std::string xmlpath = dataFile ("electrophysiology/monodomain_xml_path", "./");
    std::string xmlfile = dataFile ("electrophysiology/monodomain_xml_file", "MonodomainSolverParamList.xml");
    Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( xmlpath + xmlfile ) );
    std::string filename = list.get ("filename", "");
    if(filename == "" )
    {
    	M_rescalingVector.reset();
    }
    else
    {
        M_rescalingVector.reset(new vector_Type( M_gammafSolid -> map() ) );
		std::string fieldname = list.get ("fieldname", "");
		HeartUtility::importScalarField(M_rescalingVector, filename, fieldname, super::solver() -> solidLocalMeshPtr() );
    }
    if( super::solver() -> data().dataFluid()->dataTime()->initialTime() > 0)
    {
    	std::string prefix = multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) +  "_Electro_" + number2string ( multiscaleProblemStep - 1);
    	M_monodomain -> importSolution(dataFile,prefix,multiscaleProblemFolder,super::solver() -> data().dataFluid()->dataTime()->initialTime() );
    }
    	//    *M_activationSolver = *M_monodomain -> linearSolverPtr();
    //setting up activation solver
    vectorPtr_Type fiber(new vector_Type( super::solver() -> dFESpace().map() ) );
    HeartUtility::importFibers ( fiber,
                                 super::solver() -> solid().material() -> materialData() -> fileFiberDirections(),
                                 solidLocalMeshPtr );
    super::solver() ->solid().material() -> setFiberVector(*fiber);

    // sheet direction
    vectorPtr_Type sheet(new vector_Type( super::solver() -> dFESpace().map() ) );
    HeartUtility::importVectorField(sheet,
    									sheetFileName,
    									sheetFieldName,
    									super::solver() -> solidLocalMeshPtr() );
    super::solver() ->solid().material() -> setSheetVector(*sheet);

}

void
MultiscaleModelFSI3DActivated::buildModel()
{
    super::buildModel();
    if ( M_usingDifferentMeshes )
    {
        M_fineToCoarseInterpolant -> buildOperators();
        //if( !M_oneWayCoupling ) M_coarseToFineInterpolant -> buildOperators();
    }

    super::solver() -> solid().material() -> showMyParameters();

    M_monodomain -> setupLumpedMassMatrix();
    M_monodomain -> setupStiffnessMatrix();
    M_monodomain -> setupGlobalMatrix();

}

void
MultiscaleModelFSI3DActivated::updateModel()
{
    super::updateModel();

	if(M_preloadInTime)
	{
		if( M_monodomain -> globalSolution().at(3)-> minValue() < 0.02158)
		{
			int d = M_monodomain -> globalSolution().at(3) -> epetraVector().MyLength();
			int size =  M_monodomain -> globalSolution().at(3) -> size();
			for(int l(0); l < d; l++)
			{
				int m1 = M_monodomain -> globalSolution().at(3) -> blockMap().GID(l);
				//cout << m1 << "\t" << size << "\n";
				if( (*(M_monodomain -> globalSolution().at(3)))[m1] <= 0.02158)
				{
//					std::cout << "\nchanging at: " << m1 ;
					int m2 = super::solver() -> solid().displacementPtr() -> blockMap().GID(l + size);
					int m3 = super::solver() -> solid().displacementPtr() -> blockMap().GID(l + 2 * size);
					//std::cout << "\n" << (*emDisp0)[m1] << " has become: ";
					(*M_preloadVector)[m1] = (*(super::solver() -> solid().displacementPtr()))[m1];
					//std::cout << (*emDisp0)[m1] << "\n";
					(*M_preloadVector)[m2] = (*(super::solver() -> solid().displacementPtr()))[m2];
					(*M_preloadVector)[m3] = (*(super::solver() -> solid().displacementPtr()))[m3];
				}
				//std::cout << "\n";

			}
		}
	}
}

void
MultiscaleModelFSI3DActivated::solveModel()
{

    if ( M_nonLinearRichardsonIteration == 0 )
    {
        if (!M_oneWayCoupling)
        {
            if ( M_usingDifferentMeshes )
            {
                M_coarseToFineInterpolant -> updateRhs ( M_solidDisplacement );
                M_coarseToFineInterpolant -> interpolate();
                M_coarseToFineInterpolant -> solution ( M_displacementMonodomain );
            }
            else
            {
                M_displacementMonodomain = M_solidDisplacement;
            }

            M_monodomain -> setDisplacementPtr ( M_displacementMonodomain );
            M_monodomain -> setupStiffnessMatrix();
            M_monodomain -> setupGlobalMatrix();
        }
        // TODO: Better handling of different time steps (HeartETAMonodomainSolver uses ms for time)
        Real timeStep = base::globalData() -> dataTime() -> timeStep();
        Real tn       = base::globalData() -> dataTime() -> time() - timeStep;

        M_monodomain -> setInitialTime ( 1000.0 * tn );
        M_monodomain -> setEndTime ( 1000.0 * (tn + timeStep) );
        M_monodomain -> solveSplitting();

        switch (M_activationModelType)
        {
            case Algebraic:
            {
                // Simplistic activation model (\gammaf = a * Ca2)
                *M_gammaf = * ( M_monodomain -> globalSolution().at (3) );
                if ( M_maxCalciumLikeVariable < M_gammaf -> maxValue() )
                    M_maxCalciumLikeVariable = M_gammaf -> maxValue();
                if ( M_minCalciumLikeVariable > M_gammaf -> minValue() )
                    M_minCalciumLikeVariable = M_gammaf -> minValue();

                Real beta = -0.3;
                HeartUtility::rescaleVector ( *M_gammaf, M_minCalciumLikeVariable, M_maxCalciumLikeVariable, beta);
                break;
            }
            case SimpleODE:
            {
                // More prolonged activation model (\gammaf' = a * Ca2 + b * \gammaf)
                *M_gammaf += 1000 * timeStep * ( -0.02 * *( M_monodomain -> globalSolution().at(3) ) - 0.04 * (*M_gammaf));
                break;
            }
            case StretchDependentODE:
            {
                //ASSERT(false, "ERROR: Stretch-dependent activation is not yet implemented.");

                switch (M_activationType)
                {
                	case TransverselyIsotropic:
                		*M_gammas *= 0.0;
                		*M_gamman *= 0.0;

                		*M_gammas = 1.0 ;
                		*M_gammas /= ( 1.0 + *M_gammaf );
                		M_gammas -> sqrt();
                		*M_gammas -= 1.0;
                		*M_gamman = *M_gammas;
                		break;
                	case Orthotropic:
                		*M_gammas *= 0.0;
                		*M_gamman *= 0.0;

                		*M_gamman = *M_gammaf;
                		*M_gamman *= M_orthotropicActivationFactor;
                		*M_gammas = 1.0 ;
                		*M_gammas /= (1.0 + *M_gammaf);
                		*M_gammas /= ( 1.0 + *M_gamman);
                		*M_gammas -= 1.0;
                		break;
                }

               	//boost::shared_ptr<FLRelationshipGamma> flg (new FLRelationshipGamma);
               	boost::shared_ptr<FLRelationship> fl (new FLRelationship);

               	boost::shared_ptr<HeavisideFct> H (new HeavisideFct);

               	//boost::shared_ptr<Exp> EXP(new Exp);
               	//boost::shared_ptr<Exp2> EXP2(new Exp2);
               	//boost::shared_ptr<Psi4f> psi4f (new Psi4f);

                MatrixSmall<3,3> Id;
                Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
                Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
                Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;

              	vectorPtr_Type rhsActivation( new vector_Type( *M_gammaf ) );
               	vectorPtr_Type tmpRhsActivation( new vector_Type ( rhsActivation -> map(), Repeated ) );
  			  {


  					{
  						using namespace ExpressionAssembly;


  						BOOST_AUTO_TPL(I,      value(Id) );
  						BOOST_AUTO_TPL(Grad_u, grad( M_monodomainDisplacementETFESpace, *M_displacementMonodomain, 0) );
  						BOOST_AUTO_TPL(F,      ( Grad_u + I ) );
  						BOOST_AUTO_TPL(FmT,    minusT(F) );
  						BOOST_AUTO_TPL(J,       det(F) );
  						BOOST_AUTO_TPL(Jm23,    pow(J, -2./3));
  						BOOST_AUTO_TPL(I1,     dot(F, F));

  						// Fibres
  						BOOST_AUTO_TPL(f0,     value( M_monodomainDisplacementETFESpace, *( M_monodomain -> fiberPtr() ) ) );
  						BOOST_AUTO_TPL(f,      F * f0 );
  						BOOST_AUTO_TPL(I4f,    dot(f, f) );

  					    //BOOST_AUTO_TPL(s0,     value(M_monodomain -> ETFESpacePtr(), *(super::solver() -> solid().material() -> sheetVectorPtr() ) ) );
  						//BOOST_AUTO_TPL(I4s,    dot(F * s0, F * s0));

  						//BOOST_AUTO_TPL(I1iso,   Jm23 * I1);
  						BOOST_AUTO_TPL(I4fiso,  Jm23 * I4f);
  					    //BOOST_AUTO_TPL(I4siso,  Jm23 * I4s);

  						// shortenings
  					    BOOST_AUTO_TPL(gf,  value(M_activationETFESpace, *M_gammaf));
 				    	//BOOST_AUTO_TPL(gs,  value(M_activationETFESpace, *M_gammas));
  				    	//BOOST_AUTO_TPL(gn,  value(M_activationETFESpace, *M_gamman));

  				        //BOOST_AUTO_TPL(dI1edI1,   value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  )    );
  				        //BOOST_AUTO_TPL(dI1edI4f,  value(1.0) / (   (gf + value(1.0))  *  (gf + value(1.0))  ) -  value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  ) );
  				        //BOOST_AUTO_TPL(dI1edI4s,  value(1.0) / (   (gs + value(1.0))  *  (gs + value(1.0))  ) -  value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  ) );
  				        //BOOST_AUTO_TPL(dI4fedI4f, value(1.0) / (   (gf + value(1.0)) *   (gf + value(1.0))  )    );

  					    //BOOST_AUTO_TPL(I1eiso,   dI1edI1 * I1iso + dI1edI4f * I4fiso + dI1edI4s * I4siso);
  					    //BOOST_AUTO_TPL(I4feiso,  dI4fedI4f * I4fiso);
  						//BOOST_AUTO_TPL(I4feisom1, ( I4feiso - value(1.0) ) );

  						//initial
  						BOOST_AUTO_TPL(Grad_u_i, grad(M_monodomainDisplacementETFESpace, *M_preloadVector, 0) );
  						BOOST_AUTO_TPL(F_i,      ( Grad_u_i + I ) );
  						BOOST_AUTO_TPL(J_i,       det(F_i) );
  						BOOST_AUTO_TPL(Jm23_i,    pow(J_i, -2./3));
  						// Fibres
  						BOOST_AUTO_TPL(f_i,      F_i * f0 );
  						BOOST_AUTO_TPL(I4f_i,    dot(f_i, f_i) );
  						BOOST_AUTO_TPL(I4fiso_i,  Jm23_i * I4f_i);

  						BOOST_AUTO_TPL(dW0, value(-2.0) * I4fiso_i) ;
  						BOOST_AUTO_TPL(dW, value(-2.0) * I4fiso * pow(gf + value(1.0), -3) );

  						BOOST_AUTO_TPL(Ca,    value( M_activationETFESpace, *( M_monodomain -> globalSolution().at(3)  ) ) );
  						BOOST_AUTO_TPL(Ca2, Ca * Ca );

  						Real Ca_diastolic = -0.02155;
  						//dataFile( "solid/physics/Ca_diastolic", -0.02155 );
  						BOOST_AUTO_TPL(dCa, Ca + value(Ca_diastolic) );
  						Real active_coefficient = -2.5;
  						//dataFile( "solid/physics/active_coefficient", -2.5 );
  						BOOST_AUTO_TPL(Pa, value(active_coefficient) * eval(H, dCa) * eval(H, dCa) * eval(fl, I4fiso) + dW0 );
  						Real betaA = 0.0005;
  						//dataFile( "solid/physics/viscosity", 0.0005 );
  						BOOST_AUTO_TPL(beta, value( betaA ) );

  						BOOST_AUTO_TPL(gamma_dot, beta / ( Ca2 ) * ( Pa - dW )  );

  						{
							integrate ( elements ( M_monodomain -> localMeshPtr() ),
									M_monodomain -> feSpacePtr() -> qr() ,
									M_monodomain -> ETFESpacePtr(),
									gamma_dot * phi_i
							) >> tmpRhsActivation;
  						}

  					*rhsActivation *= 0;
  					*rhsActivation = ( (*M_activationOperator) * ( *M_gammaf ) );
  					*rhsActivation += ( ( M_monodomain -> timeStep() * *tmpRhsActivation ) );

  					M_activationSolver -> setRightHandSide(rhsActivation);
  					std::cout << "\nMax of GammaF before solving: " << M_gammaf -> maxValue()<< std::endl;
  					M_activationSolver -> solve(M_gammaf);
  					std::cout << "\nMax of GammaF after solving: " << M_gammaf -> maxValue()<< std::endl;
  					}

  					if( M_gammaf -> maxValue() > 0.0)
					{
						int d = M_gammaf -> epetraVector().MyLength();
						int size =  M_gammaf -> size();
						for(int l(0); l < d; l++)
						{
							std::cout << "\n*****************************************************";
							std::cout << "\nChanging the gamma back to zero: " ;
							std::cout << "\n*****************************************************";

							int m1 = M_gammaf -> blockMap().GID(l);
							//cout << m1 << "\t" << size << "\n";
							if( (*M_gammaf)[m1] > 0)
							{
								(*M_gammaf)[m1] = 0.0;
							}

						}
					}


  			  }
            }
        }

        // Transfer gammaf to solid mesh
        if ( M_usingDifferentMeshes )
        {
            M_fineToCoarseInterpolant -> updateRhs ( M_gammaf );
            M_fineToCoarseInterpolant -> interpolate();
            M_fineToCoarseInterpolant -> solution ( M_gammafSolid );
        }
        else
        {

            *M_gammafSolid = *M_gammaf;
        }

        if ( M_activationModelType == SimpleODE)
            *M_gammafSolid *= 0.3 / 0.415;

        *M_gammafSolid *= ( M_gammafSolid -> operator <= (0.0) );
        if(M_rescalingVector)
        {
        	*M_gammafSolid /= 0.3;
        	*M_gammafSolid *= *M_rescalingVector;
        }

        super::solver() -> solid().material() -> setGammaf ( *M_gammafSolid );
        switch (M_activationType)
        {
        	case TransverselyIsotropic:
        		*M_gammasSolid *= 0. ;
        		*M_gammanSolid *= 0. ;

        		*M_gammasSolid = 1.0 ;
        		*M_gammasSolid /= ( 1.0 + *M_gammafSolid );
        		M_gammasSolid -> sqrt();
        		*M_gammasSolid -= 1.0;
        		*M_gammanSolid = *M_gammasSolid;
        		break;
        	case Orthotropic:
        		*M_gammasSolid *= 0. ;
        		*M_gammanSolid *= 0. ;

        		*M_gammanSolid = *M_gammafSolid;
        		*M_gammanSolid *= M_orthotropicActivationFactor;
        		*M_gammasSolid = 1.0 ;
        		*M_gammasSolid /= (1.0 + *M_gammafSolid);
        		*M_gammasSolid /= ( 1.0 + *M_gammanSolid);
        		*M_gammasSolid -= 1.0;
        		break;
        }
        super::solver() -> solid().material() -> setGamman( *M_gammanSolid );
        super::solver() -> solid().material() -> setGammas( *M_gammasSolid );
        std::cout << "\nMax gammaf " << super::solver() -> solid().material() -> gammaf() -> maxValue() << std::endl;
        std::cout << "\nMin gammaf " << super::solver() -> solid().material() -> gammaf() -> minValue() << std::endl;

        std::cout << "\nMax gammas " << super::solver() -> solid().material() -> gammas() -> maxValue() << std::endl;
        std::cout << "\nMin gammas " << super::solver() -> solid().material() -> gammas() -> minValue() << std::endl;

        std::cout << "\nMax gamman " << super::solver() -> solid().material() -> gamman() -> maxValue() << std::endl;
        std::cout << "\nMin gamman " << super::solver() -> solid().material() -> gamman() -> minValue() << std::endl;

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
    M_monodomain -> exportSolution ( *M_exporterElectro, base::globalData() -> dataTime() -> time() );
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
    GetPot dataFile ( M_dataFileName );
    std::string xmlpath = dataFile ("interpolation/interpolation_xml_path", "./");
    std::string xmlfile = dataFile ("interpolation/interpolation_xml_file", "InterpolationParamList.xml");

    Teuchos::RCP< Teuchos::ParameterList > monodomainList = Teuchos::rcp ( new Teuchos::ParameterList );
    monodomainList = Teuchos::getParametersFromXmlFile (  xmlpath + xmlfile );

    Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( xmlpath + xmlfile ) );
    std::string meshName = list.get ("solid_mesh_name", "lid16.mesh");
    std::string meshPath = list.get ("solid_mesh_path", "./");
    boost::shared_ptr<mesh_Type> solidLocalMeshPtr ( new super::mesh_Type ( super::solver() -> solidLocalMesh() ) );
    MeshUtility::fillWithFullMesh (solidLocalMeshPtr, M_fullSolidMesh, meshName, meshPath);

    if (M_usingDifferentMeshes)
    {
        M_fineToCoarseInterpolant -> setup ( M_monodomain -> fullMeshPtr(), M_monodomain -> localMeshPtr(),
                                             M_fullSolidMesh, M_FSIoperator -> solidLocalMeshPtr(),  std::vector<int> (1, -1) );

        std::string fineToCoarseInterpolationType = list.get ("gammaf_interpolation_type", "RBFrescaledScalar" );
        std::string coarseToFineInterpolationType = list.get ("disp_interpolation_type", "RBFrescaledVectorial" );

        if (fineToCoarseInterpolationType != "RBFlocallyRescaledScalar" && fineToCoarseInterpolationType != "RBFhtp")
        {
            // Set radius only when using standard RBF interpolation
            M_fineToCoarseInterpolant -> setRadius ( (double) MeshUtility::MeshStatistics::computeSize ( * (M_fullSolidMesh ) ).maxH );
        }

        M_fineToCoarseInterpolant -> setupRBFData ( M_gammaf, M_gammafSolid, dataFile, monodomainList);

        if (fineToCoarseInterpolationType == "RBFscalar")
        {
            // Set RBF basis function only when using standard RBF interpolation (default = thin-plate spline)
            M_fineToCoarseInterpolant->setBasis ( list.get ("gammaf_interpolation_basis", "TPS" ) );
        }

        if ( !M_oneWayCoupling )
        {
            M_coarseToFineInterpolant -> setup ( M_fullSolidMesh, M_FSIoperator -> solidLocalMeshPtr(),
                                                 M_monodomain -> fullMeshPtr(), M_monodomain -> localMeshPtr(),
                                                 std::vector<int> (1, -1) );

            if (coarseToFineInterpolationType != "RBFlocallyRescaledVectorial" && coarseToFineInterpolationType != "RBFhtpVectorial")
            {
                // Set radius only when using standard RBF interpolation
                M_coarseToFineInterpolant -> setRadius ( (double) MeshUtility::MeshStatistics::computeSize ( * (M_fullSolidMesh ) ).maxH );
            }

            M_coarseToFineInterpolant -> setupRBFData ( M_solidDisplacement , M_monodomain -> displacementPtr(), dataFile, monodomainList);

            if (coarseToFineInterpolationType == "RBFvectorial")
            {
                // Set RBF basis function only when using standard RBF interpolation (default = thin-plate spline)
                M_coarseToFineInterpolant->setBasis ( list.get ("disp_interpolation_basis", "TPS" ) );
            }

            M_coarseToFineInterpolant -> buildOperators();
        }
    }
}



} // Namespace multiscale
} // Namespace LifeV
