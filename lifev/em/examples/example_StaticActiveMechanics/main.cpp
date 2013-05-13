#include <lifev/core/LifeV.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>

using namespace LifeV;


Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}
Real d0(const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real initialVhumanSphere(const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{

  double r = std::sqrt(pow(X-6.25,2)+pow(Y-3.75,2)+pow(Z-7.25,2));
  double auxexp=1.0-1.0/(1.0+exp(-90.0*(r-1.7)));

  return auxexp;
}


int main (int argc, char** argv)
{

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef VectorEpetra				vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef  StructuralOperator< RegionMesh<LinearTetra> >		physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >              bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif


  	//===========================================================
  	//===========================================================
  	//				SOLID MECHANICS
  	//===========================================================
  	//===========================================================


    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;
    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;
    if ( comm->MyPID() == 0 )
    {
        std::cout << "\n\ninitialization bc handler" << std::endl;
    }


    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nparameters" << std::endl;
    }

    Real rho, poisson, young, bulk, alpha, gamma, mu;
    rho     = dataFile ( "solid/physics/density", 1. );
    young   = dataFile ( "solid/physics/young",   1. );
    poisson = dataFile ( "solid/physics/poisson", 1. );
    bulk    = dataFile ( "solid/physics/bulk",    1. );
    alpha   = dataFile ( "solid/physics/alpha",   1. );
    gamma   = dataFile ( "solid/physics/gamma",   1. );
    mu   = dataFile ( "solid/physics/mu",   1. );
  //  M_gammaf  = dataFile ( "solid/physics/gammaf",  0. );

    if ( comm->MyPID() == 0 )
        {
    std::cout << "density = " << rho     << std::endl
              << "young   = " << young   << std::endl
              << "poisson = " << poisson << std::endl
              << "bulk    = " << bulk    << std::endl
              << "alpha   = " << alpha   << std::endl
              << "gamma   = " << gamma   << std::endl;
        }


    if ( comm->MyPID() == 0 )
    {
        std::cout << "\ninitialization constitutive law" << std::endl;
    }

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces" << std::endl;
    }
    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (meshPart, dOrder, 3, parameters->comm ));
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup boundary conditions" << std::endl;
    }

     //! ========================================================================
     //! Setting solid BCs from data file
	bcInterfacePtr_Type                     solidBC( new bcInterface_Type() );
	solidBC->createHandler();
	solidBC->fillHandler ( data_file_name, "solid" );



    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup structural operator" << std::endl;
    }
    //! 1. Constructor of the structuralSolver
     StructuralOperator< RegionMesh<LinearTetra> > solid;
     solid.setup (dataStructure,
                  dFESpace,
                  dETFESpace,
                  solidBC -> handler(),
                  comm);
     if ( comm->MyPID() == 0 )
     {
         std::cout << "\ninitial guess" << std::endl;
     }

     solid.setDataFromGetPot (dataFile);

     vectorPtr_Type fibers( new vector_Type( dFESpace -> map() ) );
     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nread fibers" << std::endl;
     }

     HeartUtility::importFibers(fibers, parameterList.get ("fiber_file", ""), meshPart );

     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nset fibers" << std::endl;
     }


     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nset gammaf and fibers" << std::endl;
     }
     solid.material() -> setGammaf( *( monodomain -> globalSolution().at(3) ) );

     solid.material() -> setFiberVector( * ( monodomain -> fiberPtr() ) );

     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nbuild solid system" << std::endl;
     }

     solid.buildSystem(1.0);
     vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
     vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
     vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );
     solid.initialize ( initialDisplacement );


     MPI_Barrier (MPI_COMM_WORLD);

     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nsetup solid exporter" << std::endl;
     }

      boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
      exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );

      exporter->setPostDir ( "./" );
      exporter->setMeshProcId ( monodomain -> localMeshPtr(), comm->MyPID() );

      vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
      exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
      exporter->postProcess ( 0 );



      //================================================================//
      //================================================================//
      //					SETUP COUPLING SOLVER						//
      //																//
      //================================================================//
      //================================================================//
      ExporterHDF5< RegionMesh <LinearTetra> > expGammaf;
	expGammaf.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
	expGammaf.setPrefix("gammaf");

      vectorPtr_Type gammaf( new vector_Type( monodomain -> globalSolution().at(3) -> map() ) );
      *gammaf *= 0;
  	vectorPtr_Type rhsActivation( new vector_Type( *gammaf ) );
  	*rhsActivation *= 0;




      if ( comm->MyPID() == 0 )
      {
          std::cout << "\nSolve system" << std::endl;
      }

      //==================================================================//
      //==================================================================//
      //					SETUP LINEAR SOLVER	    						//
      //						ACTIVATION									//
      //==================================================================//
      //==================================================================//

      if ( comm->MyPID() == 0 )
      {
          std::cout << "\nset up linear solver... Does it work???" << std::endl;
      }
  	typedef LinearSolver linearSolver_Type;
  	typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;
	typedef LifeV::Preconditioner basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack prec_Type;
	typedef boost::shared_ptr<prec_Type> precPtr_Type;


	prec_Type* precRawPtr;
	basePrecPtr_Type precPtr;
	precRawPtr = new prec_Type;
	precRawPtr->setDataFromGetPot(dataFile, "prec");
	precPtr.reset(precRawPtr);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nprec done!!!!" << std::endl;
    }

	Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp(
			new Teuchos::ParameterList);

	std::string xmlpath = dataFile("electrophysiology/monodomain_xml_path",
			"./");
	std::string xmlfile = dataFile("electrophysiology/monodomain_xml_file",
			"MonodomainSolverParamList.xml");

	solverParamList = Teuchos::getParametersFromXmlFile(xmlpath + xmlfile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nreading file done!!!!" << std::endl;
    }

	linearSolver_Type linearSolver;
    linearSolver.setCommunicator ( comm );
    linearSolver.setParameters ( *solverParamList );
    linearSolver.setPreconditioner ( precPtr );
	linearSolver.setOperator( monodomain -> massMatrixPtr() );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nIt does!!!!" << std::endl;
    }

  	//===========================================================
  	//				TIME LOOP
  	//===========================================================
 
   	boost::shared_ptr<FLRelationship> fl (new FLRelationship);
#define Ca    ( value( aETFESpace, *( monodomain -> globalSolution().at(3)  ) ) )
#define Gammaf 			( value( aETFESpace, *gammaf ) )

	//Activation as in the IUTAM proceedings
#define activationEquation value(-0.02) *Ca + value(-0.04)*Gammaf  



   	vectorPtr_Type tmpRhsActivation( new vector_Type ( rhsActivation -> map(), Repeated ) );

  	expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
  			monodomain -> feSpacePtr(), gammaf, UInt(0));
  	expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "rhs",
  			monodomain -> feSpacePtr(), rhsActivation, UInt(0));


  	expGammaf.postProcess(0.0);



     for( Real t(0.0); t< monodomain -> endTime(); )
	 {
		  t = t + monodomain -> timeStep();
		  k++;
		  monodomain -> solveOneSplittingStep( exp, t );

		  *tmpRhsActivation *= 0;
		  	{
		  		using namespace ExpressionAssembly;
				integrate ( elements ( monodomain -> localMeshPtr() ),
						monodomain -> feSpacePtr() -> qr() ,
						monodomain -> ETFESpacePtr(),
						activationEquation * phi_i
				) >> tmpRhsActivation;

		  	}
			*rhsActivation *= 0;
		  	*rhsActivation = ( *(monodomain -> massMatrixPtr() ) * ( *gammaf ) );
		  	*rhsActivation += ( monodomain -> timeStep() * *tmpRhsActivation );

			linearSolver.setRightHandSide(rhsActivation);
			linearSolver.solve(gammaf);

			  solid.material() -> setGammaf( *gammaf );
			  solid.iterate ( solidBC -> handler() );
			  *solidDisp = solid.displacement();

		  expGammaf.postProcess(t);

		  exporter->postProcess ( t );
      }
      exp.closeFile();
      expGammaf.closeFile();

      exporter -> closeFile();
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Active strain example: Passed!" << std::endl;
    }

#undef Gammaf
#undef deformationGradientTensor
#undef RIGHTCAUCHYGREEN
#undef firstInvariantC
#undef fiber0
#undef fiber
#undef Pa
#undef I4f
#undef beta
#undef GammaPlusOne
#undef dgGammaf
#undef activationEquation

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
