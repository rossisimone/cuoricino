#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IntracellularCalciumGoldbeter.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>
#include <lifev/structure/solver/GeneralizedActiveHolzapfelOgdenMaterial.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/interpolation/RBFhtp.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <lifev/core/mesh/MeshTransformer.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledScalar.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledScalar.hpp>
//#include <lifev/core/interpolation/RBFscalar.hpp>
#include <lifev/core/interpolation/RBFvectorial.hpp>

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <sys/stat.h>

using namespace LifeV;

void EpetraPow ( VectorEpetra& vector, const Real p )
{
	Int size = vector.epetraVector().MyLength();

	for(int j(0); j < size; j++ )
	{
	int gid = vector.blockMap().GID(j);
	vector[gid] = std::pow(vector[gid],p);
	}
}
void EpetraSqrt ( VectorEpetra& vector )
{
	Int size = vector.epetraVector().MyLength();

	for(int j(0); j < size; j++ )
	{
	int gid = vector.blockMap().GID(j);
	vector[gid] = std::sqrt(vector[gid]);
	}
}


Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}
Real d0(const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real initialStimulus(const Real& /*t*/, const Real&  X, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
	if( X == 0 ) return 1.0;
	else return  0.;
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
    typedef IntracellularCalciumGoldbeter					ionicModel_Type;
    typedef boost::shared_ptr< ionicModel_Type >  ionicModelPtr_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type >        monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
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

    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
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

  	//===========================================================
  	//===========================================================
  	//				ELECTROPHYSIOLOGY
  	//===========================================================
  	//===========================================================

    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "ParamList1.xml" ) );
    Teuchos::ParameterList parameterList2 = * ( Teuchos::getParametersFromXmlFile ( "ParamList2.xml" ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");




    if ( comm->MyPID() == 0 )
        {
            std::cout << " Done!" << endl;
        }

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // Creates a new model object representing the//
    // model from Aliev and Panfilov 1996.  The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Constructor for ionic Model with parameters ... ";
    }
    ionicModelPtr_Type  ionicModel ( new ionicModel_Type() );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // set up the monodomain solver               //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    monodomainSolverPtr_Type monodomain ( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicModel ) );
    monodomainSolverPtr_Type monodomain2 ( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicModel ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }

    bool load4restart = parameterList.get("load4restart", false);
//    ionicModel -> initialize( monodomain -> globalSolution() );



	monodomain -> setInitialConditions();

//	function_Type f = &initialStimulus;
//	monodomain -> setPotentialFromFunction(f);


	monodomain2 -> setVariablePtr( monodomain -> globalSolution().at(0), 1 );
	monodomain2 -> setVariablePtr( monodomain -> globalSolution().at(1), 0 );

	function_Type f2 = &initialStimulus;
	monodomain2 -> setPotentialFromFunction(f2);


    std::cout << "Norm Inf potential = " <<  (  *( monodomain -> globalSolution().at(0) ) ).normInf() << std::endl;

    monodomain -> setParameters ( parameterList );
    monodomain2 -> setParameters ( parameterList2 );
    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > expElectro;

    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        cout << "\nExporter setup:  " ;
    }
    }
    monodomain -> setupExporter ( expElectro, parameterList.get ("ElectroOutputFile", "ElectroOutput") );
    expElectro.setPostDir ( problemFolder );
    if ( comm->MyPID() == 0 )
    {
        cout << "\nExport at 0:  " ;
    }

    monodomain -> exportSolution ( expElectro, 0.0 );

    if ( comm->MyPID() == 0 )
    {
        cout << "\nsolve system:  " ;
    }


    //********************************************//
    // Activation time						      //
    //********************************************//
    vectorPtr_Type activationTimeVector( new vector_Type( monodomain -> potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;

    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId(monodomain -> localMeshPtr(), monodomain -> commPtr() ->MyPID());
    activationTimeExporter.addVariable(ExporterData<mesh_Type>::ScalarField, "Activation Time",
    				monodomain -> feSpacePtr(), activationTimeVector, UInt(0));
    activationTimeExporter.setPrefix("ActivationTime");
    activationTimeExporter.setPostDir(problemFolder);
  	//===========================================================
  	//===========================================================
  	//				SOLID MECHANICS
  	//===========================================================
  	//===========================================================



    if ( comm->MyPID() == 0 )
    {
        std::cout << "monodomain: passed!" << std::endl;
    }

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


    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        std::cout << "\nparameters" << std::endl;
    }
    }
    Real rho, poisson, young, bulk, alpha, gammai, mu;
    rho     = dataFile ( "solid/physics/density", 1. );
    young   = dataFile ( "solid/physics/young",   1. );
    poisson = dataFile ( "solid/physics/poisson", 1. );
    bulk    = dataFile ( "solid/physics/bulk",    1. );
    alpha   = dataFile ( "solid/physics/alpha",   1. );
    gammai   = dataFile ( "solid/physics/gamma",   1. );
    mu   = dataFile ( "solid/physics/mu",   1. );
  //  M_gammaf  = dataFile ( "solid/physics/gammaf",  0. );

    if ( comm->MyPID() == 0 )
        {
    std::cout << "density = " << rho     << std::endl
              << "young   = " << young   << std::endl
              << "poisson = " << poisson << std::endl
              << "bulk    = " << bulk    << std::endl
              << "alpha   = " << alpha   << std::endl
              << "gamma   = " << gammai   << std::endl;
        }

    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        std::cout << "\ninitialization constitutive law" << std::endl;
    }
    }
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces" << std::endl;
    }

    meshPtr_Type fullSolidMesh;
    meshPtr_Type localSolidMesh;



    fullSolidMesh = monodomain -> fullMeshPtr();
    localSolidMesh = monodomain -> localMeshPtr();

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidFESpacePtr_Type aFESpace ( new solidFESpace_Type (monodomain -> localMeshPtr(), dOrder, 1, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
    scalarETFESpacePtr_Type aETFESpace ( new scalarETFESpace_Type (monodomain -> localMeshPtr(), & (aFESpace->refFE() ), & (aFESpace->fe().geoMap() ), comm) );
    solidFESpacePtr_Type solidaFESpace ( new solidFESpace_Type (localSolidMesh, "P1", 1, comm) );


    if ( comm->MyPID() == pid )
    {
        std::cout << "\nsetup boundary conditions" << std::endl;
    }
    bcInterfacePtr_Type                     solidBC( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );


    if ( comm->MyPID() == pid )
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

     if ( comm->MyPID() == pid )
     {
         std::cout << "\ninitial guess" << std::endl;
     }

     solid.setDataFromGetPot (dataFile);

    	//===========================================================
    	//===========================================================
    	//				FIBERS
    	//===========================================================
    	//===========================================================
     for(int pid(0); pid < 4 ; pid ++){
     if ( comm->MyPID() == pid )
     {
         std::cout << "\nreading fibers ... " << std::endl;
     }
     }


     vectorPtr_Type solidFibers( new vector_Type( dFESpace -> map() ) );


     std::vector<Real> fvec(3, 0.0);
     fvec.at(0)  = parameterList.get ("fiber_X", 1.0);
     fvec.at(1)  = parameterList.get ("fiber_Y", 0.0);
     fvec.at(2)  = parameterList.get ("fiber_Z", 0.0);
     HeartUtility::setupFibers(*solidFibers, fvec);

     MPI_Barrier(MPI_COMM_WORLD);


     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nread fibers" << std::endl;
     }

     solid.material() -> setFiberVector( *solidFibers );

//     monodomain -> setupFibers();

     vectorPtr_Type gammaf( new vector_Type( ( monodomain -> globalSolution().at(0) ) -> map() ) );
     vectorPtr_Type gammas( new vector_Type( ( monodomain -> globalSolution().at(0) ) -> map() ) );
     vectorPtr_Type gamman( new vector_Type( ( monodomain -> globalSolution().at(0) ) -> map() ) );
     vectorPtr_Type solidGammaf;
     vectorPtr_Type emDisp;
     solidFESpacePtr_Type electroFiberFESpace;
     solidETFESpacePtr_Type electrodETFESpace;

    solidGammaf = gammaf;
    monodomain -> setFiberPtr( solidFibers );
    monodomain2 -> setFiberPtr( solidFibers );
    emDisp = solid.displacementPtr();
    electroFiberFESpace = dFESpace;
    electrodETFESpace = dETFESpace;



     monodomain -> exportFiberDirection(problemFolder);
     //********************************************//
     // Create the global matrix: mass + stiffness in ELECTROPHYSIOLOGY //
     //********************************************//
     if ( comm->MyPID() == 0 )
     {
         cout << "\nSetup operators:  1 = " << monodomain -> timeStep() << "\n" ;
     }

     monodomain -> setDisplacementPtr( emDisp );
     monodomain -> setupMassMatrix();
     monodomain -> setupStiffnessMatrix();
     monodomain -> setupGlobalMatrix();


     if ( comm->MyPID() == 0 )
     {
         cout << "\nSetup operators:  2 = " << monodomain -> timeStep() << "\n" ;
     }
     monodomain2 -> setDisplacementPtr( emDisp );
     monodomain2 -> setupMassMatrix();
     monodomain2 -> setupStiffnessMatrix();
     monodomain2 -> setupGlobalMatrix();

     if ( comm->MyPID() == 0 )
     {
         cout << "Done! \n" ;
     }

     //==================================================================//
     //==================================================================//
     //					SETUP Activation								//
     //==================================================================//
     //==================================================================//
     *gammaf *= 0.0;
     *solidGammaf =0.0;

     solid.material() -> setGammaf( *solidGammaf );

     vectorPtr_Type solidGammas( new vector_Type( solidGammaf -> map() ) );
     vectorPtr_Type solidGamman( new vector_Type( solidGammaf -> map() ) );

	 *solidGammas = 1.0;
	 *solidGammas /= (1.0 + *solidGammaf);
	 EpetraSqrt(*solidGammas);
	 *solidGammas -= 1.0;
	 solid.material() -> setGamman(*solidGammas);
	 *solidGamman = *solidGammas;
	 solid.material() -> setGammas(*solidGamman);

	 gamman = solidGamman;
	 gammas = solidGammas;

     //==================================================================//
     //==================================================================//
     //					SETUP INTERPOLATION								//
     //==================================================================//
     //==================================================================//

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
      exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, parameterList.get ("StructureOutputFile", "StructureOutput") ) );

      //      exporter->setPostDir ( "./" );
            exporter -> setPostDir ( problemFolder );
      exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );

      vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
      exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
      exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "solid_gammaf", solidaFESpace, solidGammaf, UInt (0) );





      //================================================================//
      //================================================================//
      //					SETUP COUPLING SOLVER						//
      //																//
      //================================================================//
      //================================================================//
      ExporterHDF5< RegionMesh <LinearTetra> > expGammaf;
	expGammaf.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
	expGammaf.setPrefix(parameterList.get ("ActivationOutputFile", "ActivationOutput"));
	expGammaf.setPostDir ( problemFolder );

    matrixPtr_Type mass(new matrix_Type( monodomain -> massMatrixPtr() -> map() ) ) ;

	boost::shared_ptr<FLRelationship> fl (new FLRelationship);

	boost::shared_ptr<HeavisideFct> H (new HeavisideFct);

	boost::shared_ptr<Psi4f> psi4f (new Psi4f);
	boost::shared_ptr<ShowValue> sv(new ShowValue);

      MatrixSmall<3,3> Id;
      Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
      Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
      Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;
      VectorSmall<3> E1;
      E1(0) = 1.;
      E1(1) = 0.;
      E1(2) = 0.;
  	{
  		using namespace ExpressionAssembly;

  		integrate(elements(monodomain -> localMeshPtr() ), monodomain -> feSpacePtr() -> qr(), monodomain -> ETFESpacePtr(),
  				monodomain -> ETFESpacePtr(), phi_i * phi_j) >> mass;

  	}
  	mass -> globalAssemble();


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
	linearSolver.setOperator( mass );
//	linearSolver.setOperator( monodomain -> massMatrixPtr() );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nIt does!!!!" << std::endl;
    }





   	vectorPtr_Type tmpRhsActivation( new vector_Type ( rhsActivation -> map(), Repeated ) );
    solidFESpacePtr_Type emDispFESpace ( new solidFESpace_Type ( monodomain -> localMeshPtr(), "P1", 3, comm) );
  	expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
  			monodomain -> feSpacePtr(), gammaf, UInt(0));
  	expGammaf.addVariable(ExporterData<mesh_Type>::VectorField, "interpolated displacement",
  			emDispFESpace, emDisp, UInt(0));
  	expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "rhs",
  			monodomain -> feSpacePtr(), rhsActivation, UInt(0));


  	expGammaf.postProcess(0.0);

    //===========================================================
  	//===========================================================
  	//				Initializing solid
  	//===========================================================
  	//===========================================================


    exporter->postProcess ( 0 );


	vectorPtr_Type emDisp0(new vector_Type( emDisp -> map() ) );
    *emDisp0 = *emDisp;

    //===========================================================
  	//===========================================================
  	//				TIME LOOP
  	//===========================================================
  	//===========================================================
      Real emdt = parameterList.get("emdt",1.0);
      int iter((emdt / monodomain -> timeStep()));
      int k(0);
      Real saveStep = parameterList.get("save_step",1.0);
      int saveIter((saveStep / monodomain -> timeStep()));
      Real meth = parameterList.get("meth",1.0);
      int subiter = parameterList.get("subiter",100);

        Real dt_min = 0.01;


		Real Ca_diastolic = dataFile( "solid/physics/Ca_diastolic", 0.02155 );




     for( Real t(0.0); t< monodomain -> endTime(); )
	 {
		  t = t + monodomain -> timeStep();
		  k++;

         LifeChrono timer;

    	 for(int j(0); j<subiter; j++) monodomain -> solveOneReactionStepFE(subiter);

          timer.stop();


		if ( comm->MyPID() == 0 )
		{
			std::cout << "\nSolve DIFFUSION step with BE!\n" << std::endl;
		}

	    (*monodomain -> rhsPtrUnique()) *= 0.0;
	      monodomain -> updateRhs();
          monodomain -> solveOneDiffusionStepBE();

          (*monodomain2 -> rhsPtrUnique()) *= 0.0;
  	      monodomain2 -> updateRhs();
            monodomain2 -> solveOneDiffusionStepBE();




		  *tmpRhsActivation *= 0;
			if ( comm->MyPID() == 0 )
			{
				std::cout << "\nASSEMBLING ACTIVATION EQUATION!\n" << std::endl;
			}


					{
						using namespace ExpressionAssembly;
//
//
//						BOOST_AUTO_TPL(I,      value(Id) );
//						BOOST_AUTO_TPL(Grad_u, grad( electrodETFESpace, *emDisp, 0) );
//						BOOST_AUTO_TPL(F,      ( Grad_u + I ) );
//						BOOST_AUTO_TPL(J,       det(F) );
//						BOOST_AUTO_TPL(Jm23,    pow(J, -2./3));
//						BOOST_AUTO_TPL(I1,     dot(F, F));
//
//						// Fibres
//						BOOST_AUTO_TPL(f0,     value( electrodETFESpace, *( monodomain -> fiberPtr() ) ) );
//						BOOST_AUTO_TPL(f,      F * f0 );
//						BOOST_AUTO_TPL(I4f,    dot(f, f) );
//
//						BOOST_AUTO_TPL(I1iso,   Jm23 * I1);
//						BOOST_AUTO_TPL(I4fiso,  Jm23 * I4f);
//						// Generalised invariants
//					    BOOST_AUTO_TPL(gf,  value(aETFESpace, *gammaf));
//					    BOOST_AUTO_TPL(gs,  value(aETFESpace, *gammas));
//					    BOOST_AUTO_TPL(gn,  value(aETFESpace, *gamman));
//
//
//
//				        BOOST_AUTO_TPL(dI1edI1,   value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  )    );
//				        BOOST_AUTO_TPL(dI1edI4f,  value(1.0) / (   (gf + value(1.0))  *  (gf + value(1.0))  ) -  value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  ) );
//				        BOOST_AUTO_TPL(dI4fedI4f, value(1.0) / (   (gf + value(1.0)) *   (gf + value(1.0))  )    );
//
//					    BOOST_AUTO_TPL(I1eiso,   dI1edI1 * I1iso + dI1edI4f * I4fiso );
//					    BOOST_AUTO_TPL(I4feiso,  dI4fedI4f * I4fiso);
//						BOOST_AUTO_TPL(I4feisom1, ( I4feiso - value(1.0) ) );
//
//						Real A = dataFile( "solid/physics/a", 4960 );
//						Real B = dataFile( "solid/physics/b_activation", 0. );
//						Real Af = dataFile( "solid/physics/af", 0. );
//						Real Bf = dataFile( "solid/physics/bf", 0. );
//
//						BOOST_AUTO_TPL(a, value( A ) );
//						BOOST_AUTO_TPL(b, value(  B ) );
//						BOOST_AUTO_TPL(af, value( Af ) );
//						BOOST_AUTO_TPL(bf, value(  Bf ) );
//
//						Real viscosity = dataFile( "solid/physics/viscosity", 0.0005 );
//						BOOST_AUTO_TPL(beta, value(viscosity ) /*/ coeff  /* Jm23 * pow( eval(EXP, ( I1iso + value(-3.0) ) ), -B )*/ );
//
////						BOOST_AUTO_TPL(dWs, a * pow(g, -1) );
//						BOOST_AUTO_TPL(gamma_dot, beta / ( Ca2 ) * ( Pa - dW )  );

						BOOST_AUTO_TPL(Ca ,   value( aETFESpace, *( monodomain -> globalSolution().at(0)  ) ) );
						BOOST_AUTO_TPL(Gammaf,	 value( aETFESpace, *gammaf )  );
					//Activation as in the IUTAM proceedings
						BOOST_AUTO_TPL(activationEquation, value(-0.02) *Ca + value(-0.04)*Gammaf  );

						integrate ( elements ( monodomain -> localMeshPtr() ),
								monodomain -> feSpacePtr() -> qr() ,
								monodomain -> ETFESpacePtr(),
								activationEquation  * phi_i
						) >> tmpRhsActivation;

					}


					*rhsActivation *= 0;
					*rhsActivation = ( *(mass) * ( *gammaf ) );
					*rhsActivation += ( ( monodomain -> timeStep() * *tmpRhsActivation ) );

					linearSolver.setRightHandSide(rhsActivation);



					if ( comm->MyPID() == 0 )
					{
						std::cout << "\nSOLVING ACTIVATION EQUATION!\n" << std::endl;
					}


					linearSolver.solve(gammaf);




			  if ( k % iter == 0)
			  {
				  solidGammaf = gammaf;

					solid.material() -> setGammaf( *solidGammaf );

						 *solidGammas = 1.0;
						 *solidGammas /= (1.0 + *solidGammaf);
						 EpetraSqrt(*solidGammas);
						 *solidGammas -= 1.0;
						 solid.material() -> setGamman(*solidGammas);
						 *solidGamman = *solidGammas;
						 solid.material() -> setGammas(*solidGamman);



					if ( comm->MyPID() == 0 )
					{
						std::cout << "\nSOLVING STATIC MECHANICS!\n" << std::endl;
					}

					if ( comm->MyPID() == 0 )
					{
						std::cout << "\n*****************************************************";
						std::cout << "\nWE ARE AT TIME: "<< t;
						std::cout << "\n*****************************************************";

					}

							solid.iterate ( solidBC -> handler() );
							*solidDisp = solid.displacement();





							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nREASSEMBLING STIFFNESS MATRIX FOR TOW WAY COUPLING!\n" << std::endl;
							}

							 monodomain -> setupStiffnessMatrix();
							 monodomain -> setupGlobalMatrix();

			  }

		  //cout << "\n\n save every " << saveIter << "iteration\n";
		  if ( k % saveIter == 0)
		  {
			  monodomain -> registerActivationTime(*activationTimeVector, t, 1.0);
			  monodomain -> exportSolution(expElectro, t);
			  expGammaf.postProcess(t);
			  exporter->postProcess ( t );
		  }




      }


      expElectro.closeFile();
      expGammaf.closeFile();

      exporter -> closeFile();

    activationTimeExporter.postProcess(0);
    activationTimeExporter.closeFile();



    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nActive strain example: Passed!" << std::endl;
    }




#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}



