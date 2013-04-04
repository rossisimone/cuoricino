////@HEADER
///*
//*******************************************************************************
//
//    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
//    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
//
//    This file is part of LifeV.
//
//    LifeV is free software; you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    LifeV is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.
//
//*******************************************************************************
//*/
////@HEADER
///**
//   \file main.cpp
//   \author Paolo Tricerri <paolo.tricerri@epfl.ch>
//
//
//   Attention: At the moment the restart works only if the solution is saved at
//   each time step
//   \date 2005-04-16
// */
//#undef HAVE_HDF5
//#ifdef TWODIM
//#error test_structure cannot be compiled in 2D
//#endif
//
//// Tell the compiler to ignore specific kind of warnings:
//#pragma GCC diagnostic ignored "-Wunused-variable"
//#pragma GCC diagnostic ignored "-Wunused-parameter"
//
//#include <Epetra_ConfigDefs.h>
//#ifdef EPETRA_MPI
//#include <mpi.h>
//#include <Epetra_MpiComm.h>
//#else
//#include <Epetra_SerialComm.h>
//#endif
//
////Tell the compiler to restore the warning previously silented
//#pragma GCC diagnostic warning "-Wunused-variable"
//#pragma GCC diagnostic warning "-Wunused-parameter"
//
//#include <lifev/core/LifeV.hpp>
//#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
//#include <lifev/core/algorithm/PreconditionerML.hpp>
//
//
////Include fils which were in the structure.cpp file
//#include <lifev/core/array/MapEpetra.hpp>
//
//#include <lifev/core/fem/TimeAdvance.hpp>
//#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
//#include <lifev/core/fem/TimeAdvanceBDF.hpp>
//
//#include <lifev/core/mesh/MeshData.hpp>
//#include <lifev/core/mesh/MeshPartitioner.hpp>
//
//#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
//
//#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
//#include <lifev/structure/solver/StructuralOperator.hpp>
//#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
//#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
//#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>
//#include <lifev/structure/solver/ExponentialMaterialNonLinear.hpp>
//
//#include <lifev/core/filter/Exporter.hpp>
//#include <lifev/core/filter/ExporterEnsight.hpp>
//#ifdef HAVE_HDF5
//#include <lifev/core/filter/ExporterHDF5.hpp>
//#endif
//#include <lifev/core/filter/ExporterEmpty.hpp>
//#include <lifev/core/filter/ExporterVTK.hpp>
//
//#include <iostream>
//#include <lifev/core/LifeV.hpp>
//
//using namespace LifeV;
//
//
//int returnValue = EXIT_FAILURE;
//enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };
//
//namespace
//{
//static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
//static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
//}
//
//std::set<UInt> parseList ( const std::string& list )
//{
//    std::string stringList = list;
//    std::set<UInt> setList;
//    if ( list == "" )
//    {
//        return setList;
//    }
//    size_t commaPos = 0;
//    while ( commaPos != std::string::npos )
//    {
//        commaPos = stringList.find ( "," );
//        setList.insert ( atoi ( stringList.substr ( 0, commaPos ).c_str() ) );
//        stringList = stringList.substr ( commaPos + 1 );
//    }
//    setList.insert ( atoi ( stringList.c_str() ) );
//    return setList;
//}
//
//
//class Structure
//{
//public:
//
//    typedef RegionMesh<LifeV::LinearTetra >                       mesh_Type;
//    typedef StructuralOperator<mesh_Type >::vector_Type             vector_Type;
//    typedef boost::shared_ptr<vector_Type>                        vectorPtr_Type;
//    typedef boost::shared_ptr< TimeAdvance< vector_Type > >       timeAdvance_type;
//
//    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
//    //Exporters Typedefs
//    typedef LifeV::Exporter<mesh_Type >                  filter_Type;
//    typedef boost::shared_ptr<filter_Type >                       filterPtr_Type;
//
//    typedef LifeV::ExporterEmpty<mesh_Type >                      emptyExporter_Type;
//    typedef boost::shared_ptr<emptyExporter_Type>                 emptyExporterPtr_Type;
//
//    typedef LifeV::ExporterEnsight<mesh_Type >                    ensightFilter_Type;
//    typedef boost::shared_ptr<ensightFilter_Type>                 ensightFilterPtr_Type;
//
//    typedef LifeV::ExporterData<mesh_Type>                        exporterData_Type;
//    typedef boost::shared_ptr<exporterData_Type>                  exporterDataPtr_Type;
//
//    typedef LifeV::ExporterVTK<mesh_Type>                         exporterVTK_Type;
//    typedef boost::shared_ptr<exporterVTK_Type>                   exporterVTKPtr_Type;
//
//#ifdef HAVE_HDF5
//    typedef LifeV::ExporterHDF5<mesh_Type >                       hdf5Filter_Type;
//    typedef boost::shared_ptr<hdf5Filter_Type>                    hdf5FilterPtr_Type;
//#endif
//
//    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >         solidFESpace_Type;
//    typedef boost::shared_ptr<solidFESpace_Type>                  solidFESpacePtr_Type;
//
//
//    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >      solidETFESpace_Type;
//    typedef boost::shared_ptr<solidETFESpace_Type>                 solidETFESpacePtr_Type;
//
//    /** @name Constructors, destructor
//     */
//    //@{
//    Structure ( int                                   argc,
//                char**                                argv,
//                boost::shared_ptr<Epetra_Comm>        structComm );
//
//    ~Structure()
//    {}
//    //@}
//
//    //@{
//    void run()
//    {
//        run3d();
//    }
//    //@}
//
//protected:
//
//private:
//
//    /**
//     * run the 2D driven cylinder simulation
//     */
//    void run2d();
//
//    /**
//     * run the 3D driven cylinder simulation
//     */
//    void run3d();
//
//private:
//    struct Private;
//    boost::shared_ptr<Private> parameters;
//    filterPtr_Type importerSolid;
//    filterPtr_Type exporterSolid;
//
//    filterPtr_Type exporterCheck;
//
//    solidFESpacePtr_Type dFESpace;
//    solidETFESpacePtr_Type dETFESpace;
//    solidFESpacePtr_Type exporterFESpace;
//
//};
//
//
//
//struct Structure::Private
//{
//    Private() :
//        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
//    {}
//    double rho, poisson, young, bulk, alpha, gamma;
//
//    std::string data_file_name;
//
//    boost::shared_ptr<Epetra_Comm>     comm;
//
//    static Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
//    {
//        return  0.;
//    }
//
//    static Real bcNonZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& i)
//    {
//        switch (i)
//        {
//            case 0:
//                return  0.0;
//                break;
//            case 1:
//                return  300000.0;
//                break;
//            case 2:
//                return 0.0;
//                break;
//
//        }
//        return 0;
//    }
//
//    static Real bcPressure (const Real& /*t*/, const Real&  x, const Real& y, const Real& /*Z*/, const ID& i)
//    {
//        Real radius = 0.5;
//        Real pressure = 5000;
//        //3000/(2*0.5*40*3.1415962);
//        switch (i)
//        {
//            case 0:
//                return  pressure * std::fabs ( ( x / radius ) );
//                break;
//            case 1:
//                return  pressure * std::fabs ( ( y / radius ) );
//                break;
//            case 2:
//                return 0.0;
//                break;
//
//        }
//        return 0;
//
//    }
//
//    static Real pressureUsingNormal (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
//    {
//
//        return -5000;
//        // if( t < 15.0 )
//        //   return  -(300000/(2*3.1415962*0.5*40))*(1/15)*t;
//        // else
//        //   return  -300000/(2*3.1415962*0.5*40);
//    }
//
//
//    static Real gf (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
//    {
//        return  (- 0.3 * (z - 0.5 ) * ( z + 0.5 ) );
//
//    }
//
//};
//
//
//
//Structure::Structure ( int                                   argc,
//                       char**                                argv,
//                       boost::shared_ptr<Epetra_Comm>        structComm) :
//    parameters ( new Private() ),
//    dFESpace(),
//    dETFESpace(),
//    exporterFESpace()
//{
//    GetPot command_line (argc, argv);
//    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
//    GetPot dataFile ( data_file_name );
//    parameters->data_file_name = data_file_name;
//
//    parameters->rho     = dataFile ( "solid/physics/density", 1. );
//    parameters->young   = dataFile ( "solid/physics/young",   1. );
//    parameters->poisson = dataFile ( "solid/physics/poisson", 1. );
//    parameters->bulk    = dataFile ( "solid/physics/bulk",    1. );
//    parameters->alpha   = dataFile ( "solid/physics/alpha",   1. );
//    parameters->gamma   = dataFile ( "solid/physics/gamma",   1. );
//
//    std::cout << "density = " << parameters->rho     << std::endl
//              << "young   = " << parameters->young   << std::endl
//              << "poisson = " << parameters->poisson << std::endl
//              << "bulk    = " << parameters->bulk    << std::endl
//              << "alpha   = " << parameters->alpha   << std::endl
//              << "gamma   = " << parameters->gamma   << std::endl;
//
//    parameters->comm = structComm;
//    int ntasks = parameters->comm->NumProc();
//
//    if (!parameters->comm->MyPID() )
//    {
//        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
//    }
//}
//
//
//
//void
//Structure::run2d()
//{
//    std::cout << "2D cylinder test case is not available yet\n";
//}
//
//
//
//void
//Structure::run3d()
//{
//
//    bool verbose = (parameters->comm->MyPID() == 0);
//
//    //! Number of boundary conditions for the velocity and mesh motion
//    boost::shared_ptr<BCHandler> BCh ( new BCHandler() );
//
//    //! dataElasticStructure
//    GetPot dataFile ( parameters->data_file_name.c_str() );
//    cout << "\n\nCIPPALIPPA 1!!!\n\n";
//    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
//    cout << "\n\nCIPPALIPPA 2!!!\n\n";
//    dataStructure->setup (dataFile);
//    cout << "\n\nCIPPALIPPA 3!!!\n\n";
//    MeshData             meshData;
//    meshData.setup (dataFile, "solid/space_discretization");
//
//    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra> ( ( parameters->comm ) ) );
//    readMesh (*fullMeshPtr, meshData);
//
//    //fullMeshPtr->showMe( );
//
//    MeshPartitioner< RegionMesh<LinearTetra> > meshPart ( fullMeshPtr, parameters->comm );
//
//    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
//
//    dFESpace.reset ( new solidFESpace_Type (meshPart, dOrder, 3, parameters->comm) );
//    dETFESpace.reset ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );
//
//    // if( dOrder.compare("P2") == 0 )
//    //     exporterFESpace.reset( new solidFESpace_Type(meshPart,"P1",3,parameters->comm) );
//
//    if (verbose)
//    {
//        std::cout << std::endl;
//    }
//
//    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "Newmark");
//
//    timeAdvance_type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );
//
//    UInt OrderDev = 2;
//
//    //! initialization of parameters of time Advance method:
//    if (timeAdvanceMethod == "Newmark")
//    {
//        timeAdvance->setup ( dataStructure->dataTimeAdvance()->coefficientsNewmark() , OrderDev);
//    }
//
//    if (timeAdvanceMethod == "BDF")
//    {
//        timeAdvance->setup (dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);
//    }
//
//    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
//    //    timeAdvance->showMe();
//
//    //! #################################################################################
//    //! BOUNDARY CONDITIONS
//    //! #################################################################################
//    vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
//    compx[0] = 0;
//    compy[0] = 1, compz[0] = 2;
//    compxy[0] = 0;
//    compxy[1] = 1;
//    compxz[0] = 0;
//    compxz[1] = 2;
//    compyz[0] = 1;
//    compyz[1] = 2;
//
//    BCFunctionBase zero (Private::bcZero);
//    BCFunctionBase nonZero (Private::bcNonZero);
//    BCFunctionBase pressure (Private::bcPressure);
//    BCFunctionBase pressureNormal (Private::pressureUsingNormal);
//
//
//    //! =================================================================================
//    //! BC for cookMembraneMedit.mesh
//    //! =================================================================================
//    //Condition for Extension
//    //This conditions have to be changed according to script used to generate the mesh
//    //cookMembraneStructure.geo or cookMembrane.geo
//    BCh->addBC ("EdgesIn",      40,  Essential, Full, zero, 3);
//    BCh->addBC ("EdgesIn",      60,  Natural, Full, zero, 3);
//    BCh->addBC ("EdgesIn",      20,  Natural, Full, zero, 3);
//    BCh->addBC ("EdgesIn",      70,  Natural, Full, zero, 3);
//    BCh->addBC ("EdgesIn",      30,  Natural, Full, zero, 3);
//    BCh->addBC ("EdgesIn",      50,  Natural, Full, zero, 3);
//
//    // BCh->addBC("EdgesIn",      40,  Essential, Component, zero, compz);
//    // BCh->addBC("EdgesIn",      20,  Essential, Component, zero, compz);
//    // BCh->addBC("EdgesIn",      30,  Essential, Component, zero, compz);
//    // BCh->addBC("EdgesIn",      70,  Essential, Component, zero, compz);
//
//    //! 1. Constructor of the structuralSolver
//    StructuralOperator< RegionMesh<LinearTetra> > solid;
//
//    //! 2. Setup of the structuralSolver
//    solid.setup (dataStructure,
//                 dFESpace,
//                 dETFESpace,
//                 BCh,
//                 parameters->comm);
//
//
//    solidFESpacePtr_Type gammaFESpace( new solidFESpace_Type (meshPart, dOrder, 1, parameters->comm) );
//    VectorEpetra gammaf( gammaFESpace -> map() );
//
//    Structure::fct_type fg = &(Private::gf);
//    gammaFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( fg ), gammaf , 0);
//
//    //! 3. Setting data from getPot
//    solid.setDataFromGetPot (dataFile);
//
//    solid.material() -> setGammaf(gammaf);
//    //! 4. Building system using TimeAdvance class
//    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
//    solid.buildSystem (timeAdvanceCoefficient);
//
//    if (verbose)
//    {
//        std::cout << "S- initialization ... ";
//    }
//
//    //Initialization of TimeAdvance
//    std::string const restart =  dataFile ( "importer/restart", "none");
//    std::vector<vectorPtr_Type> solutionStencil;
//
//    if ( restart.compare ( "none" ) )
//    {
//        //Reading fileNames - setting data for reading
//        std::string const importerType =  dataFile ( "importer/type", "ensight");
//        std::string const fileName     =  dataFile ( "importer/filename", "structure");
//        std::string const initialLoaded     =  dataFile ( "importer/initialSol", "NO_DEFAULT_VALUE");
//        //LifeV::Real initialTime        =  dataFile ( "importer/initialTime", 0.0);
//
//        //Creating the importer for the case P2 if vtk is used.
//        //For other exporters we have to go from P1 to P2
//        // if( dataStructure->order().compare("P2") == 0 )
//        // {
//        //     importerSolid.reset( new exporterVTK_Type ( dataFile, fileName ) );
//        // }
//
//        //Creating the importer
//#ifdef HAVE_HDF5
//        if ( !importerType.compare ("hdf5") )
//        {
//            importerSolid.reset ( new  hdf5Filter_Type ( dataFile, fileName) );
//        }
//        else
//#endif
//        {
//            if ( !importerType.compare ("none") )
//            {
//                importerSolid.reset ( new emptyExporter_Type ( dataFile, dFESpace->mesh(), "solid", dFESpace->map().comm().MyPID() ) );
//            }
//            else
//            {
//                importerSolid.reset ( new  ensightFilter_Type ( dataFile, fileName) );
//            }
//        }
//        importerSolid->setMeshProcId (dFESpace->mesh(), dFESpace->map().comm().MyPID() );
//
//        //Creation of Exporter to check the loaded solution (working only for HDF5)
//        // std::string expVerFile = "verificationDisplExporter";
//        // LifeV::ExporterHDF5<RegionMesh<LinearTetra> > exporter( dataFile, meshPart.meshPartition(), expVerFile, parameters->comm->MyPID());
//        // vectorPtr_Type vectVer ( new vector_Type(solid.displacement(),  LifeV::Unique ) );
//
//        // exporter.addVariable( ExporterData<mesh_Type >::VectorField, "displVer", dFESpace, vectVer, UInt(0) );
//
//        // exporter.postProcess(0.0);
//
//        //Reading the displacement field and the timesteps need by the TimeAdvance class
//        vectorPtr_Type solidDisp;
//
//        solidDisp.reset (new vector_Type (dFESpace->map(), Unique ) );
//
//        std::string iterationString;
//
//        std::cout << "size TimeAdvance:" << timeAdvance->size() << std::endl;
//
//        //Loading the stencil
//        exporterDataPtr_Type solidDataReader;
//        iterationString = initialLoaded;
//        for (UInt iterInit = 0; iterInit < timeAdvance->size(); iterInit++ )
//        {
//            *solidDisp *= 0.0;
//
//            solidDataReader.reset (new exporterData_Type (exporterData_Type::VectorField, std::string ("displacement." + iterationString), dFESpace, solidDisp, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime) );
//
//            importerSolid->readVariable ( *solidDataReader );
//
//            std::cout << "Norm of the " << iterInit + 1 << "-th solution : " << solidDisp->norm2() << std::endl;
//
//            //Exporting the just loaded solution (debug purposes)
//            // Real currentLoading(iterInit + 1.0);
//            // *vectVer = *solidDisp;
//            // exporter.postProcess( currentLoading )
//
//            solutionStencil.push_back ( solidDisp );
//
//
//            //initializing the displacement field in the StructuralSolver class with the first solution
//            if ( !iterInit )
//            {
//                solid.initialize ( solidDisp );
//            }
//
//            //Updating string name
//            int iterations = std::atoi (iterationString.c_str() );
//            iterations--;
//
//            std::ostringstream iter;
//            iter.fill ( '0' );
//            iter << std::setw (5) << ( iterations );
//            iterationString = iter.str();
//
//        }
//
//        importerSolid->closeFile();
//
//        //Putting the vector in the TimeAdvance Stencil
//        timeAdvance->setInitialCondition (solutionStencil);
//    }
//    else //Initialize with zero vectors
//    {
//
//        std::cout << "Starting from scratch" << std::endl;
//        vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
//
//        if (timeAdvanceMethod == "Newmark")
//        {
//            solutionStencil.push_back (disp);
//            solutionStencil.push_back (disp);
//            solutionStencil.push_back (disp);
//        }
//
//        if (timeAdvanceMethod == "BDF")
//        {
//            for ( UInt previousPass = 0; previousPass < dataStructure->dataTimeAdvance()->orderBDF() ; previousPass++)
//            {
//                solutionStencil.push_back (disp);
//            }
//        }
//
//        timeAdvance->setInitialCondition (solutionStencil);
//
//    }
//
//
//
//    //! =================================================================================
//    //! Temporal data and initial conditions
//    //! =================================================================================
//
//    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
//
//    timeAdvance->updateRHSContribution (dataStructure->dataTime()->timeStep() );
//
//    MPI_Barrier (MPI_COMM_WORLD);
//
//    if (verbose )
//    {
//        std::cout << "ok." << std::endl;
//    }
//
//    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterSolid;
//
//    //boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterCheck;
//
//    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
//    std::string const exportFileName = dataFile ( "exporter/nameFile", "structure");
//    std::string const exportCheckName = "checkExporter";
//
//    //If exporterVTK is used
//    // if( dataStructure->order().compare("P2") == 0 )
//    // {
//    //     exporterSolid.reset( new exporterVTK_Type ( dataFile, exportFileName ) );
//    //     exporterCheck.reset( new exporterVTK_Type ( dataFile, exportCheckName ) );
//    // }
//
//#ifdef HAVE_HDF5
//    if (exporterType.compare ("hdf5") == 0)
//    {
//        exporterSolid.reset ( new hdf5Filter_Type ( dataFile, exportFileName ) );
//        //exporterCheck.reset( new hdf5Filter_Type ( dataFile, exportCheckName ) );
//    }
//    else
//#endif
//    {
//        if (exporterType.compare ("none") == 0)
//        {
//            exporterSolid.reset ( new emptyExporter_Type ( dataFile, meshPart.meshPartition(), exportFileName, parameters->comm->MyPID() ) );
//        }
//        else
//        {
//            exporterSolid.reset ( new ensightFilter_Type ( dataFile, meshPart.meshPartition(), exportFileName, parameters->comm->MyPID() ) );
//        }
//    }
//
//    exporterSolid->setPostDir ( "./" );
//    exporterSolid->setMeshProcId ( meshPart.meshPartition(), parameters->comm->MyPID() );
//    //exporterCheck->setMeshProcId( meshPart.meshPartition(), parameters->comm->MyPID() );
//
//    // if( dataStructure->order().compare("P2") == 0 )
//    // {
//    //     exporterSolid->setDataFromGetPot( dataFile, "exporter" );
//    //     exporterCheck->setDataFromGetPot( dataFile, "exporter" );
//    // }
//
//    //discretization FESpace : the possibilities offered by the FESpace class
//    //export FESpace : P1 since Paraview does not read higher fields
//
//    //This are the vector to extract the displacement, velocity and acceleration.
//    //Their map has to be the one of the discretization FESpace
//    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(),  Unique ) );
//    vectorPtr_Type solidVel  ( new vector_Type (solid.displacement(),  Unique ) );
//    vectorPtr_Type solidAcc  ( new vector_Type (solid.displacement(),  Unique ) );
//    //Created for debug reasons. The map has to be the one of the discretization FESpace
//    //vectorPtr_Type rhsVector ( new vector_Type(solid.rhs(),  Unique ) );
//
//    vectorPtr_Type solidDispReduced;
//    vectorPtr_Type solidVelReduced;
//    vectorPtr_Type solidAccReduced;
//    vectorPtr_Type rhsVectorReduced;
//
//    // if( dataStructure->order().compare("P2") == 0 )
//    // {
//    //     solidDispReduced.reset( new vector_Type( exporterFESpace->map(),  Unique ) );
//    //     solidVelReduced.reset( new vector_Type( exporterFESpace->map(),  Unique ) );
//    //     solidAccReduced.reset( new vector_Type( exporterFESpace->map(),  Unique ) );
//    //     //rhsVectorReduced.reset( new vector_Type( exporterFESpace->map(),  Unique ) );
//    //     exporterSolid->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", exporterFESpace, solidDispReduced, UInt(0) );
//    //     exporterSolid->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     exporterFESpace, solidVelReduced,  UInt(0) );
//    //     exporterSolid->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", exporterFESpace, solidAccReduced,  UInt(0) );
//    //     //exporterCheck->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "rhs", exporterFESpace, rhsVectorReduced,  UInt(0) );
//    // }
//    // else
//    // {
//    exporterSolid->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
//    exporterSolid->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt (0) );
//    exporterSolid->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt (0) );
//    //exporterCheck->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "rhs", dFESpace, rhsVector,  UInt(0) );
//    // }
//
//    exporterSolid->postProcess ( 0 );
//    //exporterCheck->postProcess( 0 );
//
//    //!--------------------------------------------------------------------------------------------
//    //!The update of the RHS is done by the TimeAdvance class
//    //solid.updateSystem();
//    //! =================================================================================
//
//    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
//
//    //! 5. Initial data
//    Real initialTime = dataStructure->dataTime()->initialTime();
//    Real dt = dataStructure->dataTime()->timeStep();
//    Real T  = dataStructure->dataTime()->endTime();
//
//    //! =============================================================================
//    //! Temporal loop
//    //! =============================================================================
//    for (Real time = initialTime + dt; time <= T; time += dt)
//    {
//        dataStructure->dataTime()->setTime (time);
//
//        if (verbose)
//        {
//            std::cout << std::endl;
//            std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
//        }
//
//        //! 6. Updating right-hand side
//        *rhs *= 0;
//        timeAdvance->updateRHSContribution ( dt );
//        *rhs += *solid.massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;
//
//        std::cout << "Norm of the rhsNoBC: " << (*rhs).norm2() << std::endl;
//        solid.setRightHandSide ( *rhs );
//
//        //! 7. Iterate --> Calling Newton
//        solid.iterate ( BCh );
//
//        timeAdvance->shiftRight ( solid.displacement() );
//
//        *solidDisp = solid.displacement();
//        *solidVel  = timeAdvance->firstDerivative();
//        *solidAcc  = timeAdvance->secondDerivative();
//        //*rhsVector = solid.rhs();
//
//        // if( dataStructure->order().compare("P2") == 0 )
//        // {
//        //     *solidDispReduced = exporterFESpace->feToFEInterpolate( *dFESpace, *solidDisp);
//        //     *solidVelReduced = exporterFESpace->feToFEInterpolate( *dFESpace, *solidVel );
//        //     *solidAccReduced = exporterFESpace->feToFEInterpolate( *dFESpace, *solidAcc );
//        //     //*rhsVectorReduced = exporterFESpace->feToFEInterpolate( *dFESpace, *rhsVector );
//        //     exporterSolid->postProcess( time );
//        // }
//        // else
//        // {
//        exporterSolid->postProcess ( time );
//        //exporterCheck->postProcess( time );
//        // }
//
//        Real normVect;
//        normVect =  solid.displacement().norm2();
//        std::cout << "The norm 2 of the displacement field is: " << normVect << std::endl;
//
//
//        //!--------------------------------------------------------------------------------------------------
//
//        MPI_Barrier (MPI_COMM_WORLD);
//    }
//
//    exporterSolid->closeFile();
//}
//
//
//int
//main ( int argc, char** argv )
//{
//
//#ifdef HAVE_MPI
//    MPI_Init (&argc, &argv);
//    boost::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
//    if ( Comm->MyPID() == 0 )
//    {
//        cout << "% using MPI" << endl;
//    }
//#else
//    boost::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
//    cout << "% using serial Version" << endl;
//#endif
//
//    Structure structure ( argc, argv, Comm );
//    structure.run();
//
//#ifdef HAVE_MPI
//    MPI_Finalize();
//#endif
//    return returnValue ;
//}
//@HEADER
/*

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
/**
   \file main.cpp

   This test is the case of traction of a cube. It does not use the symmetry BCs
   This test uses the FESpace which is the standard in LifeV and the ETFESpace
   The FESpace is used for the BCs of Neumann type since in the ET branch there
   is not the integration over the boundary faces.

   \author Paolo Tricerri <paolo.tricerri@epfl.ch>
   \date 1861-03-17
 */

#ifdef TWODIM
#error test_structure cannot be compiled in 2D
#endif

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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinearPenalized.hpp>
#include <lifev/structure/solver/SecondOrderExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

//Includes for the Expression Template
#include <lifev/eta/fem/ETFESpace.hpp>

#include <iostream>


using namespace LifeV;

int returnValue = EXIT_FAILURE;
enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}


std::set<UInt> parseList ( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    size_t commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find ( "," );
        setList.insert ( atoi ( stringList.substr ( 0, commaPos ).c_str() ) );
        stringList = stringList.substr ( commaPos + 1 );
    }
    setList.insert ( atoi ( stringList.c_str() ) );
    return setList;
}


class Structure
{
public:

    /** @name Constructors, destructor
     */
    //@{
    Structure ( int                                   argc,
                char**                                argv,
                boost::shared_ptr<Epetra_Comm>        structComm );

    ~Structure()
    {}
    //@}

    //@{
    void run()
    {
        run3d();
    }
    void CheckResultLE (const Real& dispNorm, const Real& time);
    void CheckResultSVK (const Real& dispNorm, const Real& time);
    void CheckResultSVKPenalized (const Real& dispNorm, const Real& time);
    void CheckResultEXP (const Real& dispNorm, const Real& time);
    void CheckResultNH (const Real& dispNorm, const Real& time);
    void CheckResult2ndOrderExponential (const Real& dispNorm, const Real& time);
    void resultChanged (Real time);
    //@}

protected:

private:

    /**
     * run the 2D driven cylinder simulation
     */
    void run2d();

    /**
     * run the 3D driven cylinder simulation
     */
    void run3d();

private:
    struct Private;
    boost::shared_ptr<Private> parameters;
};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

    static Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
    {
        return  0.;
    }

    static Real bcNonZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
    {
        return  300000.;
    }

    static Real bcNonZeroSecondOrderExponential (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
    {
        return  19180.;
    }

    static Real d0 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
    {

        switch (i)
        {
            case 0:
                return  0.088002 * ( x + 0.5 );
                break;
            case 1:
                return - ( 0.02068 * 2.0 ) * ( y );
                break;
            case 2:
                return - ( 0.02068 * 2.0 ) * ( z );
                break;
            default:
                ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
                return 0.;
                break;
        }

    }

    static Real gf (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
    {
        return  ( - 0.1 );
    }

};



Structure::Structure ( int                                   argc,
                       char**                                argv,
                       boost::shared_ptr<Epetra_Comm>        structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->rho     = dataFile ( "solid/physics/density", 1. );
    parameters->young   = dataFile ( "solid/physics/young",   1. );
    parameters->poisson = dataFile ( "solid/physics/poisson", 1. );
    parameters->bulk    = dataFile ( "solid/physics/bulk",    1. );
    parameters->alpha   = dataFile ( "solid/physics/alpha",   1. );
    parameters->gamma   = dataFile ( "solid/physics/gamma",   1. );

    std::cout << "density = " << parameters->rho     << std::endl
              << "young   = " << parameters->young   << std::endl
              << "poisson = " << parameters->poisson << std::endl
              << "bulk    = " << parameters->bulk    << std::endl
              << "alpha   = " << parameters->alpha   << std::endl
              << "gamma   = " << parameters->gamma   << std::endl;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID() )
    {
        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}



void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}



void
Structure::run3d()
{
    typedef StructuralOperator< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                              vectorPtr_Type;
    typedef boost::shared_ptr< TimeAdvance< vector_Type > >             timeAdvance_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;


    bool verbose = (parameters->comm->MyPID() == 0);

    //! Number of boundary conditions for the velocity and mesh motion
    boost::shared_ptr<BCHandler> BCh ( new BCHandler() );

    //! dataElasticStructure
    GetPot dataFile ( parameters->data_file_name.c_str() );

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra> (  parameters->comm  ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< RegionMesh<LinearTetra> > meshPart ( fullMeshPtr, parameters->comm );

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");

    //Mainly used for BCs assembling (Neumann type)
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (meshPart, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    if (verbose)
    {
        std::cout << std::endl;
    }

    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "Newmark");

    timeAdvance_Type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );

    UInt OrderDev = 2;

    //! initialization of parameters of time Advance method:
    if (timeAdvanceMethod == "Newmark")
    {
        timeAdvance->setup ( dataStructure->dataTimeAdvance()->coefficientsNewmark() , OrderDev);
    }

    if (timeAdvanceMethod == "BDF")
    {
        timeAdvance->setup (dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);
    }

    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
    //timeAdvance->showMe();

    //! #################################################################################
    //! BOUNDARY CONDITIONS
    //! #################################################################################
    vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
    compx[0] = 0;
    compy[0] = 1, compz[0] = 2;
    compxy[0] = 0;
    compxy[1] = 1;
    compxz[0] = 0;
    compxz[1] = 2;
    compyz[0] = 1;
    compyz[1] = 2;

    BCFunctionBase zero (Private::bcZero);
    BCFunctionBase nonZero;

    if ( dataStructure->solidType().compare ("secondOrderExponential") )
    {
        nonZero.setFunction (Private::bcNonZero);
    }
    else
    {
        nonZero.setFunction (Private::bcNonZeroSecondOrderExponential);
    }
    //! =================================================================================
    //! BC for StructuredCube4_test_structuralsolver.mesh
    //! =================================================================================
//    BCh->addBC ("lateral",      1,  Natural,   Full, zero, 3);
//    BCh->addBC ("top",      3,  Natural,   Full, zero, 3);
//    BCh->addBC ("bottom",     2,  Essential, Full, zero, 3);


//    BCh->addBC ("EdgesIn",      20,  Natural,   Full, zero, 3);


//    BCh->addBC ("EdgesIn",      20,  Natural,   Component, zero, compx);
//    BCh->addBC ("EdgesIn",      40,  Essential, Component, zero,    compx);
//
//    BCh->addBC ("SymmetryX",    7,  Essential, Component, zero,    compy);
//    BCh->addBC ("SymmetryY",    3,  Essential, Component, zero,    compz);
//    BCh->addBC ("edgeone",      1000,  EssentialVertices, Full, zero, 3);
//    BCh->addBC ("edgetwo",      500,  EssentialVertices, Component, zero,    compxz);
//    BCh->addBC ("edgetwo",      300,  EssentialVertices, Component, zero,    compxy);
//    BCh->addBC ("edgetwo",      800,  EssentialVertices, Component, zero,    compyz);


    	//My boundary conditions for sym CUBE
    BCh->addBC ("EdgesIn",      100,  Essential, Component, zero,    compx);
    BCh->addBC ("SymmetryX",    200,  Essential, Component, zero,    compz);
    BCh->addBC ("SymmetryY",    300,  Essential, Component, zero,    compy);
    BCh->addBC ("edgetwo",      10,  EssentialVertices, Component, zero,    compxz);
    BCh->addBC ("edgetwo",      20,  EssentialVertices, Component, zero,    compxy);
    BCh->addBC ("edgetwo",      30,  EssentialVertices, Component, zero,    compyz);
//    BCh->addBC ("edgetwo",      10,  EssentialEdges, Component, zero,    compxz);
//    BCh->addBC ("edgetwo",      20,  EssentialEdges, Component, zero,    compxy);
//    BCh->addBC ("edgetwo",      30,  EssentialEdges, Component, zero,    compyz);
    //! =================================================================================


    //! 1. Constructor of the structuralSolver
    StructuralOperator< RegionMesh<LinearTetra> > solid;

//    //! 2. Setup of the structuralSolver
//    solid.setup (dataStructure,
//                 dFESpace,
//                 dETFESpace,
//                 BCh,
//                 parameters->comm);

    //! 2. Setup of the structuralSolver
    solid.setup (dataStructure,
                 dFESpace,
                 dETFESpace,
                 BCh,
                 parameters->comm);


    solidFESpacePtr_Type gammaFESpace( new solidFESpace_Type (meshPart, dOrder, 1, parameters->comm) );
    VectorEpetra gammaf( gammaFESpace -> map() );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    fct_type fg = &(Private::gf);
    gammaFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( fg ), gammaf , 0);
    solid.material() -> setGammaf(gammaf);
    //! 3. Setting data from getPot
    solid.setDataFromGetPot (dataFile);


    Real fx=1.0;
    Real fy = 0.0;
    Real fz = 0.0;
    solid.material() -> setupFiberVector(fx, fy, fz);



    //! 4. Building system using TimeAdvance class
    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
    solid.buildSystem (timeAdvanceCoefficient);


    //dataStructure->showMe();
    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    //! 5. Initial data
    Real dt = dataStructure->dataTime()->timeStep();
    // Real T  = dataStructure->dataTime()->endTime();

    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type vel (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type acc (new vector_Type (solid.displacement(), Unique) );

    if (verbose)
    {
        std::cout << "S- initialization ... ";
    }

    std::vector<vectorPtr_Type> uv0;

    if (timeAdvanceMethod == "Newmark")
    {
        uv0.push_back (disp);
        uv0.push_back (vel);
        uv0.push_back (acc);
    }

    vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );

    if ( !dataStructure->solidType().compare ("secondOrderExponential") )
    {
        dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type> ( Private::d0 ), *initialDisplacement, 0.0 );
    }

    if (timeAdvanceMethod == "BDF")
    {
        Real tZero = dataStructure->dataTime()->initialTime();

        for ( UInt previousPass = 0; previousPass < timeAdvance->size() ; previousPass++)
        {
            Real previousTimeStep = tZero - previousPass * dt;
            std::cout << "BDF " << previousTimeStep << "\n";
            if ( !dataStructure->solidType().compare ("secondOrderExponential") )
            {
                uv0.push_back (initialDisplacement);
            }
            else
            {
                uv0.push_back (disp);
            }
        }
    }

    timeAdvance->setInitialCondition (uv0);

    timeAdvance->setTimeStep ( dt );

    timeAdvance->updateRHSContribution ( dt );

    if ( !dataStructure->solidType().compare ("secondOrderExponential") )
    {
        solid.initialize ( initialDisplacement );
    }
    else
    {
        solid.initialize ( disp );
    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose )
    {
        std::cout << "ok." << std::endl;
    }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }

        else
        {
            exporter.reset ( new ExporterEnsight<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( meshPart.meshPartition(), parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), Repeated /*,  exporter->mapType()*/ ) );
    vectorPtr_Type solidVel  ( new vector_Type (solid.displacement(), Repeated/*,  exporter->mapType()*/ ) );
    vectorPtr_Type solidAcc  ( new vector_Type (solid.displacement(), Repeated/*,  exporter->mapType()*/ ) );

//    vectorPtr_Type solidgamma  ( new vector_Type ( *( solid.material() -> gammaf() ), Unique ) );
//    vectorPtr_Type vecgamma  ( new vector_Type ( gammaf , Unique ) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt (0) );
//    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "gamma_function", gammaFESpace, vecgamma,  UInt (0) );
//    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "gamma_material", gammaFESpace, solidgamma,  UInt (0) );


    exporter->postProcess ( 0 );


    /*
     //!--------------------------------------------------------------------------------------------
     //! MATLAB FILE WITH DISPLACEMENT OF A CHOSEN POINT
     //!--------------------------------------------------------------------------------------------
     cout.precision(16);
     ofstream file_comp( "Displacement_components_NL.m" );
     if ( !file_comp )
     {
       std::cout <<" Unable to open file! You need to specify the output folder in the data file " << std::endl;
     }

     int IDPoint = 73; // StructuredCube4
     //int IDPoint = 401; // StructuredCube8
     //int IDPoint = 2593; // StructuredCube16

     //int IDPoint = 74;// cube4x4.mesh
     //int IDPoint = 315;// cube8x8.mesh
     //int IDPoint = 1526;// cube16x16.mesh

     file_comp << " % TEST NONLINEAR ELASTICITY" << endl;
     file_comp << " % Displacement components of ID node  " << IDPoint << " :" << endl;
     file_comp << " % Each row is a time step" << endl;
     file_comp << " % First column = comp x, second = comp y and third = comp z. " << endl;
     file_comp <<  endl;
     file_comp << " SolidDisp_NL = [ " ;

     for ( UInt k = IDPoint - 1; k <= solid.displacement().size() - 1; k = k + solid.displacement().size()/nDimensions )
     {
     file_comp<< solid.displacement()[ k ] << " ";
     }

     file_comp<< endl;
     */
     //!--------------------------------------------------------------------------------------------
     //!The update of the RHS is done by the TimeAdvance class
     //solid.updateSystem();
     //! =================================================================================


     Real normVect;
     normVect =  solid.displacement().norm2();
     std::cout << "The norm 2 of the displacement field is: " << normVect << std::endl;

     //! =============================================================================
     //! Temporal loop
     //! =============================================================================
     //    for (Real time = dt; time <= T; time += dt)
//     for (dataStructure->dataTime()->setTime ( dt ) ; dataStructure->dataTime()->canAdvance( ); dataStructure->dataTime()->updateTime( ) )
//     {
//         //      dataStructure->dataTime()->setTime(time);
//
//         if (verbose)
//         {
//             std::cout << std::endl;
//             std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
//         }
//
//         if (verbose)
//         {
//             std::cout << std::endl;
//             std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
//         }
//
//         //! 6. Updating right-hand side
//         *rhs *= 0;
//         timeAdvance->updateRHSContribution ( dt );
//         *rhs += *solid.massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;
//         solid.setRightHandSide ( *rhs );

         //! 7. Iterate --> Calling Newton
         solid.iterate ( BCh );

         timeAdvance->shiftRight ( solid.displacement() );

         *solidDisp = solid.displacement();
         *solidVel  = timeAdvance->firstDerivative();
         *solidAcc  = timeAdvance->secondDerivative();

         exporter->postProcess ( dataStructure->dataTime()->time() );

         /* This part lets to save the displacement at one point of the mesh and to check the result
            w.r.t. manufactured solution.
            //!--------------------------------------------------------------------------------------------------
            //! MATLAB FILE WITH DISPLACEMENT OF A CHOOSEN POINT
            //!--------------------------------------------------------------------------------------------------
            cout <<"*******  DISPLACEMENT COMPONENTS of ID node "<< IDPoint << " *******"<< std::endl;
            for ( UInt k = IDPoint - 1 ; k <= solid.displacement().size() - 1; k = k + solid.displacement().size()/nDimensions )
            {
            file_comp<< solid.displacement()[ k ] << " ";
            cout.precision(16);
            cout <<"*********************************************************"<< std::endl;
            cout <<" solid.disp()[ "<< k <<" ] = "<<  solid.displacement()[ k ]  << std::endl;
            cout <<"*********************************************************"<< std::endl;
            }
            file_comp<< endl;
         */

//         Real normVect;
         normVect =  solid.displacement().norm2();
         std::cout << "The norm 2 of the displacement field is: " << normVect << std::endl;

         ///////// CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
//         if (!dataStructure->solidType().compare ("linearVenantKirchhoff") )
//         {
//             CheckResultLE (normVect, dataStructure->dataTime()->time() );
//         }
//         else if (!dataStructure->solidType().compare ("nonLinearVenantKirchhoff") )
//         {
//             CheckResultSVK (normVect, dataStructure->dataTime()->time() );
//         }
//         else if (!dataStructure->solidType().compare ("nonLinearVenantKirchhoffPenalized") )
//         {
//             CheckResultSVKPenalized (normVect, dataStructure->dataTime()->time() );
//         }
//         else if (!dataStructure->solidType().compare ("exponential") )
//         {
//             CheckResultEXP (normVect, dataStructure->dataTime()->time() );
//         }
//         else if (!dataStructure->solidType().compare ("secondOrderExponential") )
//         {
//             CheckResult2ndOrderExponential (normVect, dataStructure->dataTime()->time() );
//         }
//         else
//         {
//             CheckResultNH (normVect, dataStructure->dataTime()->time() );
//         }

         ///////// END OF CHECK


         //!--------------------------------------------------------------------------------------------------
        MPI_Barrier (MPI_COMM_WORLD);
    //}
}



void Structure::CheckResultLE (const Real& dispNorm, const Real& time)
{
    if ( time == 0.1  && std::fabs (dispNorm - 0.276527) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.2  && std::fabs (dispNorm - 0.276536) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.3  && std::fabs (dispNorm - 0.276529) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.4  && std::fabs (dispNorm - 0.276531) <= 1e-5 )
    {
        this->resultChanged (time);
    }
}

void Structure::CheckResultSVK (const Real& dispNorm, const Real& time)
{
    if ( time == 0.1  && std::fabs (dispNorm - 0.263348) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.2  && std::fabs (dispNorm - 0.263350) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.3  && std::fabs (dispNorm - 0.263350) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.4  && std::fabs (dispNorm - 0.263351) <= 1e-5 )
    {
        this->resultChanged (time);
    }
}
void Structure::CheckResultSVKPenalized (const Real& dispNorm, const Real& time)
{
    if ( time == 0.1  && std::fabs (dispNorm - 0.254316) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.2  && std::fabs (dispNorm - 0.254322) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.3  && std::fabs (dispNorm - 0.254317) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.4  && std::fabs (dispNorm - 0.254318) <= 1e-5 )
    {
        this->resultChanged (time);
    }
}
void Structure::CheckResultEXP (const Real& dispNorm, const Real& time)
{
    if ( time == 0.1  && std::fabs (dispNorm - 0.284844) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.2  && std::fabs (dispNorm - 0.284853) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.3  && std::fabs (dispNorm - 0.284846) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.4  && std::fabs (dispNorm - 0.284848) <= 1e-5 )
    {
        this->resultChanged (time);
    }
}

void Structure::CheckResultNH (const Real& dispNorm, const Real& time)
{
    if ( time == 0.1  && std::fabs (dispNorm - 0.286120) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.2  && std::fabs (dispNorm - 0.286129) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.3  && std::fabs (dispNorm - 0.286122) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.4  && std::fabs (dispNorm - 0.286123) <= 1e-5 )
    {
        this->resultChanged (time);
    }
}


void Structure::CheckResult2ndOrderExponential (const Real& dispNorm, const Real& time)
{
    if ( time == 0.1  && std::fabs (dispNorm - 0.561523) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.2  && std::fabs (dispNorm - 0.561496) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.3  && std::fabs (dispNorm - 0.561517) <= 1e-5 )
    {
        this->resultChanged (time);
    }
    if ( time == 0.4  && std::fabs (dispNorm - 0.561512) <= 1e-5 )
    {
        this->resultChanged (time);
    }
}



void Structure::resultChanged (Real time)
{
    std::cout << "Correct value at time: " << time << std::endl;
    returnValue = EXIT_SUCCESS;
}



int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }
#else
    boost::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif

    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
