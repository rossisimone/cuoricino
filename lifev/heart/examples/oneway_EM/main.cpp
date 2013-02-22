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
    @brief main for the one way electromechanical coupling

    @date 09−2012
    @author Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>

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

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/heart/solver/HeartUpdatedMonodomainSolver.hpp>
#include <lifev/heart/solver/NeoHookeanActivatedMaterial.hpp>
#include <lifev/heart/solver/HeartIonicSolver.hpp>
#include <lifev/heart/solver/HeartAlievPanfilov.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/ExponentialMaterialNonLinear.hpp>


#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

using namespace LifeV;

int returnValue = EXIT_SUCCESS;

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


class ElectroMech

{
public:

    /** @name Constructors, destructor
     */
    //@{
    ElectroMech ( int                                   argc,
                  char**                                argv,
                  boost::shared_ptr<Epetra_Comm>        structComm );

    ~ElectroMech()
    {}
    //@}

    //@{
    void run();
    //! To compute the displacements
    void computeDispl (HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type& dep,
                       HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >& electricModel,
                       HeartMonodomainData& dataMonodomain);

    //! To compute the righthand side of the electrical system
    void computeRhs ( HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type& rhs,
                      HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >& electricModel,
                      boost::shared_ptr< HeartIonicSolver< RegionMesh<LinearTetra> > > ionicModel,
                      HeartMonodomainData& dataMonodomain );

    //@}


private:
    UInt ion_model;
    UInt nbeq;
    boost::shared_ptr<HeartFunctors> M_heart_fct;
    struct Private;
    boost::shared_ptr<Private> parameters;
};


struct ElectroMech::Private
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

    static Real zero_scalar ( const Real& /* t */,
                              const Real& /* x */,
                              const Real& /* y */,
                              const Real& /* z */,
                              const ID& /* i */ )
    {
        return 0.;
    }


    static Real leftfront_init ( const Real& /* t */,
                                 const Real& x,
                                 const Real& /* y*/,
                                 const Real& /* z*/,
                                 const ID&   /* i */ )
    {
        Real const oldi (1.0 - 1.0 / (1.0 + exp (-90.0 * (x - 0.15) ) ) );

        return -84.0 + 100.0 * oldi;
    }

    static Real minus84_scalar ( const Real& /* t */,
                                 const Real& /* x */,
                                 const Real& /* y */,
                                 const Real& /* z */,
                                 const ID& /* i */ )
    {
        return  -84.0;
    }

    Real scrollwaveInitEllisoid ( const Real& /* t */,
                                  const Real& x,
                                  const Real& y,
                                  const Real& z,
                                  const ID& /* i */ )
    {
        double r = std::sqrt ( (x - 0.0) * (x - 0.0) / 3.0625
                               + (y - 0.0) * (y - 0.0) / 3.0625
                               + (z + 0.15) * (z + 0.15) / 21.16); // elli
        Real const oldi (1.0 - 1.0 / (1.0 + exp (-90.0 * (r - 1) ) ) ); // elli
        return -84.0 + 110.0 * oldi;
    }
    Real scrollwaveInitHeart ( const Real& /* t */,
                               const Real& x,
                               const Real& y,
                               const Real& z,
                               const ID& /* i */ )
    {
        double r1 = std::sqrt ( (x - 45.0) * (x - 45.0) / 1600 //40^2
                                + (y - 115.0) * (y - 115.0) / 625 //25^2
                                + (z - 70.0) * (z - 70.0) / 2025); //45^2
        double r2 = std::sqrt ( (x - 40.0) * (x - 40.0) / 3600 //60^2
                                + (y - 60.0) * (y - 60.0) / 3025 //55^2
                                + (z - 70.0) * (z - 70.0) / 1600); //40^2


        Real const oldi (2.0 - 1.0 / (1.0 + exp (-90.0 * (r1 - 1) ) )
                         - 1.0 / (1.0 + exp (-90.0 * (r2 - 1) ) ) ); // heart
        return -84.0 + 110.0 * oldi;
    }


    Real scrollwaveInitCanineHeart ( const Real& /* t */,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const ID& /* i */ )
    {
        double r1 = std::sqrt ( (x - 24.0) * (x - 24.0) / 61 //9^2  /left
                                + (y - 40.0) * (y - 40.0) / 525 //15^2
                                + (z - 44.0) * (z - 44.0) / 90); //14^2
        double r2 = std::sqrt ( (x - 50.0) * (x - 50.0) / 240 //13^2  /right
                                + (y - 39.0) * (y - 39.0) / 240 //13^2
                                + (z - 36.0) * (z - 36.0) / 250); //100^2


        Real const oldi (2.0 - 1.0 / (1.0 + exp (-90.0 * (r1 - 1) ) )
                         - 1.0 / (1.0 + exp (-90.0 * (r2 - 1) ) ) ); // heart
        return -84.0 + 110.0 * oldi;
    }

    Real d0elli (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
    {
        double rr = std::sqrt ( (x + 3.0) * (x + 3.0) + (y - 0.0) * (y - 0.0) + (z + 80.0) * (z + 80.0) ); //elli
        if (t == 0)
        {
            switch (i)
            {
                case 1:
                    return 10.0 * (1.0 - 1.0 / (1.0 + exp (-90.0 * (rr - 35.0) ) ) ); //-z*(z - 5.)*x/50.; // -z*(z - 80.)*x/5000.;//
                    return 0;
                    break;
                case 2:
                    return 10.0 * (1.0 - 1.0 / (1.0 + exp (-90.0 * (rr - 35.0) ) ) ); //-z*(z - 5.)*y/50.; // -z*(z - 80.)*y/5000.;//
                    return 0;
                    break;
                case 3:
                    return 0.0;
                    break;
                default:
                    ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
                    break;
            }
        }

        return 0.0;
    }

    Real d0heart (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
    {
        double rr = std::sqrt ( (x - 110.0) * (x - 110.0) + (y - 80.0) * (y - 80.0) + (z - 70.0) * (z - 70.0) ); //heart
        if (t == 0)
        {
            switch (i)
            {
                case 1:
                    return 10.0 * (1.0 - 1.0 / (1.0 + exp (-90.0 * (rr - 45.0) ) ) ); //-z*(z - 5.)*x/50.; // -z*(z - 80.)*x/5000.;//
                    return 0;
                    break;
                case 2:
                    return 10.0 * (1.0 - 1.0 / (1.0 + exp (-90.0 * (rr - 45.0) ) ) ); //-z*(z - 5.)*y/50.; // -z*(z - 80.)*y/5000.;//
                    return 0;
                    break;
                case 3:
                    return 0.0;
                    break;
                default:
                    ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
                    break;
            }
        }

        return 0.0;
    }

static Real d0CanineHeart(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    double rr = std::sqrt((x-1.0)*(x-1.0)+(y-1.5)*(y-1.5)+(z-0.9)*(z-0.9));//canineheart
    if (t == 0)
    {
          return 10.0*(1.0-1.0/(1.0+exp(-90.0*(rr-0.5))));//-z*(z - 5.)*x/50.; // -z*(z - 80.)*x/5000.;//
    }

        return 0.0;
    }



    static Real w0 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
    {
        return 0.0;
    }


};


//! Identifiers for heart boundaries
const int EPICARDIUM    = 40;
const int ENDOCARDIUM   = 60;
const int TRUNC_SEC     = 50;


// ===================================================
// Constructors & Destructor
// ===================================================

ElectroMech::ElectroMech ( int argc, char** argv,
                           boost::shared_ptr<Epetra_Comm> structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //! Pointer to access functors
    M_heart_fct.reset (new HeartFunctors ( dataFile) );
    ion_model = dataFile ("electric/physics/ion_model", 1);
    std::cout << "mpi initialization ... " << std::endl;
    M_heart_fct->M_comm.reset (new Epetra_MpiComm ( MPI_COMM_WORLD ) );

    if (!M_heart_fct->M_comm->MyPID() )
    {
        std::cout << "My PID = " << M_heart_fct->M_comm->MyPID() << std::endl;
    }
    parameters.reset ( new Private() );
    parameters->data_file_name = data_file_name;
    // density = 0 to consider the pseudo static mechanics
    parameters->rho     = dataFile ( "solid/physics/density", 0. );
    parameters->young   = dataFile ( "solid/physics/young",   1. );
    parameters->poisson = dataFile ( "solid/physics/poisson", 1. );
    parameters->bulk    = dataFile ( "solid/physics/bulk",    1. );
    parameters->alpha   = dataFile ( "solid/physics/alpha",   1. );
    parameters->gamma   = dataFile ( "solid/physics/gamma",   1. );
    std::cout << "density = " << parameters->rho << std::endl
              << "young   = " << parameters->young << std::endl
              << "poisson = " << parameters->poisson << std::endl;
    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();
    if (!parameters->comm->MyPID() )
    {
        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}


// ===================================================
// Methods
// ===================================================
void
ElectroMech::run()
{
  typedef HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
  typedef HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::matrix_Type  matrix_Type;
  typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
  typedef boost::shared_ptr< TimeAdvance< vector_Type > >       timeAdvance_type;

  LifeChrono chronoinitialsettings;
  LifeChrono chronototaliterations;
  chronoinitialsettings.start();
  Real normu;
  Real meanu;
  Real minu;
  //Real maxu;
  Real maxdispl;
  Real mindispl;
  //! Construction of data classes

  HeartMonodomainData _dataFunctors(M_heart_fct);
  HeartIonicData _dataIonic(M_heart_fct->M_dataFile);
  bool verbose = (M_heart_fct->M_comm->MyPID() == 0);

  //! Number of boundary conditions for the velocity and mesh motion
  //boost::shared_ptr<BCHandler> BChE( new BCHandler() );
  BCHandler BChE;
  BCFunctionBase uZero(Private::zero_scalar);
  BChE.addBC( "Endo",   	ENDOCARDIUM,	Natural,	Full,	uZero,  1 );
  BChE.addBC( "Epi",   	        EPICARDIUM, 	Natural,   	Full,   uZero, 	1 );
  BChE.addBC( "Trunc",    	TRUNC_SEC,  	Natural, 	Full,   uZero, 	1 );

  boost::shared_ptr<BCHandler> BChS( new BCHandler() );
  BCFunctionBase dZero(Private::bcZero);

  //! dataElasticStructure
  GetPot dataFile( parameters->data_file_name.c_str() );

  boost::shared_ptr<StructuralConstitutiveLawData> dataStructure(new StructuralConstitutiveLawData( ));
  dataStructure->setup(dataFile);

  MeshData             meshData;
  meshData.setup(dataFile, "solid/space_discretization");

  boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr(new RegionMesh<LinearTetra>( parameters->comm ) );
  readMesh(*fullMeshPtr, meshData);

  MeshPartitioner< RegionMesh<LinearTetra> > meshPart( fullMeshPtr, parameters->comm );

  BCFunctionBase fixed(dZero);

  BChS->addBC( "Endo",  	ENDOCARDIUM,	Natural,	Full,	fixed,  3 );
  BChS->addBC( "Epi",   	EPICARDIUM, 	Natural,   	Full,   dZero, 	3 );
  BChS->addBC( "Trunc",  TRUNC_SEC,  	Essential, 	Full,   dZero, 	3 );

  std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");
  std::string uOrder =  M_heart_fct->M_dataFile( "electric/space_discretization/u_order", "P1");
  std::string wOrder =  M_heart_fct->M_dataFile( "electric/space_discretization/w_order", "P1");

  typedef FESpace< RegionMesh<LinearTetra>, MapEpetra > solidFESpace_type;
  typedef boost::shared_ptr<solidFESpace_type> solidFESpace_ptrtype;
  solidFESpace_ptrtype dFESpace( new solidFESpace_type(meshPart,dOrder,3,parameters->comm) );
  if (verbose) std::cout << std::endl;

  std::string timeAdvanceMethod =  dataFile( "solid/time_discretization/method", "BDF");
  timeAdvance_type  timeAdvance( TimeAdvanceFactory::instance().createObject( timeAdvanceMethod ) );
  UInt OrderDev = 2;

    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "BDF");
    timeAdvance_type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );
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
    timeAdvance->showMe();

    const ReferenceFE*    refFE_w;
    const QuadratureRule* qR_w;
    const QuadratureRule* bdQr_w;

    const ReferenceFE*    refFE_u;
    const QuadratureRule* qR_u;
    const QuadratureRule* bdQr_u;

    refFE_u = &feTetraP1;
    qR_u    = &quadRuleTetra15pt;
    bdQr_u  = &quadRuleTria3pt;

    refFE_w = &feTetraP1;
    qR_w    = &quadRuleTetra4pt;
    bdQr_w  = &quadRuleTria3pt;

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra > uFESpace_type;
    typedef boost::shared_ptr<uFESpace_type> uFESpace_ptrtype;

    uFESpace_ptrtype uFESpacePtr ( new uFESpace_type (meshPart, *refFE_u, *qR_u,  *bdQr_u,
                                                      1,
                                                      M_heart_fct->M_comm) );


    UInt totalUDof  = uFESpacePtr->map().map (Unique)->NumGlobalElements();

    if (verbose)
    {
        std::cout << "Total Potential DOF = " << totalUDof << std::endl;
    }


    if (verbose)
    {
        std::cout << "Calling the electric model constructor ... " << std::endl << std::flush;
    }

    HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> > electricModel (_dataFunctors, *uFESpacePtr, BChE, M_heart_fct->M_comm);
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Calling the ionic model constructor ... ";
    }

    boost::shared_ptr< HeartIonicSolver< RegionMesh<LinearTetra> > > ionicModel;
    if (ion_model == 1)
    {
        if (verbose)
        {
            std::cout << "Ion Model = Rogers-McCulloch" << std::endl << std::flush;
        }
        ionicModel.reset (new RogersMcCulloch< RegionMesh<LinearTetra> > (_dataIonic,
                                                                          *meshPart.meshPartition(),
                                                                          *uFESpacePtr,
                                                                          *M_heart_fct->M_comm) );
    }
    /*else if (ion_model==2)
      {
        if (verbose) std::cout<<"Ion Model = Luo-Rudy"<<std::endl<<std::flush;
        ionicModel.reset(new HeartLuoRudy< RegionMesh<LinearTetra> >(_dataIonic,
                                 *meshPart.meshPartition(),
                                 *uFESpacePtr,
                                 *M_heart_fct->M_comm));
                                 }*/
    else if (ion_model == 3)
    {
        if (verbose)
        {
            std::cout << "Ion Model = Mitchell-Schaeffer" << std::endl << std::flush;
        }
        ionicModel.reset (new MitchellSchaeffer< RegionMesh<LinearTetra> > (_dataIonic,
                                                                            *meshPart.meshPartition(),
                                                                            *uFESpacePtr,
                                                                            *M_heart_fct->M_comm) );
    }

    else if (ion_model == 4)
    {
        if (verbose)
        {
            std::cout << "Ion Model = Bueno-Orovio minimal human model" << std::endl << std::flush;
        }
        ionicModel.reset (new MinimalModel< RegionMesh<LinearTetra> > (_dataIonic,
                                                                       *meshPart.meshPartition(),
                                                                       *uFESpacePtr,
                                                                       *M_heart_fct->M_comm) );
    }
    else if (ion_model == 5)
    {
        if (verbose)
        {
            std::cout << "Ion Model = Courtemanche-Ramirez-Nattel model for the atrium" << std::endl << std::flush;
        }
        ionicModel.reset (new CourtemancheRamirezNattel< RegionMesh<LinearTetra> > (_dataIonic,
                                                                                    *meshPart.meshPartition(),
                                                                                    *uFESpacePtr,
                                                                                    *M_heart_fct->M_comm) );
    }
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    //! Constructor of the structuralSolver
    StructuralOperator< RegionMesh<LinearTetra> > solid;


    electricModel.setup ( M_heart_fct->M_dataFile );
    std::cout << "electrical setup ok" << std::endl;

    //! Building time-independent part of the system
    electricModel.buildSystem( );
    std::cout << "buildsystem electrical ok" << std::endl;


    //! Setup of the structuralSolver
    solid.setup (dataStructure,
                 dFESpace,
                 BChS,
                 parameters->comm);
    std::cout << "mechanical setup ok" << std::endl;
    solid.setDataFromGetPot (dataFile);

    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
    solid.buildSystem (timeAdvanceCoefficient);

    std::cout << "buildsystem mechanical ok" << std::endl;


    vectorPtr_Type rhsSolid (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type vel (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type acc (new vector_Type (solid.displacement(), Unique) );

    vector_Type displElech (electricModel.disp(), Repeated);
    vector_Type displMech (electricModel.disp(), Repeated);
    vector_Type displ (electricModel.disp(), Repeated);


    dFESpace->interpolate (Private::d0CanineHeart, *disp, 0.0);
    dFESpace->interpolate (Private::w0, *vel , 0.0);


    if (verbose)
    {
        std::cout << "S- initialization ... ";
    }

    electricModel.initialize (Private::minus84_scalar);
    ionicModel->initialize();
    electricModel.resetPreconditioner();
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    //  MapEpetra fullMap(electricModel.getMap());
    //  vector_Type rhs(fullMap);
    vector_Type rhs (electricModel.getMap() );


    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    Real dt = dataStructure->dataTime()->timeStep();
    Real T  = dataStructure->dataTime()->endTime();

    std::vector<vectorPtr_Type> uv0;


    if (timeAdvanceMethod == "Newmark")
    {
        uv0.push_back (disp);
        uv0.push_back (vel);
        uv0.push_back (acc);
    }

    if (timeAdvanceMethod == "BDF")
    {
        for ( UInt previousPass = 0; previousPass < dataStructure->dataTimeAdvance()->orderBDF() ; previousPass++)
        {
            Real previousTimeStep = -previousPass * dt;
            std::cout << "BDF " << previousTimeStep << "\n";
            uv0.push_back (disp);
        }
    }

    timeAdvance->setInitialCondition (uv0);
    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
    timeAdvance->updateRHSContribution (dataStructure->dataTime()->timeStep() );

    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose )
    {
        std::cout << "ok." << std::endl;
    }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    std::string const exporterType =  M_heart_fct->M_dataFile ( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( M_heart_fct->M_dataFile, "electro_mech" ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( M_heart_fct->M_dataFile, meshPart.meshPartition(), "ElectroMech", M_heart_fct->M_comm->MyPID() ) );
        }
        else
        {
            exporter.reset ( new ExporterEnsight<RegionMesh<LinearTetra> > ( M_heart_fct->M_dataFile, meshPart.meshPartition(), "ElectroMech", M_heart_fct->M_comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( meshPart.meshPartition(), M_heart_fct->M_comm->MyPID() );
    vectorPtr_Type Uptr ( new vector_Type (electricModel.solutionTransmembranePotential(), exporter->mapType() ) );
    vectorPtr_Type Dptr ( new vector_Type (electricModel.disp(), exporter->mapType()  ) );
    vectorPtr_Type Fptr ( new vector_Type (electricModel.fiberVector(), exporter->mapType() ) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField,  "potential", uFESpacePtr, Uptr,
                            UInt (0) );
    exporter->postProcess ( 0 );



    //! =============================================================================
    //! Temporal loop
    //! =============================================================================

    LifeChrono chrono;
    chronototaliterations.start();
    for (Real time = dt; time <= T; time += dt)
    {
        dataStructure->dataTime()->setTime (time);
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
        }

        chrono.start();
        MPI_Barrier (MPI_COMM_WORLD);
        ionicModel->solveIonicModel ( electricModel.solutionTransmembranePotential(), _dataFunctors.timeStep() );
        rhs *= 0.;
        computeRhs ( rhs, electricModel, ionicModel, _dataFunctors );

        //! time-dependent updatePDESystem
        electricModel.updatePDESystem ( rhs );
        electricModel.PDEiterate ( BChE );

        *rhsSolid *= 0;
        timeAdvance->updateRHSContribution ( dt );
        *rhsSolid += *solid.massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;
        solid.setRightHandSide ( *rhsSolid );

        //! 7. Iterate --> Calling Newton
        solid.iterate ( BChS );

        timeAdvance->shiftRight ( solid.displacement() );

        *disp = solid.displacement();
        *vel  = timeAdvance->velocity();
        *acc  = timeAdvance->acceleration();


    displMech=solid.displacement();
    computeDispl(displElech, electricModel, _dataFunctors );
    //if (time < 22 )
    //  {
    //	displ=0.0*displElech;
    //  }
    //else
    //  {        //! Solving the solid mechanics
    displ=2*(0.65*displMech+0.25*displElech);
	//}

    electricModel.moveMesh( displ );
    normu=electricModel.solutionTransmembranePotential().norm2();
    electricModel.solutionTransmembranePotential().epetraVector().MeanValue(&meanu);
    electricModel.solutionTransmembranePotential().epetraVector().MaxValue(&minu);
    displ.epetraVector().MaxValue(&maxdispl);
    displ.epetraVector().MinValue(&mindispl);
    if (verbose)
      {
	  std::cout << "norm u " << normu << std::endl;
	  std::cout << "mean u " << meanu << std::endl;
	  std::cout << "max u " << minu << std::endl<<std::flush;
	  std::cout << "max disp " << maxdispl << std::endl;
	  std::cout << "min disp " << mindispl << std::endl;
      }

    *Uptr = electricModel.solutionTransmembranePotential();
    exporter->postProcess( time );
    MPI_Barrier(MPI_COMM_WORLD);
    chrono.stop();
    if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    chronototaliterations.stop();
  }// end of time loop.

        *Uptr = electricModel.solutionTransmembranePotential();
        exporter->postProcess ( time );
        MPI_Barrier (MPI_COMM_WORLD);
        chrono.stop();
        if (verbose)
        {
            std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
        }
        chronototaliterations.stop();
    }// end of time loop.

    if (verbose)
    {
        std::cout << "Total iterations time " << chronototaliterations.diff() << " s." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total initial settings time " << chronoinitialsettings.diff() << " s." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total execution time " << chronoinitialsettings.diff() + chronototaliterations.diff() << " s." << std::endl;
    }
}


void ElectroMech::computeDispl (HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type& displElech, HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >& electricModel, HeartMonodomainData& data)
{
    bool verbose = (M_heart_fct->M_comm->MyPID() == 0);
    typedef HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    if (verbose)
    {
        std::cout << "  f-  Computing displacements ... " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    vector_Type uVecRep (electricModel.solutionTransmembranePotential(), Repeated);
    UInt nbNode = electricModel.potentialFESpace().mesh()->pointList.size();
    UInt dim = electricModel.potentialFESpace().dof().numTotalDof();
    Real aux;
    Vector xx (3), dep (3), cen (3);
    if (verbose)
    {
        std::cout << "  tt-  Vector OK... " << "\n" << std::flush;
    }
    //cen(0)=0; cen(1)=0.75; cen(2)=1.75; //elli
    //cen(0)=70; cen(1)=65; cen(2)=70; //heart
    //cen(0)=40; cen(1)=40; cen(2)=40; //canine
    //cen(0)=2; cen(1)=2.9; cen(2)=3; //canine2
    //cen(0)=-30; cen(1)=0; cen(2)=-30; //biventricular
    cen (0) = 2;
    cen (1) = 2.9;
    cen (2) = 3.0; //heartJean
    for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
    {
        Int id = electricModel.potentialFESpace().mesh()->pointList[ iNode ].id();
        aux = data.beta() * (uVecRep[id] + 84.0) / (184.0 + uVecRep[id]) + 0.001 * uVecRep[id]; //slab
        //aux = data.beta()*(uVecRep[id]+150.0)/(250.0+uVecRep[id]);//heart and elli
        for ( ID j = 1; j <= 3; ++j )
        {
            if ( displElech.blockMap().LID (id + dim * (j - 1) ) >= 0 )
            {
                xx (j - 1) = electricModel.potentialFESpace().mesh()->pointList[ iNode ].coordinate ( j - 1);
                //dep(j-1) = -3*data.beta()*aux*xx(j-1) + 0.25*aux; //slab
                //dep(j-1) = -data.beta()*aux*(xx(j-1)-cen(j-1)) + 5*aux;//heart
                dep (j - 1) = -data.beta() * aux * (xx (j - 1) - cen (j - 1) ) + 0.04 * aux; //elli
                displElech[id + (j - 1) *dim] = dep (j - 1);
            }
        }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}



void ElectroMech::computeRhs ( HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type& rhs,
                               HeartUpdatedMonodomainSolver<  RegionMesh<LinearTetra> >& electricModel,
                               boost::shared_ptr< HeartIonicSolver<  RegionMesh<LinearTetra> > > ionicModel,
                               HeartMonodomainData& data)

{
    typedef HeartUpdatedMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    bool verbose = (M_heart_fct->M_comm->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type uVecRep (electricModel.solutionTransmembranePotential(), Repeated);
    ionicModel->updateRepeated();
    VectorElemental elvec_Iapp ( electricModel.potentialFESpace().fe().nbFEDof(), 2 ),
                    elvec_u ( electricModel.potentialFESpace().fe().nbFEDof(), 1 ),
                    elvec_Iion ( electricModel.potentialFESpace().fe().nbFEDof(), 1 );

    for (UInt iVol = 0; iVol < electricModel.potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        electricModel.potentialFESpace().fe().updateJacQuadPt ( electricModel.potentialFESpace().mesh()->volumeList ( iVol ) );
        elvec_Iapp.zero();
        elvec_u.zero();
        elvec_Iion.zero();

        UInt eleIDu = electricModel.potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) electricModel.potentialFESpace().fe().nbFEDof();

        //! Filling local elvec_u with potential values in the nodes
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int  ig = electricModel.potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u.vec() [ iNode ] = uVecRep[ig];
        }

        ionicModel->updateElementSolution (eleIDu);
        ionicModel->computeIonicCurrent (data.membraneCapacitance(), elvec_Iion, elvec_u, electricModel.potentialFESpace() );

        //! Computing the current source of the righthand side, repeated
        source (M_heart_fct->stimulus(),
                elvec_Iapp,
                electricModel.potentialFESpace().fe(),
                data.time(),
                0);
        source (M_heart_fct->stimulus(),
                elvec_Iapp,
                electricModel.potentialFESpace().fe(),
                data.time(),
                1);

        //! Assembling the righthand side
        for ( UInt i = 0 ; i < electricModel.potentialFESpace().fe().nbFEDof() ; i++ )
        {
            Int  ig = electricModel.potentialFESpace().dof().localToGlobalMap ( eleIDu, i );
            rhs.sumIntoGlobalValues (ig, (data.conductivityRatio() * elvec_Iapp.vec() [i] +
                                          elvec_Iapp.vec() [i + nbNode]) /
                                     (1 + data.conductivityRatio() ) + data.volumeSurfaceRatio() * elvec_Iion.vec() [i] );
        }
    }
    rhs.globalAssemble();
    Real coeff = data.volumeSurfaceRatio() * data.membraneCapacitance() / data.timeStep();
    vector_Type tmpvec (electricModel.solutionTransmembranePotential() );
    tmpvec *= coeff;
    rhs += electricModel.massMatrix() * tmpvec;
    MPI_Barrier (MPI_COMM_WORLD);
    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
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

    ElectroMech electroMech ( argc, argv, Comm );
    electroMech.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
