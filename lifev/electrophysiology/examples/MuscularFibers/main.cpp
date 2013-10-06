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
    @brief Generation muscular fiber architecture on the ventricles by solving a laplacian problem.

    @date 10−2013
    @author Simone Palamara <palamara.simone@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>
    @mantainer Simone Palamara <palamara.simone@gmail.com>
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



#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

// ---------------------------------------------------------------
// In order to use the ETA framework, a special version of the
// FESpace structure must be used. It is called ETFESpace and
// has basically the same role as the FESpace.
// ---------------------------------------------------------------

#include <lifev/eta/fem/ETFESpace.hpp>


// ---------------------------------------------------------------
// The most important file to include is the Integrate.hpp file
// which contains all the definitions required to perform the
// different integrations.
// ---------------------------------------------------------------

#include <lifev/eta/expression/Integrate.hpp>


// ---------------------------------------------------------------
// Finally, we include shared pointer from boost since we use
// them explicitly in this tutorial.
// ---------------------------------------------------------------

#include <boost/shared_ptr.hpp>

#include <lifev/core/util/HeartUtility.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <lifev/electrophysiology/solver/ElectroMonodomainSolver.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/fem/GradientRecovery.hpp>
#include <sys/stat.h>
using std::cout;
using std::endl;
using namespace LifeV;

namespace
{
    typedef LifeV::Real real_t;

    static real_t bcDirichletZero(const real_t&, const real_t&, const real_t&, const real_t&, const LifeV::ID&)
    {
       return 0.;
    }

    static real_t bcDirichletOne(const real_t&, const real_t&, const real_t&, const real_t&, const LifeV::ID&)
    {
       return 1;
    }

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

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef  ElectroMonodomainSolver< RegionMesh<LinearTetra> >		physicalSolver_Type;

    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr< vector_Type > vectorPtr_Type;

    typedef FESpace< mesh_Type, MapEpetra >    feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>    feSpacePtr_Type;
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >    ETfeSpace1D_Type;
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 >    ETfeSpace3D_Type;
    typedef boost::shared_ptr<ETfeSpace1D_Type>    ETfeSpacePtr1D_Type;
    typedef boost::shared_ptr<ETfeSpace3D_Type>    ETfeSpacePtr3D_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) > function_Type;




    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList List = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }
    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( Comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }
    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    std::string meshName = List.get ("mesh_name", "lid16.mesh");
    std::string meshPath = List.get ("mesh_path", "./");


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Read the mesh...";
    }

    meshPtr_Type fullMeshPtr (new mesh_Type ( Comm ) );
    MeshUtility::fillWithMesh(fullMeshPtr,0,meshName,meshPath);


    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!"<<std::endl;
    }

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    const string dataSection = ("problem");
    const string name = ("dirichlet");
    GetPot dataFile (data_file_name);

    if (Comm->MyPID() == 0)
    {
        std::cout << " -- Building ETFESpaces ... " << std::flush;
    }

    ETfeSpacePtr1D_Type uSpace ( new ETfeSpace1D_Type (fullMeshPtr, &feTetraP1, Comm) );
    feSpacePtr_Type uFESpace ( new feSpace_Type (fullMeshPtr, "P1", 1, Comm) );
    feSpacePtr_Type sFESpace( new feSpace_Type (fullMeshPtr, "P1", 3, Comm) );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!"<<std::endl;
    }


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "-- Reading bc ..";
    }

    boost::shared_ptr<LifeV::BCHandler> bcs(new LifeV::BCHandler());

    LifeV::BCFunctionBase zero(bcDirichletZero);
    LifeV::BCFunctionBase one(bcDirichletOne);
    std::vector<LifeV::ID> compx(1, 0), compy(1, 1), compz(1, 2);
    bcs->addBC("boundaryDirichletEndo", 0, LifeV::Essential, LifeV::Full, zero,1);
    bcs->addBC("boundaryDirichletEpi", 1, LifeV::Essential, LifeV::Full, one,1);
    bcs->addBC("boundaryNeumannBase", 2, LifeV::Natural, LifeV::Full, zero,1);
    bcs->bcUpdate ( *uFESpace->mesh(),uFESpace->feBd(), uFESpace->dof() );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!"<<std::endl;
    }



    if (Comm->MyPID() == 0)
    {
        std::cout << " -- Defining the matrix and rhs ... " << std::flush;
    }

    matrixPtr_Type systemMatrix (new matrix_Type ( uSpace->map() ) );
    *systemMatrix*=0;
    vectorPtr_Type Rhs(new vector_Type ( uSpace->map() ) );
    *Rhs*=0;
    Rhs -> globalAssemble();

  {
        using namespace ExpressionAssembly;

        // ---------------------------------------------------------------
        // We can now proceed with assembly. The next instruction
        // assembles the laplace operator.
        //
        // The first argument of the integrate function indicates that the
        // integration is done on the elements of the mesh located in the
        // ETFESpace defined earlier.
        //
        // The second argument is simply the quadrature rule to be used.
        //
        // The third argument is the finite element space of the test
        // functions.
        //
        // The fourth argument is the finite element space of the trial
        // functions (those used to represent the solution).
        //
        // The last argument is the expression to be integrated, i.e.
        // that represents the weak formulation of the problem. The
        // keyword phi_i stands for a generic test function and phi_j
        // a generic trial function. The function grad applied to them
        // indicates that the gradient is considered and the dot function
        // indicates a dot product between the two gradients. The
        // expression to be integrated is then the dot product between
        // the gradient of the test function and the gradient of the trial
        // function. This corresponds to the left hand side of the weak
        // formulation of the Laplace problem.
        //
        // Finally, the operator >> indicates that the result of the
        // integration must be added to the systemMatrix.
        // ---------------------------------------------------------------

        integrate (  elements (uSpace->mesh() ),
                     quadRuleTetra4pt,
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
                  )
                >> systemMatrix;
    }


    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!"<<std::endl;
    }

    if (Comm->MyPID() == 0)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }

    systemMatrix->globalAssemble();

    if (Comm->MyPID() == 0)
    {
        std::cout << " done ! " << std::endl;
    }

    bcManage ( *systemMatrix, *Rhs, *uSpace->mesh(), uSpace->dof(), *bcs, uFESpace->feBd(), 1.0, 0.0 );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating Preconditioner...\n";
    }

    //Preconditioner
    typedef LifeV::Preconditioner             basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack       prec_Type;
    typedef boost::shared_ptr<prec_Type>      precPtr_Type;

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating linear solver...\n";
    }

    //linear solver
    Teuchos::RCP< Teuchos::ParameterList > belosList3 = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList3 = Teuchos::getParametersFromXmlFile ( "ParamList.xml" );

    LinearSolver linearSolver;
    linearSolver.setCommunicator(Comm);
//    .setCommunicator ( Comm );
    linearSolver.setParameters ( *belosList3 );
    linearSolver.setPreconditioner ( precPtr );

    vectorPtr_Type solution( new vector_Type( uFESpace -> map() ) );

    vectorPtr_Type sx(new vector_Type( uSpace -> map() ) );
    vectorPtr_Type sy(new vector_Type( uSpace -> map() ) );
    vectorPtr_Type sz(new vector_Type( uSpace -> map() ) );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Exporter...\n";
    }

    ExporterHDF5< RegionMesh <LinearTetra> > exporter;
      exporter.setMeshProcId (fullMeshPtr , Comm->MyPID() );
      exporter.setPrefix ("laplacian");
      exporter.setPostDir(problemFolder);
      exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential", uFESpace,
                             solution, UInt (0) );
      exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "sx", uFESpace,
                             sx, UInt (0) );
      exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "sy", uFESpace,
                             sy, UInt (0) );
      exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "sz", uFESpace,
                             sz, UInt (0) );



      if ( Comm->MyPID() == 0 )
      {
          std::cout << "Solve System...\n";
      }
    	linearSolver.setOperator(systemMatrix);
      linearSolver.setRightHandSide(Rhs);
      linearSolver.solve(solution);



	  std::cout << "\nRecovering gradient ... ";
      *sx = GradientRecovery::ZZGradient(uSpace, *solution, 0);
      *sy = GradientRecovery::ZZGradient(uSpace, *solution, 1);
      *sz = GradientRecovery::ZZGradient(uSpace, *solution, 2);
	  std::cout << " Done! \n";

	  exporter.postProcess ( 0 );





      vectorPtr_Type rbSheet( new vector_Type( sFESpace -> map() ) );
//      vectorPtr_Type fx(new vector_Type( uSpace -> map() ) );
//      vectorPtr_Type fy(new vector_Type( uSpace -> map() ) );
//      vectorPtr_Type fz(new vector_Type( uSpace -> map() ) );
      vectorPtr_Type rbFiber( new vector_Type( sFESpace -> map() ) );
      vectorPtr_Type cl( new vector_Type( sFESpace -> map() ) );
      vectorPtr_Type projection( new vector_Type( sFESpace -> map() ) );

		int n = (*rbSheet).epetraVector().MyLength();
		int d = n / 3;

	    if ( Comm->MyPID() == 0 )
	    {
	        std::cout << "Assembling sheets...\n";
	    }
		for ( int l (0); l < d; l++)
		{
			int i = (*rbSheet).blockMap().GID (l);
			int j = (*rbSheet).blockMap().GID (l + d);
			int k = (*rbSheet).blockMap().GID (l + 2 * d);

			(*rbSheet) [i] = (*sx) [i];
			(*rbSheet) [j] = (*sy) [i];
			(*rbSheet) [k] = (*sz) [i];
		}

	    if ( Comm->MyPID() == 0 )
	    {
	        std::cout << "Normalizing sheets...\n";
	    }

		HeartUtility::normalize(*rbSheet);

		//human biventricular
		//Real cx = -0.655333989927163;
		//Real cy = 0.699986220544843;
		//Real cz = -0.283825038876930;

		//human ventricular
		Real cx = -21.3;
		Real cy = 44.33;
		Real cz = 33.36;
		Real norm=std::sqrt(cx*cx+cy*cy+cz*cz);
		cx=cx/norm;
		cy=cy/norm;
		cz=cz/norm;
		//idealized
		//Real cx = 0.0;
		//Real cy = 0.0;
		//Real cz = 1.0;

		//idealized biventricular
		//Real cx = -0.585264869348945;
		//Real cy = 0.724137970380166;
		//Real cz = -0.364813969797835;


	    if ( Comm->MyPID() == 0 )
	    {
	        std::cout << "Creating fibers and projection...\n";
	    }
		for ( int l (0); l < d; l++)
		{
			int i = (*rbSheet).blockMap().GID (l);
			int j = (*rbSheet).blockMap().GID (l + d);
			int k = (*rbSheet).blockMap().GID (l + 2 * d);

			(*rbFiber) [i] = (*rbSheet) [j];
			(*rbFiber) [j] = -(*rbSheet) [i];
			(*rbFiber) [k] = 0.0;

			Real cdot = cx * (*rbSheet) [i] + cy * (*rbSheet) [j] + cz * (*rbSheet) [k];

//			if( cdot > 0 )
//			{
//				(*rbSheet) [i] = -(*rbSheet) [i];
//				(*rbSheet) [j] = -(*rbSheet) [j];
//				(*rbSheet) [k] = -(*rbSheet) [k];
//
//			}
//
//			cdot = cx * (*rbSheet) [i] + cy * (*rbSheet) [j] + cz * (*rbSheet) [k];

			(*projection) [i] = cx - cdot *  (*rbSheet) [i];
			(*projection) [j] = cy - cdot *  (*rbSheet) [j];
			(*projection) [k] = cz - cdot *  (*rbSheet) [k];

			(*cl) [i] = cx;
			(*cl) [j] = cy;
			(*cl) [k] = cz;


			(*rbFiber) [i] = (*rbSheet) [j];
			(*rbFiber) [j] = -(*rbSheet) [i];
			(*rbFiber) [k] = 0.0;


			Real scalarp;
			scalarp = (*projection) [i] * (*rbFiber) [i]
                      + (*projection) [j] * (*rbFiber) [j]
                      + (*projection) [k] * (*rbFiber) [k];

			if( scalarp == 1.0 )
			{
				(*rbFiber) [i] = 0;
				(*rbFiber) [j] = -(*rbSheet) [k];
				(*rbFiber) [k] = (*rbSheet) [j];
			}

		}

	    if ( Comm->MyPID() == 0 )
	    {
	        std::cout << "Normalizing fibers and projection...\n";
	    }





        ExporterHDF5< mesh_Type > exporter2;
        exporter2.setMeshProcId ( fullMeshPtr, Comm -> MyPID() );
        exporter2.setPrefix ("sheets");
        exporter2.setPostDir(problemFolder);
        exporter2.addVariable ( ExporterData<mesh_Type>::VectorField,  "sheet", sFESpace, rbSheet, UInt (0) );
//        exporter2.addVariable ( ExporterData<mesh_Type>::ScalarField,  "fx", uFESpace,
//                               fx, UInt (0) );
//        exporter2.addVariable ( ExporterData<mesh_Type>::ScalarField,  "fy", uFESpace,
//                               fy, UInt (0) );
//        exporter2.addVariable ( ExporterData<mesh_Type>::ScalarField,  "fz", uFESpace,
//                               fz, UInt (0) );
        exporter2.addVariable ( ExporterData<mesh_Type>::VectorField,  "fiber", sFESpace, rbFiber, UInt (0) );
        exporter2.addVariable ( ExporterData<mesh_Type>::VectorField,  "centerline", sFESpace, cl, UInt (0) );
        exporter2.addVariable ( ExporterData<mesh_Type>::VectorField,  "projection", sFESpace, projection, UInt (0) );
        exporter2.postProcess (0);


		HeartUtility::normalize(*rbFiber);
		HeartUtility::normalize(*projection);
        exporter2.postProcess (1);



	    if ( Comm->MyPID() == 0 )
	    {
	        std::cout << "Find out angles...\n";
	    }
		for ( int l (0); l < d; l++)
		{
			int i = (*rbSheet).blockMap().GID (l);
			int j = (*rbSheet).blockMap().GID (l + d);
			int k = (*rbSheet).blockMap().GID (l + 2 * d);

//			Real scalarp;
//			scalarp = (*projection) [i] * (*rbFiber) [i]
//                      + (*projection) [j] * (*rbFiber) [j]
//                      + (*projection) [k] * (*rbFiber) [k];


			//fiber =  sheet x projection
			(*rbFiber) [i] = (*rbSheet) [j]  * (*projection) [k] - (*rbSheet) [k]  * (*projection) [j];
			(*rbFiber) [j] = (*rbSheet) [k]  * (*projection) [i] - (*rbSheet) [i]  * (*projection) [k];
			(*rbFiber) [k] = (*rbSheet) [i]  * (*projection) [j] - (*rbSheet) [j]  * (*projection) [i];


//			if(scalarp < 0)
//			{
//				(*rbFiber) [i] = -(*rbFiber) [i];
//				(*rbFiber) [j] = -(*rbFiber) [j];
//				(*rbFiber) [k] = -(*rbFiber) [k];
//
//			}

		}

		HeartUtility::normalize(*rbFiber);

        exporter2.postProcess (2);
        Real epi_angle = List.get ("epi_angle", -60.0);
        Real endo_angle = List.get ("endo_angle", 60.0);

	    if ( Comm->MyPID() == 0 )
	    {
	        std::cout << "Find out angles...\n";
	    }
		for ( int l (0); l < d; l++)
		{
			int i = (*rbSheet).blockMap().GID (l);
			int j = (*rbSheet).blockMap().GID (l + d);
			int k = (*rbSheet).blockMap().GID (l + 2 * d);

			Real scalarp;
			scalarp = (*projection) [i] * (*rbFiber) [i]
                      + (*projection) [j] * (*rbFiber) [j]
                      + (*projection) [k] * (*rbFiber) [k];


//			Real epi_angle =  -50.0;
//			Real endo_angle =  50.0;
			Real p = 3.14159265358979;
			Real teta1 = p * epi_angle / 180; //67.5 degrees
			Real teta2 = p * endo_angle / 180;
			Real m = (teta1 - teta2 );
			Real q = teta2;// - m * (*solution)[i];
			Real teta;

			teta = m * (*solution)[i] + q;
			(*solution)[i] = teta;

			Real s01 = (*rbSheet)[i];
			Real s02 = (*rbSheet)[j];
			Real s03 = (*rbSheet)[k];
			Real f01 = (*rbFiber)[i];
			Real f02 = (*rbFiber)[j];
			Real f03 = (*rbFiber)[k];

			Real R11 = std::cos(teta) + s01 * s01 * ( 1 - std::cos(teta) );
			Real R12 = s01 * s02 *  ( 1 - std::cos(teta) ) - s03 * std::sin(teta);
			Real R13 = s01 * s03 *  ( 1 - std::cos(teta) ) + s02 * std::sin(teta);

			Real R21 = s02 * s01 *  ( 1 - std::cos(teta) ) + s03 * std::sin(teta);
			Real R22 = std::cos(teta) + s02 * s02 * ( 1 - std::cos(teta) );
			Real R23 = s02 * s03 *  ( 1 - std::cos(teta) ) - s01 * std::sin(teta);

			Real R31 = s03 * s01 *  ( 1 - std::cos(teta) ) - s02 * std::sin(teta);
			Real R32 = s03 * s02 *  ( 1 - std::cos(teta) ) + s01 * std::sin(teta);
			Real R33 = std::cos(teta) + s03 * s03 * ( 1 - std::cos(teta) );

			Real ca = std::cos(teta);
			Real sa = std::sin(teta);
//			(*rbFiber) [i] = R11 * (*rbFiber) [i] + R12 * (*rbFiber) [j] + R13 * (*rbFiber) [k];
//			(*rbFiber) [j] = R21 * (*rbFiber) [i] + R22 * (*rbFiber) [j] + R23 * (*rbFiber) [k];
//			(*rbFiber) [k] = R31 * (*rbFiber) [i] + R32 * (*rbFiber) [j] + R33 * (*rbFiber) [k];
//
//
//			(*rbFiber) [i]=s01*scalarp+(f01*(s02*s02+s03*s03)-s01*(s02*f02+s03*f03))*ca+(-s03*f02+s02*f03)*sa;
//			(*rbFiber) [j]=s02*scalarp+(f02*(s01*s01+s03*s03)-s02*(s01*f01+s03*f03))*ca+(s03*f01-s01*f03)*sa;
//			(*rbFiber) [k]=s03*scalarp+(f03*(s01*s01+s02*s02)-s03*(s01*f01+s02*f02))*ca+(-s02*f01+s01*f02)*sa;


			Real W11 = 0.0;
			Real W12 = -s03;
			Real W13 = s02;
			Real W21 = s03;
			Real W22 = 0.0;
			Real W23 = -s01;
			Real W31 = -s02;
			Real W32 = s01;
			Real W33 = 0.0;
			Real sa2 = 2.0*std::sin(0.5*teta)*std::sin(0.5*teta);
//
			R11 = 1.0 + sa * W11 + sa2 * ( s01 * s01 - 1.0 );
			R12 = 0.0 + sa * W12 + sa2 * ( s01 * s02 );
			R13 = 0.0 + sa * W13 + sa2 * ( s01 * s03 );
			R21 = 0.0 + sa * W21 + sa2 * ( s02 * s01 );
			R22 = 1.0 + sa * W22 + sa2 * ( s02 * s02 - 1.0 );
			R23 = 0.0 + sa * W23 + sa2 * ( s02 * s03 );
			R31 = 0.0 + sa * W31 + sa2 * ( s03 * s01 );
			R32 = 0.0 + sa * W32 + sa2 * ( s03 * s02 );
			R33 = 1.0 + sa * W33 + sa2  * ( s03 * s03 - 1.0 );

			(*rbFiber) [i] = R11 * f01 + R12 * f02 + R13 * f03;
			(*rbFiber) [j] = R21 * f01 + R22 * f02 + R23 * f03;
			(*rbFiber) [k] = R31 * f01 + R32 * f02 + R33 * f03;



		}

		//HeartUtility::normalize(*rbFiber);

        exporter2.postProcess (3);
        exporter2.closeFile();

        exporter.postProcess (1);
        exporter.closeFile();

        ExporterHDF5< mesh_Type > exporterFibers;
        exporterFibers.setMeshProcId ( fullMeshPtr, Comm -> MyPID() );
        exporterFibers.setPostDir(problemFolder);
        std::string outputFiberFileName = List.get ("output_fiber_filename", "FiberDirection");
        exporterFibers.setPrefix (outputFiberFileName);
        exporterFibers.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", sFESpace, rbFiber, UInt (0) );
        exporterFibers.postProcess (0);
        exporterFibers.closeFile();
        ExporterHDF5< mesh_Type > exporterSheets;
        exporterSheets.setMeshProcId ( fullMeshPtr, Comm -> MyPID() );
        exporterSheets.setPostDir(problemFolder);
        std::string outputSheetsFileName = List.get ("output_sheets_filename", "SheetsDirection");
        exporterSheets.setPrefix (outputSheetsFileName);
        exporterSheets.addVariable ( ExporterData<mesh_Type>::VectorField,  "sheets", sFESpace, rbSheet, UInt (0) );
        exporterSheets.postProcess (0);
        exporterSheets.closeFile();


    // ---------------------------------------------------------------
    // We finalize the MPI session if MPI was used
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif



}
