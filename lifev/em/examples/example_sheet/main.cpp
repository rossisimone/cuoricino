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
    @brief Tutorial introducing the expression assembly

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 28-06-2012

    In this first tutorial, we assemble the matrix
    associated to a scalar laplacian problem. The basics
    of the ETA module are explained and are pushed further
    in the next tutorials.

    ETA stands Expression Template Assembly, in reference
    to the metaprogramming technique used.

 */

// ---------------------------------------------------------------
// We include here the MPI headers for the parallel computations.
// The specific "pragma" instructions are used to avoid warning
// coming from the MPI library, that are not useful to us.
// ---------------------------------------------------------------

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


// ---------------------------------------------------------------
// We include then the required headers from LifeV. First of all,
// the definition file and mesh related files. We also include
// the MatrixEpetra since this is the kind of object that we want
// to assemble.
// ---------------------------------------------------------------

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

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <sys/stat.h>
#include <lifev/core/fem/GradientRecovery.hpp>
// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. For clarity, we also
// make two typedefs for the mesh type and matrix type.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef boost::shared_ptr< mesh_Type > meshPtr_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef boost::shared_ptr< matrix_Type > matrixPtr_Type;
typedef VectorEpetra vector_Type;
typedef boost::shared_ptr< vector_Type > vectorPtr_Type;
typedef BCHandler                                          bc_Type;
typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
typedef  StructuralOperator< RegionMesh<LinearTetra> >		physicalSolver_Type;
typedef BCInterface3D< bc_Type, physicalSolver_Type >              bcInterface_Type;
typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;


// ---------------------------------------------------------------
// We start the programm by the definition of the communicator
// (as usual) depending on whether MPI is available or not. We
// also define a boolean to allow only one process to display
// messages.
// ---------------------------------------------------------------

Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real bcOne (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  1.0;
}

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);


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
    // ---------------------------------------------------------------
    // The next step is to build the mesh. We use here a structured
    // cartesian mesh over the square domain (-1,1)x(-1,1)x(-1,1).
    // The mesh is the partitioned for the parallel computations and
    // the original mesh is deleted.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }
    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");

    meshPtr_Type meshPart (new mesh_Type ( Comm ) );
    MeshUtility::fillWithMesh(meshPart,meshName,meshPath);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We define now the ETFESpace that we need for the assembly.
    // Remark that we use a shared pointer because other structures
    // will require this ETFESpace to be alive. We can also observe
    // that the ETFESpace has more template parameters than the
    // classical FESpace (this is the main difference). The 3
    // indicates that the problem is in 3D while the 1 indicate that
    // the unknown is scalar.
    //
    // After having constructed the ETFESpace, we display the number
    // of degrees of freedom of the problem.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building ETFESpaces ... " << std::flush;
    }

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPart, "P1", 1, Comm) );
    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > sSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPart, &feTetraP1, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > sFESpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPart, "P1", 3, Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << uSpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // The matrix is then defined using the map of the FE space.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Defining the matrix ... " << std::flush;
    }

    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uSpace->map() ) );

    *systemMatrix *= 0.0;



    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We start now the assembly of the matrix.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembling the Laplace matrix ... " << std::flush;
    }


    // ---------------------------------------------------------------
    // To use the ETA framework, it is mandatory to use a special
    // namespace, called ExpressionAssembly. This namespace is useful
    // to avoid collisions with keywords used for the assembly. A
    // special scope is opened to keep only that part of the code
    // in the ExpressionAssembly namespace.
    // ---------------------------------------------------------------

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

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }




    // ---------------------------------------------------------------
    // As we are already done with the assembly of the matrix, we
    // finalize it to be able to work on it, e.g. to solve a linear
    // system.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }

    systemMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating Boundary Conditions...\n";
    }
    //-----------------------
    //  BOundary Conditions
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    bcInterfacePtr_Type                     BC( new bcInterface_Type() );
    BC->createHandler();
    BC->fillHandler ( data_file_name, "problem" );
    BC -> handler() -> bcUpdate ( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

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


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating rhs...\n";
    }
    //Create right hand side
    vectorPtr_Type rhs(new vector_Type( uSpace -> map() ) );
    *rhs *= 0.0;
    rhs -> globalAssemble();


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Applying BC...\n";
    }
    bcManage ( *systemMatrix, *rhs, *uSpace->mesh(), uSpace->dof(), *BC -> handler(), uFESpace->feBd(), 1.0, 0.0 );




    linearSolver.setOperator(systemMatrix);
    vectorPtr_Type solution( new vector_Type( uFESpace -> map() ) );

    vectorPtr_Type sx(new vector_Type( uSpace -> map() ) );
    vectorPtr_Type sy(new vector_Type( uSpace -> map() ) );
    vectorPtr_Type sz(new vector_Type( uSpace -> map() ) );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Exporter...\n";
    }

    ExporterHDF5< RegionMesh <LinearTetra> > exporter;
      exporter.setMeshProcId ( meshPart, Comm->MyPID() );
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
      linearSolver.setRightHandSide(rhs);
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
		Real cx = -0.655333989927163;
		Real cy = 0.699986220544843;
		Real cz = -0.283825038876930;

		//idealized
//		Real cx = 0.0;
//		Real cy = 0.0;
//		Real cz = 1.0;


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
        exporter2.setMeshProcId ( meshPart, Comm -> MyPID() );
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
        Real epi_angle = parameterList.get ("epi_angle", -60.0);
        Real endo_angle = parameterList.get ("endo_angle", 60.0);

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
        exporterFibers.setMeshProcId ( meshPart, Comm -> MyPID() );
        exporterFibers.setPostDir(problemFolder);
        std::string outputFiberFileName = parameterList.get ("output_fiber_filename", "FiberDirection");
        exporterFibers.setPrefix (outputFiberFileName);
        exporterFibers.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", sFESpace, rbFiber, UInt (0) );
        exporterFibers.postProcess (0);
        exporterFibers.closeFile();
        ExporterHDF5< mesh_Type > exporterSheets;
        exporterSheets.setMeshProcId ( meshPart, Comm -> MyPID() );
        exporterSheets.setPostDir(problemFolder);
        std::string outputSheetsFileName = parameterList.get ("output_sheets_filename", "SheetsDirection");
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


