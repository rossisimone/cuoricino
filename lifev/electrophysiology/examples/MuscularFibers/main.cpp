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
    @brief 0D test with the FitzHugh Nagumo model.

    @date 01−2013
    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @contributor
    @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
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
    bcs->addBC("boundaryDirichletEndo", 10, LifeV::Essential, LifeV::Full, zero,1);
    bcs->addBC("boundaryDirichletEpi", 20, LifeV::Essential, LifeV::Full, one,1);
    //bcs->addBC("boundaryNeumannBase", 99, LifeV::Natural, LifeV::Full, one,1);
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

    vectorPtr_Type solution(new vector_Type ( uSpace->map() ) );
    *solution*=0;
    ExporterHDF5< RegionMesh <LinearTetra> > exporter;
      exporter.setMeshProcId ( fullMeshPtr, Comm->MyPID() );
      exporter.setPrefix ("laplacian");
      exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "laplacian", uFESpace,
                             solution, UInt (0) );


      if ( Comm->MyPID() == 0 )
      {
          std::cout << "Solve System...\n";
      }


      linearSolver.setRightHandSide(Rhs);
        linearSolver.setOperator(systemMatrix);
      linearSolver.solve(solution);


      exporter.postProcess ( 1 );
      exporter.closeFile();
    //Compute the potential at t0
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
