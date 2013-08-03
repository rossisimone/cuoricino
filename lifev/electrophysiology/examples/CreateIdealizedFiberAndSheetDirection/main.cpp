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

#include <lifev/core/array/VectorEpetra.hpp>


// ---------------------------------------------------------------
// Finally, we include shared pointer from boost since we use
// them explicitly in this tutorial.
// ---------------------------------------------------------------

#include <boost/shared_ptr.hpp>

#include <lifev/core/util/HeartUtility.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>
// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. For clarity, we also
// make two typedefs for the mesh type and matrix type.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef boost::shared_ptr< mesh_Type > meshPtr_Type;
typedef VectorEpetra vector_Type;
typedef boost::shared_ptr< vector_Type > vectorPtr_Type;
typedef boost::function < Real (const Real& /*t*/,
                                const Real &   x,
                                const Real &   y,
                                const Real& z,
                                const ID&   i ) >   function_Type;
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

Real sheets (const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& i)
{
	Real Z0 = Z - 6.4;
//
//	Real  varphi = std::asin(Z0 / R);
	Real ae1 = 2.8 + 0.75;
	Real be1 = 2.8 + 0.75 ;
	Real ce1 = 0.8*8.0 + 0.3;

	Real N = std::sqrt(X*X/ae1/ae1/ae1/ae1+Y*Y/be1/be1/be1/be1+Z0*Z0/ce1/ce1/ce1/ce1);
	Real compx = X / ae1 / ae1 / N;
	Real compy = Y / be1 / be1 / N;
	Real compz = Z0 / ce1 / ce1 /N;


	if(i==0) return compx;
//		return a*std::sqrt( ( sigma*sigma - 1.0 ) * ( 1.0 -  tau * tau) ) * std::cos(phi);
	else if(i ==1) return compy;
	//	return a*std::sqrt( ( sigma*sigma - 1.0 ) * ( 1.0 -  tau * tau) ) * std::sin(phi);
	else if(i==2) return compz;
	//	return a * sigma * tau;
	else
		return 0.0;
}

Real fibers (const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& i)
{


	//Compute sheets
	Real Z0 = Z - 6.4;
	Real R = std::sqrt(X*X+Y*Y+Z0*Z0);
	Real phi   = std::acos( Z0 / R );

	Real p = 3.14159265358979;


	Real R0 = std::sqrt(X*X+Y*Y);

	Real  varphi = std::asin(Z0 / R);
	Real ae1 = 2.8 + 0.75;
	Real be1 = 2.8 + 0.75 ;
	Real ce1 = 0.8*8.0 + 0.3;
	Real Re = std::sqrt( 1.0 / ( ( std::cos(phi) * std::cos(varphi) / ae1 ) * ( std::cos(phi) * std::cos(varphi) / ae1 )
							   + ( std::sin(phi) * std::cos(varphi) / be1 ) * ( std::sin(phi) * std::cos(varphi) / be1 )
							   + ( std::sin(varphi) / ce1 ) * ( std::sin(varphi) / ce1 ) ) );
	Real dR = Re - R;
	//compute fibers
	Real N = std::sqrt(X*X/ae1/ae1/ae1/ae1+Y*Y/be1/be1/be1/be1+Z0*Z0/ce1/ce1/ce1/ce1);
	Real s01 = X / ae1 / ae1 / N;
	Real s02 = Y / be1 / be1 / N;
	Real s03 = Z0 / ce1 / ce1 /N;


	Real f001;
	Real f002;
	Real f003;

	f001 = s02;
	f002 = -s01;
	f003 = 0;

//	Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
//	Real epi_angle = parameterList.get ("epi_angle", -60.0);
//	Real endo_angle = parameterList.get ("endo_angle", 60.0);
	Real epi_angle =  -60.0;
	Real endo_angle =  60.0;
//	Real teta = - 1.2566; //72 degrees
	Real dR1 = -0.87177357881155; //on the epi
	Real teta1 = p * epi_angle / 180; //67.5 degrees
	Real dR2 = 0.741440817851517; //on the endo
	Real teta2 = p * endo_angle / 180;
	Real m = (teta1 - teta2 ) / (dR1 - dR2);
	Real q = teta1 - m * dR1;
	Real teta;

	teta = m * dR + q;
	if( Z > 6.4 ) teta = teta * ( 1 - 0.5 / 1.6 * Z0 );


	Real R11 = std::cos(teta) + s01 * s01 * ( 1 - std::cos(teta) );
	Real R12 = s01 * s02 *  ( 1 - std::cos(teta) ) - s03 * std::sin(teta);
	Real R13 = s01 * s03 *  ( 1 - std::cos(teta) ) + s02 * std::sin(teta);

	Real R21 = s02 * s01 *  ( 1 - std::cos(teta) ) + s03 * std::sin(teta);
	Real R22 = std::cos(teta) + s02 * s02 * ( 1 - std::cos(teta) );
	Real R23 = s02 * s03 *  ( 1 - std::cos(teta) ) - s01 * std::sin(teta);

	Real R31 = s03 * s01 *  ( 1 - std::cos(teta) ) - s02 * std::sin(teta);
	Real R32 = s03 * s02 *  ( 1 - std::cos(teta) ) + s01 * std::sin(teta);
	Real R33 = std::cos(teta) + s03 * s03 * ( 1 - std::cos(teta) );


	Real compz(0.0);
	Real compx(0.0);
	Real compy(0.0);
	compx = R11 * f001 + R12 * f002 + R13 * f003;
	compy = R21 * f001 + R22 * f002 + R23 * f003;
	compz = R31 * f001 + R32 * f002 + R33 * f003;

	//compx = phi;

	if(i==0) return -compx;
	else if(i ==1) return -compy;
	else if(i==2) return -compz;
	else
		return 0.0;
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







    // ---------------------------------------------------------------
    // We start now the assembly of the matrix.
    // ---------------------------------------------------------------



    // ---------------------------------------------------------------
    // To use the ETA framework, it is mandatory to use a special
    // namespace, called ExpressionAssembly. This namespace is useful
    // to avoid collisions with keywords used for the assembly. A
    // special scope is opened to keep only that part of the code
    // in the ExpressionAssembly namespace.
    // ---------------------------------------------------------------

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





    // ---------------------------------------------------------------
    // As we are already done with the assembly of the matrix, we
    // finalize it to be able to work on it, e.g. to solve a linear
    // system.
    // ---------------------------------------------------------------

    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);


        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > vSpace;
        ( new FESpace< mesh_Type, MapEpetra > ( meshPart, "P1", 3, Comm ) );

        boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > coarse
        ( new FESpace< mesh_Type, MapEpetra > ( meshPart, "P1", 3, Comm ) );
        vectorPtr_Type sheet( new vector_Type( coarse -> map() ) );
        function_Type s = &sheets;
        coarse -> interpolate(static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type>(s),
        				*sheet, 0.0);
        HeartUtility::normalize(*sheet);

        vectorPtr_Type fiber( new vector_Type( coarse -> map() ) );
        function_Type f = &fibers;
        coarse -> interpolate(static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type>(f),
        				*fiber, 0.0);

        HeartUtility::normalize(*fiber);


//        int n1 = sheet -> epetraVector().MyLength();
//        int d1 = n1 / 3;
//        int i (0);
//        int j (0);
//        int k (0);
//
//        for ( int l (0); l < d1; l++)
//        {
//            i = sheet->blockMap().GID (l);
//            j = sheet->blockMap().GID (l + d1);
//            k = sheet->blockMap().GID (l + 2 * d1);
//            Real norm = std::sqrt((*sheet)[i]*(*sheet)[i]+(*sheet)[j]*(*sheet)[j]+(*sheet)[k]*(*sheet)[k]);
//            if(norm!=0)
//            {
//            (*sheet)[i] = (*sheet)[i]/norm;
//            (*sheet)[j] = (*sheet)[j]/norm;
//            (*sheet)[k] = (*sheet)[k]/norm;
//
//            }
//            else
//            {
//                (*sheet)[i] = 0.0;
//                (*sheet)[j] = 0.0;
//                (*sheet)[k] = 1.0;
//            }
//        }
        ExporterHDF5< RegionMesh <LinearTetra> > exporter1;
          exporter1.setMeshProcId ( meshPart, Comm->MyPID() );
          exporter1.setPrefix (parameterList.get ("fiber_output_file", "fibers"));
          exporter1.addVariable ( ExporterData<mesh_Type>::VectorField,  parameterList.get ("fiber_output_field", "fibers"), coarse,
                                 fiber, UInt (0) );
          exporter1.postProcess ( 0 );
          exporter1.closeFile();

        ExporterHDF5< RegionMesh <LinearTetra> > exporter2;
          exporter2.setMeshProcId ( meshPart, Comm->MyPID() );
          exporter2.setPrefix (parameterList.get ("sheet_output_file", "sheets"));
          exporter2.addVariable ( ExporterData<mesh_Type>::VectorField,  parameterList.get ("sheet_output_field", "sheets"), coarse,
                                 sheet, UInt (0) );
          exporter2.postProcess ( 0 );
          exporter2.closeFile();


    // ---------------------------------------------------------------
    // We finalize the MPI session if MPI was used
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif




}

