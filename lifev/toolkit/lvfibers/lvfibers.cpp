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
    @brief Generation muscular fibers and sheets

    @author Simone Rossi <simone.rossi@epfl.ch>
    @maintainer Simone Palamara <palamara.simone@gmail.com>
    @date 31-01-2014


    Generation of the muscular fibers and sheets on a generic
    geometry representing the left or right ventricle, generated
    according to geometrical rules based on anatomical knowledge.
    For more details about the method see [S.Rossi et al,European Journal of
    Mechanics A/Solids (2013), http://dx.doi.org/10.1016/j.euromechsol.2013.10.009]

 */

//#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <boost/shared_ptr.hpp>

#include <Teuchos_ParameterList.hpp>

#include <sys/stat.h>

#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <lifev/core/fem/GradientRecovery.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>



// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. Moreover,
// to make the code more readable, we also make typedefs for the mesh type,
// matrix type, vector type, boundary condition
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra>                        mesh_Type;
typedef boost::shared_ptr< mesh_Type >                 meshPtr_Type;
typedef MatrixEpetra<Real>                             matrix_Type;
typedef boost::shared_ptr< matrix_Type >               matrixPtr_Type;
typedef VectorEpetra                                   vector_Type;
typedef boost::shared_ptr< vector_Type >               vectorPtr_Type;
typedef BCHandler                                      bc_Type;
typedef boost::shared_ptr< bc_Type >                   bcPtr_Type;
typedef StructuralOperator< RegionMesh<LinearTetra> >  physicalSolver_Type;
typedef BCInterface3D< bc_Type, physicalSolver_Type >  bcInterface_Type;
typedef boost::shared_ptr< bcInterface_Type >          bcInterfacePtr_Type;
typedef MeshUtility::MeshTransformer<mesh_Type>        meshTransformer_Type;



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
    std::string problemFolder = commandLine.follow ( "work", 2, "-w", "--work" );

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
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    const string fibercreation_datafile_name = commandLine.follow ("FiberGenerationParamList.xml", 2, "-f", "--file");
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( fibercreation_datafile_name ) );

    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    std::string meshName = commandLine.follow ("default", 2, "-m", "--model") + ".mesh";

    meshPtr_Type meshPart (new mesh_Type ( Comm ) );
    MeshUtility::fillWithMesh (meshPart, meshName, problemFolder);

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPart, "P1", 1, Comm) );
    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > sSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPart, &feTetraP1, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > sFESpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPart, "P1", 3, Comm) );

    //***************************************************//
    // Assembly of the stiffness matrix related to the   //
    // laplacian problem used for the computation of the //
    // sheetlet direction by using the ETA.              //
    //***************************************************//

    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uSpace->map() ) );

    *systemMatrix *= 0.0;
    {
        using namespace ExpressionAssembly;


        integrate (  elements (uSpace->mesh() ),
                     quadRuleTetra4pt,
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
        )
        >> systemMatrix;
    }

    systemMatrix->globalAssemble();

    //****************************************************//
    //  Boundary conditions read from the getPot object   //
    //  related to the file FiberGenerationPreconditioner //
    //****************************************************//

    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("FiberGenerationPreconditioner", 2, "-p", "--prec");

    GetPot dataFile (data_file_name);
    bcInterfacePtr_Type                     BC ( new bcInterface_Type() );
    BC->createHandler();
    BC->fillHandler ( data_file_name, "problem" );
    BC->handler()->bcUpdate( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

    //**********************************************************//
    //  Preconditioner information read from the getPot object  //
    //  related to the file FiberGenerationPreconditioner       //
    //**********************************************************//

    typedef LifeV::Preconditioner             basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack       prec_Type;
    typedef boost::shared_ptr<prec_Type>      precPtr_Type;

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    //*********************************************************//
    //  Linear solver information read from the parameterList  //
    //  related to the xml file FiberGenerationParamList.xml   //
    //*********************************************************//

    LinearSolver linearSolver;
    linearSolver.setCommunicator (Comm);
    linearSolver.setParameters ( parameterList );
    linearSolver.setPreconditioner ( precPtr );

    //*************************************************************//
    //  Create right hand side and correction taking into account  //
    //  the boundary condition by using the bcManage method        //
    //*************************************************************//

    vectorPtr_Type rhs (new vector_Type ( uSpace -> map() ) );
    *rhs *= 0.0;
    rhs -> globalAssemble();

    bcManage ( *systemMatrix, *rhs, *uSpace->mesh(), uSpace->dof(), *BC -> handler(), uFESpace->feBd(), 1.0, 0.0 );

    linearSolver.setOperator (systemMatrix);

    //******************************************//
    //  Computing solution of the linear system //
    //  related to the laplacian operator  	//
    //******************************************//

    vectorPtr_Type solution ( new vector_Type ( uFESpace -> map() ) );

    vectorPtr_Type sx (new vector_Type ( uSpace -> map() ) );
    vectorPtr_Type sy (new vector_Type ( uSpace -> map() ) );
    vectorPtr_Type sz (new vector_Type ( uSpace -> map() ) );

    linearSolver.setRightHandSide (rhs);
    linearSolver.solve (solution);

    //******************************************//
    //  Computing the sheetlet direction as     //
    //  the gradient of the laplacian solution, //
    //  by using the method GradientRecovery    //
    //******************************************//

    *sx = GradientRecovery::ZZGradient (uSpace, *solution, 0);
    *sy = GradientRecovery::ZZGradient (uSpace, *solution, 1);
    *sz = GradientRecovery::ZZGradient (uSpace, *solution, 2);

    vectorPtr_Type rbSheet ( new vector_Type ( sFESpace -> map() ) );
    vectorPtr_Type rbFiber ( new vector_Type ( sFESpace -> map() ) );
    vectorPtr_Type cl ( new vector_Type ( sFESpace -> map() ) );
    vectorPtr_Type projection ( new vector_Type ( sFESpace -> map() ) );

    int n = (*rbSheet).epetraVector().MyLength();
    int d = n / 3;

    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

        (*rbSheet) [i] = (*sx) [i];
        (*rbSheet) [j] = (*sy) [i];
        (*rbSheet) [k] = (*sz) [i];
    }


    HeartUtility::normalize (*rbSheet);

    //**************************************************//
    //  Computing the centerline projection     	//
    //  on the plane orthogonal to the sheet direction  //
    //**************************************************//

    Real cx = parameterList.get ("centerline_x", 0.0);
    Real cy = parameterList.get ("centerline_y", 0.0);
    Real cz = parameterList.get ("centerline_z", 1.0);

    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

        (*rbFiber) [i] = (*rbSheet) [j];
        (*rbFiber) [j] = - (*rbSheet) [i];
        (*rbFiber) [k] = 0.0;

        Real cdot = cx * (*rbSheet) [i] + cy * (*rbSheet) [j] + cz * (*rbSheet) [k];

        (*projection) [i] = cx - cdot *  (*rbSheet) [i];
        (*projection) [j] = cy - cdot *  (*rbSheet) [j];
        (*projection) [k] = cz - cdot *  (*rbSheet) [k];

        (*cl) [i] = cx;
        (*cl) [j] = cy;
        (*cl) [k] = cz;


        (*rbFiber) [i] = (*rbSheet) [j];
        (*rbFiber) [j] = - (*rbSheet) [i];
        (*rbFiber) [k] = 0.0;


        Real scalarp;
        scalarp = (*projection) [i] * (*rbFiber) [i]
                   + (*projection) [j] * (*rbFiber) [j]
                   + (*projection) [k] * (*rbFiber) [k];

        if ( scalarp == 1.0 )
        {
            (*rbFiber) [i] = 0;
            (*rbFiber) [j] = - (*rbSheet) [k];
            (*rbFiber) [k] = (*rbSheet) [j];
        }

    }

    HeartUtility::normalize (*rbFiber);
    HeartUtility::normalize (*projection);

    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

	    //******************************************//
	    //  Definition of the fiber field as the    //
	    //	cross product between the sheet and the //
	    //	projection field.                       //
	    //******************************************//

        (*rbFiber) [i] = (*rbSheet) [j]  * (*projection) [k] - (*rbSheet) [k]  * (*projection) [j];
        (*rbFiber) [j] = (*rbSheet) [k]  * (*projection) [i] - (*rbSheet) [i]  * (*projection) [k];
        (*rbFiber) [k] = (*rbSheet) [i]  * (*projection) [j] - (*rbSheet) [j]  * (*projection) [i];
    }

    HeartUtility::normalize (*rbFiber);


    //**************************************************//
    //  Fixed angle on the endocardium and epicardium. 	//
    //  Angle varies uniformly from the endocardium to  //
    //  the epicardium					//
    //**************************************************//

    Real epi_angle = parameterList.get ("epi_angle", -60.0);
    Real endo_angle = parameterList.get ("endo_angle", 60.0);

    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

        Real scalarp;
        scalarp = (*projection) [i] * (*rbFiber) [i]
                   + (*projection) [j] * (*rbFiber) [j]
                   + (*projection) [k] * (*rbFiber) [k];

	    //*************************************************//
	    //  Creation of the rotation matrix applied to the //
	    //  fiber field defined as the cross product       //
	    //  between the sheet and projection field         //
	    //*************************************************//


        Real p = 3.14159265358979;
        Real teta1 = p * epi_angle / 180;
        Real teta2 = p * endo_angle / 180;
        Real m = (teta1 - teta2 );
        Real q = teta2;
        Real teta;

        teta = m * (*solution) [i] + q;
        (*solution) [i] = teta;

        Real s01 = (*rbSheet) [i];
        Real s02 = (*rbSheet) [j];
        Real s03 = (*rbSheet) [k];
        Real f01 = (*rbFiber) [i];
        Real f02 = (*rbFiber) [j];
        Real f03 = (*rbFiber) [k];

        Real R11 = std::cos (teta) + s01 * s01 * ( 1 - std::cos (teta) );
        Real R12 = s01 * s02 *  ( 1 - std::cos (teta) ) - s03 * std::sin (teta);
        Real R13 = s01 * s03 *  ( 1 - std::cos (teta) ) + s02 * std::sin (teta);

        Real R21 = s02 * s01 *  ( 1 - std::cos (teta) ) + s03 * std::sin (teta);
        Real R22 = std::cos (teta) + s02 * s02 * ( 1 - std::cos (teta) );
        Real R23 = s02 * s03 *  ( 1 - std::cos (teta) ) - s01 * std::sin (teta);

        Real R31 = s03 * s01 *  ( 1 - std::cos (teta) ) - s02 * std::sin (teta);
        Real R32 = s03 * s02 *  ( 1 - std::cos (teta) ) + s01 * std::sin (teta);
        Real R33 = std::cos (teta) + s03 * s03 * ( 1 - std::cos (teta) );

        Real ca = std::cos (teta);
        Real sa = std::sin (teta);

        Real W11 = 0.0;
        Real W12 = -s03;
        Real W13 = s02;
        Real W21 = s03;
        Real W22 = 0.0;
        Real W23 = -s01;
        Real W31 = -s02;
        Real W32 = s01;
        Real W33 = 0.0;
        Real sa2 = 2.0 * std::sin (0.5 * teta) * std::sin (0.5 * teta);
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
    //**************************************//
    //  Export the fibers and sheets field  //
    //**************************************//

    ExporterHDF5< mesh_Type > exporterFibers;
    exporterFibers.setMeshProcId ( meshPart, Comm -> MyPID() );
    exporterFibers.setPostDir (problemFolder);
    std::string outputFiberFileName = parameterList.get ("output_fiber_filename", "FiberDirection");
    exporterFibers.setPrefix (outputFiberFileName);
    exporterFibers.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", sFESpace, rbFiber, UInt (0) );
    exporterFibers.postProcess (0);
    exporterFibers.closeFile();
    ExporterHDF5< mesh_Type > exporterSheets;
    exporterSheets.setMeshProcId ( meshPart, Comm -> MyPID() );
    exporterSheets.setPostDir (problemFolder);
    std::string outputSheetsFileName = parameterList.get ("output_sheets_filename", "SheetsDirection");
    exporterSheets.setPrefix (outputSheetsFileName);
    exporterSheets.addVariable ( ExporterData<mesh_Type>::VectorField,  "sheets", sFESpace, rbSheet, UInt (0) );
    exporterSheets.postProcess (0);
    exporterSheets.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


}


