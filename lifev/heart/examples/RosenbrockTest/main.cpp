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
    @brief 0D test with the Negroni Lascano model of 1996.

    @date 01−2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
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



#include <fstream>
#include <string>

#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/heart/solver/HeartETAMonodomainSolver.hpp>
#include <lifev/heart/solver/HeartIonicSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/heart/solver/IonicModels/IonicFitzHughNagumo.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/fem/RosenbrockTransformed.hpp>



using namespace std;
using namespace LifeV;

typedef MatrixEpetra<double>                    matrix_Type;
typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;
typedef VectorEpetra                            vector_Type;
typedef boost::shared_ptr<vector_Type >         vectorPtr_Type;
typedef MapEpetra						 	    map_Type;

void Playing(boost::shared_ptr<Epetra_Comm>  Comm);
matrix_Type getIdentity(const map_Type& map, UInt n);

Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }


    cout<<"Creating RosenbrockTransformed..."<<endl;
    RosenbrockTransformed<4> ros;

    cout<<"Done!"<<endl;

    Playing(Comm);

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}

void Playing(boost::shared_ptr<Epetra_Comm>  Comm)
{
    cout<<"Creating map ..."<<endl;
    MapEpetra mappa(4, Comm);

    cout<<"\n\n\n		Playing with vectors\n\n\n"<<endl;

    cout<<"Creating vector ..."<<endl;
    VectorEpetra vett1( mappa );			// (mappa, Unique) ???
    VectorEpetra vett2( mappa );

    cout<<"Creating matrix ..."<<endl;
    //MatrixEpetra<Real> matrice( mappa, vett1.size() );
    MatrixEpetra<Real> matrice( vett1.map(), vett1.size() );

    Int d = vett1.epetraVector().MyLength();
    const Int* j = vett1.blockMap().MyGlobalElements();
    const Int* i = vett2.blockMap().MyGlobalElements();


    for ( int k (0); k < d; k++)
    {
    	cout<<"Componente "<<k<<" di vett1 ha indice "<<j[k]<<endl;
    	cout<<"Componente "<<k<<" di vett2 ha indice "<<i[k]<<endl;

    	vett1( j[k] ) = k*k;
    	vett2( i[k] ) = k;

    	cout<<"vett1("<<j[k]<<") = "<<vett1(j[k])<<endl;
    	cout<<"vett2("<<i[k]<<") = "<<vett2(i[k])<<endl;

    }

    cout<<"\n\n\n		Playing with matrix\n\n\n"<<endl;

    matrixPtr_Type M_jac;

    int* GlobalID = new int[4];
    GlobalID[0] = 0;
    GlobalID[1] = 1;
    GlobalID[2] = 2;
    GlobalID[3] = 3;

    int* ElementsPerRow = new int[4];
    ElementsPerRow[0] = 4;
    ElementsPerRow[1] = 4;
    ElementsPerRow[2] = 4;
    ElementsPerRow[3] = 4;

    int* Indices   = new int[4];

    double* Values0 = new double[4];
    double* Values1 = new double[4];
    double* Values2 = new double[4];
    double* Values3 = new double[4];

    Indices[0] = 0;
    Indices[1] = 1;
    Indices[2] = 2;
    Indices[3] = 3;

    Values0[0] = 1;
    Values0[1] = 2;
    Values0[2] = 3;
    Values0[3] = 4;

    Values1[0] = 5;
    Values1[1] = 6;
    Values1[2] = 7;
    Values1[3] = 8;

    Values2[0] = 9;
    Values2[1] = 10;
    Values2[2] = 11;
    Values2[3] = 12;

    Values3[0] = 13;
    Values3[1] = 14;
    Values3[2] = 15;
    Values3[3] = 16;

    M_jac.reset (new matrix_Type (vett1.map(), ElementsPerRow) );

    M_jac->matrixPtr()->InsertGlobalValues (GlobalID[0], 4, Values0, Indices);
    M_jac->matrixPtr()->InsertGlobalValues (GlobalID[1], 4, Values1, Indices);
    M_jac->matrixPtr()->InsertGlobalValues (GlobalID[2], 4, Values2, Indices);
    M_jac->matrixPtr()->InsertGlobalValues (GlobalID[3], 4, Values3, Indices);

    M_jac->globalAssemble();

    delete Indices;
    delete Values0;
    delete Values1;
    delete Values2;
    delete Values3;
    delete ElementsPerRow;
    delete GlobalID;

    cout<<"Norm inf = "<<M_jac->normInf()<<endl;
    cout<<"Norm 1   = "<<M_jac->norm1()<<endl;

    vett2 = (*M_jac)*vett1;
    cout<<"vett1 = "<<endl;
    vett1.showMe();
    cout<<"M_jac * vett1 = \n"<<vett2.epetraVector()<<endl;


    cout<<"Setting M_jac = I...";
    //M_jac.reset(new matrix_Type(getIdentity(vett1.map(), d)));
    *M_jac = getIdentity(vett1.map(), d);
    cout<<"Done\n";
    cout<<"Norm inf = "<<matrice.normInf()<<endl;
    cout<<"Norm 1   = "<<matrice.norm1()<<endl;

    vett2 = (*M_jac)*vett1;
    cout<<"vett1 = "<<endl;
    vett1.showMe();
    cout<<"M_jac * vett1 = \n";
    vett2.showMe();

    cout<<"substracting M_jac-=3*M_jac...";
    *M_jac -= (*M_jac)*3.0;
    cout<<"Done\nMultiplying by vett1...";
    vett2 = (*M_jac)*vett1;
    cout<<"Done\nresult = "<<endl;
    vett2.showMe();


    cout<<"\n\n\n		Playing with linear solver\n\n\n"<<endl;

    Teuchos::RCP< Teuchos::ParameterList > List = Teuchos::rcp ( new Teuchos::ParameterList );
    List = Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" );

    LinearSolver solver;
    solver.setCommunicator(Comm);
    solver.setParameters(*List);

    vectorPtr_Type sol;
    sol.reset( new vector_Type( vett1.map() ) );
    *sol *= 0;

    vectorPtr_Type rhs;
    rhs.reset( new vector_Type( vett1.map() ) );
    *rhs = vett2;

    solver.setOperator(M_jac);
    solver.setRightHandSide( rhs );
    solver.solve(sol);

    cout<<"Solution of M_jac^{-1}*(M_jac*vett1) = "<<endl<<sol->epetraVector()<<endl;

    cout<<"vett1.size = "<<vett1.size()<<endl;

}

matrix_Type getIdentity(const map_Type& map, UInt n)
{

	matrix_Type I(map, n, false);
	int* Indices = new int[n];
	double* Values = new double[n];

	for( int i=0; i<n; i++)
		Indices[i] = i;

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
			Values[j] = (double)(i==j);

		I.matrixPtr()->InsertGlobalValues (i, n, Values, Indices);
	}

	I.globalAssemble();

	delete Indices;
	delete Values;

	return I;
}












