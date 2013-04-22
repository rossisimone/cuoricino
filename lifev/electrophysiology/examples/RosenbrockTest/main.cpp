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



#include <fstream>
#include <string>

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/VectorLU.hpp>
#include <lifev/core/array/MatrixLU.hpp>

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
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroIonicSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/electrophysiology/solver/IonicModels/IonicFitzHughNagumo.hpp>
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

/*    const int n=4;
    RosenbrockTransformed<n> ros;
    typedef vector<vector<Real> > matrix;
    typedef vector<Real> vector;

    matrix A(n,vector(n,0.0));
    matrix I(ros.getIdentity(n));
    matrix P(I),Q(I),L(I),U(I), B1(I), B2(I);
    A[0][0] = 1.;	A[0][1] = 2.;	A[0][2] = 3.;	A[0][3] = 4.;
    A[1][0] = 9.;	A[1][1] = 8.;	A[1][2] = 7.;	A[1][3] = 6.;
    A[2][0] = 20.;	A[2][1] = 12.;	A[2][2] = 13.;	A[2][3] = 45.;
    A[3][0] = 89.;	A[3][1] = 23.;	A[3][2] = 23.;	A[3][3] = 23.;

    matrix AA(A);

    ros.LU(A, P, Q, L, U, n);
    B1 = ros.mult(P,AA,n);
    B1 = ros.mult(B1,Q,n);
    B2 = ros.mult(L,U,n);

    ros.disp(P,"P", n);
    ros.disp(Q,"Q", n);
    ros.disp(L,"L", n);
    ros.disp(U,"U", n);
    ros.disp(B1,"B1", n);
    ros.disp(B2,"B2", n);

    vector v(n,0.0);
    v[0] = 1.;	v[1] = 1.;	v[2] = 1.;	v[3] = 1.;
    vector w(ros.mult(AA,v,n));
    cout<<"\n\nv = \n"; for(int i=0; i<n; i++) cout<<"    "<<v[i]<<endl;
    cout<<"\n\nw = \n"; for(int i=0; i<n; i++) cout<<"    "<<w[i]<<endl;
    v = ros.mult(P,w,n);
    cout<<"\n\nv =Pw = \n"; for(int i=0; i<n; i++) cout<<"    "<<v[i]<<endl;
    v = ros.solvel(L,v,n);
    cout<<"\n\nv = L-1v = \n"; for(int i=0; i<n; i++) cout<<"    "<<v[i]<<endl;
    v = ros.solveu(U,v,n);
    cout<<"\n\nv = U-1v\n"; for(int i=0; i<n; i++) cout<<"    "<<v[i]<<endl;
    v = ros.mult(Q,v,n);
    cout<<"\n\nv = Qv\n"; for(int i=0; i<n; i++) cout<<"    "<<v[i]<<endl;
*/
    //Playing(Comm);

    int n=4;
    MatrixLU A(n);
    cout<<"I = \n";
    A.disp();
    A[0][0] = 1.;	A[0][1] = 2.;	A[0][2] = 3.;	A[0][3] = 4.;
    A[1][0] = 9.;	A[1][1] = 8.;	A[1][2] = 7.;	A[1][3] = 6.;
    A[2][0] = 20.;	A[2][1] = 12.;	A[2][2] = 13.;	A[2][3] = 45.;
    A[3][0] = 89.;	A[3][1] = 23.;	A[3][2] = 23.;	A[3][3] = 23.;
    cout<<"A = \n";
    A.disp();

    MatrixLU B(n,n,1.0);
    MatrixLU C(n);
    C = A+B;
    cout<<"A+B = \n";
    C.disp();
    C-=B;
    cout<<"A+B-B = \n";
    C.disp();

    B*=0.5;
    cout<<"B*0.5 = \n";
    B.disp();
    B = B*2.0;
    cout<<"B = \n";
    B.disp();

    MatrixLU P(n),Q(n),L(n),U(n);
    A.LU(P,Q,L,U);
    VectorLU v(n,1.0);
    VectorLU w(A*v);

    cout<<"v = "; v.disp();
    cout<<"\nw = A*v = "; w.disp();

    v = P*w;
    v = L.solveL(v);
    v = U.solveU(v);
    v = Q*v;

    cout<<"\nv = A^(-1)*w = "; v.disp();

    A = A/2.0;

    cout<<"\nA/2 = /n";
    A.disp();
    A/= 0.5;
    cout<<"A = \n";
    A.disp();



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












