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
//
///*!
//  @file
//  @brief Mother class for ionic models
//
//  @date 01-2013
//  @author Simone Rossi <simone.rossi@epfl.ch>
//
//  @contributors
//  @mantainer Simone Rossi <simone.rossi@epfl.ch>
//  @last update 02-2012
// */
//


	#include <lifev/heart/solver/IonicModels/HeartIonicModel.hpp>

namespace LifeV
{
// ===================================================
//! Constructors
// ===================================================
HeartIonicModel::HeartIonicModel() :
    M_numberOfEquations (0)
{
}

HeartIonicModel::HeartIonicModel ( int n ) :
    M_numberOfEquations (n)
{
}



HeartIonicModel::HeartIonicModel ( const HeartIonicModel& Ionic ) :
    M_numberOfEquations ( Ionic.Size() )
{
}

// ===================================================
//! Methods
// ===================================================
HeartIonicModel& HeartIonicModel::operator = ( const HeartIonicModel& Ionic )
{
    M_numberOfEquations = Ionic.M_numberOfEquations;
    return      *this;
}

vector< vector<Real> > HeartIonicModel::getJac (const vector<Real>& v, Real h)
{
	vector< vector<Real> > J( M_numberOfEquations, vector<Real>(M_numberOfEquations,0.0) );
	vector<Real> f1(M_numberOfEquations,0.0);
	vector<Real> f2(M_numberOfEquations,0.0);
	vector<Real> y1(M_numberOfEquations,0.0);
	vector<Real> y2(M_numberOfEquations,0.0);

	for(int i=0; i<M_numberOfEquations; i++)
	{
		for(int j=0; j<M_numberOfEquations; j++)
		{
			y1[j] = v[j] + ((double)(i==j))*h;
			y2[j] = v[j] - ((double)(i==j))*h;
		}
		this->computeRhs(y1, 0.0, f1);
		this->computeRhs(y2, 0.0, f2);

		for(int j=0; j<M_numberOfEquations; j++)
			J[j][i] = (f1[j]-f2[j])/(2.0*h);
	}

	return J;
}

MatrixEpetra<Real> HeartIonicModel::getJac (const vector_Type& v, Real h)
{
	matrix_Type J(v.map(), M_numberOfEquations, false);
	vector<Real> f1(M_numberOfEquations,0.0);
	vector<Real> f2(M_numberOfEquations,0.0);
	vector< vector<Real> > df ( M_numberOfEquations, vector<Real> (M_numberOfEquations,0.0) );
	vector<Real> y1(M_numberOfEquations,0.0);
	vector<Real> y2(M_numberOfEquations,0.0);
	const Int* k = v.blockMap().MyGlobalElements();

	int* Indices = new int[M_numberOfEquations];
	double* Values =  new double[M_numberOfEquations];
	for(int i=0; i<M_numberOfEquations; i++)
		Indices[i] = i;

	for(int i=0; i<M_numberOfEquations; i++)
	{
		for(int j=0; j<M_numberOfEquations; j++)
		{
			y1[j] = v[ k[j] ] + ((double)(i==j))*h;
			y2[j] = v[ k[j] ] - ((double)(i==j))*h;
		}
		this->computeRhs(y1, 0.0, f1);
		this->computeRhs(y2, 0.0, f2);

		for(int j=0; j<M_numberOfEquations; j++)
			df[j][i] = (f1[j]-f2[j])/(2.0*h);
	}


	for(int i=0; i<M_numberOfEquations; i++)
	{
		for(int j=0; j<M_numberOfEquations; j++)
			Values[j] = df[i][j];
		J.matrixPtr()->InsertGlobalValues (i, M_numberOfEquations, Values, Indices);
	}

	J.globalAssemble();

	delete[] Indices;
	delete[] Values;

	return J;
}

void HeartIonicModel::computeRhs ( const vector_Type& v, const Real& Iapp, vector_Type& rhs)
{
	std::vector<Real> vec(M_numberOfEquations, 0.0);
	std::vector<Real> f(M_numberOfEquations, 0.0);
	const Int* j = v.blockMap().MyGlobalElements();
	const Int* i = rhs.blockMap().MyGlobalElements();

	for(int k=0; k<M_numberOfEquations; k++)
		vec[k] = v[j[k]];

	this->computeRhs(vec, Iapp, f);

	for(int k=0; k<M_numberOfEquations; k++)
			rhs[i[k]] = f[k];

	//delete j;
	//delete i;
}

void HeartIonicModel::computeRhs (   const std::vector<vectorPtr_Type>& v,
                                     std::vector<vectorPtr_Type>& rhs )
{

    int nodes = ( * (v.at (1) ) ).epetraVector().MyLength();


    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );
    std::vector<Real>   localRhs ( M_numberOfEquations - 1, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = ( * (v.at (1) ) ).blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }

        computeRhs ( localVec, localRhs );

        for ( int i = 1; i < M_numberOfEquations; i++ )
        {
            ( * ( rhs.at (i) ) ) [j] =  localRhs.at (i - 1);
        }

    }

}


void HeartIonicModel::computeRhs (   const std::vector<vectorPtr_Type>& v,
                                     const VectorEpetra& Iapp,
                                     std::vector<vectorPtr_Type>& rhs )
{

    int nodes = Iapp.epetraVector().MyLength();


    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );
    std::vector<Real>   localRhs ( M_numberOfEquations, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = Iapp.blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }

        computeRhs ( localVec, Iapp[j], localRhs );

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            ( * ( rhs.at (i) ) ) [j] =  localRhs.at (i);
        }

    }

}

void HeartIonicModel::computePotentialRhsICI (   const std::vector<vectorPtr_Type>& v,
                                                 const VectorEpetra& Iapp,
                                                 std::vector<vectorPtr_Type>& rhs,
                                                 matrix_Type&                    massMatrix  )
{
    int nodes = ( * (v.at (0) ) ).epetraVector().MyLength();


    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = ( * (v.at (0) ) ).blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }

        ( * ( rhs.at (0) ) ) [j] =  computeLocalPotentialRhs ( localVec, Iapp[j] );

    }

    ( * ( rhs.at (0) ) ) = massMatrix * ( * ( rhs.at (0) ) );

}


void HeartIonicModel::computePotentialRhsSVI (   const std::vector<vectorPtr_Type>& v,
                                                 const VectorEpetra& Iapp,
                                                 std::vector<vectorPtr_Type>& rhs,
                                                 FESpace<mesh_Type, MapEpetra>& uFESpace )
{

    std::vector<Real> U (M_numberOfEquations, 0.0);
    Real I (0.0);
    ( * ( rhs.at (0) ) ) *= 0.0;

    std::vector<vectorPtr_Type>      URepPtr;
    for ( int k = 0; k < M_numberOfEquations; k++ )
    {
        URepPtr.push_back ( * ( new vectorPtr_Type ( new VectorEpetra (  * ( v.at (k) )     , Repeated ) ) ) );
    }

    VectorEpetra    IappRep ( Iapp       , Repeated );


    std::vector<elvecPtr_Type>      elvecPtr;
    for ( int k = 0; k < M_numberOfEquations; k++ )
    {
        elvecPtr.push_back ( * ( new elvecPtr_Type ( new VectorElemental (  uFESpace.fe().nbFEDof(), 1  ) ) ) );
    }

    VectorElemental elvec_Iapp ( uFESpace.fe().nbFEDof(), 1 );
    VectorElemental elvec_Iion ( uFESpace.fe().nbFEDof(), 1 );

    for (UInt iVol = 0; iVol < uFESpace.mesh()->numVolumes(); ++iVol)
    {

        uFESpace.fe().updateJacQuadPt ( uFESpace.mesh()->volumeList ( iVol ) );


        for ( int k = 0; k < M_numberOfEquations; k++ )
        {
            ( * ( elvecPtr.at (k) ) ).zero();
        }
        elvec_Iapp.zero();
        elvec_Iion.zero();

        UInt eleIDu = uFESpace.fe().currentLocalId();
        UInt nbNode = ( UInt ) uFESpace.fe().nbFEDof();

        //! Filling local elvec_u with potential values in the nodes
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {

            Int  ig = uFESpace.dof().localToGlobalMap ( eleIDu, iNode );

            for ( int k = 0; k < M_numberOfEquations; k++ )
            {
                ( * ( elvecPtr.at (k) ) ).vec() [iNode] = ( * ( URepPtr.at (k) ) ) [ig];
            }

            elvec_Iapp.vec() [ iNode ] = IappRep[ig];

        }

        //compute the local vector
        for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
        {

            for ( int k = 0; k < M_numberOfEquations; k++ )
            {
                U.at (k) = 0;
            }
            I = 0;

            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {

                for ( int k = 0; k < M_numberOfEquations; k++ )
                {
                    U.at (k) +=  ( * ( elvecPtr.at (k) ) ) (i) *  uFESpace.fe().phi ( i, ig );
                }

                I += elvec_Iapp (i) * uFESpace.fe().phi ( i, ig );

            }

            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {

                elvec_Iion ( i ) += computeLocalPotentialRhs (U, I) * uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );

            }

        }

        //assembly
        for ( UInt i = 0 ; i < uFESpace.fe().nbFEDof(); i++ )
        {
            Int  ig = uFESpace.dof().localToGlobalMap ( eleIDu, i );
            ( * ( rhs.at (0) ) ).sumIntoGlobalValues (ig,  elvec_Iion.vec() [i] );
        }
    }
}
}


