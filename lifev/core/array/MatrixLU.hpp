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
    @brief Implemetation of a small square matrix which performs LU decomposition

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 01-04-2013
 */

#ifndef MATRIXLU_HPP_
#define MATRIXLU_HPP_

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorLU.hpp>

namespace LifeV
{

class MatrixLU
{
	typedef std::vector<std::vector <Real> > Matrix;
	typedef std::vector<Real> Vector;

public:
	MatrixLU(UInt n);
	MatrixLU(UInt n, UInt m, Real r = 0.0);
	MatrixLU(Matrix A);
	virtual ~MatrixLU(){};

	void setIdentity(UInt n = 0);
	void setValues(Real r);
	void setZero();
	void setCol(UInt j, const VectorLU& v);
	void setLine(UInt i, const VectorLU& v);

	VectorLU getLine(UInt i) const;
	VectorLU getCol(UInt j) const;
	UInt size(UInt n = 1) const;

	MatrixLU triU() const;
	MatrixLU triL() const;
	VectorLU solveL(const VectorLU& b) const;
	VectorLU solveU(const VectorLU& b) const;

	void LU(MatrixLU& P, MatrixLU& Q, MatrixLU& L, MatrixLU& U) const;
	void Pivot(MatrixLU& Pk, MatrixLU& Qk, UInt k) const;
	void MGauss(MatrixLU& Mk, UInt k) const;

	MatrixLU  operator+ (const MatrixLU& B) const;
	MatrixLU& operator+=(const MatrixLU& B);
	MatrixLU  operator- (const MatrixLU& B) const;
	MatrixLU& operator-=(const MatrixLU& B);
	MatrixLU  operator* (const MatrixLU& B) const;
	MatrixLU& operator*=(const MatrixLU& B);
	MatrixLU  operator* (const Real r) const;
	MatrixLU& operator*=(const Real r);
	MatrixLU  operator/ (const Real r) const;
	MatrixLU& operator/=(const Real r);
	Vector 	  operator* (const Vector& x) const;
	VectorLU  operator* (const VectorLU& x) const;
	Vector&   operator[](UInt i);
	Real&	  operator()(UInt i, UInt j);
	const Vector&   operator[](UInt i) const;
	const Real&	  operator()(UInt i, UInt j) const;

	template<UInt s>
	VectorLU  timesVectorSmall(const VectorSmall<s>& v) const;

	void disp();

protected:
	UInt M_n;
	UInt M_m;
	Matrix M_A;
};


MatrixLU::MatrixLU(UInt n)
: M_n(n), M_m(n)
{
	setIdentity(n);
}

MatrixLU::MatrixLU(UInt n, UInt m, Real r)
: M_n(n), M_m(m)
{
	setValues(r);
}

MatrixLU::MatrixLU(std::vector<std::vector<Real> > A)
: M_A(A), M_n(A.size()), M_m(A[0].size())
{
}

void MatrixLU::setIdentity(UInt n)
{
	if(n==0)
	{
		n = M_n;
		M_m = M_n;
	}

	M_A = std::vector< std::vector<Real> >(n,std::vector<Real>(n,0.0));

	for(Int i = 0; i<n; i++)
		M_A[i][i] = 1.0;
}

void MatrixLU::setValues(Real r)
{
	M_A = std::vector< std::vector<Real> >(M_n, std::vector<Real>(M_m,r));
}

void MatrixLU::setZero()
{
	setValues(0.0);
}

void MatrixLU::setCol(UInt j, const VectorLU& v)
{
	if(v.size() != M_n)
	{
		std::cout<<"ERROR : the size of the column does not agree\n";
		return;
	}

	for(UInt i=0; i<M_n; i++)
		M_A[i][j] = v[i];
}

void MatrixLU::setLine(UInt i, const VectorLU& v)
{
	if(v.size() != M_m)
	{
		std::cout<<"ERROR : the size of the line does not agree\n";
		return;
	}

	for(UInt j=0; j<M_m; j++)
		M_A[i][j] = v[j];
}

VectorLU MatrixLU::getLine(UInt i) const
{
	VectorLU v(M_m);

	for(UInt j = 0; j<M_m; j++)
		v[j] = M_A[i][j];

	return v;
}

VectorLU MatrixLU::getCol(UInt j) const
{
	VectorLU v(M_n);

	for(UInt i = 0; i<M_n; i++)
		v[i] = M_A[i][j];

	return v;
}

UInt MatrixLU::size(UInt n) const
{
	if(n==1)
		return M_n;
	else
		return M_m;
}

MatrixLU MatrixLU::triU() const
{
	MatrixLU C(M_n, M_m);

	for(UInt i=0; i<M_n; i++)
		for(UInt j=i; j<M_m; j++)
			C(i,j) = M_A[i][j];

	return C;
}

MatrixLU MatrixLU::triL() const
{
	MatrixLU C(M_n,M_m);

	for(UInt i=0; i<M_n; i++)
		for(UInt j=0; j<=i; j++)
			C(i,j) = M_A[i][j];

	return C;
}

VectorLU MatrixLU::solveL(const VectorLU& b) const
{
	VectorLU x(M_n);
	Real tmp;

	x[0] = b[0]/M_A[0][0];

	for(UInt i=1; i<M_n; i++)
	{
		tmp = 0.;
		for(UInt k = 0; k<=i-1; k++)
			tmp += M_A[i][k]*x[k];
		x[i] = (b[i]-tmp)/M_A[i][i];
	}

	return x;
}

VectorLU MatrixLU::solveU(const VectorLU& b) const
{
	VectorLU x(M_n);
	Real tmp;

	x[M_n-1] = b[M_n-1]/M_A[M_n-1][M_n-1];

	for(UInt i=2; i<=M_n; i++)
	{
		tmp = 0.0;
		for(UInt k = M_n-i+1; k<=M_n-1; k++)
			tmp += M_A[M_n-i][k]*x[k];
		x[M_n-i] = (b[M_n-i]-tmp)/M_A[M_n-i][M_n-i];
	}

	return x;
}

void MatrixLU::LU(MatrixLU& P, MatrixLU& Q, MatrixLU& L, MatrixLU& U) const
{
	MatrixLU A(*this);
	MatrixLU I(M_n);
	MatrixLU Minv(M_n), Pk(M_n), Qk(M_n), Mk(M_n), Mkinv(M_n);
	P = I;
	Q = I;

	for(UInt k = 1; k<M_n; k++)
	{
		Pk = I;
		Qk = I;
		A.Pivot(Pk, Qk, k);
		A = Pk*A;
		A = A*Qk;
		Mk = I;
		A.MGauss(Mk, k);
		Mkinv = I*2.0 - Mk;
		A = Mk*A;
		P = Pk*P;
		Q = Q*Qk;
		Minv = Minv*Pk;
		Minv = Minv*Mkinv;
	}

	U = A.triU();
	L = P*Minv;
}

void MatrixLU::Pivot(MatrixLU& Pk, MatrixLU& Qk, UInt k) const
{
	UInt i = k-1;
	UInt j = k-1;
	Real max = abs(M_A[i][j]);

	for(UInt l = k-1; l<M_n; l++)
	{
		for(UInt c = k-1; c<M_n; c++)
		{
			if( abs(M_A[l][c]) > max)
			{
				max = std::abs(M_A[l][c]);
				i = l;
				j = c;
			}
		}
	}

	Pk[i][i] = 0.;	Pk[k-1][k-1] = 0.;	Pk[k-1][i] = 1.;	Pk[i][k-1] = 1.;
	Qk[j][j] = 0.;	Qk[k-1][k-1] = 0.;	Qk[k-1][j] = 1.;	Qk[j][k-1] = 1.;
}

void MatrixLU::MGauss(MatrixLU& Mk, UInt k) const
{
	for(UInt i = k; i<M_n; i++)
		Mk[i][k-1] = -M_A[i][k-1]/M_A[k-1][k-1];
}

MatrixLU MatrixLU::operator+(const MatrixLU& B) const
{
	return MatrixLU(*this) += B;
}

MatrixLU& MatrixLU::operator+=(const MatrixLU& B)
{
	for(UInt i = 0; i<M_n; i++)
		for(UInt j = 0; j<M_m; j++)
			M_A[i][j] += B[i][j];

	return *this;
}

MatrixLU MatrixLU::operator-(const MatrixLU& B) const
{
	return MatrixLU(*this) -= B;
}

MatrixLU& MatrixLU::operator-=(const MatrixLU& B)
{
	for(UInt i = 0; i<M_n; i++)
		for(UInt j = 0; j<M_m; j++)
			M_A[i][j] -= B[i][j];

	return *this;
}

MatrixLU  MatrixLU::operator*(const MatrixLU& B) const
{
	return MatrixLU(*this) *= B;
}

MatrixLU& MatrixLU::operator*=(const MatrixLU& B)
{
	MatrixLU C(M_n, B.size(2));

	for(UInt i = 0; i<M_n; i++)
		for(UInt j = 0; j<B.size(2); j++)
			for(UInt k = 0; k<M_m; k++)
				C[i][j] += M_A[i][k]*B[k][j];
	*this = C;

	return *this;
}

MatrixLU MatrixLU::operator*(const Real r) const
{
	return MatrixLU(*this) *= r;
}

MatrixLU& MatrixLU::operator*=(const Real r)
{
	for(UInt i = 0; i<M_n; i++)
		for(UInt j = 0; j<M_m; j++)
			M_A[i][j] *= r;

	return *this;
}

MatrixLU MatrixLU::operator/(const Real r) const
{
	return MatrixLU(*this)*=(1.0/r);
}

MatrixLU& MatrixLU::operator/=(const Real r)
{
	*this *= (1.0/r);

	return *this;
}

std::vector<Real> MatrixLU::operator*(const std::vector<Real>& x) const
{
	std::vector<Real> b(M_n, 0.0);

	for(Int i = 0; i<M_n; i++)
		for(Int j = 0; j<M_m; j++)
			b[i] += M_A[i][j]*x[j];

	return b;
}

VectorLU MatrixLU::operator*(const VectorLU& x) const
{
	VectorLU b(M_n);

	for(Int i = 0; i<M_n; i++)
		for(Int j = 0; j<M_m; j++)
			b[i] += M_A[i][j]*x[j];

	return b;
}

std::vector<Real>& MatrixLU::operator[](UInt i)
{
	return M_A[i];
}

Real& MatrixLU::operator()(UInt i, UInt j)
{
	return M_A[i][j];
}

const std::vector<Real>& MatrixLU::operator[](UInt i) const
{
	return M_A[i];
}

const Real& MatrixLU::operator()(UInt i, UInt j) const
{
	return M_A[i][j];
}

template<UInt s>
VectorLU MatrixLU::timesVectorSmall(const VectorSmall<s>& v) const
{
	VectorLU b(M_n);

	for(UInt i = 0; i<M_n; i++)
		for(UInt j = 0; j<M_m; j++)
			b[i] += M_A[i][j]*v[j];

	return b;
}

void MatrixLU::disp()
{
	for(UInt i = 0; i<M_n; i++)
	{
		for(UInt j = 0; j<M_m; j++)
			cout<<"   "<<M_A[i][j];
		cout<<"\n";
	}
}


}// namespace LifeV


#endif /* MATRIXLU_HPP_ */
