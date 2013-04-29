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

#include <lifev/core/array/MatrixStandard.hpp>

namespace LifeV
{

MatrixStandard::MatrixStandard(UInt n)
: std::vector<std::vector <Real> >(n, std::vector<Real>(n,0.0))
{
	setIdentity(n);
}

MatrixStandard::MatrixStandard(UInt n, UInt m, Real r)
: std::vector<std::vector <Real> >(n, std::vector<Real>(m, r))
{
}

MatrixStandard::MatrixStandard(std::vector<std::vector<Real> > A)
: std::vector<std::vector <Real> >(A)
{
}

void MatrixStandard::setIdentity(UInt n)
{
	*this *= 0.0;

	for(Int i = 0; i<n; i++)
		this->at(i).at(i) = 1.0;
}

void MatrixStandard::setValues(Real r)
{
	for (UInt i=0; i<this->size(); i++)
		for(UInt j=0; j<this->at(0).size(); j++)
			this->at(i).at(j) = r;
}

void MatrixStandard::setCol(UInt j, const VectorStandard& v)
{
	for(UInt i=0; i<this->size(); i++)
		this->at(i).at(j) = v[i];
}

void MatrixStandard::setLine(UInt i, const VectorStandard& v)
{
	for(UInt j=0; j<this->at(0).size(); j++)
		this->at(i).at(j) = v[j];
}

VectorStandard MatrixStandard::getLine(UInt i) const
{
	VectorStandard v(this->at(0).size(), 0.0);
	getLine(i, v);

	return v;
}

VectorStandard MatrixStandard::getCol(UInt j) const
{
	VectorStandard v(this->size(), 0.0);
	getCol(j, v);

	return v;
}

void MatrixStandard::getLine(UInt i, VectorStandard& v) const
{
	for(UInt j = 0; j<this->at(0).size(); j++)
		v[j] = this->at(i).at(j);
}

void MatrixStandard::getCol(UInt j, VectorStandard& v) const
{
	for(UInt i = 0; i<this->size(); i++)
		v[i] = this->at(i).at(j);
}

MatrixStandard MatrixStandard::triU() const
{
	MatrixStandard C(this->size(), this->at(0).size());
	triU(C);

	return C;
}

MatrixStandard MatrixStandard::triL() const
{
	MatrixStandard C(this->size(), this->at(0).size());
	triL(C);

	return C;
}

void MatrixStandard::triU(MatrixStandard& U) const
{
	for(UInt i=0; i<this->size(); i++)
		for(UInt j=i; j<this->at(0).size(); j++)
			U[i][j] = this->at(i).at(j);
}

void MatrixStandard::triL(MatrixStandard& L) const
{
	for(UInt i=0; i<this->size(); i++)
		for(UInt j=0; j<=i; j++)
			L[i][j] = this->at(i).at(j);
}

MatrixStandard MatrixStandard::operator+(const MatrixStandard& B) const
{
	return MatrixStandard(*this) += B;
}

MatrixStandard& MatrixStandard::operator+=(const MatrixStandard& B)
{
	for(UInt i = 0; i<this->size(); i++)
		for(UInt j = 0; j<this->at(0).size(); j++)
			this->at(i).at(j) += B[i][j];

	return *this;
}

MatrixStandard MatrixStandard::operator-(const MatrixStandard& B) const
{
	return MatrixStandard(*this) -= B;
}

MatrixStandard& MatrixStandard::operator-=(const MatrixStandard& B)
{
	for(UInt i = 0; i<this->size(); i++)
		for(UInt j = 0; j<this->at(0).size(); j++)
			this->at(i).at(j) -= B[i][j];

	return *this;
}

MatrixStandard& MatrixStandard::operator-=(const std::vector<std::vector<Real> >& B)
{
	for(UInt i = 0; i<this->size(); i++)
		for(UInt j = 0; j<this->at(0).size(); j++)
			this->at(i).at(j) -= B[i][j];

	return *this;
}

MatrixStandard  MatrixStandard::operator*(const MatrixStandard& B) const
{
	MatrixStandard C(this->size(), B.at(0).size());

	for(UInt i = 0; i<this->size(); i++)
		for(UInt j = 0; j<B.at(0).size(); j++)
			for(UInt k = 0; k<B.size(); k++)
				C[i][j] += (this->at(i).at(k))*B[k][j];
	return C;
}

void MatrixStandard::times(const MatrixStandard& B, MatrixStandard& C) const
{
	for(UInt i = 0; i<this->size(); i++)
	{
		for(UInt j = 0; j<B.at(0).size(); j++)
		{
			C[i][j] = this->at(i).at(0)*B[0][j];
			for(UInt k = 1; k<B.size(); k++)
				C[i][j] += (this->at(i).at(k))*B[k][j];
		}
	}
}

MatrixStandard MatrixStandard::operator*(const Real r) const
{
	return MatrixStandard(*this) *= r;
}

MatrixStandard& MatrixStandard::operator*=(const Real r)
{
	for(UInt i = 0; i<this->size(); i++)
		for(UInt j = 0; j<this->at(0).size(); j++)
			this->at(i).at(j) *= r;

	return *this;
}

MatrixStandard MatrixStandard::operator/(const Real r) const
{
	return MatrixStandard(*this)*=(1.0/r);
}

MatrixStandard& MatrixStandard::operator/=(const Real r)
{
	*this *= (1.0/r);

	return *this;
}

std::vector<Real> MatrixStandard::operator*(const std::vector<Real>& x) const
{
	std::vector<Real> b(this->size(), 0.0);

	for(UInt i = 0; i<this->size(); i++)
		for(UInt j = 0; j<this->at(0).size(); j++)
			b[i] += (this->at(i).at(j))*x[j];

	return b;
}

VectorStandard MatrixStandard::operator*(const VectorStandard& x) const
{
	VectorStandard b(this->size(), 0.0);

	for(UInt i = 0; i<this->size(); i++)
		for(UInt j = 0; j<x.size(); j++)
			b[i] += (this->at(i).at(j))*x[j];

	return b;
}

void MatrixStandard::times(const VectorStandard& x, VectorStandard& b) const
{
	for(UInt i = 0; i<this->size(); i++)
	{
		b[i] = (this->at(i).at(0))*x[0];
		for(UInt j = 1; j<x.size(); j++)
			b[i] += (this->at(i).at(j))*x[j];
	}
}

std::vector<Real>& MatrixStandard::operator[](UInt i)
{
	return this->at(i);
}

const std::vector<Real>& MatrixStandard::operator[](UInt i) const
{
	return this->at(i);
}

void MatrixStandard::disp()
{
	for(UInt i = 0; i<this->size(); i++)
	{
		for(UInt j = 0; j<this->at(0).size(); j++)
			cout<<"   "<<this->at(i).at(j);
		cout<<"\n";
	}
}

} //namespace LifeV
