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

#ifndef MatrixStandard_HPP_
#define MatrixStandard_HPP_

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

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorStandard.hpp>
#include <lifev/core/array/VectorSmall.hpp>


namespace LifeV
{

class MatrixStandard : public std::vector< std::vector <Real> >
{
    typedef std::vector<std::vector <Real> > Matrix;
    typedef std::vector<Real> Vector;

public:
    MatrixStandard() {};
    MatrixStandard (UInt n);
    MatrixStandard (UInt n, UInt m, Real r = 0.0);
    MatrixStandard (Matrix A);
    virtual ~MatrixStandard() {};

    void setIdentity (UInt n = 0);
    void setValues (Real r);
    void setCol (UInt j, const VectorStandard& v);
    void setLine (UInt i, const VectorStandard& v);

    VectorStandard getLine (UInt i) const;
    VectorStandard getCol (UInt j) const;
    void getLine (UInt i, VectorStandard& v) const;
    void getCol (UInt j, VectorStandard& v) const;

    MatrixStandard triU() const;
    MatrixStandard triL() const;
    void triU (MatrixStandard& U) const;
    void triL (MatrixStandard& L) const;

    MatrixStandard  operator+ (const MatrixStandard& B) const;
    MatrixStandard& operator+= (const MatrixStandard& B);
    MatrixStandard  operator- (const MatrixStandard& B) const;
    MatrixStandard& operator-= (const MatrixStandard& B);
    MatrixStandard& operator-= (const Matrix& B);
    MatrixStandard  operator* (const MatrixStandard& B) const;
    MatrixStandard  operator* (const Real r) const;
    MatrixStandard& operator*= (const Real r);
    MatrixStandard  operator/ (const Real r) const;
    MatrixStandard& operator/= (const Real r);
    Vector    operator* (const Vector& x) const;
    VectorStandard  operator* (const VectorStandard& x) const;
    Vector&   operator[] (UInt i);
    const Vector&   operator[] (UInt i) const;

    void times (const MatrixStandard& B, MatrixStandard& C) const;
    void times (const VectorStandard& x, VectorStandard& b) const;
    template<UInt s>
    void  timesVectorSmall (const VectorSmall<s>& v, VectorStandard& b) const;

    void disp();
};

template<UInt s>
void MatrixStandard::timesVectorSmall (const VectorSmall<s>& v, VectorStandard& b) const
{
    for (UInt i = 0; i < this->size(); i++)
    {
        b[i] = (this->at (i).at (0) ) * v[0];
        for (UInt j = 1; j < this->at (0).size(); j++)
        {
            b[i] += (this->at (i).at (j) ) * v[j];
        }
    }
}

}// namespace LifeV


#endif /* MatrixStandard_HPP_ */
