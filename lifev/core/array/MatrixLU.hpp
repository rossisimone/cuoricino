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
#include <lifev/core/array/VectorSmall.hpp>


namespace LifeV
{

class MatrixLU
{
    typedef std::vector<std::vector <Real> > Matrix;
    typedef std::vector<Real> Vector;

public:
    MatrixLU (UInt n);
    MatrixLU (UInt n, UInt m, Real r = 0.0);
    MatrixLU (Matrix A);
    virtual ~MatrixLU() {};

    void setIdentity (UInt n = 0);
    void setValues (Real r);
    void setZero();
    void setCol (UInt j, const VectorLU& v);
    void setLine (UInt i, const VectorLU& v);

    VectorLU getLine (UInt i) const;
    VectorLU getCol (UInt j) const;
    UInt size (UInt n = 1) const;

    MatrixLU triU() const;
    MatrixLU triL() const;
    VectorLU solveL (const VectorLU& b) const;
    VectorLU solveU (const VectorLU& b) const;

    void LU (MatrixLU& P, MatrixLU& Q, MatrixLU& L, MatrixLU& U) const;
    void Pivot (MatrixLU& Pk, MatrixLU& Qk, UInt k) const;
    void MGauss (MatrixLU& Mk, UInt k) const;

    MatrixLU  operator+ (const MatrixLU& B) const;
    MatrixLU& operator+= (const MatrixLU& B);
    MatrixLU  operator- (const MatrixLU& B) const;
    MatrixLU& operator-= (const MatrixLU& B);
    MatrixLU  operator* (const MatrixLU& B) const;
    MatrixLU& operator*= (const MatrixLU& B);
    MatrixLU  operator* (const Real r) const;
    MatrixLU& operator*= (const Real r);
    MatrixLU  operator/ (const Real r) const;
    MatrixLU& operator/= (const Real r);
    Vector    operator* (const Vector& x) const;
    VectorLU  operator* (const VectorLU& x) const;
    Vector&   operator[] (UInt i);
    Real&     operator() (UInt i, UInt j);
    const Vector&   operator[] (UInt i) const;
    const Real&   operator() (UInt i, UInt j) const;

    template<UInt s>
    VectorLU  timesVectorSmall (const VectorSmall<s>& v) const;

    void disp();

protected:
    UInt M_n;
    UInt M_m;
    Matrix M_A;
};

template<UInt s>
VectorLU MatrixLU::timesVectorSmall (const VectorSmall<s>& v) const
{
    VectorLU b (M_n);

    for (UInt i = 0; i < M_n; i++)
        for (UInt j = 0; j < M_m; j++)
        {
            b[i] += M_A[i][j] * v[j];
        }

    return b;
}

}// namespace LifeV


#endif /* MATRIXLU_HPP_ */
