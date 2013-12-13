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
    @brief Implemetation of a small vector used with MatrixLU

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 01-04-2013
 */

#ifndef VECTORLU_HPP_
#define VECTORLU_HPP_

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


namespace LifeV
{

class VectorLU
{
    typedef std::vector< std::vector<Real> > Matrix;
    typedef std::vector<Real> Vector;

public:
    VectorLU (UInt n, Real r = 0.0);
    VectorLU (const Vector& v);
    virtual ~VectorLU() {};

    Real norm2() const;
    UInt size() const;      //togli

    Vector& getVector();
    void disp();

    Real& operator[] (UInt i);
    Real& operator() (UInt i);
    const Real& operator[] (UInt i) const;
    const Real& operator() (UInt i) const;
    VectorLU operator+ (const VectorLU& w) const;
    VectorLU& operator+= (const VectorLU& w);
    VectorLU operator- (const VectorLU& w) const;
    VectorLU& operator-= (const VectorLU& w);
    VectorLU operator* (const Real r) const;
    VectorLU& operator*= (const Real r);
    VectorLU operator/ (const Real r) const;
    VectorLU& operator/= (const Real r);


private:
    UInt M_n; //togli
    Vector M_v;
};

}// namespace LifeV



#endif /* VECTORLU_HPP_ */
