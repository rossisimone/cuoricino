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
    @brief Implementation of LU for MatrixStandard

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 22-04-2013
 */


#include <lifev/core/array/MatrixStandard.hpp>

#ifndef MATRIXSTANDARDSOLVER_HPP_
#define MATRIXSTANDARDSOLVER_HPP_

namespace LifeV
{

class MatrixStandardSolver
{
public:
    MatrixStandardSolver() {};
    virtual ~MatrixStandardSolver() {};

    void LU (MatrixStandard& A, MatrixStandard& P, MatrixStandard& Q, MatrixStandard& L, MatrixStandard& U) const;
    void LU (MatrixStandard& A, MatrixStandard& P, MatrixStandard& Q, MatrixStandard& L, MatrixStandard& U,
             const MatrixStandard& I) const;
    void Pivot (const MatrixStandard& A, MatrixStandard& Pk, MatrixStandard& Qk, UInt k) const;
    void solveU (const MatrixStandard& A, VectorStandard& b) const;
    void solveL (const MatrixStandard& A, VectorStandard& b) const;

};

} //namespace LifeV


#endif /* MATRIXSTANDARDSOLVER_HPP_ */
