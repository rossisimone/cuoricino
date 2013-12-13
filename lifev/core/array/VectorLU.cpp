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

#include <lifev/core/array/VectorLU.hpp>

namespace LifeV
{

VectorLU::VectorLU (UInt n, Real r)
    : M_n (n), M_v (n, r)
{
}

VectorLU::VectorLU (const std::vector<Real>& v)
    : M_n (v.size() ), M_v (M_n, 0.0)
{
    for (Int i = 0; i < M_n; i++)
    {
        M_v[i] = v[i];
    }
}

Real VectorLU::norm2() const
{
    Real norm (0.0);

    for (UInt i = 0; i < M_n; i++)
    {
        norm += M_v[i] * M_v[i];
    }

    return std::sqrt (norm);
}

UInt VectorLU::size() const
{
    return M_n;
}

std::vector<Real>& VectorLU::getVector()
{
    return M_v;
}

Real& VectorLU::operator[] (UInt i)
{
    return M_v[i];
}

Real& VectorLU::operator() (UInt i)
{
    return M_v[i];
}

const Real& VectorLU::operator[] (UInt i) const
{
    return M_v[i];
}

const Real& VectorLU::operator() (UInt i) const
{
    return M_v[i];
}

VectorLU VectorLU::operator+ (const VectorLU& w) const
{
    return VectorLU (*this) += w;
}

VectorLU& VectorLU::operator+= (const VectorLU& w)
{
    for (Int i = 0; i < M_n; i++)
    {
        M_v[i] += w[i];
    }

    return *this;
}

VectorLU VectorLU::operator- (const VectorLU& w) const
{
    return VectorLU (*this) -= w;
}

VectorLU& VectorLU::operator-= (const VectorLU& w)
{
    for (Int i = 0; i < M_n; i++)
    {
        M_v[i] -= w[i];
    }

    return *this;
}

VectorLU VectorLU::operator* (const Real r) const
{
    return VectorLU (*this) *= r;
}

VectorLU& VectorLU::operator*= (const Real r)
{
    for (Int i = 0; i < M_n; i++)
    {
        M_v[i] *= r;
    }

    return *this;
}

VectorLU VectorLU::operator/ (const Real r) const
{
    return (*this) * (1.0 / r);
}

VectorLU& VectorLU::operator/= (const Real r)
{
    *this *= (1.0 / r);

    return *this;
}

void VectorLU::disp()
{
    for (UInt i = 0; i < M_n; i++)
    {
        cout << "   " << M_v[i];
    }
}


} //namespace LifeV
