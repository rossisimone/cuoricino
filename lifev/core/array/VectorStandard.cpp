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
    @brief Implemetation of a small vector used with MatrixStandard

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 01-04-2013
 */

#include <lifev/core/array/VectorStandard.hpp>

namespace LifeV
{

VectorStandard::VectorStandard (UInt size, Real r)
    : std::vector<Real> (size, r)
{
}

VectorStandard::VectorStandard (std::vector<Real> vector)
    : std::vector<Real> (vector)
{
}


Real VectorStandard::norm2() const
{
    Real norm (0.0);

    for (UInt i = 0; i < this->size(); i++)
    {
        norm += (this->at (i) ) * (this->at (i) );
    }

    return std::sqrt (norm);
}

VectorStandard VectorStandard::operator+ (const VectorStandard& w) const
{
    return VectorStandard (*this) += w;
}

VectorStandard& VectorStandard::operator+= (const VectorStandard& w)
{
    for (Int i = 0; i < this->size(); i++)
    {
        this->at (i) += w[i];
    }

    return *this;
}

VectorStandard VectorStandard::operator- (const VectorStandard& w) const
{
    return VectorStandard (*this) -= w;
}

VectorStandard& VectorStandard::operator-= (const VectorStandard& w)
{
    for (Int i = 0; i < this->size(); i++)
    {
        this->at (i) -= w[i];
    }

    return *this;
}

VectorStandard VectorStandard::operator* (const Real r) const
{
    return VectorStandard (*this) *= r;
}

VectorStandard& VectorStandard::operator*= (const Real r)
{
    for (Int i = 0; i < this->size(); i++)
    {
        this->at (i) *= r;
    }

    return *this;
}

VectorStandard VectorStandard::operator/ (const Real r) const
{
    return VectorStandard (*this) /= r;
}

VectorStandard& VectorStandard::operator/= (const Real r)
{
    for (Int i = 0; i < this->size(); i++)
    {
        this->at (i) /= r;
    }

    return *this;
}

void VectorStandard::disp() const
{
    for (UInt i = 0; i < this->size(); i++)
    {
        cout << "   " << this->at (i);
    }
}


} //namespace LifeV
