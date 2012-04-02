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
    @brief This file implements the standard selector for internal entities

    @author
    @contributor Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date

    A more detailed description of the file (if necessary)
 */

#include <lifev/core/mesh/InternalEntitySelector.hpp>

namespace LifeV
{
// LF REGIONMESH
// This class was meant to separate internal from boundary flags. With the
// new way of selecting boundary entities this will be useless!
const markerID_Type InternalEntitySelector::defMarkFlag(markerID_Type(100000));

// ===================================================
// Constructors & Destructor
// ===================================================

InternalEntitySelector::InternalEntitySelector():
M_watermarkFlag( defMarkFlag )
{}

InternalEntitySelector::InternalEntitySelector(const markerID_Type & w):
M_watermarkFlag( w )
{}

// ===================================================
// Operators
// ===================================================

bool
InternalEntitySelector::operator()(markerID_Type const & test) const
{
    return (test==markerID_Type(0) || test > M_watermarkFlag );
}


} // Namespace LifeV
