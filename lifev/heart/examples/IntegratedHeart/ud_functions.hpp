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
 *  @file
 *  @brief Deprecated file for the Integrated Heart example
 *
 *  @date 2012-09-25
 *  @author Toni Lassila <toni.lassila@epfl.ch>
 *          Paolo Crosetto <crosetto@iacspc70.epfl.ch>

 *  @maintainer Toni Lassila <toni.lassila@epfl.ch>
 *
*/

#ifndef UDF_HPP
#define UDF_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>

namespace LifeV
{

class aortaVelIn
{
public:
    static Real S_timestep;
};

Real f (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real u1 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real fZero (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real inPressure (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
Real outPressure (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

// Initial velocity
Real u0 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
Real p0 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
Real E (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);
Real hydrostatic (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);

Real hydro (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i);
Real u2 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real flow_mitral (const Real& t);
Real u_mitral (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);


// Initial displacement and velocity
Real d0 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real w0 (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real uInterpolated (const Real& time, const Real& x, const Real& y, const Real& z, const ID& i);

Real aortaPhisPress (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i);
Real vinit (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i);

Real u2normal (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);

Real fluxFunction (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i);

Real squareSinusoidalFluxFunction (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);

Real u_test (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
}



#endif
