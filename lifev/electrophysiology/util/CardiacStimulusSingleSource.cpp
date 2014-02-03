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
 @brief Class for applying cardiac stimulus represented by
	a current at a single point with a given time and duration

 @date 02-2014
 @author Simone Palamara <palamara.simone@gmail.com>

 @last update 02-2014
 */

#include <lifev/electrophysiology/util/CardiacStimulusSingleSource.hpp>

namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
CardiacStimulusSingleSource::CardiacStimulusSingleSource() :
    M_radius ( 0 ),
    M_totalCurrent ( 0 ),
    M_pacingSite_X ( 0 ),
    M_pacingSite_Y ( 0 ),
    M_pacingSite_Z ( 0 ),
    M_stimulusValue ( 0 ),
    M_startingActivationTime ( 0 ),
    M_endingActivationTime ( 0 ),
{

}


// ===================================================
//! Methods
// ===================================================
Real CardiacStimulusPMJ::appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/ )
{

    Real current = 0.0;
    const Real volumeOfBall = (4. / 3.) * M_PI * M_radius * M_radius * M_radius;
    Real distance = std::sqrt ( (x - M_pacingSite_X) * (x - M_pacingSite_X) + (y - M_pacingSite_Y) * (y - M_pacingSite_Y) + (z - M_pacingSite_Z) * (z - M_pacingSite_Z) );
    if (distance <= M_radius && t >= (M_startingActivationTime) && t <= (M_endingActivationTime) )
    {
       current += M_totalCurrent / volumeOfBall;
    }
    return current;

}

