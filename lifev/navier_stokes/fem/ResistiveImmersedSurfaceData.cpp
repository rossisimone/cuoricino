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
    @brief This file contains the data for resistive immersed surface imposition

    @author Toni Lassila <toni.lassila@epfl.ch>
    @contributor Aymen Laadhari <aymen.laadhari@epfl.ch>

    @date 27-03-2014
 */

#include <lifev/navier_stokes/fem/ResistiveImmersedSurfaceData.hpp>

namespace LifeV {

ResistiveImmersedSurfaceData::ResistiveImmersedSurfaceData() {}

void
ResistiveImmersedSurfaceData::setup( const GetPot& dataFile)
{
    M_resistance = dataFile( "immersedSurface/resistance", 0.0 );
    M_nofSurfaces = dataFile( "immersedSurface/numberOfSurfaces", 0 );

    for (int i=1; i <= M_nofSurfaces; i++)
    {
        //dataFile( "immersedSurface/surface" + i + "/", 0 );
    }
}

Real
ResistiveImmersedSurfaceData::eval( )
{
    return 0.;
}

}
