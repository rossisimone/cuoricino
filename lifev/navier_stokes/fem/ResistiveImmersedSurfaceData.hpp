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

#ifndef RESISTIVEIMMERSEDSURFACEDATA_H
#define RESISTIVEIMMERSEDSURFACEDATA_H

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>

namespace LifeV {


class SmoothHeavisideFct
{
public:
    typedef Real return_Type;

    return_Type operator()(const Real& value)
    {
        return 0.0;
    }

    SmoothHeavisideFct(){}
    SmoothHeavisideFct(const SmoothHeavisideFct&){}
    ~SmoothHeavisideFct(){}
};


class SmoothDeltaFct
{
public:
    typedef Real return_Type;

    return_Type operator()(const Real& value)
    {
        return 0.0;
    }

    SmoothDeltaFct(){}
    SmoothDeltaFct(const SmoothDeltaFct&){}
    ~SmoothDeltaFct(){}
};

class ResistiveImmersedSurfaceData
{
public:
    //! @name Public Types
    //@{

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    ResistiveImmersedSurfaceData();

    //! virtual destructor
    virtual ~ResistiveImmersedSurfaceData() {}

    //@}

    //! @name Methods
    //@{

    //! Setup that initializes the class from a data file
    void setup( const GetPot& dataFile );

    //! Evaluates the smooth level set term of the immersed surface(s)
    Real eval( );

    //! Set methods

    //! Get methods
    inline Real resistance()
    {
        return M_resistance;
    }

    inline Real psiFunction( const Real& t, const Real& x, const Real& y, const Real& z , const ID& i)
    {

    }

    inline Real phiFunction( const Real& t, const Real& x, const Real& y, const Real& z , const ID& i)
    {

    }

    inline boost::shared_ptr<SmoothHeavisideFct> heavisideFunctor()
    {

    }

    inline boost::shared_ptr<SmoothDeltaFct> diracFunctor()
    {

    }

    //@}

private:
    //private members

    Real M_resistance;
    Real M_epsilon;
    Real M_unitsFactor;
    UInt M_nofSurfaces;

};

}

#endif /* RESISTIVEIMMERSEDSURFACEDATA_H */
