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
 @brief Base class for applying cardiac stimulus

 @date 11-2013
 @author Toni Lassila <toni.lassila@epfl.ch>

 @last update 11-2013
 */


#ifndef CARDIACSTIMULUS_HPP_
#define CARDIACSTIMULUS_HPP_

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{

class CardiacStimulus
{

public:

    //! @name Type definitions
    //@{
    typedef VectorEpetra                    vector_Type;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    CardiacStimulus();

    //! Destructor
    virtual ~CardiacStimulus() {};

    //@}

    //! @name Get Methods
    //@{

    //@}

    //! @name Set Methods
    //@{


    //@}

    //! @name Copy Methods
    //@{

    //@}

    //! @name Methods
    //@{
    inline virtual Real appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
    {
        return 0.0;
    }

    //@}

private:

};


} // namespace LifeV

#endif /* CARDIACSTIMULUS_HPP_ */
