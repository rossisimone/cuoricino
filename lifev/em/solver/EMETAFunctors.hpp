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
    @brief  Heart Functors for the Luo-Rudy Kinetics

    @date 04−2010
    @author

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */


#ifndef _EMFUNCTORS_H_
#define _EMFUNCTORS_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/fem/FESpace.hpp>



namespace LifeV
{


class FLRelationship
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<1>& I4f)
    {
    	Real i4f = I4f[0];
    	Real Force = this->operator()(i4f);
    	return Force;
    }

    return_Type operator() (const Real& I4f)
    {

    	if(I4f > 0.87277 && I4f < 1.334)
    	{
			Real d0 = -4.333618335582119e3;
			Real d1 = 2.570395355352195e3;
			Real e1 = -2.051827278991976e3;
			Real d2 = 1.329536116891330e3;
			Real e2 = 0.302216784558222e3;
			Real d3 = 0.104943770305116e3;
			Real e3 = 0.218375174229422e3;
			Real l0 = 1.95;

			Real Force = d0/2 + d1 * std::sin(I4f * l0)
							  + e1 * std::cos(I4f * l0)
							  + d2 * std::sin(2 * I4f * l0)
							  + e2 * std::cos(2 * I4f * l0)
							  + d3 * std::sin(3 * I4f * l0)
							  + e3 * std::cos(3 * I4f * l0);
			return Force;
    	}
    	else
    		return 0.0;
    }

    FLRelationship() {}
    FLRelationship (const FLRelationship&) {}
    ~FLRelationship() {}
};



};



#endif
