////@HEADER
///*
//*******************************************************************************
//
//    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
//    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
//
//    This file is part of LifeV.
//
//    LifeV is free software; you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    LifeV is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.
//
//*******************************************************************************
//*/
////@HEADER
//
///*!
//  @file
//  @brief Mother class for ionic models
//
//  @date 01-2013
//  @author Simone Rossi <simone.rossi@epfl.ch>
//
//  @contributors
//  @mantainer Simone Rossi <simone.rossi@epfl.ch>
//  @last update 02-2012
// */
//
//
#ifndef _HEARTIONICMODEL_H_
#define _HEARTIONICMODEL_H_

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV{
class HeartIonicModel{

public:
    //! @name Type definitions
    //@{

	typedef VectorEpetra vector_Type;

	//@}

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param Epetra communicator
     */
	HeartIonicModel();

	HeartIonicModel( int n );

	HeartIonicModel( const HeartIonicModel &Ionic );

    //! Destructor
    virtual ~HeartIonicModel() {};

    //@}

    //! @name Methods
    //@{
    inline const short int& Size() const { return M_numberOfEquations; }

    //virtual void updateRepeated( )=0;

    //! Solves the ODE model
    //virtual void solveODEModel( const vector_Type& Calcium,
    //                              const Real timeStep )=0;

    //@}

    //! @name Overloads
    //@{

   HeartIonicModel& operator=( const HeartIonicModel &Ionic );

    //@}


    //@}

protected:

    short int  M_numberOfEquations;


};


// ===================================================
//! Constructors
// ===================================================
HeartIonicModel::HeartIonicModel():
		M_numberOfEquations(0)
{
}

HeartIonicModel::HeartIonicModel( int n ):
		M_numberOfEquations(n)
{
}



HeartIonicModel::HeartIonicModel( const HeartIonicModel &Ionic ):
		M_numberOfEquations( Ionic.Size() )
{
}

// ===================================================
//! Methods
// ===================================================
HeartIonicModel& HeartIonicModel::operator =( const HeartIonicModel &Ionic )
{
	return 		*this;
}

}

#endif
