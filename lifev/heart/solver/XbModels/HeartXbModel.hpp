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
//  @brief Mother class for mean-field corssbridge models (Xb stands for crossbridge)
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
#ifndef _HEARTXBMODEL_H_
#define _HEARTXBMODEL_H_

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV{
class HeartXbModel{

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
	HeartXbModel();

	HeartXbModel( const HeartXbModel &Xb );

    //! Destructor
    virtual ~HeartXbModel() {};

    //@}

    //! @name Methods
    //@{

    //virtual void updateRepeated( )=0;

    //! Solves the xb model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                              const Real timeStep )=0;

    //@}

    //! @name Overloads
    //@{

    HeartXbModel& operator=( const HeartXbModel &Xb );

    //@}


    //@}



};


// ===================================================
//! Constructors
// ===================================================
HeartXbModel::HeartXbModel()
{
}



HeartXbModel::HeartXbModel( const HeartXbModel &Xb )
{
}

// ===================================================
//! Methods
// ===================================================
HeartXbModel& HeartXbModel::operator =( const HeartXbModel &Xb )
{
	return 		*this;
}

}

#endif
