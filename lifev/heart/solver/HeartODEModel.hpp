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
//  @brief Mother class for mean-field corssbridge models (ODE stands for crossbridge)
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
#ifndef _HEARTODEMODEL_H_
#define _HEARTODEMODEL_H_

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV{
class HeartODEModel{

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
	HeartODEModel();

	HeartODEModel( const HeartODEModel &ODE );

    //! Destructor
    virtual ~HeartODEModel() {};

    //@}

    //! @name Methods
    //@{

    //virtual void updateRepeated( )=0;

    //! Solves the ODE model
    //virtual void solveODEModel( const vector_Type& Calcium,
    //                              const Real timeStep )=0;

    //@}

    //! @name Overloads
    //@{

    HeartODEModel& operator=( const HeartODEModel &ODE );

    //@}


    //@}



};


// ===================================================
//! Constructors
// ===================================================
HeartODEModel::HeartODEModel()
{
}



HeartODEModel::HeartODEModel( const HeartODEModel &ODE )
{
}

// ===================================================
//! Methods
// ===================================================
HeartODEModel& HeartODEModel::operator =( const HeartODEModel &ODE )
{
	return 		*this;
}

}

#endif
