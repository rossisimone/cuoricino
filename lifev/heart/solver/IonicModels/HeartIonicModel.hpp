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

#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/MapEpetra.hpp>


namespace LifeV
{
class HeartIonicModel
{

public:
    //! @name Type definitions
    //@{

    typedef VectorEpetra                            vector_Type;
    typedef boost::shared_ptr<VectorEpetra>         vectorPtr_Type;
    typedef boost::shared_ptr<VectorElemental>  elvecPtr_Type;
    typedef RegionMesh<LinearTetra>                 mesh_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    //@}

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param Epetra communicator
     */
    HeartIonicModel();

    HeartIonicModel ( int n );

    HeartIonicModel ( const HeartIonicModel& Ionic );

    //! Destructor
    virtual ~HeartIonicModel() {};

    //@}

    //! @name Methods
    //@{
    inline const short int& Size() const
    {
        return M_numberOfEquations;
    }

    //virtual void updateRepeated( )=0;

    //! Solves the ODE model
    //virtual void solveODEModel( const vector_Type& Calcium,
    //                              const Real timeStep )=0;

    //@}

    //! @name Overloads
    //@{

    HeartIonicModel& operator= ( const HeartIonicModel& Ionic );

    virtual matrix_Type getJac(const vector_Type& v, Real h=1.0e-8);

    virtual vector< vector<Real> > getJac(const vector<Real>& v, Real h=1.0e-8);

    virtual void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs) = 0;

    virtual void computeRhs ( const std::vector<Real>& v, const Real& Iapp, std::vector<Real>& rhs) = 0;

    virtual void computeRhs ( const vector_Type& v, const Real& Iapp, vector_Type& rhs);

    //Compute the rhs on a mesh/ 3D case
    virtual void computeRhs ( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );

    virtual void computeRhs ( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Iapp, std::vector<vectorPtr_Type>& rhs );

    // compute the rhs with state variable interpolation
    virtual Real computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp) = 0;

    virtual void computePotentialRhsICI ( const std::vector<vectorPtr_Type>& v,
                                          const VectorEpetra&                 Iapp,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          matrix_Type&                        massMatrix );

    virtual void computePotentialRhsSVI ( const std::vector<vectorPtr_Type>& v,
                                          const VectorEpetra&                 Iapp,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          FESpace<mesh_Type, MapEpetra>&  uFESpace );
    //@}

    virtual void showMe() = 0;
    //@}

protected:

    short int  M_numberOfEquations;


};


}

#endif
