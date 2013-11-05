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
#ifndef _ELECTROIONICMODEL_H_
#define _ELECTROIONICMODEL_H_

#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/electrophysiology/util/CardiacStimulus.hpp>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

namespace LifeV
{
class ElectroIonicModel
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
	typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
	typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
	typedef boost::function<
			Real(const Real& t, const Real& x, const Real& y, const Real& z,
					const ID& i)> function_Type;
    //@}

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param Epetra communicator
     */
    ElectroIonicModel();

    ElectroIonicModel ( int n );

    ElectroIonicModel ( int n, int g );

    ElectroIonicModel ( const ElectroIonicModel& Ionic );

    //! Destructor
    virtual ~ElectroIonicModel() {};

    //@}

    //! @name Methods
    //@{
    inline const short int Size() const
    {
        return M_numberOfEquations;
    }
    inline const short int numberOfGatingVariables() const
    {
        return M_numberOfGatingVariables;
    }
    inline const Real membraneCapacitance() const
    {
    	return M_membraneCapacitance;
    }
    inline const Real appliedCurrent() const
    {
    	return M_appliedCurrent;
    }
    inline vectorPtr_Type appliedCurrentPtr()
    {
    	return M_appliedCurrentPtr;
    }

    inline const std::vector<Real> restingConditions() const
    {
        return M_restingConditions;
    }
    inline const function_Type pacaingProtocol() const
    {
    	return M_pacingProtocol;
    }

    inline void setMembraneCapacitance( const Real p )
    {
    	M_membraneCapacitance = p;
    }
    inline void setAppliedCurrent(const Real p)
    {
    	M_appliedCurrent = p;
    }

	inline void setAppliedCurrentPtr(const vectorPtr_Type p) {
		this->M_appliedCurrentPtr = p;
	}
	inline void setAppliedCurrent(const vector_Type& p) {
		M_appliedCurrentPtr.reset( new vector_Type( p ) );
	}
	inline void setAppliedCurrentFromFunction(function_Type& f, feSpacePtr_Type feSpacePtr,
	                                          Real time = 0.0) {

	    feSpacePtr -> interpolate(
	                    static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type>(f),
	                    *M_appliedCurrentPtr, time);
	}
	inline void setAppliedCurrentFromCardiacStimulus( CardiacStimulus& stimulus, feSpacePtr_Type feSpacePtr, Real time = 0.0) {

	    // boost::ref() is needed here because otherwise a copy of the base object is reinstantiated
	    function_Type f = boost::bind(&CardiacStimulus::appliedCurrent, boost::ref(stimulus), _1, _2, _3, _4, _5 );

	    feSpacePtr -> interpolate(
	                    static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type>(f),
	                    *M_appliedCurrentPtr, time);
	}

	inline void setPacingProtocol( function_Type pacingProtocol )
	{
		M_pacingProtocol = pacingProtocol;
	}


    //virtual void updateRepeated( )=0;

    //! Solves the ODE model
    //virtual void solveODEModel( const vector_Type& Calcium,
    //                              const Real timeStep )=0;

    //@}

    //! @name Overloads
    //@{

    ElectroIonicModel& operator= ( const ElectroIonicModel& Ionic );

    virtual matrix_Type getJac(const vector_Type& v, Real h=1.0e-8);

    virtual vector< vector<Real> > getJac(const vector<Real>& v, Real h=1.0e-8);

    virtual void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs ) = 0;

    virtual void computeNonGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs )
    {
    	if(M_numberOfGatingVariables == 0 ) computeGatingRhs(v, rhs);
    };


    virtual void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs) = 0;

    virtual void computeGatingVariablesWithRushLarsen ( std::vector<Real>& /*v*/, const Real /*dt*/ ) {}

    //Compute the rhs on a mesh/ 3D case
    virtual void computeGatingRhs ( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );

    //Compute the rhs on a mesh/ 3D case
    virtual void computeNonGatingRhs ( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );

    virtual void computeRhs ( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );

    virtual void computeGatingVariablesWithRushLarsen ( std::vector<vectorPtr_Type>& v, const Real dt );
    // compute the rhs with state variable interpolation
    virtual Real computeLocalPotentialRhs ( const std::vector<Real>& v ) = 0;

    virtual void computePotentialRhsICI ( const std::vector<vectorPtr_Type>& v,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          matrix_Type&                        massMatrix );

    virtual void computePotentialRhsSVI ( const std::vector<vectorPtr_Type>& v,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          FESpace<mesh_Type, MapEpetra>&  uFESpace );

    virtual void computePotentialRhsSVI ( const std::vector<vectorPtr_Type>& v,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          FESpace<mesh_Type, MapEpetra>&  uFESpace,
                                          vector_Type& disp,
                                          boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > dispFESPace);
    //@}

    virtual void showMe() = 0;

    //initialize with given conditions
    virtual void initialize( std::vector<Real>& v );

    //initialize with given conditions
    virtual void initialize( std::vector<vectorPtr_Type>& v );
    //@}

protected:

    short int  M_numberOfEquations;
    short int  M_numberOfGatingVariables;
    std::vector<Real> M_restingConditions;
    Real M_membraneCapacitance;
    Real M_appliedCurrent;
    vectorPtr_Type M_appliedCurrentPtr;
    function_Type M_pacingProtocol;


};


}

#endif
