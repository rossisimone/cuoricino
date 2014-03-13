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
  @brief Ionic model of Hodgkin and Huxley
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICHODGKINHUXLEY_H_
#define _IONICHODGKINHUXLEY_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.

class IonicHodgkinHuxley : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroIonicModel                         super;
    typedef boost::shared_ptr<VectorEpetra>         vectorPtr_Type;
    typedef boost::shared_ptr<VectorElemental>  elvecPtr_Type;
    typedef RegionMesh<LinearTetra>                 mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicHodgkinHuxley();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicHodgkinHuxley ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicHodgkinHuxley object
     */
    IonicHodgkinHuxley ( const IonicHodgkinHuxley& model );
    //! Destructor
    virtual ~IonicHodgkinHuxley() {}

    //@}

    //! @name Overloads
    //@{

    IonicHodgkinHuxley& operator= ( const IonicHodgkinHuxley& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& gNa()             const
    {
        return M_gNa;
    }

    inline const Real& gK()             const
    {
        return M_gK;
    }

    inline const Real& gL()             const
    {
        return M_gL;
    }


    inline const Real& vL()             const
    {
        return M_vL;
    }

    inline const Real& vK()             const
    {
        return M_vK;
    }

    inline const Real& vNa()             const
    {
        return M_vNa;
    }




    inline void setgNa           ( const Real& p )
    {
        this->M_gNa        = p;
    }

    inline void setgK           ( const Real& p )
    {
        this->M_gK        = p;
    }

    inline void setgL           ( const Real& p )
    {
        this->M_gL        = p;
    }

    inline void setvL           ( const Real& p )
    {
        this->M_vL        = p;
    }

    inline void setvNa          ( const Real& p )
    {
        this->M_vNa        = p;
    }

    inline void setvK           ( const Real& p )
    {
        this->M_vK         = p;
    }


    //inline const short int& Size() const { return M_numberOfEquations; }
    //@}


    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt );


    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_gNa;
    Real M_gK;
    Real M_gL;
    Real M_vNa;
    Real M_vK;
    Real M_vL;


    //! Xb states == equivalent to the number of equations
    //short int M_numberOfEquations;

    //@}

}; // class IonicMinimalModel



}

#endif
