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
  @brief Mean-field cross-bridge model by Negroni and Lascano 1996.
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _XBNEGRONILASCANO96_H_
#define _XBNEGRONILASCANO96_H_

#include <lifev/heart/solver/XbModels/HeartXbModel.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! XbModel - This class implements a mean field model.


class XbNegroniLascano96 : public virtual HeartXbModel
{

public:
    //! @name Type definitions
    //@{
    typedef HeartXbModel super;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    XbNegroniLascano96();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    XbNegroniLascano96 ( Teuchos::ParameterList& parameterList );

    /*!
     * @param XbNegroniLascano96 object
     */
    XbNegroniLascano96 ( const XbNegroniLascano96& Xb );
    //! Destructor
    virtual ~XbNegroniLascano96() {}

    //@}

    //! @name Overloads
    //@{

    XbNegroniLascano96& operator= ( const XbNegroniLascano96& Xb );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& Alpha1() const
    {
        return M_alpha1;
    }
    inline const Real& Alpha2() const
    {
        return M_alpha2;
    }
    inline const Real& Alpha3()     const
    {
        return M_alpha3;
    }
    inline const Real& Alpha4()     const
    {
        return M_alpha4;
    }
    inline const Real& Alpha5()     const
    {
        return M_alpha5;
    }
    inline const Real& Beta1()      const
    {
        return M_beta1;
    }
    inline const Real& Beta2()      const
    {
        return M_beta2;
    }
    inline const Real& Beta3()      const
    {
        return M_beta3;
    }

    inline void setAlpha1   ( const Real& alpha1 )
    {
        this->M_alpha1 = alpha1;
    }
    inline void setAlpha2   ( const Real& alpha2 )
    {
        this->M_alpha2 = alpha2;
    }
    inline void setAlpha3   ( const Real& alpha3 )
    {
        this->M_alpha3 = alpha3;
    }
    inline void setAlpha4   ( const Real& alpha4 )
    {
        this->M_alpha4 = alpha4;
    }
    inline void setAlpha5   ( const Real& alpha5 )
    {
        this->M_alpha5 = alpha5;
    }
    inline void setBeta1    ( const Real& beta1  )
    {
        this->M_beta1  =  beta1;
    }
    inline void setBeta2    ( const Real& beta2  )
    {
        this->M_beta2  =  beta2;
    }
    inline void setBeta3    ( const Real& beta3  )
    {
        this->M_beta3  =  beta3;
    }


    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeRhs ( const std::vector<Real>& v, const Real& Ca, const Real& vel,  std::vector<Real>& rhs);
    //Compute the rhs on a mesh/ 3D case
    void computeRhs ( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Ca, const VectorEpetra& vel, std::vector<vectorPtr_Type>& rhs );

    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_alpha1;
    Real M_alpha2;
    Real M_alpha3;
    Real M_alpha4;
    Real M_alpha5;
    Real M_beta1;
    Real M_beta2;
    Real M_beta3;



    //! Xb states == equivalent to the number of equations
    //short int M_numberOfStates;

    //@}

}; // class XbNegroniLascano96

// ===================================================
//! Constructors
// ===================================================
XbNegroniLascano96::XbNegroniLascano96()    :
    super    ( 3   ),
    M_alpha1 ( 0.0 ),
    M_alpha2 ( 0.0 ),
    M_alpha3 ( 0.0 ),
    M_alpha4 ( 0.0 ),
    M_alpha5 ( 0.0 ),
    M_beta1  ( 0.0 ),
    M_beta2  ( 0.0 ),
    M_beta3  ( 0.0 )
{
}

XbNegroniLascano96::XbNegroniLascano96 ( Teuchos::ParameterList& parameterList   )   :
    super    ( 3   )
{
    M_alpha1 =  parameterList.get ("alpha1", 0.0);
    M_alpha2 =  parameterList.get ("alpha2", 0.0);
    M_alpha3 =  parameterList.get ("alpha3", 0.0);
    M_alpha4 =  parameterList.get ("alpha4", 0.0);
    M_alpha5 =  parameterList.get ("alpha5", 0.0);
    M_beta1  =  parameterList.get ("beta1",  0.0);
    M_beta2  =  parameterList.get ("beta1",  0.0);
    M_beta3  =  parameterList.get ("beta1",  0.0);
}

XbNegroniLascano96::XbNegroniLascano96 ( const XbNegroniLascano96& Xb )
{
    M_alpha1 =  Xb.M_alpha1;
    M_alpha2 =  Xb.M_alpha2;
    M_alpha3 =  Xb.M_alpha3;
    M_alpha4 =  Xb.M_alpha4;
    M_alpha5 =  Xb.M_alpha5;
    M_beta1  =  Xb.M_beta1;
    M_beta2  =  Xb.M_beta2;
    M_beta3  =  Xb.M_beta3;

    M_numberOfEquations = Xb.M_numberOfEquations;
}

// ===================================================
//! Operator
// ===================================================
XbNegroniLascano96& XbNegroniLascano96::operator= ( const XbNegroniLascano96& Xb )
{
    M_alpha1 =  Xb.M_alpha1;
    M_alpha2 =  Xb.M_alpha2;
    M_alpha3 =  Xb.M_alpha3;
    M_alpha4 =  Xb.M_alpha4;
    M_alpha5 =  Xb.M_alpha5;
    M_beta1  =  Xb.M_beta1;
    M_beta2  =  Xb.M_beta2;
    M_beta3  =  Xb.M_beta3;

    M_numberOfEquations = Xb.M_numberOfEquations;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
// v(0) = TCa
// v(1) = TCas
// v(2) = Ts
void XbNegroniLascano96::computeRhs (    const   std::vector<Real>&  v,
                                         const   Real&           Ca,
                                         const   Real&           vel,
                                         std::vector<Real>& rhs )
{

    Real Q1 = M_alpha1 * Ca * ( 1.0 - v[0] - v[1] - v[2] ) - M_beta1 * v[0];

    Real TCaEff (1.0);
    Real Q2 = M_alpha2 * v[0] * TCaEff - M_beta2 * v[1];
    Real Q3 = M_alpha3 * v[1] - M_beta3 * Ca * v[2];
    Real param (1.0);
    Real Q4 = M_alpha4 * v[1] + M_alpha5 * ( param * vel ) * ( param * vel ) * v[2];
    Real Q5 = M_alpha5 * ( param * vel ) * ( param * vel ) * v[1];

    rhs[0] = Q1 - Q2;
    rhs[1] = Q2 - Q3 - Q5;
    rhs[2] = Q3 - Q4;


}

void XbNegroniLascano96::computeRhs (    const std::vector<vectorPtr_Type>& v,
                                         const VectorEpetra& Ca,
                                         const VectorEpetra& vel,
                                         std::vector<vectorPtr_Type>& rhs )
{

    int nodes = Ca.epetraVector().MyLength();


    std::vector<Real>   localVec ( 3, 0.0 );
    std::vector<Real>   localRhs ( 3, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = Ca.blockMap().GID (k);

        localVec.at (0) = ( * ( v.at (0) ) ) [j];
        localVec.at (1) = ( * ( v.at (1) ) ) [j];
        localVec.at (2) = ( * ( v.at (2) ) ) [j];

        computeRhs ( localVec, Ca[j], vel[j], localRhs );

        ( * ( rhs.at (0) ) ) [j] =  localRhs.at (0);
        ( * ( rhs.at (1) ) ) [j] =  localRhs.at (1);
        ( * ( rhs.at (2) ) ) [j] =  localRhs.at (2);

    }

}


void XbNegroniLascano96::showMe()
{
    std::cout << "\n\n\t\tXbNegroniLascano96 Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "alpha1: " << this->Alpha1() << std::endl;
    std::cout << "alpha2: " << this->Alpha2() << std::endl;
    std::cout << "alpha3: " << this->Alpha3() << std::endl;
    std::cout << "alpha4: " << this->Alpha4() << std::endl;
    std::cout << "alpha5: " << this->Alpha5() << std::endl;
    std::cout << "beta1: "  << this->Beta1()  << std::endl;
    std::cout << "beta2: "  << this->Beta2()  << std::endl;
    std::cout << "beta3: "  << this->Beta3()  << std::endl;
    std::cout << "\n\t\t End of XbNegroniLascano96 Informations\n\n\n";

}


}

#endif
