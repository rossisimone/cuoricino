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
  @brief Mean-field cross-bridge model by FitzHugh-Nagumo.
  @date 01-2013
  @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

  @contributors
  @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICFITZHUGHNAGUMO_H_
#define _IONICFITZHUGHNAGUMO_H_

#include <lifev/heart/solver/IonicModels/HeartIonicModel.hpp>


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! XbModel - This class implements a mean field model.


class IonicFitzHughNagumo : public virtual HeartIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef HeartIonicModel super;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
    typedef RegionMesh<LinearTetra> mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicFitzHughNagumo();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicFitzHughNagumo ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicFitzHughNagumo object
     */
    IonicFitzHughNagumo ( const IonicFitzHughNagumo& model );
    //! Destructor
    virtual ~IonicFitzHughNagumo() {}

    //@}

    //! @name Overloads
    //@{

    IonicFitzHughNagumo& operator= ( const IonicFitzHughNagumo& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& G()        const
    {
        return M_G;
    }
    inline const Real& Vth()      const
    {
        return M_Vth;
    }
    inline const Real& Vp()    const
    {
        return M_Vp;
    }
    inline const Real& Eta1()    const
    {
        return M_Eta1;
    }
    inline const Real& Eta2()    const
    {
        return M_Eta2;
    }
    inline const Real& Eta3()    const
    {
        return M_Eta3;
    }
    inline const Real& Eta()    const
    {
        return M_Eta;
    }
    inline const Real& Gamma()    const
    {
        return M_Gamma;
    }

    inline void setG      ( const Real& G )
    {
        this->M_G       = G;
    }
    inline void setVth    ( const Real& Vth )
    {
        this->M_Vth     = Vth;
    }
    inline void setVp     ( const Real& Vp )
    {
        this->M_Vp      = Vp;
        this->M_Eta     = M_Eta2 / M_Vp;
    }
    inline void setEta1   ( const Real& Eta1 )
    {
        this->M_Eta1    = Eta1;
    }
    inline void setEta2   ( const Real& Eta2 )
    {
        this->M_Eta2    = Eta2;
        this->M_Eta     = M_Eta2 / M_Vp;
        this->M_Gamma   = M_Eta2 * M_Eta3;
    }
    inline void setEta3   ( const Real& Eta3 )
    {
        this->M_Eta3    = Eta3;
        this->M_Gamma   = M_Eta2 * M_Eta3;
    }

    //Compute the rhs on a single node or for the 0D case
    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeRhs ( const std::vector<Real>& v, const Real& Iapp, std::vector<Real>& rhs);


    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp);

    //compute the Jacobian
    void computeJ (Real& a, Real& b, Real& c, Real& d, const std::vector<Real>& v);

    //! Display information about the model
    void showMe();


private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_G;
    Real M_Vth;
    Real M_Vp;
    Real M_Eta1;
    Real M_Eta2;
    Real M_Eta3;
    Real M_Eta;
    Real M_Gamma;


    //@}

}; // class IonicFitzHughNagumo

// ===================================================
//! Constructors
// ===================================================
IonicFitzHughNagumo::IonicFitzHughNagumo()    :
    super   ( 2   ),
    M_G     ( 1.5 ),
    M_Vth   ( 13. ),
    M_Vp    ( 100. ),
    M_Eta1  ( 4.4 ),
    M_Eta2  ( 0.012 ),
    M_Eta3  ( 1.)
{
    M_Eta = M_Eta2 / M_Vp;
    M_Gamma = M_Eta2 * M_Eta3;
}

IonicFitzHughNagumo::IonicFitzHughNagumo ( Teuchos::ParameterList& parameterList   )   :
    super       ( 2 )
{
    M_G     =   parameterList.get ("G", 1.5);
    M_Vth   =   parameterList.get ("Vth", 13.);
    M_Vp    =   parameterList.get ("Vp", 100.);
    M_Eta1  =   parameterList.get ("Eta1", 4.4);
    M_Eta2  =   parameterList.get ("Eta2", 0.012);
    M_Eta3  =   parameterList.get ("Eta3", 1.);
    M_Eta   =   M_Eta2 / M_Vp;
    M_Gamma =   M_Eta2 * M_Eta3;
}

IonicFitzHughNagumo::IonicFitzHughNagumo ( const IonicFitzHughNagumo& model )
{
    M_G     =   model.M_G;
    M_Vth   =   model.M_Vth;
    M_Vp    =   model.M_Vp;
    M_Eta1  =   model.M_Eta1;
    M_Eta2  =   model.M_Eta2;
    M_Eta3  =   model.M_Eta3;
    M_Eta   =   model.M_Eta;
    M_Gamma =   model.M_Gamma;

    M_numberOfEquations = model.M_numberOfEquations;
}

// ===================================================
//! Operator
// ===================================================
IonicFitzHughNagumo& IonicFitzHughNagumo::operator= ( const IonicFitzHughNagumo& model )
{
    M_G     =   model.M_G;
    M_Vth   =   model.M_Vth;
    M_Vp    =   model.M_Vp;
    M_Eta1  =   model.M_Eta1;
    M_Eta2  =   model.M_Eta2;
    M_Eta3  =   model.M_Eta3;
    M_Eta   =   model.M_Eta;
    M_Gamma =   model.M_Gamma;

    M_numberOfEquations = model.M_numberOfEquations;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
//Only gating variables
void IonicFitzHughNagumo::computeRhs (    const   std::vector<Real>&  v,
                                          std::vector<Real>& rhs )
{

    //Real dr = - ( M_epsilon + M_mu1 * v[1] / ( M_mu2 + v[0] ) ) * ( v[1] + M_k * v[0] * ( v[0] - M_a  - 1.0 ) );
    Real dr = M_Eta * v[0] - M_Gamma * v[1] ;

    rhs[0] = dr;

}

//Potential and gating variables
void IonicFitzHughNagumo::computeRhs (    const   std::vector<Real>&  v,
                                          const   Real&           Iapp,
                                          std::vector<Real>& rhs )
{
    Real dV = - ( M_G * v[0] * ( 1.0 - v[0] / M_Vth ) * ( 1.0 - v[0] / M_Vp ) + M_Eta1 * v[0] * v[1] ) + Iapp;
    Real dr = M_Eta * v[0] - M_Gamma * v[1] ;

    rhs[0] = dV;
    rhs[1] = dr;

}


Real IonicFitzHughNagumo::computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp)
{
    return ( - ( M_G * v[0] * ( 1.0 - v[0] / M_Vth ) * ( 1.0 - v[0] / M_Vp ) + M_Eta1 * v[0] * v[1] ) + Iapp );
}

void IonicFitzHughNagumo::computeJ (Real& a, Real& b, Real& c, Real& d, const std::vector<Real>& v)
{
    a = - ( M_G / ( M_Vth * M_Vp ) ) * ( M_Vth * ( M_Vp - 2.0 * v[0] ) + v[0] * ( 3.0 * v[0] - 2.0 * M_Vp ) ) - M_Eta1 * v[1];
    b = -M_Eta1 * v[0];
    c = M_Eta;
    d = -M_Gamma;
}


void IonicFitzHughNagumo::showMe()
{
    std::cout << "\n\n\t\tIonicFitzHughNagumo Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "G    : " << this->G() << std::endl;
    std::cout << "Vth  : " << this->Vth() << std::endl;
    std::cout << "Vp   : " << this->Vp() << std::endl;
    std::cout << "Eta1 : " << this->Eta1() << std::endl;
    std::cout << "Eta2 : " << this->Eta2() << std::endl;
    std::cout << "Eta3 : " << this->Eta3() << std::endl;
    std::cout << "Eta  : " << this->Eta() << std::endl;
    std::cout << "Gamma: " << this->Gamma() << std::endl;

    std::cout << "\n\t\t End of IonicFitzHughNagumo Informations\n\n\n";

}


}

#endif
