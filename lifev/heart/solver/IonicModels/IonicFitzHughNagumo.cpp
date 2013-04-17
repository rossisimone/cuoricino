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

#include <lifev/heart/solver/IonicModels/IonicFitzHughNagumo.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
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
void IonicFitzHughNagumo::computeRhs ( const VectorSmall<2>& v, const Real& Iapp, VectorSmall<2>& rhs)
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

MatrixSmall<2,2> IonicFitzHughNagumo::computeJ (const Real& t, const std::vector<Real>& v)
{
	MatrixSmall<2,2> J;
    J(0,0) = - ( M_G / ( M_Vth * M_Vp ) ) * ( M_Vth * ( M_Vp - 2.0 * v[0] ) + v[0] * ( 3.0 * v[0] - 2.0 * M_Vp ) ) - M_Eta1 * v[1];
    J(0,1) = -M_Eta1 * v[0];
    J(1,0) = M_Eta;
    J(1,1) = -M_Gamma;

    return J;
}

MatrixSmall<2,2> IonicFitzHughNagumo::computeJ (const Real& t, const VectorSmall<2>& v)
{
	MatrixSmall<2,2> J;
    J(0,0) = - ( M_G / ( M_Vth * M_Vp ) ) * ( M_Vth * ( M_Vp - 2.0 * v[0] ) + v[0] * ( 3.0 * v[0] - 2.0 * M_Vp ) ) - M_Eta1 * v[1];
    J(0,1) = -M_Eta1 * v[0];
    J(1,0) = M_Eta;
    J(1,1) = -M_Gamma;

    return J;
}
MatrixSmall<2,2> IonicFitzHughNagumo::computeNumJ (const Real& t, const VectorSmall<2>& v, const Real& h)
{
	Real Iapp = 0.0;

	VectorSmall<2> y1(v);			y1(0) = y1(0) + h;
	VectorSmall<2> y2(v);			y2(0) = y2(0) - h;
	VectorSmall<2> f1,f2;

	this->computeRhs(y1, Iapp, f1);
	this->computeRhs(y2, Iapp, f2);

	VectorSmall<2> df1 = (f1-f2)/(2.0*h);

	y1 = v;			y1(1) = y1(1) + h;
	y2 = v;			y2(1) = y2(1) - h;

	this->computeRhs(y1, Iapp, f1);
	this->computeRhs(y2, Iapp, f2);

	VectorSmall<2> df2 = (f1-f2)/(2.0*h);

	MatrixSmall<2,2> J;
	J(0,0) = df1(0);	J(0,1) = df2(0);
	J(1,0) = df1(1);	J(1,1) = df2(1);

	return J;
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

