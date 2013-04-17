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

#include <lifev/heart/solver/IonicModels/IonicAlievPanfilov.hpp>


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! XbModel - This class implements a mean field model.
// ===================================================
//! Constructors
// ===================================================
IonicAlievPanfilov::IonicAlievPanfilov()    :
    super       ( 2   ),
    M_mu1       ( 0.12 ),
    M_mu2       ( 0.3 ),
    M_k         ( 8.0 ),
    M_a         ( 0.1 ),
    M_epsilon   ( 0.01 )
{
}

IonicAlievPanfilov::IonicAlievPanfilov ( Teuchos::ParameterList& parameterList   )   :
    super       ( 2 )
{
    M_mu1       =  parameterList.get ("mu1", 0.12);
    M_mu2       =  parameterList.get ("mu2", 0.3);
    M_k         =  parameterList.get ("k", 8.0);
    M_a         =  parameterList.get ("a", 0.1);
    M_epsilon   =  parameterList.get ("epsilon", 0.01);
}

IonicAlievPanfilov::IonicAlievPanfilov ( const IonicAlievPanfilov& model )
{

    M_mu1       =  model.M_mu1;
    M_mu2       =  model.M_mu2;
    M_k         =  model.M_k;
    M_a         =  model.M_a;
    M_epsilon   =  model.M_epsilon;

    M_numberOfEquations = model.M_numberOfEquations;
}

// ===================================================
//! Operator
// ===================================================
IonicAlievPanfilov& IonicAlievPanfilov::operator= ( const IonicAlievPanfilov& model )
{
    M_mu1       =  model.M_mu1;
    M_mu2       =  model.M_mu2;
    M_k         =  model.M_k;
    M_a         =  model.M_a;
    M_epsilon   =  model.M_epsilon;

    M_numberOfEquations = model.M_numberOfEquations;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
//Only gating variables
void IonicAlievPanfilov::computeRhs (    const   std::vector<Real>&  v,
                                         std::vector<Real>& rhs )
{

    Real dr = - ( M_epsilon + M_mu1 * v[1] / ( M_mu2 + v[0] ) ) * ( v[1] + M_k * v[0] * ( v[0] - M_a  - 1.0 ) );

    rhs[0] = dr;

}

//Potential and gating variables
void IonicAlievPanfilov::computeRhs (    const   std::vector<Real>&  v,
                                         const   Real&           Iapp,
                                         std::vector<Real>& rhs )
{

    Real dr = - ( M_epsilon + M_mu1 * v[1] / ( M_mu2 + v[0] ) ) * ( v[1] + M_k * v[0] * ( v[0] - M_a  - 1.0 ) );
    Real dV = - M_k * v[0] * ( v[0] - M_a ) * ( v[0] - 1.0) - v[0] * v[1] + Iapp;

    rhs[0] = dV;
    rhs[1] = dr;

}

//Only gating variables
//void IonicAlievPanfilov::computeRhs(  const std::vector<vectorPtr_Type>& v,
//                                          std::vector<vectorPtr_Type>& rhs )
//{
//
//  int nodes = ( *(v.at(1) ) ).epetraVector().MyLength();
//
//
//  std::vector<Real>   localVec( 2, 0.0 );
//  std::vector<Real>   localRhs( 1, 0.0 );
//
//  int j(0);
//
//  for( int k = 0; k < nodes; k++ ){
//
//      j = ( *(v.at(1) ) ).blockMap().GID(k);
//
//      localVec.at(0) = ( *( v.at(0) ) )[j];
//      localVec.at(1) = ( *( v.at(1) ) )[j];
//
//      computeRhs( localVec, localRhs );
//
//      ( *( rhs.at(1) ) )[j] =  localRhs.at(0);
//
//  }
//
//}
//
////Potential and gating variables
//void IonicAlievPanfilov::computeRhs(  const std::vector<vectorPtr_Type>& v,
//                                          const VectorEpetra& Iapp,
//                                          std::vector<vectorPtr_Type>& rhs )
//{
//
//  int nodes = Iapp.epetraVector().MyLength();
//
//
//  std::vector<Real>   localVec( 2, 0.0 );
//  std::vector<Real>   localRhs( 2, 0.0 );
//
//  int j(0);
//
//  for( int k = 0; k < nodes; k++ ){
//
//      j = Iapp.blockMap().GID(k);
//
//      localVec.at(0) = ( *( v.at(0) ) )[j];
//      localVec.at(1) = ( *( v.at(1) ) )[j];
//
//      computeRhs( localVec, Iapp[j], localRhs );
//
//      ( *( rhs.at(0) ) )[j] =  localRhs.at(0);
//      ( *( rhs.at(1) ) )[j] =  localRhs.at(1);
//
//  }
//
//}



Real IonicAlievPanfilov::computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp)
{
    return ( - M_k * v[0] * ( v[0] - M_a ) * ( v[0] - 1.0) - v[0] * v[1] + Iapp );
}



//void IonicAlievPanfilov::computePotentialRhs(     const std::vector<vectorPtr_Type>& v,
//                                                  const VectorEpetra& Iapp,
//                                                  std::vector<vectorPtr_Type>& rhs,
//                                                  FESpace<mesh_Type, MapEpetra>& uFESpace )
//{
//
//  std::vector<Real> U(2, 0.0);
//  Real V(0.0);
//  Real r(0.0);
//  Real I(0.0);
//
//  VectorEpetra    VRep( *( v.at(0) )  , Repeated);
//  VectorEpetra    rRep( *( v.at(1) )  , Repeated);
//  VectorEpetra    IappRep( Iapp       , Repeated );
//
//  VectorElemental elvec_V( uFESpace.fe().nbFEDof(), 1 );
//  VectorElemental elvec_Iapp( uFESpace.fe().nbFEDof(), 1 );
//  VectorElemental elvec_r( uFESpace.fe().nbFEDof(), 1 );
//  VectorElemental elvec_Iion( uFESpace.fe().nbFEDof(), 1 );
//
//  for (UInt iVol=0; iVol< uFESpace.mesh()->numVolumes(); ++iVol){
//
//      uFESpace.fe().updateJacQuadPt( uFESpace.mesh()->volumeList( iVol ) );
//
//      elvec_V.zero();
//      elvec_Iapp.zero();
//      elvec_r.zero();
//      elvec_Iion.zero();
//
//      UInt eleIDu = uFESpace.fe().currentLocalId();
//      UInt nbNode = ( UInt ) uFESpace.fe().nbFEDof();
//
//      //! Filling local elvec_u with potential values in the nodes
//      for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ ){
//
//          Int  ig = uFESpace.dof().localToGlobalMap( eleIDu, iNode );
//          elvec_V.vec()[ iNode ] = VRep[ig];
//          elvec_Iapp.vec()[ iNode ] = IappRep[ig];
//          elvec_r.vec()[ iNode ] = rRep[ig];
//
//      }
//
//      //compute the local vector
//      for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt();ig++ ){
//
//          V = 0;
//          r = 0;
//          I = 0;
//
//          for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ ){
//
//              V += elvec_V(i) * uFESpace.fe().phi( i, ig );
//              r += elvec_r(i) * uFESpace.fe().phi( i, ig );
//              I += elvec_Iapp(i) * uFESpace.fe().phi( i, ig );
//
//          }
//
//          for ( UInt i = 0;i < uFESpace.fe().nbFEDof();i++ ){
//
//              elvec_Iion( i ) += ( - M_k * V * ( V - M_a ) * ( V - 1.0 ) - r * V + I) *
//                                  uFESpace.fe().phi( i, ig ) * uFESpace.fe().weightDet( ig );
//
//          }
//
//      }
//
//      //assembly
//      for ( UInt i = 0 ; i < uFESpace.fe().nbFEDof() ; i++ ) {
//          Int  ig = uFESpace.dof().localToGlobalMap( eleIDu, i );
//          ( *( rhs.at(0) ) ).sumIntoGlobalValues (ig,  elvec_Iion.vec()[i] );
//      }
//  }
//}

void IonicAlievPanfilov::showMe()
{
    std::cout << "\n\n\t\tIonicAlievPanfilov Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "mu1: " << this->Mu1() << std::endl;
    std::cout << "mu2: " << this->Mu2() << std::endl;
    std::cout << "k: " << this->K() << std::endl;
    std::cout << "a: " << this->A() << std::endl;
    std::cout << "epsilon: " << this->Epsilon() << std::endl;

    std::cout << "\n\t\t End of IonicAlievPanfilov Informations\n\n\n";

}


}

