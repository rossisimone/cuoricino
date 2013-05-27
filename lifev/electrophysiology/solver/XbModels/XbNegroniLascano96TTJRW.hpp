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

  @contributors Simone Rossi <simone.rossi@epfl.ch>
  @mantainer Luis Miguel de Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
  @last update 01-2013
 */


#ifndef _XBNEGRONILASCANO96TTJRW_H_
#define _XBNEGRONILASCANO96TTJRW_H_

#include <lifev/electrophysiology/solver/XbModels/ElectroXbModel.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <cmath>
#include <string>

namespace LifeV
{
//! XbModel - This class implements a mean field model.


class XbNegroniLascano96TTJRW : public virtual ElectroXbModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroXbModel super;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    XbNegroniLascano96TTJRW();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    XbNegroniLascano96TTJRW ( Teuchos::ParameterList& parameterList );

    /*!
     * @param XbNegroniLascano96TTJRW object
     */
    XbNegroniLascano96TTJRW ( const XbNegroniLascano96TTJRW& Xb );

    //! Destructor
    virtual ~XbNegroniLascano96TTJRW() {}

    //@}

    //! @name Overloads
    //@{

    XbNegroniLascano96TTJRW& operator= ( const XbNegroniLascano96TTJRW& Xb );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& Y1() const
    {
        return M_Y1;
    }
    inline const Real& Y2() const
    {
        return M_Y2;
    }
    inline const Real& Y3() const
    {
        return M_Y3;
    }
    inline const Real& Y4() const
    {
        return M_Y4;
    }
    inline const Real& Y5() const
    {
        return M_Yd;
    }
    inline const Real& Z1() const
    {
        return M_Z1;
    }
    inline const Real& Z2() const
    {
        return M_Z2;
    }
    inline const Real& Z3() const
    {
        return M_Z3;
    }
    inline const Real& Tt() const
    {
        return M_Tt;
    }
    inline const Real& B()  const
    {
        return M_B;
    }
    inline const Real& Hc() const
    {
        return M_Hc;
    }
    inline const Real& La() const
    {
        return M_La;
    }
    inline const Real& R()  const
    {
        return M_R;
    }
    inline const std::string& contrType() const
    {
        return M_contrType;
    }
    inline const Real& length() const
    {
       return M_L;
    }
    inline const Real& stepL() const
    {
       return M_dL;
    }
    inline const Real& tChangeL() const
    {
      return M_tChangeL;
    }
    inline const Real& dtChangeL() const
    {
      return M_dtChangeL;
    }


    inline void setY1   ( const Real& y1 )
    {
        this->M_Y1 = y1;
    }
    inline void setY2   ( const Real& y2 )
    {
        this->M_Y2 = y2;
    }
    inline void setY3   ( const Real& y3 )
    {
        this->M_Y3 = y3;
    }
    inline void setY4   ( const Real& y4 )
    {
        this->M_Y4 = y4;
    }
    inline void setY5   ( const Real& y5 )
    {
        this->M_Yd = y5;
    }
    inline void setZ1    ( const Real& z1 )
    {
        this->M_Z1  =  z1;
    }
    inline void setZ2    ( const Real& z2 )
    {
        this->M_Z2  =  z2;
    }
    inline void setZ3    ( const Real& z3 )
    {
        this->M_Z3  =  z3;
    }
    inline void setTt    ( const Real& tt )
    {
        this->M_Tt  =  tt;
    }
    inline void setB    ( const Real& b )
    {
        this->M_B  =  b;
    }
    inline void setHc    ( const Real& hc )
    {
        this->M_Hc  =  hc;
    }
    inline void setLa    ( const Real& la )
    {
        this->M_La  =  la;
    }
    inline void setR    ( const Real& r )
    {
        this->M_R  =  r;
    }
    inline void setContrType   ( const std::string& contrType )
    {
        this->M_contrType  =  contrType;
    }
    inline void setLength   ( const Real& length )
    {
        this->M_L  =  length;
    }
    inline void setdLength   ( const Real& dLength )
    {
        this->M_dL  =  dLength;
    }
    inline void setTChangeL    ( const Real& tChangeL )
    {
        this->M_tChangeL  =  tChangeL;
    }
    inline void setDtChangeL    ( const Real& dtChangeL )
    {
        this->M_dtChangeL  =  dtChangeL;
    }

    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeRhs ( const std::vector<Real>& v, const Real& Ca, const Real& L, const Real& vel, std::vector<Real>& rhs);

    //Compute the rhs on a mesh/ 3D case
//    void computeRhs ( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Ca, const Real& t, const VectorEpetra& vel, std::vector<vectorPtr_Type>& rhs );

    void computeBackwardEuler( std::vector<Real>&  v, const Real& Ca, const Real& vel, const Real dt );

    void computeVelocity( const Real& X, const Real& L, Real& vel );

    void computeX( const Real& dt, const Real& L, Real& X );

    void computeHalfSarcomereLength( const Real& t, Real& L );

    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_Y1;
    Real M_Y2;
    Real M_Y3;
    Real M_Y4;
    Real M_Yd;
    Real M_Z1;
    Real M_Z2;
    Real M_Z3;
    Real M_Tt;
    Real M_B;
    Real M_Hc;
    Real M_La;
    Real M_R;

    //! Choice of the contraction type and parameters
    std::string M_contrType;
    Real M_L;
    Real M_dL;
    Real M_tChangeL;
    Real M_dtChangeL;

    //! Xb states == equivalent to the number of equations
    //short int M_numberOfStates;

    //@}

}; // class XbNegroniLascano96TTJRW

// ===================================================
//! Constructors
// ===================================================
XbNegroniLascano96TTJRW::XbNegroniLascano96TTJRW()    :
    super       (  3  ),
    M_Y1        ( 39.0e-3 ),
    M_Y2        ( 1.3e-3 ),
    M_Y3        ( 30.0e-3 ),
    M_Y4        ( 40.0e-3 ),
    M_Yd        ( 9.0e3 ),
    M_Z1        ( 30.0e-3 ),
    M_Z2        ( 1.3e-3 ),
    M_Z3        ( 1560.0e-2 ),
    M_Tt        ( 70.0 ),
    M_B         ( 1200.0e-3 ),
    M_Hc        ( 0.005 ),
    M_La        ( 1.17 ),
    M_R         ( 20.0 ),
    M_contrType ("isometric"),
    M_L         ( 1.05 ),
    M_dL        ( -0.006 ),
    M_tChangeL  ( 3000.0 ),
    M_dtChangeL ( 0.2 )
{
}

XbNegroniLascano96TTJRW::XbNegroniLascano96TTJRW ( Teuchos::ParameterList& parameterList )   :
    super ( 3 )
{
    M_Y1        = parameterList.get ("y1", 39.0e-3);
    M_Y2        = parameterList.get ("y2", 1.3e-3);
    M_Y3        = parameterList.get ("y3", 30.0e-3);
    M_Y4        = parameterList.get ("y4", 40.0e-3);
    M_Yd        = parameterList.get ("yd", 9.0e3);
    M_Z1        = parameterList.get ("z1", 30.0e-3);
    M_Z2        = parameterList.get ("z2", 1.3e-2);
    M_Z3        = parameterList.get ("z3", 1560.0e-3);
    M_Tt        = parameterList.get ("tt", 70.0);
    M_B         = parameterList.get ("b", 1200.0e-3);
    M_Hc        = parameterList.get ("hc", 0.005);
    M_La        = parameterList.get ("la", 1.17);
    M_R         = parameterList.get ("r", 20.0 );
    M_contrType = parameterList.get ("contrType", "isometric");
    M_L         = parameterList.get ("inLength", 1.05);
    M_dL        = parameterList.get ("stepL", -0.006);
    M_tChangeL  = parameterList.get ("tChangeL", 3000.0);
    M_dtChangeL = parameterList.get ("dtChangeL", 0.2);
}

XbNegroniLascano96TTJRW::XbNegroniLascano96TTJRW ( const XbNegroniLascano96TTJRW& Xb )
{
    M_Y1        = Xb.M_Y1;
    M_Y2        = Xb.M_Y2;
    M_Y3        = Xb.M_Y3;
    M_Y4        = Xb.M_Y4;
    M_Yd        = Xb.M_Yd;
    M_Z1        = Xb.M_Z1;
    M_Z2        = Xb.M_Z2;
    M_Z3        = Xb.M_Z3;
    M_Tt        = Xb.M_Tt;
    M_B         = Xb.M_B;
    M_Hc        = Xb.M_Hc;
    M_La        = Xb.M_La;
    M_R         = Xb.M_R;
    M_contrType = Xb.M_contrType;
    M_L         = Xb.M_L;
    M_dL        = Xb.M_dL;
    M_tChangeL  = Xb.M_tChangeL;
    M_dtChangeL = Xb.M_dtChangeL;

    M_numberOfEquations = Xb.M_numberOfEquations;
}

// ===================================================
//! Operator
// ===================================================
XbNegroniLascano96TTJRW& XbNegroniLascano96TTJRW::operator= ( const XbNegroniLascano96TTJRW& Xb )
{
	 M_Y1        = Xb.M_Y1;
	 M_Y2        = Xb.M_Y2;
	 M_Y3        = Xb.M_Y3;
	 M_Y4        = Xb.M_Y4;
	 M_Yd        = Xb.M_Yd;
	 M_Z1        = Xb.M_Z1;
	 M_Z2        = Xb.M_Z2;
	 M_Z3        = Xb.M_Z3;
	 M_Tt        = Xb.M_Tt;
	 M_B         = Xb.M_B;
	 M_Hc        = Xb.M_Hc;
	 M_La        = Xb.M_La;
	 M_R         = Xb.M_R;
	 M_contrType = Xb.M_contrType;
	 M_L         = Xb.M_L;
	 M_dL        = Xb.M_dL;
	 M_tChangeL  = Xb.M_tChangeL;
	 M_dtChangeL = Xb.M_dtChangeL;

     M_numberOfEquations = Xb.M_numberOfEquations;

     return *this;
}


// ===================================================
//! Methods
// ===================================================
// v(0) = TCa
// v(1) = TCas
// v(2) = Ts
void XbNegroniLascano96TTJRW::computeRhs ( const std::vector<Real>&  v,
			                               const Real& Ca,
			                               const Real& L,
                                           const Real& vel,
                                           std::vector<Real>& rhs )
{

	Real TCaEff = v[0] * std::exp( - M_R * ( L- M_La ) * ( L - M_La ) );

    Real Qb  = M_Y1 * Ca * ( M_Tt - v[0] - v[1] - v[2] ) - M_Z1 * v[0];
    Real Qa  = M_Y2 * TCaEff - M_Z2 * v[1];
    Real Qr  = M_Y3 * v[1] - M_Z3 * Ca * v[2];
    Real Qd  = M_Y4 * v[2];
    Real Qd1 = M_Yd * vel * vel * v[2];
    Real Qd2 = M_Yd * vel * vel * v[1];


    rhs[0] = Qb - Qa;
    rhs[1] = Qa - Qr - Qd2;
    rhs[2] = Qr - Qd - Qd1;

}

void XbNegroniLascano96TTJRW::computeBackwardEuler( std::vector<Real>&  v,
                                               const Real& Ca,
                                               const Real& vel,
                                               const Real dt )
{
	Real TCaEff = v[0] * std::exp( - M_R * ( M_L - M_La ) * ( M_L - M_La ) );

	v[0] = ( v[0] / dt + M_Y1 * Ca * ( M_Tt - v[1] - v[2] ) + M_Z2 * v[1] )
				/ ( 1.0 / dt + M_Y1 * Ca + M_Z1 + M_Y2 * std::exp( - M_R * ( M_L - M_La ) * ( M_L - M_La ) ) );
    v[1] = ( v[1] / dt + M_Y2 * TCaEff + M_Z3 * v[2] * Ca )
				/ ( 1.0 / dt + M_Z2 + M_Y3 + M_Yd * vel * vel );
	v[2] = ( v[2] / dt + M_Y3 * v[1] )
				/ ( 1.0 / dt + M_Z3 * Ca + M_Y4 + M_Yd * vel * vel );
}

void XbNegroniLascano96TTJRW::computeVelocity( const Real& X, const Real& L, Real& vel )
{
    vel = M_B * ( L - X - M_Hc );
}

void XbNegroniLascano96TTJRW::computeX( const Real& dt, const Real&L, Real& X )
{
	X = ( X / dt + M_B * ( L - M_Hc ) ) / ( 1.0 / dt + M_B );
}


void XbNegroniLascano96TTJRW::computeHalfSarcomereLength( const Real& t, Real& L )
{
	if ( M_contrType == "noIsometric")
	{
		if ( t >= M_tChangeL && t <= M_tChangeL + M_dtChangeL )
		  L = L + M_dL;
	}
	else if ( M_contrType == "Stiffness" )
	{
		  L = 0.98 + M_dL * std::sin( 2 * 3.14 * t * 5e-3 );
	}
	else
		L = L;
}

//void XbNegroniLascano96TTJRW::computeRhs (    const std::vector<vectorPtr_Type>& v,
//                                         const VectorEpetra& Ca,
//                                         const VectorEpetra& vel,
//                                         std::vector<vectorPtr_Type>& rhs )
//{
//
//    int nodes = Ca.epetraVector().MyLength();
//
//
//    std::vector<Real>   localVec ( 3, 0.0 );
//    std::vector<Real>   localRhs ( 3, 0.0 );
//
//    int j (0);
//
//   for ( int k = 0; k < nodes; k++ )
//    {
//
//        j = Ca.blockMap().GID (k);
//
//        localVec.at (0) = ( * ( v.at (0) ) ) [j];
//        localVec.at (1) = ( * ( v.at (1) ) ) [j];
//        localVec.at (2) = ( * ( v.at (2) ) ) [j];
//
//        computeRhs ( localVec, Ca[j], vel[j], localRhs );
//
//        ( * ( rhs.at (0) ) ) [j] =  localRhs.at (0);
//        ( * ( rhs.at (1) ) ) [j] =  localRhs.at (1);
//        ( * ( rhs.at (2) ) ) [j] =  localRhs.at (2);
//
//    }
//
//}


void XbNegroniLascano96TTJRW::showMe()
{
    std::cout << "\n\n\t\tXbNegroniLascano96TTJRW Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "y1: " << this->Y1() << std::endl;
    std::cout << "y2: " << this->Y2() << std::endl;
    std::cout << "y3: " << this->Y3() << std::endl;
    std::cout << "y4: " << this->Y4() << std::endl;
    std::cout << "yd: " << this->Y5() << std::endl;
    std::cout << "z1: " << this->Z1() << std::endl;
    std::cout << "z2: " << this->Z2() << std::endl;
    std::cout << "z3: " << this->Z3() << std::endl;
    std::cout << "tt: " << this->Tt() << std::endl;
    std::cout << "b: " << this->B() << std::endl;
    std::cout << "hc: " << this->Hc() << std::endl;
    std::cout << "la: " << this->La() << std::endl;
    std::cout << "contrType: "  << this->contrType() << std::endl;
    std::cout << "length: "  << this->length() << std::endl;
    std::cout << "stepLength: "  << this->stepL() << std::endl;
    std::cout << "tChangeL: "  << this->tChangeL() << std::endl;
    std::cout << "dtChangeL: "  << this->dtChangeL() << std::endl;
    std::cout << "\n\t\t End of XbNegroniLascano96TTJRW Informations\n\n\n";

}


}

#endif
