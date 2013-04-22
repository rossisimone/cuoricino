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


#ifndef _XBNEGRONILASCANO96_H_
#define _XBNEGRONILASCANO96_H_

#include <lifev/heart/solver/XbModels/HeartXbModel.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <cmath>

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
    inline const Real& Y1() const
    {
        return M_Y1;
    }
    inline const Real& Y2() const
    {
        return M_Y2;
    }
    inline const Real& Y3()     const
    {
        return M_Y3;
    }
    inline const Real& Y4()     const
    {
        return M_Y4;
    }
    inline const Real& Y5()     const
    {
        return M_Yd;
    }
    inline const Real& Z1()      const
    {
        return M_Z1;
    }
    inline const Real& Z2()      const
    {
        return M_Z2;
    }
    inline const Real& Z3()      const
    {
        return M_Z3;
    }
    inline const Real& Tt()      const
    {
        return M_Tt;
    }
    inline const Real& B()      const
    {
        return M_B;
    }
    inline const Real& Hc()      const
    {
        return M_Hc;
    }
    inline const Real& La()      const
    {
        return M_La;
    }
    inline const Real& R()      const
    {
        return M_R;
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

    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeRhs ( const std::vector<Real>& v, const Real& Ca, const Real& vel,  std::vector<Real>& rhs);

    //Compute the rhs on a mesh/ 3D case
    void computeRhs ( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Ca, const VectorEpetra& vel, std::vector<vectorPtr_Type>& rhs );

    std::vector<Real> computeBackwardEuler( const std::vector<Real>&  v, const Real& Ca, const Real& vel, const Real dt );

    void computeVelocity( const Real& dt, Real& X, Real& vel );

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


    //! Xb states == equivalent to the number of equations
    //short int M_numberOfStates;

    //@}

}; // class XbNegroniLascano96

// ===================================================
//! Constructors
// ===================================================
XbNegroniLascano96::XbNegroniLascano96()    :
    super (  3  ),
    M_Y1  ( 39.0 ),
    M_Y2  ( 1.3 ),
    M_Y3  ( 30.0 ),
    M_Y4  ( 40.0 ),
    M_Yd  ( 9.0 ),
    M_Z1  ( 30.0 ),
    M_Z2  ( 1.3 ),
    M_Z3  ( 1560.0 ),
    M_Tt  ( 70.0 ),
    M_B   ( 1200.0 ),
    M_Hc  ( 0.005 ),
    M_La  ( 1.17 ),
    M_R   ( 20.0 )
{
}

XbNegroniLascano96::XbNegroniLascano96 ( Teuchos::ParameterList& parameterList )   :
    super ( 3 )
{
    M_Y1 = parameterList.get ("y1", 39.0);
    M_Y2 = parameterList.get ("y2", 1.3);
    M_Y3 = parameterList.get ("y3", 30.0);
    M_Y4 = parameterList.get ("y4", 40.0);
    M_Yd = parameterList.get ("yd", 9.0);
    M_Z1 = parameterList.get ("z1", 30.0);
    M_Z2 = parameterList.get ("z2", 1.3);
    M_Z3 = parameterList.get ("z3", 1560.0);
    M_Tt = parameterList.get ("tt", 70.0);
    M_B  = parameterList.get ("b", 70.0);
    M_Hc = parameterList.get ("hc", 70.0);
    M_La = parameterList.get ("la", 1.17);
    M_R  = parameterList.get ("r", 20.0 );

}

XbNegroniLascano96::XbNegroniLascano96 ( const XbNegroniLascano96& Xb )
{
    M_Y1 = Xb.M_Y1;
    M_Y2 = Xb.M_Y2;
    M_Y3 = Xb.M_Y3;
    M_Y4 = Xb.M_Y4;
    M_Yd = Xb.M_Yd;
    M_Z1 = Xb.M_Z1;
    M_Z2 = Xb.M_Z2;
    M_Z3 = Xb.M_Z3;
    M_Tt = Xb.M_Tt;
    M_B  = Xb.M_B;
    M_Hc = Xb.M_Hc;
    M_La = Xb.M_La;
    M_R  = Xb.M_R;

    M_numberOfEquations = Xb.M_numberOfEquations;
}

// ===================================================
//! Operator
// ===================================================
XbNegroniLascano96& XbNegroniLascano96::operator= ( const XbNegroniLascano96& Xb )
{
    M_Y1 = Xb.M_Y1;
    M_Y2 = Xb.M_Y2;
    M_Y3 = Xb.M_Y3;
    M_Y4 = Xb.M_Y4;
    M_Yd = Xb.M_Yd;
    M_Z1 = Xb.M_Z1;
    M_Z2 = Xb.M_Z2;
    M_Z3 = Xb.M_Z3;
    M_Tt = Xb.M_Tt;
    M_B  = Xb.M_B;
    M_Hc = Xb.M_Hc;
    M_La = Xb.M_La;
    M_R  = Xb.M_R;

    M_numberOfEquations = Xb.M_numberOfEquations;

    return *this;
}


// ===================================================
//! Methods
// ===================================================
// v(0) = TCa
// v(1) = TCas
// v(2) = Ts
void XbNegroniLascano96::computeRhs (    const   std::vector<Real>&  v,
                                         const   Real& Ca,
                                         const   Real& vel,
                                         std::vector<Real>& rhs )
{
	Real length ( 1.05 );
    Real TCaEff = v[0] * std::exp( - M_R * ( length - M_La ) * ( length - M_La ) );

    Real Qb  = M_Y1 * Ca * ( M_Tt - v[0] - v[1] - v[2] ) - M_Z1 * v[0];
    Real Qa  = M_Y2 * TCaEff - M_Z2 * v[1];
    Real Qr  = M_Y3 * v[1] - M_Z3 * Ca * v[2];
    Real Qd  = M_Y4 * v[3];
    Real Qd1 = M_Yd * vel * vel * v[2];
    Real Qd2 = M_Yd * vel * vel * v[1];



    rhs[0] = Qb - Qa;
    rhs[1] = Qa - Qr - Qd2;
    rhs[2] = Qr - Qd - Qd1;


}

std::vector<Real> XbNegroniLascano96::computeBackwardEuler( const std::vector<Real>&  v,
                                                            const Real& Ca,
                                                            const Real& vel,
                                                            const Real dt )
{
    std::vector<Real> BErhs(3);

	Real length ( 1.05 );
    Real TCaEff = v[0] * std::exp( - M_R * ( length - M_La ) * ( length - M_La ) );

	BErhs[0] = ( v[0] / dt + M_Y1 * Ca * ( M_Tt - v[1] - v[2] ) - M_Z2 * v[1] )
				/ ( 1.0 / dt + M_Y1 * Ca + M_Z1 + M_Y2 * std::exp( - M_R * ( length - M_La ) * ( length - M_La ) ) );
    BErhs[1] = ( v[1] / dt + M_Y2 * TCaEff - M_Z3 * v[2] * Ca )
				/ ( 1.0 / dt + M_Z3 + M_Y3 + M_Yd * vel * vel );
	BErhs[2] = ( v[2] / dt + M_Y3 * v[1] )
				/ ( 1.0 / dt + M_Z3 * Ca + M_Y4 + M_Yd * vel * vel );

	return BErhs;
}

void XbNegroniLascano96::computeVelocity( const Real& dt, Real& X, Real& vel )
{
    Real length ( 1.05 );

    X   = ( X / dt + M_B * ( length - M_Hc ) ) / ( 1.0 / dt + M_B );
    vel = M_B * ( length - X - M_Hc );
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
    std::cout << "y1: " << this->Y1() << std::endl;
    std::cout << "y2: " << this->Y2() << std::endl;
    std::cout << "y3: " << this->Y3() << std::endl;
    std::cout << "y4: " << this->Y4() << std::endl;
    std::cout << "yd: " << this->Y5() << std::endl;
    std::cout << "z1: " << this->Z1() << std::endl;
    std::cout << "z2: " << this->Z2() << std::endl;
    std::cout << "z3: " << this->Z3() << std::endl;
    std::cout << "tt: " << this->Tt() << std::endl;
    std::cout << "la: " << this->La() << std::endl;
    std::cout << "r: "  << this->R() << std::endl;
    std::cout << "\n\t\t End of XbNegroniLascano96 Informations\n\n\n";

}


}

#endif
