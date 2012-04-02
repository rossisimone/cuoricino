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
    @brief  A class for static boundary finite element

    @author Jean-Frederic Gerbeau
            Vincent Martin
    @date 00-09-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef _StaticBDFE_H
#define _StaticBDFE_H

#include <vector>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/numeric/ublas/vector.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/array/RNM.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{
/*!
  \class CurrentBoundaryFEBase
  \brief A class for static boundary finite element
  \author J.-F. Gerbeau & V. Martin
  \date 09/2002

  This class has two purposes:
  \par
  (1) it is a base class for standard boundary element (see CurrentBoundaryFE.h)
  \par
  (2) it is used by refHybridFE as static boundary for a reference element

*/

class CurrentBoundaryFEBase
{

public:
    //! @name Public typedefs
    //@{
    typedef boost::numeric::ublas::vector<Real> Vector;
    typedef boost::numeric::ublas::zero_vector<Real> ZeroVector;
    //@}

    //! @name Constructor & Destructor
    //@{

    CurrentBoundaryFEBase( const ReferenceFE& refFE, const GeometricMap& geoMap );

    CurrentBoundaryFEBase( const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr );

    //! new optionnal argument for hybrid hdiv fe invArea
    CurrentBoundaryFEBase( const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr,
                const Real* refcoor, UInt currentid, Real invarea = 1. );

    virtual ~CurrentBoundaryFEBase();

    //@}


    //! @name Methods
    //@{

    /*!
      return the coordinate (x,y,z)= F(xi,eta), (F: geo mappping)
      where (xi,eta) are the coor in the ref element
      and   (x,y,z) the coor in the current element
      (if the code is compiled in 2D mode then z=0 and eta is disgarded)
    */
    void coorMap( Real& x, Real& y, Real& z,
                  const Real & xi, const Real & eta ) const;
    /*!
      return (x,y,z) = the global coordinates of the quadrature point ig
      in the current element. \warning this function is almost obsolete since if
      you call the function updateMeasQuadPt rather than updateMeas
      (for example), the coordinates of the quadrature points have already
      been computed and may be obtained via quadPt(ig,icoor). This is usually
      much less expensive since it avoids many calls to coorQuadPt
    */
    void coorQuadPt( Real& x, Real& y, Real& z, const int ig ) const
    {
        ASSERT_PRE( M_hasQR, "Needs a quadrature rule" )
        coorMap( x, y, z, qr.quadPointCoor( ig, 0 ), qr.quadPointCoor( ig, 1 ) );
    }
    /*!
      Return the measure of the current element
      \warning either updateMeas(...) or updateMeasNormal(...)
      must have been called before
    */
    Real measure() const;
    /*!
      compute the integral over the current boundary element
      \warning either updateMeas(...) or updateMeasNormal(...)
      must have been called before
    */
    template <typename FunctorType>
    Real integral( const FunctorType & f ) const;

    /*!
      compute the integral of f . n over the current boundary element
      \warning  updateMeasNormal(...) must have been called before
    */
    template <typename FunctorType>
    Real integral_n( const FunctorType & f ) const;

    //@}


    //! @name Get Methods
    //@{

#ifdef TEST_PRE
    /*!
      return true if a quadrature rule has been given
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    const bool& hasQR() const
    {
        return M_hasQR;
    }
    /*!
      return true if the measure has been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    const bool& hasMeas() const
    {
        return M_hasMeasure;
    }
    /*!
      return true if the first derivatives have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    const bool& hasFirstDeriv() const
    {
        return M_hasFirstDerivative;
    }
    /*!
      return true if the tangents have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    const bool& hasTangent() const
    {
        return M_hasTangent;
    }
    /*!
      return true if the normal has been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    const bool& hasNormal() const
    {
        return M_hasNormal;
    }
    /*!
      return true if the coordinate of the quadrature points have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    const bool& hasQuadPtCoor() const
    {
        return M_hasQuadPtCoor;
    }
#endif
    /*!
      return the id of the current element (updated with the update* functions)
    */
    const UInt& currentId() const
    {
        return M_currentID;
    }

    //@}


    //! Number of geometrical nodes
    const UInt& nbGeoNode() const {return M_nbGeoNode; }

    //! Number of finite element node
    const UInt& nbNode() const {return M_nbNode; }

    //! Number of coordinates
    const UInt& nbCoor() const {return M_nbCoor; }

    //! Number of quadrature points
    const UInt& nbQuadPt() const {return M_nbQuadPt; }

    //! The point that define the geometry
    const Real& point(const UInt& i, const UInt& coor) const
    {
        return M_point(int(i),int(coor));
    }

    //! Values of the basis functions on quadrature points
    const Real& phi(const UInt& i, const UInt& iQuadPt) const
    {
    	return M_phi(int(i),int(iQuadPt));
    }

    //! Values of the derivatives of the basis functions on quadrature points on the reference finite element
    const Real& dPhiRef(const UInt& i, const UInt& dxj, const UInt& iQuadPt) const
    {
        return M_dPhiRef(int(i),int(dxj),int(iQuadPt));
    }

    //! Values of the derivatives of the basis functions on quadrature points on the current finite element
    const Real& dPhi(const UInt& i, const UInt& dxj, const UInt& iQuadPt) const
    {
        return M_dPhi(int(i),int(dxj),int(iQuadPt));
    }

    //! Values of the geometric basis functions on quadrature points
    const Real& phiGeo(const UInt& i, const UInt& iQuadPt) const
    {
        return M_phiGeo(int(i),int(iQuadPt));

    }

    //! Values of the derivatives of the geometric basis functions on quadrature points
    const Real& dPhiGeo(const UInt& i, const UInt& dxj, const UInt& iQuadPt) const
    {
        return M_dPhiGeo(int(i),int(dxj),int(iQuadPt));
    }

    //! Values of the weight times the measure on the quadrature points
    const Real& weightMeas(const UInt& iQuadPt) const
    {
        return M_weightMeas(int(iQuadPt));
    }

    //! Values of the measures on the quadrature points
    const Real& meas(const UInt& iQuadPt) const
    {
        return M_meas(int(iQuadPt));
    }

    //! Values of the normal on the quadrature points
    const Real& normal(const UInt& coor, const UInt& iQuadPt) const
    {
        return M_normal(int(coor),int(iQuadPt));
    }

    //! Values of the tangents on the quadrature points
    const Real& tangent(const UInt& i, const UInt& coor, const UInt& iQuadPt) const
    {
        return M_tangent(int(i),int(coor),int(iQuadPt));
    }

    //! Metric tensor on the quadrature points
    const Real& metric(const UInt& iCoor, const UInt& jCoor, const UInt& iQuadPt) const
    {
        return M_metric(int(iCoor),int(jCoor),int(iQuadPt));
    }

    //! Coordinates of the quadrature points on the current element
    const Real& quadPt(const UInt& i, const UInt& coor) const
    {
        return M_quadPt(int(i),int(coor));
    }



protected:
    void computeMeasure();
    void computeMeasureNormal();
    void computeQuadPointCoordinate();

#ifdef TEST_PRE
    bool M_hasQR;
    bool M_hasMeasure;
    bool M_hasFirstDerivative;
    bool M_hasTangent;
    bool M_hasNormal;
    bool M_hasQuadPtCoor;
#endif

    const UInt M_nbGeoNode;
    const UInt M_nbNode;
    const UInt M_nbCoor;
    const UInt M_nbQuadPt;

    KNM<Real>  M_point;

public:
    //! The reference finite element
    const ReferenceFE& refFE;

    //! The geometical mapping
    const GeometricMap&   geoMap;

    //! The quadrature rule
    const QuadratureRule& qr;
protected:

    KNM<Real>  M_phi;
    KNMK<Real> M_dPhiRef;
    KNMK<Real> M_dPhi;
    KNM<Real>  M_phiGeo;
    KNMK<Real> M_dPhiGeo;
    KN<Real>   M_weightMeas;
    KN<Real>   M_meas;
    KNM<Real>  M_normal;
    KNMK<Real> M_tangent;
    KNMK<Real> M_metric;
    KNM<Real>  M_quadPt;

public:
    //! Inverse of the area
    const Real invArea;

protected:

    UInt M_currentID;

};


template <typename FunctorType>
Real
CurrentBoundaryFEBase::
integral( const FunctorType & f ) const
{
    ASSERT_PRE( M_hasMeasure, "integral needs measure. Call an update function" );
    Real integ( 0.0 );
    Real x, y, z;
    for ( UInt iQuadPt(0); iQuadPt < M_nbQuadPt; ++iQuadPt )
    {
        coorQuadPt( x, y, z, iQuadPt );
        integ += f( x, y, z ) * M_weightMeas( iQuadPt );
    }
    return integ;
}

template <typename FunctorType>
Real
CurrentBoundaryFEBase::
integral_n( const FunctorType & f ) const
{
    ASSERT_PRE( M_hasNormal, "integral_n needs measure and normal. Call the appropriate update function" );
    Real integ( 0.0 );
    Real x, y, z;
    std::vector<Real> ret( M_nbCoor + 1 );
    Real tmp;
    for ( UInt iQuadPt(0); iQuadPt < M_nbQuadPt; ++iQuadPt )
    {
        coorQuadPt( x, y, z, iQuadPt );
        f( x, y, z, &(ret[0]) );
        tmp = 0;
        for ( UInt d(0); d <= M_nbCoor; ++d )
        {
            tmp += ret[ d ] * normal( d, iQuadPt );
        }
        integ += tmp * M_weightMeas( iQuadPt );
    }
    return integ;
}


}
#endif
