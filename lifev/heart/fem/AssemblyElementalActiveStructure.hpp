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
    @brief Local assembly of the jacobian terms, active neo-Hookean elasticity

    @contributor Ricardo Ruiz <ricardo.ruiz@epfl.ch>
 */


#ifndef _ELEMOPERSTRUCTURE_H_INCLUDED
#define _ELEMOPERSTRUCTURE_H_INCLUDED

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/CurrentBoundaryFE.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/DOF.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{
  //! @name Public typedefs
  //@{
  typedef boost::numeric::ublas::matrix<Real> Matrix;
  typedef boost::numeric::ublas::vector<Real> Vector;
  typedef boost::numeric::ublas::zero_matrix<Real> ZeroMatrix;
  //@}

  /*! /namespace AssemblyElementalActiveStructure

   */
  namespace AssemblyElementalActiveStructure
  {

    //! METHODS FOR ACTIVE NEO-HOOKEAN MODEL
    /*!
      @param Gammaf The activation scalar field
      @param I_1E Activated first invariant I_1E=(1+gammaf)I_1+g(Gammaf)I_4f
      @param coef The coefficient of the vector
      @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
      @param Fk The deformation gradient that depends on the local displacement uk_loc
      @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
      @param I_4f The fourth invariant of the right Cauchy-Green tensor C
      @param elvec The elementary vector of the current volume
      @param fe The current finite element
    */
    //! The following 6 methods are exactly as in the non-active case, but multiplied by (1+Gammaf)
    void source_P1iso_NH_Act(Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Gammaf, const KN<Real> Jk, const KN<Real> Ic_isok, VectorElemental& elvec, const CurrentFE& fe);

    void stiff_Jac_P1iso_NH_1term_Act( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Gammaf, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    void stiff_Jac_P1iso_NH_2term_Act( Real coef, const KNMK<Real> CofFk, const KN<Real> Gammaf, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

    void stiff_Jac_P1iso_NH_3term_Act( Real coef, const KN<Real> Gammaf, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    void stiff_Jac_P1iso_NH_4term_Act( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Gammaf, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    void stiff_Jac_P1iso_NH_5term_Act( Real coef, const KNMK<Real> CofFk, const KN<Real> Gammaf, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

    //!The following 9 methods contain the pure active components and the anisotropic activation

    //! 6. Jacobian matrix S1 : Int { coef * g(Gammaf) * J^(-2/3) (\nabla \delta [fo \tomes fo] : \nabla \v)}
    void stiff_Jac_NH_S1term_Act( Real coef, const KN<Real> Gammaf, const KNM<Real> fo, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    //! 7. Jacobian matrix S2 : Int { -2/3* coef * g(Gammaf) * J^(-5/3) *( CofF : \nabla \delta ) ( F [fo\otimes fo]: \nabla \v ) }
    void stiff_Jac_NH_S2term_Act( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Gammaf, const KNM<Real> fo, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    //! 8. Jacobian matrix S3: Int { 1/3 * coef * J^(-2-2/3) * I_4f * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
    void stiff_Jac_NH_S3term_Act( Real coef, const KNMK<Real> CofFk,  const KN<Real> Jk, const KN<Real> Gammaf, const KN<Real> I_4f, MatrixElemental& elmat, const CurrentFE& fe);

    //! 9. Jacobian matrix S41 : Int { -2/3* coef * g(Gammaf) * J^(-5/3) *( F [fo\otimes fo] : \nabla \delta ) ( CofF: \nabla \v ) }
    void stiff_Jac_NH_S41term_Act( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real>   Gammaf, const KNM<Real> fo, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    //! 10. Jacobian matrix S42: Int { 2/9 * coef * g(Gammaf) * I_4f * J^{-8/3} ( CofF : \nabla \delta ) ( CofF : \nabla \v ) }
    void stiff_Jac_NH_S42term_Act( Real coef, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Gammaf, const KN<Real> I_4f, MatrixElemental& elmat, const CurrentFE& fe);

    //! 11. Jacobian matrix S5: Int { 1/3 * coef * J^(-2-2/3) * I_1E * gammaf * (CofF [\nabla \delta]^t CofF [fo \otimes fo]) : \nabla \v }
    void stiff_Jac_NH_S5term_Act( Real coef, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Gammaf, const KNM<Real> fo, const KN<Real> I_1E, MatrixElemental& elmat, const CurrentFE& fe );

    //! 12. Jacobian matrix S61 : Int { 2/9* coef * gammaf * I_1E * J^(-8/3) *( CofF : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
    void stiff_Jac_NH_S61term_Act( Real coef, const KNMK<Real> CofFk, const KN<Real> Gammaf, const KNM<Real> fo, const KN<Real> I_1E, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    //! 13. Jacobian matrix S62 : Int { -2/3* coef * gammaf *(1+gammaf) * J^(-5/3) *( F : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
    void stiff_Jac_NH_S62term_Act( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real>   Gammaf, const KNM<Real> fo, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

    //! 14. Jacobian matrix S63 : Int { -2/3* coef * gammaf *(1+gammaf) * J^(-5/3) *( F : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
    void stiff_Jac_NH_S63term_Act( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Gammaf, const KNM<Real> fo, const KN<Real> Jk , MatrixElemental& elmat, const CurrentFE& fe );
  }
  //! End namespace AssemblyElementalActiveStructure

} //! End namespace LifeV
#endif
