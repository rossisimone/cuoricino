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
    @brief Local assembly of the active Neo-Hookean elasticity problem

    @author Ricardo Ruiz <ricardo.ruiz@epfl.ch>
 */

#ifndef ELEMOPERSTRUCTURE_CPP
#define ELEMOPERSTRUCTURE_CPP 1

#include <lifev/structure/fem/AssemblyElementalActiveStructure.hpp>

namespace LifeV
{

namespace AssemblyElementalActiveStructure
{



//! ***********************************************************************************************
//! METHODS for ACTIVE NEO-HOOKEAN
//! ***********************************************************************************************
//! Stiffness vector isochoric part ---------------------------------------------------------------

//! Gammaf is the activation parameter

//! Source term source_P1iso_NH: Int { coef * (1+Gammaf) *( J^(-2/3) * (F : \nabla v) - 1/3 * (Ic_iso / J) (CofF : \nabla v) ) }
void source_P1iso_NH_Act ( Real      coef,
                           const KNMK<Real> CofFk,
                           const KNMK<Real> Fk,
                           const KN<Real>   Gammaf,
                           const KN<Real>   Jk,
                           const KN<Real>   Ic_isok ,
                           VectorElemental& elvec,
                           const CurrentFE& fe )
{
    Real s1, s2;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        VectorElemental::vector_view vec =  elvec.block ( icoor );
        for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
        {
            s1 = 0.0;
            s2 = 0.0;
            for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
            {
                for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                {
                    s1 +=  (1.0 + Gammaf (ig) ) * pow ( Jk ( ig ), (-2.0 / 3.0) ) * Fk ( icoor,  k, ig ) *
                           fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                    s2 +=  1.0 / 3.0 * (1.0 + Gammaf (ig) ) * ( Ic_isok ( ig ) * ( 1 / Jk (ig) ) ) *
                           CofFk ( icoor, k, ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                }
            }
            vec ( i ) += (s1 - s2) * coef;
        }
    }
}
//! -----------------------------------------------------------------------------------------------






//! Jacobian matrix isochoric part ----------------------------------------------------------------

//! 1. Jacobian matrix : Int { -2/3 * coef * (1+Gammaf) * J^(-5/3) *( CofF : \nabla \delta ) ( F : \nabla \v ) }
void stiff_Jac_P1iso_NH_1term_Act ( Real coef,
                                    const KNMK<Real> CofFk,
                                    const KNMK<Real> Fk,
                                    const KN<Real>   Gammaf,
                                    const KN<Real> Jk ,
                                    MatrixElemental& elmat,
                                    const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                            {
                                s += (1.0 + Gammaf (ig) ) * pow ( Jk (ig), -5. / 3. ) *
                                     Fk ( jcoor, l, ig ) * fe.phiDer ( j, l, ig ) *
                                     CofFk ( icoor, k, ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



//! 2. Stiffness matrix: Int { 2/9 * coef * (1+Gammaf)* ( Ic_iso / J^2 )( CofF : \nabla \delta ) ( CofF : \nabla \v ) }
void stiff_Jac_P1iso_NH_2term_Act ( Real coef,
                                    const KNMK<Real> CofFk,
                                    const KN<Real>   Gammaf,
                                    const KN<Real> Jk ,
                                    const KN<Real> Ic_isok,
                                    MatrixElemental& elmat,
                                    const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                            {
                                s += (1.0 + Gammaf (ig) ) * ( 1 / (Jk (ig) * Jk (ig) ) ) *  Ic_isok (ig) *
                                     CofFk ( jcoor, l, ig ) * fe.phiDer ( j, l, ig ) *
                                     CofFk ( icoor, k, ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

//! 3. Stiffness matrix : Int { coef * (1+Gammaf) * J^(-2/3) (\nabla \delta : \nabla \v)}
void stiff_Jac_P1iso_NH_3term_Act ( Real         coef,
                                    const KN<Real>   Gammaf,
                                    const KN<Real>   Jk,
                                    MatrixElemental& elmat,
                                    const CurrentFE& fe )
{
    Real s;

    //! assembling diagonal block
    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
    {
        for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
        {
            s = 0.0;
            for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
            {
                for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                {
                    s += (1.0 + Gammaf (ig) ) * pow ( Jk (ig), -2. / 3.) * fe.phiDer ( i, k, ig ) *
                         fe.phiDer ( j, k, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        //! copy of diagonal block
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }
}



//! 4. Stiffness matrix : Int { -2/3 * coef * (1+Gammaf) * J^(-5/3) ( F : \nabla \delta ) ( CofF : \nabla \v ) }
void stiff_Jac_P1iso_NH_4term_Act ( Real coef,
                                    const KNMK<Real> CofFk,
                                    const KNMK<Real> Fk,
                                    const KN<Real>   Gammaf,
                                    const KN<Real> Jk ,
                                    MatrixElemental& elmat,
                                    const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                            {
                                s += (1.0 + Gammaf (ig) ) * pow ( Jk (ig), -5. / 3. ) *
                                     Fk ( icoor, k, ig )  * fe.phiDer ( i, k, ig ) *
                                     CofFk ( jcoor, l, ig ) * fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}



//! 5. Stiffness matrix : Int { 1/3 * coef * (1+Gammaf) * J^(-2) * Ic_iso * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
void stiff_Jac_P1iso_NH_5term_Act ( Real coef,
                                    const KNMK<Real> CofFk,
                                    const KN<Real>   Gammaf,
                                    const KN<Real> Jk ,
                                    const KN<Real> Ic_isok,
                                    MatrixElemental& elmat,
                                    const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                            {
                                s += (1.0 + Gammaf (ig) ) * ( 1 / ( Jk (ig) * Jk (ig) ) ) * Ic_isok (ig) *
                                     CofFk ( icoor , l , ig ) * fe.phiDer ( j, l, ig ) *
                                     CofFk ( jcoor , k , ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += s * coef;
                }
            }
        }
    }
}


//! -----------------------------------------------------------------------------------------------
//! New terms corresponding to the pure active part, Gammaf: activation, fo: fibers direction
//! -----------------------------------------------------------------------------------------------

//! 6. Jacobian matrix S1 : Int { coef * g(Gammaf) * J^(-2/3) (\nabla \delta [fo \tomes fo] : \nabla \v)}
void stiff_Jac_NH_S1term_Act ( Real       coef,
                               const KN<Real>   Gammaf,
                               const KNM<Real>   fo,
                               const KN<Real>   Jk,
                               MatrixElemental& elmat,
                               const CurrentFE& fe )
{

    Real s;

    //
    // blocks (icoor,jcoor) of elmat
    //
    MatrixElemental::matrix_type mat_tmp ( fe.nbFEDof(), fe.nbFEDof() );

    for ( UInt i = 0; i < fe.nbFEDof(); ++i )
    {
        for ( UInt j = 0; j < fe.nbFEDof(); ++j )
        {
            s = 0.0;
            for ( Int l = 0; l < static_cast<Int> (fe.nbCoor() ); ++l )
            {
                for ( Int k = 0; k < static_cast<Int> (fe.nbCoor() ); ++k )
                {
                    for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                        s += -Gammaf (ig) * (1.0 + (Gammaf (ig) + 2.0) * pow ( (1.0 + Gammaf (ig) ), -2.) ) * pow ( Jk (ig), -2. / 3. ) *
                             fo ( l , ig ) * fo ( k , ig ) * fe.phiDer ( i, k, ig ) *  fe.phiDer ( j, l, ig ) * fe.weightDet ( ig );
                }
            }
            mat_tmp ( i, j ) = coef * s;
        }
    }

    for ( UInt icoor = 0; icoor < fe.nbCoor(); ++icoor )
    {
        MatrixElemental::matrix_view mat = elmat.block ( icoor, icoor );
        mat += mat_tmp;
    }

}

//! 7. Jacobian matrix S2 : Int { -2/3* coef * g(Gammaf) * J^(-5/3) *( CofF : \nabla \delta ) ( F [fo\otimes fo]: \nabla \v ) }
void stiff_Jac_NH_S2term_Act ( Real coef,
                               const KNMK<Real> CofFk,
                               const KNMK<Real> Fk,
                               const KN<Real>   Gammaf,
                               const KNM<Real>   fo,
                               const KN<Real> Jk ,
                               MatrixElemental& elmat,
                               const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int m = 0; m < static_cast<Int> (nDimensions); ++m )
                            {
                                for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                                {
                                    s += -Gammaf (ig) * (1.0 + (Gammaf (ig) + 2.0) * pow ( (1.0 + Gammaf (ig) ), -2.) ) * pow ( Jk (ig), -5. / 3. ) *
                                         Fk ( jcoor, m, ig ) * fo (m, ig) * fo (l, ig) * fe.phiDer ( j, l, ig ) *
                                         CofFk ( icoor, k, ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                                }
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

//! 8. Jacobian matrix S3: Int { 1/3 * coef * J^(-2-2/3) * I_4f * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
void stiff_Jac_NH_S3term_Act ( Real coef,
                               const KNMK<Real> CofFk,
                               const KN<Real> Jk ,
                               const KN<Real>   Gammaf,
                               const KN<Real> I_4f,
                               MatrixElemental& elmat,
                               const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                            {
                                s += -Gammaf (ig) * (1.0 + (Gammaf (ig) + 2.0) * pow ( (1.0 + Gammaf (ig) ), -2.) ) * pow ( Jk (ig), -8. / 3. ) * I_4f (ig) *
                                     CofFk ( icoor , l , ig ) * fe.phiDer ( j, l, ig ) *
                                     CofFk ( jcoor , k , ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += s * coef;
                }
            }
        }
    }
}

//! 9. Jacobian matrix S41 : Int { -2/3* coef * g(Gammaf) * J^(-5/3) *( F [fo\otimes fo] : \nabla \delta ) ( CofF: \nabla \v ) }
void stiff_Jac_NH_S41term_Act ( Real coef,
                                const KNMK<Real> CofFk,
                                const KNMK<Real> Fk,
                                const KN<Real>   Gammaf,
                                const KNM<Real>   fo,
                                const KN<Real> Jk ,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int m = 0; m < static_cast<Int> (nDimensions); ++m )
                            {
                                for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                                {
                                    s += -Gammaf (ig) * (1.0 + (Gammaf (ig) + 2.0) * pow ( (1.0 + Gammaf (ig) ), -2.) ) * pow ( Jk (ig), -5. / 3. ) *
                                         CofFk ( jcoor, l, ig ) * fe.phiDer ( j, l, ig ) *
                                         Fk ( icoor, m, ig ) * fo (m, ig) * fo (k, ig) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                                }
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}


//! 10. Jacobian matrix S42: Int { 2/9 * coef * g(Gammaf) * I_4f * J^{-8/3} ( CofF : \nabla \delta ) ( CofF : \nabla \v ) }
void stiff_Jac_NH_S42term_Act ( Real coef,
                                const KNMK<Real> CofFk,
                                const KN<Real> Jk ,
                                const KN<Real>   Gammaf,
                                const KN<Real> I_4f,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                            {
                                s += -Gammaf (ig) * (1.0 + (Gammaf (ig) + 2.0) * pow ( (1.0 + Gammaf (ig) ), -2.) ) *
                                     pow ( Jk (ig), -8. / 3. ) * I_4f (ig) *
                                     CofFk ( jcoor, l, ig ) * fe.phiDer ( j, l, ig ) *
                                     CofFk ( icoor, k, ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}


//! 11. Jacobian matrix S5: Int { 1/3 * coef * J^(-2-2/3) * I_1E * gammaf * (CofF [\nabla \delta]^t CofF [fo \otimes fo]) : \nabla \v }
void stiff_Jac_NH_S5term_Act ( Real coef,
                               const KNMK<Real> CofFk,
                               const KN<Real> Jk ,
                               const KN<Real>   Gammaf,
                               const KNM<Real>   fo,
                               const KN<Real> I_1E,
                               MatrixElemental& elmat,
                               const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int m = 0; m < static_cast<Int> (nDimensions); ++m )
                            {
                                for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                                {
                                    s += Gammaf (ig) * pow ( Jk (ig), -8. / 3. ) * I_1E (ig) *
                                         CofFk ( icoor , l , ig ) * fe.phiDer ( j, l, ig ) *
                                         CofFk ( jcoor , m , ig ) * fo (m, ig) * fo (k, ig) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                                }
                            }
                        }
                    }
                    mat ( i, j ) += s * coef;
                }
            }
        }
    }
}


//! 12. Jacobian matrix S61 : Int { 2/9* coef * gammaf * I_1E * J^(-8/3) *( CofF : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
void stiff_Jac_NH_S61term_Act ( Real coef,
                                const KNMK<Real> CofFk,
                                const KN<Real>   Gammaf,
                                const KNM<Real>   fo,
                                const KN<Real> I_1E,
                                const KN<Real> Jk ,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int m = 0; m < static_cast<Int> (nDimensions); ++m )
                            {
                                for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                                {
                                    s += Gammaf (ig) * pow ( Jk (ig), -8. / 3. ) * I_1E (ig) *
                                         CofFk ( jcoor, m, ig ) * fo (m, ig) * fo (l, ig) * fe.phiDer ( j, l, ig ) *
                                         CofFk ( icoor, k, ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                                }
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}

//! 13. Jacobian matrix S62 : Int { -2/3* coef * gammaf *(1+gammaf) * J^(-5/3) *( F : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
void stiff_Jac_NH_S62term_Act ( Real coef,
                                const KNMK<Real> CofFk,
                                const KNMK<Real> Fk,
                                const KN<Real>   Gammaf,
                                const KNM<Real>   fo,
                                const KN<Real> Jk ,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int m = 0; m < static_cast<Int> (nDimensions); ++m )
                            {
                                for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                                {
                                    s += Gammaf (ig) * (1.0 + Gammaf (ig) ) * pow ( Jk (ig), -5. / 3. ) *
                                         CofFk ( jcoor, m, ig ) * fo (m, ig) * fo (l, ig) * fe.phiDer ( j, l, ig ) *
                                         Fk ( icoor, k, ig ) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                                }
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}


//! 14. Jacobian matrix S63 : Int { -2/3* coef * gammaf *(1+gammaf) * J^(-5/3) *( F [fo\otimes fo] : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
void stiff_Jac_NH_S63term_Act ( Real coef,
                                const KNMK<Real> CofFk,
                                const KNMK<Real> Fk,
                                const KN<Real>   Gammaf,
                                const KNM<Real>   fo,
                                const KN<Real> Jk ,
                                MatrixElemental& elmat,
                                const CurrentFE& fe )
{
    Real s;

    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); ++icoor )
    {
        for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); ++jcoor )
        {
            MatrixElemental::matrix_view mat = elmat.block ( icoor, jcoor );
            for ( Int i = 0; i < static_cast<Int> (fe.nbFEDof() ); ++i )
            {
                for ( Int j = 0; j < static_cast<Int> (fe.nbFEDof() ); ++j )
                {
                    s = 0.0;
                    for ( Int l = 0; l < static_cast<Int> (nDimensions); ++l )
                    {
                        for ( Int k = 0; k < static_cast<Int> (nDimensions); ++k )
                        {
                            for ( Int m = 0; m < static_cast<Int> (nDimensions); ++m )
                            {
                                for ( Int n = 0; m < static_cast<Int> (nDimensions); ++n )
                                {
                                    for ( Int ig = 0; ig < static_cast<Int> (fe.nbQuadPt() ); ++ig )
                                    {
                                        s += -pow (Gammaf (ig), 2.) * (1.0 + (Gammaf (ig) + 2.0) * pow ( (1.0 + Gammaf (ig) ), -2.) ) * pow ( Jk (ig), -5. / 3. ) *
                                             CofFk ( jcoor, m, ig ) * fo (m, ig) * fo (l, ig) * fe.phiDer ( j, l, ig ) *
                                             Fk ( icoor, n, ig ) * fo (n, ig) * fo (k, ig) * fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                                    }
                                }
                            }
                        }
                    }
                    mat ( i, j ) += coef * s;
                }
            }
        }
    }
}






//! ***********************************************************************************************
//! END OF ACTIVE NEO-HOOKEAN MODEL
//! ***********************************************************************************************

} //! End namespace AssemblyElementalActiveStructure

} //! End namespace LifeV

#endif
