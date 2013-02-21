//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
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
 *  @file
 *  @brief This file contains the definition for the Neo-Hookean activated material
 *  @brief The activation is assumed to be a time-space dependent scalar field acting
 *  @brief    along the fibers direction
 *  @version 1.0
 *  @date 29-02-2012
 *  @author Ricardo Ruiz-Baier

 *  @maintainer  Ricardo Ruiz-Baier  <ricardo.ruiz@epfl.ch>
  */

#ifndef _NEOHOOKEANACTIVATEDMATERIAL_H_
#define _NEOHOOKEANACTIVATEDMATERIAL_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/structure/solver/StructuralMaterial.hpp>

namespace LifeV
{

template <typename Mesh>
class NeoHookeanActivatedMaterial :
    public StructuralMaterial<Mesh>
{
    //!@name Type definitions
    //@{

public:
    typedef StructuralMaterial<Mesh>                 super;

    typedef VenantKirchhoffElasticData               data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;
    typedef typename super::dataPtr_Type    dataPtr_Type;
    typedef typename super::displayerPtr_Type    displayerPtr_Type;

    //! Definition of local scalar, vectors and tensors
    typedef KN<Real>                     KN_Type;
    typedef boost::shared_ptr<KN_Type>           KNPtr_Type;

    typedef KNM<Real>                    KNM_Type;
    typedef boost::shared_ptr<KNM_Type>          KNMPtr_Type;

    typedef KNMK<Real>                   KNMK_Type;
    typedef boost::shared_ptr<KNMK_Type>         KNMKPtr_Type;

    //@}



    //! @name Constructor &  Destructor
    //@{

    NeoHookeanActivatedMaterial();

    virtual  ~NeoHookeanActivatedMaterial();

    //@}



    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralMaterial
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    void setup ( const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& /*dataMaterial*/ );


    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix ( const vector_Type& disp,
                                const dataPtr_Type& dataMaterial,
                                const displayerPtr_Type& displayer);


    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms ( matrixPtr_Type& jacobian,
                                        const vector_Type& /*disp*/,
                                        const dataPtr_Type& /*dataMaterial*/,
                                        const displayerPtr_Type& /*displayer*/);


    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness ( const vector_Type& sol, Real factor, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );


    //! Computes the new Stiffness vector for Neo-Hookean and Exponential materials in StructuralSolver given a certain displacement field.
    //! This function is used both in StructuralSolver::evalResidual and in StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeVector ( const vector_Type& sol,
                         Real factor,
                         const dataPtr_Type& dataMaterial,
                         const displayerPtr_Type& displayer );


    //! Computes the deformation gradient F, the cofactor matrix Cof(F), the determinant of F (J = det(F)), the trace of right Cauchy-Green tensor tr(C)
    //! This function is used in StructuralMaterial::computeStiffness
    /*!
      \param dk_loc: the elemental displacement
    */
    void computeKinematicsVariables ( const VectorElemental& dk_loc);


    //! Computes the deformation Gradient F, the cofactor of F Cof(F), the determinant of F J = det(F), the trace of C Tr(C).
    /*!
      \param dk_loc: local displacement vector
    */
    void computeStress ( const vector_Type& sol);


    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMe ( std::string const& fileNameVectStiff,
                  std::string const& fileNameJacobain);

    //@}



    //! @name Get Methods
    //@{

    //! Get the Stiffness matrix
    matrixPtr_Type const stiffMatrix() const { }

    //! Get the stiffness vector
    vectorPtr_Type const stiffVector() const
    {
        return M_stiff;
    }

    void Apply ( const vector_Type& sol, vector_Type& res);

    //@}



protected:

    //! Local stress vector
    boost::scoped_ptr<VectorElemental>              M_elvecK;

    //! Elementary matrices
    boost::scoped_ptr<MatrixElemental>              M_elmatK;

    //! Vector: stiffness non-linear
    vectorPtr_Type                      M_stiff;

    //! First Piola-Kirchhoff stress tensor
    vectorPtr_Type                          M_FirstPiolaKStress;

    //! Local tensors initialization
    KNMPtr_Type                     M_FirstPiolaKStressEle;
    KNMKPtr_Type                    M_Fk;
    KNMKPtr_Type                    M_CofFk;
    KNMKPtr_Type                    M_Identity;
    KNPtr_Type                      M_Jack;
    KNPtr_Type                      M_trCisok;
    KNPtr_Type                      M_trCk;

    //! Local fibers vector
    KNMPtr_Type                                         M_fiberVector;

    KNPtr_Type                      M_Gammaf;
    KNPtr_Type                      M_I_4f;
    KNPtr_Type                      M_I_1E;

};





template <typename Mesh>
NeoHookeanActivatedMaterial<Mesh>::NeoHookeanActivatedMaterial() :
    super           ( ),
    M_elvecK            ( ),
    M_elmatK                    ( ),
    M_stiff                 ( ),
    M_FirstPiolaKStress     ( )
{
}





template <typename Mesh>
NeoHookeanActivatedMaterial<Mesh>::~NeoHookeanActivatedMaterial()
{}





template <typename Mesh>
void
NeoHookeanActivatedMaterial<Mesh>::setup ( const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                           const boost::shared_ptr<const MapEpetra>&            monolithicMap,
                                           const UInt                                           offset,
                                           const dataPtr_Type& dataMaterial,
                                           const displayerPtr_Type& displayer )
{
    this->M_displayer = displayer;
    this->M_dataMaterial  = dataMaterial;

    std::cout << "setting up the activated material" << std::endl;

    this->M_FESpace                     = dFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;
    this->M_dataMaterial                = dataMaterial;
    this->M_displayer                   = displayer;
    M_stiff.reset                   ( new vector_Type (*this->M_localMap) );


    M_FirstPiolaKStress.reset       ( new vector_Type (*this->M_localMap) );
    M_elvecK.reset          ( new VectorElemental (this->M_FESpace->fe().nbFEDof(), nDimensions) );
    this->M_elmatK.reset                ( new MatrixElemental ( this->M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );

    //! Local tensors initilization
    M_FirstPiolaKStressEle.reset    ( new KNM_Type ( nDimensions, nDimensions ) );
    M_Fk.reset                  ( new KNMK_Type ( nDimensions, nDimensions, dFESpace->fe().nbQuadPt() ) );
    M_CofFk.reset           ( new KNMK_Type ( nDimensions, nDimensions, dFESpace->fe().nbQuadPt() ) );
    M_Identity.reset                ( new KNMK_Type ( nDimensions, nDimensions, dFESpace->fe().nbQuadPt() ) );
    M_Jack.reset            ( new KN_Type ( dFESpace->fe().nbQuadPt() ) );
    M_trCisok.reset         ( new KN_Type ( dFESpace->fe().nbQuadPt() ) );
    M_trCk.reset            ( new KN_Type ( dFESpace->fe().nbQuadPt() ) );
    M_I_4f.reset            ( new KN_Type ( dFESpace->fe().nbQuadPt() ) );
    M_I_1E.reset            ( new KN_Type ( dFESpace->fe().nbQuadPt() ) );
}





template <typename Mesh>
void NeoHookeanActivatedMaterial<Mesh>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/)
{
    //! Empty method for active neo-hookean material
}





template <typename Mesh>
void NeoHookeanActivatedMaterial<Mesh>::updateJacobianMatrix ( const vector_Type&       disp,
                                                               const dataPtr_Type&      dataMaterial,
                                                               const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n*********************************\n  ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, dataMaterial, displayer);
    displayer->leaderPrint (" \n*********************************\n  ");
    std::cout << std::endl;

    this->M_jacobian->globalAssemble();
}





template <typename Mesh>
void NeoHookeanActivatedMaterial<Mesh>::updateNonLinearJacobianTerms ( matrixPtr_Type&       jacobian,
                                                                       const vector_Type&    disp,
                                                                       const dataPtr_Type&   dataMaterial,
                                                                       const displayerPtr_Type&  displayer )
{
    displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Activated Neo-Hookean)");

    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc (this->M_FESpace->fe().nbFEDof(), nDimensions);

    vector_Type dRep (disp, Repeated);

    //! Number of displacement components
    UInt nc = nDimensions;

    //! Nonlinear part of jacobian
    //! loop on volumes: assembling source term
    for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); ++i )
    {
        this->M_FESpace->fe().updateFirstDerivQuadPt ( this->M_FESpace->mesh()->volumeList ( i ) );
        this->M_elmatK->zero();

        UInt marker = this->M_FESpace->mesh()->volumeList ( i ).marker();
        Real mu     = dataMaterial->mu (marker);
        Real bulk   = dataMaterial->bulk (marker);

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof(); iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * this->M_FESpace->dim() + this->M_offset;
                dk_loc[iloc + iComp * this->M_FESpace->fe().nbFEDof()] = dRep[ig];
            }
        }

        //! Computes F, Cof(F), J = det(F), Tr(C)
        //! Passes Gammaf, fo
        //! and computes I_4f, I_1E
        computeKinematicsVariables ( dk_loc );

        //! Stiffness for non-linear terms of the active Neo-Hookean model
        /*!
         The results of the integrals are stored at each step into elmatK, until to build K matrix of the bilinear form
        */

        //! VOLUMETRIC PART (Same as in the passive mechanics)
        //! 1. Stiffness matrix: int { 1/2 * bulk * ( 2 - 1/J + 1/J^2 ) * ( CofF : \nabla \delta ) (CofF : \nabla v) }
        AssemblyElementalStructure::stiff_Jac_Pvol_1term ( 0.5 * bulk, (*M_CofFk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 2. Stiffness matrix: int { 1/2 * bulk * ( 1/J- 1 - log(J)/J^2 ) * ( CofF [\nabla \delta]^t CofF ) : \nabla v }
        AssemblyElementalStructure::stiff_Jac_Pvol_2term ( 0.5 * bulk, (*M_CofFk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! ISOCHORIC PART
        //! 1. Stiffness matrix : int { -2/3 * mu * (1+Gammaf) * J^(-5/3) *( CofF : \nabla \delta ) ( F : \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_P1iso_NH_1term_Act ( (-2.0 / 3.0) * mu, (*M_CofFk), (*M_Fk), (*M_Gammaf), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 2. Stiffness matrix : int { 2/9 * mu * (1+Gammaf) * ( Ic_iso / J^2 )( CofF : \nabla \delta ) ( CofF : \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_P1iso_NH_2term_Act ( (2.0 / 9.0) * mu, (*M_CofFk), (*M_Gammaf), (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );

        //! 3. Stiffness matrix : int { mu * (1+Gammaf) * J^(-2/3) (\nabla \delta : \nabla \v)}
        AssemblyElementalActiveStructure::stiff_Jac_P1iso_NH_3term_Act ( mu, (*M_Gammaf), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 4. Stiffness matrix : int { -2/3 * mu * (1+Gammaf) * J^(-5/3) ( F : \nabla \delta ) ( CofF : \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_P1iso_NH_4term_Act ( (-2.0 / 3.0) * mu, (*M_CofFk), (*M_Fk), (*M_Gammaf), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 5. Stiffness matrix : int { 1/3 * mu * (1+Gammaf) * J^(-2) * Ic_iso * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
        AssemblyElementalActiveStructure::stiff_Jac_P1iso_NH_5term_Act ( (1.0 / 3.0) * mu, (*M_CofFk), (*M_Gammaf), (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );

        //! 6. Jacobian matrix S1 : Int { coef * g(Gammaf) * J^(-2/3) (\nabla \delta [fo \tomes fo] : \nabla \v)}
        AssemblyElementalActiveStructure::stiff_Jac_NH_S1term_Act ( mu, (*M_Gammaf), (*M_fiberVector), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 7. Jacobian matrix S2 : Int { -2/3* coef * g(Gammaf) * J^(-5/3) *( CofF : \nabla \delta ) ( F [fo\otimes fo]: \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S2term_Act ( (-2.0 / 3.0) * mu, (*M_CofFk),  (*M_Fk), (*M_Gammaf), (*fo), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 8. Jacobian matrix S3: Int { 1/3 * coef * J^(-2-2/3) * I_4f * (CofF [\nabla \delta]^t CofF ) : \nabla \v }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S3term_Act ( (1.0 / 3.0) * mu, (*M_CofFk), (*M_Jack), (*M_Gammaf), (*I_4f), *this->M_elmatK, this->M_FESpace->fe() );

        //! 9. Jacobian matrix S41 : Int { -2/3* coef * g(Gammaf) * J^(-5/3) *( F [fo\otimes fo] : \nabla \delta ) ( CofF: \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S41term_Act (- (2.0 / 3.0) * mu, (*M_CofFk), (*M_Fk), (*M_Gammaf), (*M_fiberVector), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 10. Jacobian matrix S42: Int { 2/9 * coef * g(Gammaf) * I_4f * J^{-8/3} ( CofF : \nabla \delta ) ( CofF : \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S42term_Act ( (2.0 / 9.0) * mu, (*M_CofFk), (*M_Jack), (*M_Gammaf), (*I_4f), *this->M_elmatK, this->M_FESpace->fe() );

        //! 11. Jacobian matrix S5: Int { 1/3 * coef * J^(-2-2/3) * I_1E * gammaf * (CofF [\nabla \delta]^t CofF [fo \otimes fo]) : \nabla \v }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S5term_Act ( (1.0 / 3.0) * mu, (*M_CofFk), (*M_Jack), (*M_Gammaf), (*M_fiberVector), (*I_1E), *this->M_elmatK, this->M_FESpace->fe() );

        //! 12. Jacobian matrix S61 : Int { 2/9* coef * gammaf * I_1E * J^(-8/3) *( CofF : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S61term_Act ( (2.0 / 9.0) * mu, (*M_CofFk), (*M_Gammaf), (*M_fiberVector), (*I_1E), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 13. Jacobian matrix S62 : Int { -2/3* coef * gammaf *(1+gammaf) * J^(-5/3) *( F : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S62term_Act (- (2.0 / 3.0) * mu, (*M_CofFk), (*M_Fk), (*M_Gammaf), (*M_fiberVector), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! 14. Jacobian matrix S63 : Int { -2/3* coef * gammaf *(1+gammaf) * J^(-5/3) *( F [fo\otimes fo] : \nabla \delta ) ( CofF [fo\otimes fo]: \nabla \v ) }
        AssemblyElementalActiveStructure::stiff_Jac_NH_S63term_Act ( (1.0 / 3.0) * mu, (*M_CofFk), (*M_Fk), (*M_Gammaf), (*M_fiberVector), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

        //! assembling
        for ( UInt ic = 0; ic < nc; ++ic )
        {
            for ( UInt jc = 0; jc < nc; jc++ )
            {
                assembleMatrix ( *jacobian, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic * totalDof, this->M_offset +  jc * totalDof );
            }
        }
    }
}

template <typename Mesh>
void NeoHookeanActivatedMaterial<Mesh>::Apply ( const vector_Type& sol, vector_Type& res )
{
    computeStiffness (sol, 0., this->M_dataMaterial, this->M_displayer);
    res += *M_stiff;
}


template <typename Mesh>
void NeoHookeanActivatedMaterial<Mesh>::computeStiffness ( const vector_Type&       sol,
                                                           Real                     /*factor*/,
                                                           const dataPtr_Type&      dataMaterial,
                                                           const displayerPtr_Type& displayer )
{
    this->M_stiff.reset (new vector_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n******************************************************************\n  ");
    displayer->leaderPrint (" Non-Linear S-  Computing the stiffness vector - Active neo-Hookean "     );
    displayer->leaderPrint (" \n******************************************************************\n  ");

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    UInt dim        = this->M_FESpace->dim();

    VectorElemental dk_loc ( this->M_FESpace->fe().nbFEDof(), nDimensions );
    //vector_Type disp(sol);

    vector_Type dRep (sol, Repeated);

    for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); i++ )
    {
        this->M_FESpace->fe().updateFirstDerivQuadPt ( this->M_FESpace->mesh()->volumeList ( i ) );

        UInt marker = this->M_FESpace->mesh()->volumeList ( i ).marker();

        //! Getting the proper coefficients
        Real mu     = dataMaterial->mu (marker);
        Real bulk   = dataMaterial->bulk (marker);

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * dim + this->M_offset;
                dk_loc[ iloc + iComp * this->M_FESpace->fe().nbFEDof() ] = dRep[ig];
            }
        }

        this->M_elvecK->zero();

        computeKinematicsVariables ( dk_loc );

        //! Stiffness for non-linear terms of the Neo-Hookean model
        /*!
        The results of the integrals are stored at each step into elvecK, until to build K matrix of the bilinear form
        */
        //! Volumetric part
        /*!
        Source term Pvol: int { bulk /2* (J1^2 - J1  + log(J1) ) * 1/J1 * (CofF1 : \nabla v) }
        */
        AssemblyElementalStructure::source_Pvol ( 0.5 * bulk, (*M_CofFk), (*M_Jack), *this->M_elvecK,  this->M_FESpace->fe() );

        //! Isochoric part
        /*!
        Source term Active P1iso_NH_Act
        */
        AssemblyElementalActiveStructure::source_P1iso_NH_Act ( mu, (*M_CofFk) , (*M_Fk), (*M_Gammaf), (*M_Jack),  (*M_trCisok) , *this->M_elvecK,  this->M_FESpace->fe() );

        for ( UInt ic = 0; ic < nDimensions; ++ic )
        {
            /*!
            M_elvecK is assemble into *vec_stiff vector that is recall
            from updateSystem(matrix_ptrtype& mat_stiff, vector_ptr_type& vec_stiff)
                */
            assembleVector ( *this->M_stiff, *this->M_elvecK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
        }
    }

    this->M_stiff->globalAssemble();
}





template <typename Mesh>
void NeoHookeanActivatedMaterial<Mesh>::computeKinematicsVariables ( const VectorElemental& dk_loc, const vector_Type& rawFiberVector)
{
    Real s;
    UInt eleID = this->M_FESpace->fe().currentLocalId();
    VectorElemental aux0 (dFESpace->fe().nbQuadPt() ) );
    VectorElemental aux1 (dFESpace->fe().nbQuadPt() ) );
    VectorElemental aux2 (dFESpace->fe().nbQuadPt() ) );

    //!Definition of the local value of the fibers direction
    for ( Int ig = 0; ig < static_cast<Int> (this->M_FESpace->fe().nbQuadPt() ); ig++ )
{
    aux0[ig] = aux1[ig] = aux2[ig] = 0;
        for ( Int i = 0; i < static_cast<Int> (this->M_FESpace->fe().nbFEDof() ); i++ )
        {
            aux0[ig] += rawFiberVector[this->M_FESpace->dof().localToGlobalMap ( eleID, i)] * this->M_FESpace->fe().phi ( i, ig );
            aux1[ig] += rawFiberVector[this->M_FESpace->dof().localToGlobalMap ( eleID, i) + this->M_FESpace->dim()] * this->M_FESpace->fe().phi ( i, ig );
            aux2[ig] += rawFiberVector[this->M_FESpace->dof().localToGlobalMap ( eleID, i) + 2 * this->M_FESpace->dim()] * this->M_FESpace->fe().phi ( i, ig );
        }
        Real norm = sqrt (aux0[ig] * aux0[ig] + aux1[ig] * aux1[ig] + aux2[ig] * aux2[ig]);
        (*M_fiberVector) (0, ig) = aux0[ig] / norm;
        (*M_fiberVector) (1, ig) = aux1[ig] / norm;
        (*M_fiberVector) (2, ig) = aux2[ig] / norm;
    }



    //! loop on quadrature points (ig)
    for ( Int ig = 0; ig < static_cast<Int> (this->M_FESpace->fe().nbQuadPt() ); ig++ )
{
    //! loop on space coordinates (icoor)
    for ( Int icoor = 0; icoor < static_cast<Int> (nDimensions); icoor++ )
        {
            //! loop on space coordinates (jcoor)
            for ( Int jcoor = 0; jcoor < static_cast<Int> (nDimensions); jcoor++ )
            {
                s = 0.0;
                for ( Int i = 0; i < static_cast<Int> (this->M_FESpace->fe().nbFEDof() ); i++ )
                {
                    s += this->M_FESpace->fe().phiDer ( i, jcoor, ig ) * dk_loc[ i + icoor * this->M_FESpace->fe().nbFEDof() ]; // \grad u^k at a quadrature point
                }
                //! gradient of displacement
                (*M_Fk) ( icoor , jcoor , ig ) = s;
            }
        }
    }

    //! loop on quadrature points (ig)
    for ( Int ig = 0; ig < static_cast<Int> (this->M_FESpace->fe().nbQuadPt() ); ig++ )
{
    for ( Int  icoor = 0; icoor < static_cast<Int> (nDimensions); icoor++ )
        {
            //! deformation gradient Fk
            (*M_Fk) ( icoor , icoor , ig ) +=  1.0;
        }
    }

    Real a, b, c, d, e, f, g, h, i;

    for ( Int ig = 0; ig < static_cast<Int> (this->M_FESpace->fe().nbQuadPt() ); ig++)
{
    a = (*M_Fk) ( 0 , 0 , ig );
        b = (*M_Fk) ( 0 , 1 , ig );
        c = (*M_Fk) ( 0 , 2 , ig );
        d = (*M_Fk) ( 1 , 0 , ig );
        e = (*M_Fk) ( 1 , 1 , ig );
        f = (*M_Fk) ( 1 , 2 , ig );
        g = (*M_Fk) ( 2 , 0 , ig );
        h = (*M_Fk) ( 2 , 1 , ig );
        i = (*M_Fk) ( 2 , 2 , ig );

        //! determinant of deformation gradient Fk
        (*M_Jack) (ig) = a * ( e * i - f * h ) - b * ( d * i - f * g ) + c * ( d * h - e * g );

        (*M_CofFk) ( 0 , 0 , ig ) =   ( e * i - f * h );
        (*M_CofFk) ( 0 , 1 , ig ) = - ( d * i - g * f );
        (*M_CofFk) ( 0 , 2 , ig ) =   ( d * h - e * g );
        (*M_CofFk) ( 1 , 0 , ig ) = - ( b * i - c * h );
        (*M_CofFk) ( 1 , 1 , ig ) =   ( a * i - c * g );
        (*M_CofFk) ( 1 , 2 , ig ) = - ( a * h - g * b );
        (*M_CofFk) ( 2 , 0 , ig ) =   ( b * f - c * e );
        (*M_CofFk) ( 2 , 1 , ig ) = - ( a * f - c * d );
        (*M_CofFk) ( 2 , 2 , ig ) =   ( a * e - d * b );
    }

    for ( Int ig = 0; ig <  static_cast<Int> (this->M_FESpace->fe().nbQuadPt() )  ; ig++ )
{
    if ( (*M_Jack) (ig) < 0)
        {
            std::cout << "negative jacobian !!!!!! ERROR in computeKinematicsVariables" <<  std::endl;
        }
    }

    //! loop on quadrature points (ig)
    for ( Int ig = 0; ig < static_cast<Int> (this->M_FESpace->fe().nbQuadPt() ); ig++ )
{
    s = 0.0;
    for ( Int i = 0; i < static_cast<Int> (nDimensions); i++)
        {
            for ( Int j = 0; j < static_cast<Int> (nDimensions); j++ )
            {
                //! trace of  C1 = (F1k^t F1k)
                s +=  (*M_Fk) ( i , j , ig ) * (*M_Fk) ( i , j , ig );
            }
        }
        //! trace of C
        (*M_trCk) ( ig ) = s;
    }

    for ( Int ig = 0; ig <  static_cast<Int> (this->M_FESpace->fe().nbQuadPt() ); ig++ )
{
    //! trace of deviatoric C
    (*M_trCisok) ( ig ) =  pow ( (*M_Jack) ( ig ), -2. / 3.) * (*M_trCk) ( ig );
    }
}





template <typename Mesh>
void NeoHookeanActivatedMaterial<Mesh>::showMe ( std::string const& fileNameStiff,
                                                 std::string const& fileNameJacobian)
{
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}





template <typename Mesh>
inline StructuralMaterial<Mesh>* createNeoHookeanActivatedMaterial()
{
    return new NeoHookeanActivatedMaterial<Mesh >();
}
namespace
{
static bool registerNH = StructuralMaterial<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct ( "neoHookeanActivated", &createNeoHookeanActivatedMaterial<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __NEOHOOKENANACTIVATEDMATERIAL_H */
