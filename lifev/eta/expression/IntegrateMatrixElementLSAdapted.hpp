//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the definition of the IntegrateMatrixElementLSAdapted class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_MATRIX_ELEMENT_LS_ADAPTED_HPP
#define INTEGRATE_MATRIX_ELEMENT_LS_ADAPTED_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/level_set/fem/LevelSetQRAdapter.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETMatrixElemental.hpp>

#include <boost/shared_ptr.hpp>



namespace LifeV
{

namespace ExpressionAssembly
{


//! The class to actually perform the loop over the elements to assemble a matrix
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is used to store the data required for the assembly of a matrix and
  perform that assembly with a loop over the elements, and then, for each elements,
  using the Evaluation corresponding to the Expression (This convertion is done
  within a typedef).

  The speciality of this class with respect to the LifeV::ExpressionAssembly::IntegrateMatrixElement class
  is that the quadrature can be adapted.
 */
template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
class IntegrateMatrixElementLSAdapted
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
            TestSpaceType::S_fieldDim,
            SolutionSpaceType::S_fieldDim,
            3 >::evaluation_Type  evaluation_Type;

    //! Type of the adapter
    typedef LevelSetQRAdapter<LSFESpaceType, LSVectorType> QRAdapter_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateMatrixElementLSAdapted (const boost::shared_ptr<MeshType>& mesh,
                                     const QRAdapter_Type& QRAdapter,
                                     const boost::shared_ptr<TestSpaceType>& testSpace,
                                     const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                                     const ExpressionType& expression);

    //! Copy constructor
    IntegrateMatrixElementLSAdapted (const IntegrateMatrixElementLSAdapted
                                     < MeshType, TestSpaceType, SolutionSpaceType,
                                     ExpressionType, LSFESpaceType, LSVectorType > & integrator);

    //! Destructor
    ~IntegrateMatrixElementLSAdapted();

    //@}


    //! @name Operators
    //@{

    //! Operator wrapping the addTo method
    template <typename MatrixType>
    inline void operator>> (MatrixType& mat)
    {
        addTo (mat);
    }

    template <typename MatrixType>
    inline void operator>> (boost::shared_ptr<MatrixType> mat)
    {
        addTo (*mat);
    }


    //@}


    //! @name Methods
    //@{

    //! Ouput method
    void check (std::ostream& out = std::cout);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.
     */
    template <typename MatrixType>
    void addTo (MatrixType& mat);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.

      Specialized for the case where the matrix is passed as a shared_ptr
     */
    template <typename MatrixType>
    inline void addTo (boost::shared_ptr<MatrixType> mat)
    {
        ASSERT (mat != 0, " Cannot assemble with an empty matrix");
        addTo (*mat);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    IntegrateMatrixElementLSAdapted();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Quadrature rule adapter to be used
    QRAdapter_Type M_QRAdapter;

    // Shared pointer on the Spaces
    boost::shared_ptr<TestSpaceType> M_testSpace;
    boost::shared_ptr<SolutionSpaceType> M_solutionSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    // Duplicated CurrentFE
    ETCurrentFE<3, 1>* M_globalCFE_unadapted;
    ETCurrentFE<3, 1>* M_globalCFE_adapted;

    ETCurrentFE<3, TestSpaceType::S_fieldDim>* M_testCFE_unadapted;
    ETCurrentFE<3, TestSpaceType::S_fieldDim>* M_testCFE_adapted;

    ETCurrentFE<3, SolutionSpaceType::S_fieldDim>* M_solutionCFE_unadapted;
    ETCurrentFE<3, SolutionSpaceType::S_fieldDim>* M_solutionCFE_adapted;

    ETMatrixElemental M_elementalMatrix;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateMatrixElementLSAdapted<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
IntegrateMatrixElementLSAdapted (const boost::shared_ptr<MeshType>& mesh,
                                 const QRAdapter_Type& QRAdapter,
                                 const boost::shared_ptr<TestSpaceType>& testSpace,
                                 const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                                 const ExpressionType& expression)
    :   M_mesh (mesh),
        M_QRAdapter (QRAdapter),
        M_testSpace (testSpace),
        M_solutionSpace (solutionSpace),
        M_evaluation (expression),

        M_globalCFE_unadapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), QRAdapter.standardQR() ) ),
        M_globalCFE_adapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), QRAdapter.standardQR() ) ),

        M_testCFE_unadapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (testSpace->refFE(), testSpace->geoMap(), QRAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (testSpace->refFE(), testSpace->geoMap(), QRAdapter.standardQR() ) ),

        M_solutionCFE_unadapted (new ETCurrentFE<3, SolutionSpaceType::S_fieldDim> (solutionSpace->refFE(), testSpace->geoMap(), QRAdapter.standardQR() ) ),
        M_solutionCFE_adapted (new ETCurrentFE<3, SolutionSpaceType::S_fieldDim> (solutionSpace->refFE(), testSpace->geoMap(), QRAdapter.standardQR() ) ),

        M_elementalMatrix (TestSpaceType::S_fieldDim * testSpace->refFE().nbDof(),
                           SolutionSpaceType::S_fieldDim * solutionSpace->refFE().nbDof() )
{
    M_evaluation.setQuadrature ( M_QRAdapter.standardQR() );
    M_evaluation.setGlobalCFE ( M_globalCFE_unadapted );
    M_evaluation.setTestCFE ( M_testCFE_unadapted );
    M_evaluation.setSolutionCFE ( M_solutionCFE_unadapted );
}

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateMatrixElementLSAdapted<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
IntegrateMatrixElementLSAdapted (const IntegrateMatrixElementLSAdapted < MeshType, TestSpaceType, SolutionSpaceType,
                                 ExpressionType, LSFESpaceType, LSVectorType > & integrator)
    :   M_mesh (integrator.M_mesh),
        M_QRAdapter (integrator.M_QRAdapter),
        M_testSpace (integrator.M_testSpace),
        M_solutionSpace (integrator.M_solutionSpace),
        M_evaluation (integrator.M_evaluation),

        M_globalCFE_unadapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) ),
        M_globalCFE_adapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) ),

        M_testCFE_unadapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (M_testSpace->refFE(), M_testSpace->geoMap(), M_QRAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (M_testSpace->refFE(), M_testSpace->geoMap(), M_QRAdapter.standardQR() ) ),

        M_solutionCFE_unadapted (new ETCurrentFE<3, SolutionSpaceType::S_fieldDim> (M_solutionSpace->refFE(), M_solutionSpace->geoMap(), M_QRAdapter.standardQR() ) ),
        M_solutionCFE_adapted (new ETCurrentFE<3, SolutionSpaceType::S_fieldDim> (M_solutionSpace->refFE(), M_solutionSpace->geoMap(), M_QRAdapter.standardQR() ) ),


        M_elementalMatrix (integrator.M_elementalMatrix)
{
    M_evaluation.setQuadrature (M_QRAdapter.standardQR() );
    M_evaluation.setGlobalCFE ( M_globalCFE_unadapted );
    M_evaluation.setTestCFE ( M_testCFE_unadapted );
    M_evaluation.setSolutionCFE ( M_solutionCFE_unadapted );
}

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateMatrixElementLSAdapted<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
~IntegrateMatrixElementLSAdapted()
{
    delete M_globalCFE_unadapted;
    delete M_globalCFE_adapted;
    delete M_testCFE_unadapted;
    delete M_testCFE_adapted;
    delete M_solutionCFE_unadapted;
    delete M_solutionCFE_adapted;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
void
IntegrateMatrixElementLSAdapted<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
    out << " Elemental matrix : " << std::endl;
    M_elementalMatrix.showMe (out);
    out << std::endl;
}

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
template <typename MatrixType>
void
IntegrateMatrixElementLSAdapted<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
addTo (MatrixType& mat)
{
    const UInt nbElements (M_mesh->numElements() );
    const UInt nbQuadPt_unadapted (M_QRAdapter.standardQR().nbQuadPt() );
    const UInt nbTestDof (M_testSpace->refFE().nbDof() );
    const UInt nbSolutionDof (M_solutionSpace->refFE().nbDof() );

    bool isPreviousAdapted (true);

    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Zeros out the matrix
        M_elementalMatrix.zero();

        // Update the adapter
        M_QRAdapter.update (iElement);


        if ( M_QRAdapter.isAdaptedElement() )
        {
            // Reset the QR in the CurrentFEs
            M_globalCFE_adapted->setQuadratureRule ( M_QRAdapter.adaptedQR() );
            M_testCFE_adapted->setQuadratureRule ( M_QRAdapter.adaptedQR() );
            M_solutionCFE_adapted->setQuadratureRule ( M_QRAdapter.adaptedQR() );

            // Reset the Evaluation
            M_evaluation.setGlobalCFE ( M_globalCFE_adapted );
            M_evaluation.setTestCFE ( M_testCFE_adapted );
            M_evaluation.setSolutionCFE ( M_solutionCFE_adapted );

            M_evaluation.setQuadrature ( M_QRAdapter.adaptedQR() );

            // Update the currentFEs
            M_globalCFE_adapted->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
            M_testCFE_adapted->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);
            M_solutionCFE_adapted->update (M_mesh->element (iElement), evaluation_Type::S_solutionUpdateFlag);

            // Update the evaluation
            M_evaluation.update (iElement);

            // Loop on the blocks

            for (UInt iblock (0); iblock < TestSpaceType::S_fieldDim; ++iblock)
            {
                for (UInt jblock (0); jblock < SolutionSpaceType::S_fieldDim; ++jblock)
                {

                    // Set the row global indices in the local matrix
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        M_elementalMatrix.setRowIndex
                        (i + iblock * nbTestDof,
                         M_testSpace->dof().localToGlobalMap (iElement, i) + iblock * M_testSpace->dof().numTotalDof() );
                    }

                    // Set the column global indices in the local matrix
                    for (UInt j (0); j < nbSolutionDof; ++j)
                    {
                        M_elementalMatrix.setColumnIndex
                        (j + jblock * nbSolutionDof,
                         M_solutionSpace->dof().localToGlobalMap (iElement, j) + jblock * M_solutionSpace->dof().numTotalDof() );
                    }

                    for (UInt iQuadPt (0); iQuadPt < M_QRAdapter.adaptedQR().nbQuadPt(); ++iQuadPt)
                    {
                        for (UInt i (0); i < nbTestDof; ++i)
                        {
                            for (UInt j (0); j < nbSolutionDof; ++j)
                            {
                                M_elementalMatrix.element (i + iblock * nbTestDof, j + jblock * nbSolutionDof) +=
                                    M_evaluation.value_qij (iQuadPt, i + iblock * nbTestDof, j + jblock * nbSolutionDof)
                                    * M_globalCFE_adapted->wDet (iQuadPt);

                            }
                        }
                    }
                }
            }

            isPreviousAdapted = true;
        }
        else
        {
            if (isPreviousAdapted)
            {
                // Reset the Evaluation
                M_evaluation.setGlobalCFE ( M_globalCFE_unadapted );
                M_evaluation.setTestCFE ( M_testCFE_unadapted );
                M_evaluation.setSolutionCFE ( M_solutionCFE_unadapted );

                M_evaluation.setQuadrature ( M_QRAdapter.standardQR() );

            }

            // Update the currentFEs
            M_globalCFE_unadapted->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
            M_testCFE_unadapted->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);
            M_solutionCFE_unadapted->update (M_mesh->element (iElement), evaluation_Type::S_solutionUpdateFlag);

            // Update the evaluation
            M_evaluation.update (iElement);

            // Loop on the blocks

            for (UInt iblock (0); iblock < TestSpaceType::S_fieldDim; ++iblock)
            {
                for (UInt jblock (0); jblock < SolutionSpaceType::S_fieldDim; ++jblock)
                {

                    // Set the row global indices in the local matrix
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        M_elementalMatrix.setRowIndex
                        (i + iblock * nbTestDof,
                         M_testSpace->dof().localToGlobalMap (iElement, i) + iblock * M_testSpace->dof().numTotalDof() );
                    }

                    // Set the column global indices in the local matrix
                    for (UInt j (0); j < nbSolutionDof; ++j)
                    {
                        M_elementalMatrix.setColumnIndex
                        (j + jblock * nbSolutionDof,
                         M_solutionSpace->dof().localToGlobalMap (iElement, j) + jblock * M_solutionSpace->dof().numTotalDof() );
                    }

                    for (UInt iQuadPt (0); iQuadPt < nbQuadPt_unadapted; ++iQuadPt)
                    {
                        for (UInt i (0); i < nbTestDof; ++i)
                        {
                            for (UInt j (0); j < nbSolutionDof; ++j)
                            {
                                M_elementalMatrix.element (i + iblock * nbTestDof, j + jblock * nbSolutionDof) +=
                                    M_evaluation.value_qij (iQuadPt, i + iblock * nbTestDof, j + jblock * nbSolutionDof)
                                    * M_globalCFE_unadapted->wDet (iQuadPt);

                            }
                        }
                    }
                }
            }


            isPreviousAdapted = false;
        }


        M_elementalMatrix.pushToGlobal (mat);
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif // INTEGRATE_MATRIX_ELEMENT_LS_ADAPTED_HPP
