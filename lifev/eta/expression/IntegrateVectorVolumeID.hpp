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
     @brief This file contains the definition of the IntegrateVectorElement class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_VECTOR_VOLUMEID_HPP
#define INTEGRATE_VECTOR_VOLUMEID_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>
#include <lifev/eta/fem/QRAdapterBase.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETVectorElemental.hpp>

#include <boost/shared_ptr.hpp>



namespace LifeV
{

namespace ExpressionAssembly
{

//! The class to actually perform the loop over the elements to assemble a vector
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is used to store the data required for the assembly of a vector and
  perform that assembly with a loop over the elements, and then, for each elements,
  using the Evaluation corresponding to the Expression (this convertion is done
  within a typedef).
 */
template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename QRAdapterType >
class IntegrateVectorVolumeID
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
            TestSpaceType::field_dim,
            0,
            3 >::evaluation_Type evaluation_Type;

    typedef typename MeshType::element_Type element_Type;


    typedef boost::shared_ptr<std::vector<element_Type*> > vectorVolumesPtr_Type;
    typedef boost::shared_ptr<std::vector<UInt> > vectorIndexesPtr_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateVectorVolumeID (const vectorVolumesPtr_Type volumeList,
                             const vectorIndexesPtr_Type indexList,
                             const QRAdapterType& qrAdapter,
                             const boost::shared_ptr<TestSpaceType>& testSpace,
                             const ExpressionType& expression);

    //! Copy constructor
    IntegrateVectorVolumeID ( const IntegrateVectorVolumeID < MeshType, TestSpaceType, ExpressionType, QRAdapterType>& integrator);

    //! Destructor
    ~IntegrateVectorVolumeID();

    //@}


    //! @name Operator
    //@{

    //! Operator wrapping the addTo method
    template <typename Vector>
    inline void operator>> (Vector& vector)
    {
        addTo (vector);
    }

    //! Operator wrapping the addTo method (for shared_ptr)
    template <typename Vector>
    inline void operator>> (boost::shared_ptr<Vector> vector)
    {
        addTo (*vector);
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
      performed: update the values, update the local vector,
      sum over the quadrature nodes, assemble in the global
      vector.
     */
    template <typename Vector>
    void addTo (Vector& vec);

    //@}

private:

    //! @name Private Methods
    //@{

    // No default constructor
    IntegrateVectorVolumeID();

    //@}

    //List of volumes with a marker
    vectorVolumesPtr_Type M_volumeList;
    vectorIndexesPtr_Type M_indexList;

    // Quadrature to be used
    QRAdapterType M_qrAdapter;

    // Shared pointer on the Space
    boost::shared_ptr<TestSpaceType> M_testSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    ETCurrentFE<3, 1>* M_globalCFE_std;
    ETCurrentFE<3, 1>* M_globalCFE_adapted;

    ETCurrentFE<3, TestSpaceType::field_dim>* M_testCFE_std;
    ETCurrentFE<3, TestSpaceType::field_dim>* M_testCFE_adapted;

    ETVectorElemental M_elementalVector;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateVectorVolumeID < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
IntegrateVectorVolumeID (const vectorVolumesPtr_Type volumeList,
                         const vectorIndexesPtr_Type indexList,
                         const QRAdapterType& qrAdapter,
                         const boost::shared_ptr<TestSpaceType>& testSpace,
                         const ExpressionType& expression)
    :   M_volumeList ( volumeList ),
        M_indexList ( indexList ),
        M_qrAdapter (qrAdapter),
        M_testSpace (testSpace),
        M_evaluation (expression),

        M_globalCFE_std (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) ),
        M_globalCFE_adapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) ),

        M_testCFE_std (new ETCurrentFE<3, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<3, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),

        M_elementalVector (TestSpaceType::field_dim * testSpace->refFE().nbDof() )
{
    M_evaluation.setQuadrature (qrAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_std);
    M_evaluation.setTestCFE (M_testCFE_std);
}


template <typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateVectorVolumeID <MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
IntegrateVectorVolumeID ( const IntegrateVectorVolumeID <MeshType,  TestSpaceType, ExpressionType, QRAdapterType>& integrator)
    :   M_volumeList (integrator.M_volumeList),
        M_indexList (integrator.M_indexList),
        M_qrAdapter (integrator.M_qrAdapter),
        M_testSpace (integrator.M_testSpace),
        M_evaluation (integrator.M_evaluation),

        M_globalCFE_std (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() ) ),
        M_globalCFE_adapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() ) ),

        M_testCFE_std (new ETCurrentFE<3, TestSpaceType::field_dim> (M_testSpace->refFE(), M_testSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<3, TestSpaceType::field_dim> (M_testSpace->refFE(), M_testSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),

        M_elementalVector (integrator.M_elementalVector)
{
    M_evaluation.setQuadrature (integrator.M_qrAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_std);
    M_evaluation.setTestCFE (M_testCFE_std);
}


template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateVectorVolumeID < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
~IntegrateVectorVolumeID()
{
    delete M_globalCFE_std;
    delete M_globalCFE_adapted;
    delete M_testCFE_std;
    delete M_testCFE_adapted;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
void
IntegrateVectorVolumeID <MeshType,  TestSpaceType, ExpressionType, QRAdapterType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
    out << " Elemental vector : " << std::endl;
    M_elementalVector.showMe (out);
    out << std::endl;
}

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
template <typename Vector>
void
IntegrateVectorVolumeID <MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
addTo (Vector& vec)
{
    // Defaulted to true for security
    bool isPreviousAdapted (true);

    //number of volumes
    UInt nbElements ( (*M_volumeList).size() );
    UInt nbIndexes ( (*M_indexList).size() );

    ASSERT ( nbElements == nbIndexes, "The number of indexes is different from the number of volumes!!!");

    UInt nbQuadPt_std (M_qrAdapter.standardQR().nbQuadPt() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );

    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Zeros out the elemental vector
        M_elementalVector.zero();

        // Update the quadrature rule adapter
        M_qrAdapter.update ( (*M_indexList) [iElement] );


        if (M_qrAdapter.isAdaptedElement() )
        {
            // Reset the quadrature in the different structures
            M_evaluation.setQuadrature ( M_qrAdapter.adaptedQR() );
            M_globalCFE_adapted -> setQuadratureRule ( M_qrAdapter.adaptedQR() );
            M_testCFE_adapted -> setQuadratureRule ( M_qrAdapter.adaptedQR() );

            // Reset the CurrentFEs in the evaluation
            M_evaluation.setGlobalCFE ( M_globalCFE_adapted );
            M_evaluation.setTestCFE ( M_testCFE_adapted );

            // Update with the correct element
            M_evaluation.update ( (*M_indexList) [iElement] );

            // Update the currentFEs
            M_globalCFE_adapted->update (* ( (*M_volumeList) [iElement]), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
            M_testCFE_adapted->update (* ( (*M_volumeList) [iElement]), evaluation_Type::S_testUpdateFlag);


            // Assembly
            for (UInt iblock (0); iblock < TestSpaceType::field_dim; ++iblock)
            {
                // Set the row global indices in the local vector
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    M_elementalVector.setRowIndex
                    (i + iblock * nbTestDof,
                     M_testSpace->dof().localToGlobalMap ( (*M_indexList) [iElement], i) + iblock * M_testSpace->dof().numTotalDof() );
                }

                // Make the assembly
                for (UInt iQuadPt (0); iQuadPt < M_qrAdapter.adaptedQR().nbQuadPt(); ++iQuadPt)
                {
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        M_elementalVector.element (i + iblock * nbTestDof) +=
                            M_evaluation.value_qi (iQuadPt, i + iblock * nbTestDof)
                            * M_globalCFE_adapted->wDet (iQuadPt);

                    }
                }
            }

            // Finally, set the flag
            isPreviousAdapted = true;
        }
        else
        {

            // Check if the last one was adapted
            if (isPreviousAdapted)
            {
                M_evaluation.setQuadrature ( M_qrAdapter.standardQR() );
                M_evaluation.setGlobalCFE ( M_globalCFE_std );
                M_evaluation.setTestCFE ( M_testCFE_std );

                isPreviousAdapted = false;
            }


            // Update the currentFEs
            M_globalCFE_std->update (* ( (*M_volumeList) [iElement]), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
            M_testCFE_std->update (* ( (*M_volumeList) [iElement]), evaluation_Type::S_testUpdateFlag);

            // Update the evaluation
            M_evaluation.update ( (*M_indexList) [iElement] );

            // Loop on the blocks
            for (UInt iblock (0); iblock < TestSpaceType::field_dim; ++iblock)
            {
                // Set the row global indices in the local vector
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    M_elementalVector.setRowIndex
                    (i + iblock * nbTestDof,
                     M_testSpace->dof().localToGlobalMap ( (*M_indexList) [iElement], i) + iblock * M_testSpace->dof().numTotalDof() );
                }

                // Make the assembly
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt_std; ++iQuadPt)
                {
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        M_elementalVector.element (i + iblock * nbTestDof) +=
                            M_evaluation.value_qi (iQuadPt, i + iblock * nbTestDof) *
                            M_globalCFE_std->wDet (iQuadPt);

                    }
                }
            }

        }
        M_elementalVector.pushToGlobal (vec);
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
