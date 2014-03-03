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
     @brief This file contains the definition of the integrate function.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/OpenMPParameters.hpp>

#include <lifev/eta/expression/RequestLoopElement.hpp>
#include <lifev/eta/expression/RequestLoopVolumeID.hpp>
#include <lifev/eta/expression/RequestLoopFaceID.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/QRAdapterBase.hpp>
#include <lifev/eta/fem/QRAdapterNeverAdapt.hpp>

#include <lifev/eta/expression/IntegrateMatrixElement.hpp>
#include <lifev/eta/expression/IntegrateVectorElement.hpp>
#include <lifev/eta/expression/IntegrateValueElement.hpp>


//Integration over portions of the domain
#include <lifev/eta/expression/IntegrateMatrixVolumeID.hpp>
#include <lifev/eta/expression/IntegrateVectorVolumeID.hpp>
#include <lifev/eta/expression/IntegrateVectorFaceID.hpp>
#include <lifev/eta/expression/IntegrateMatrixFaceID.hpp>

#include <lifev/eta/expression/IntegrateValueElementLSAdapted.hpp>
#include <lifev/eta/expression/IntegrateVectorElementLSAdapted.hpp>
#include <lifev/eta/expression/IntegrateMatrixElementLSAdapted.hpp>

#include <lifev/eta/expression/IntegrateMatrixFaceIDLSAdapted.hpp>
#include <lifev/eta/expression/IntegrateVectorFaceIDLSAdapted.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

/*!
  \namespace ExpressionAssembly

  Namespace for the assembly via expressions

 */
namespace ExpressionAssembly
{

//! Integrate function for matricial expressions
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a matrix
  with a loop on the elements.

  This function is repeated 4 times:
  versions with and without QR adapter
  versions with and without Offset

 */
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const UInt offsetUp = 0,
            const UInt offsetLeft = 0);
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const UInt offsetUp,
            const UInt offsetLeft)
{
    return IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>
           (request.mesh(), qrAdapterBase.implementation(), testSpace, solutionSpace, expression, offsetUp, offsetLeft);
}

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const UInt offsetUp = 0,
            const UInt offsetLeft = 0);
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const UInt offsetUp,
            const UInt offsetLeft)
{
    return IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
           (request.mesh(), QRAdapterNeverAdapt (quadrature), testSpace, solutionSpace, expression, offsetUp, offsetLeft);
}

//! Integrate function for matricial expressions (multi-threaded path)
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a matrix
  with a loop on the elements.

  This is an overload of the integrate function for matrices, which
  uses multiple threads to do assembly

  This function is repeated 4 times:
  versions with and without QR adapter
  versions with and without Offset

 */
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const OpenMPParameters& ompParams,
            const UInt offsetUp = 0,
            const UInt offsetLeft = 0);
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const OpenMPParameters& ompParams,
            const UInt offsetUp,
            const UInt offsetLeft)
{
    return IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
           (request.mesh(), qrAdapterBase.implementation(), testSpace, solutionSpace, expression,
            ompParams, offsetUp, offsetLeft);
}
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const OpenMPParameters& ompParams,
            const UInt offsetUp = 0,
            const UInt offsetLeft = 0);
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const OpenMPParameters& ompParams,
            const UInt offsetUp,
            const UInt offsetLeft)
            {
    return IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
           (request.mesh(), QRAdapterNeverAdapt (quadrature), testSpace, solutionSpace, expression,
            ompParams, offsetUp, offsetLeft);
}



//! Integrate function for vectorial expressions
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a vector
  with a loop on the elements.

  This function is repeated 4 times:
  versions with and without QR adapter
  versions with and without Offset

 */
template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression,
            const UInt offset = 0);
template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression,
            const UInt offset)
{
    return IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType, QRAdapterType>
           (request.mesh(), qrAdapterBase.implementation(), testSpace, expression, offset);
}

template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression,
            const UInt offset = 0);
template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType,QRAdapterNeverAdapt>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression,
            const UInt offset)
{
    return IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType, QRAdapterNeverAdapt>
           (request.mesh(), QRAdapterNeverAdapt (quadrature), testSpace, expression, offset);
}

//! Integrate function for benchmark expressions
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a benchmark
  with a loop on the elements.

  This function is repeated 2 times:
  versions with and without QR adapter

 */
template < typename MeshType, typename ExpressionType, typename QRAdapterType>
IntegrateValueElement<MeshType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const ExpressionType& expression)
{
    return IntegrateValueElement<MeshType, ExpressionType, QRAdapterType>
           (request.mesh(), qrAdapterBase.implementation(), expression);
}

template < typename MeshType, typename ExpressionType>
IntegrateValueElement<MeshType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const ExpressionType& expression)
{
    return IntegrateValueElement<MeshType, ExpressionType, QRAdapterNeverAdapt>
           (request.mesh(), QRAdapterNeverAdapt (quadrature), expression);
}

// =============================================================
// Methods to integrate over a portion of the mesh
// =============================================================
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixVolumeID<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopVolumeID<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression)
{
    return IntegrateMatrixVolumeID<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
           (request.volumeList(), request.indexList(), QRAdapterNeverAdapt (quadrature), testSpace, solutionSpace, expression);
}

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixVolumeID<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopVolumeID<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapter,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression)
{
    return IntegrateMatrixVolumeID<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>
           (request.volumeList(), request.indexList(), qrAdapter.implementation(), testSpace, solutionSpace, expression);
}

template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorVolumeID<MeshType, TestSpaceType, ExpressionType, QRAdapterNeverAdapt>
integrate ( const RequestLoopVolumeID<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression)
{
    return IntegrateVectorVolumeID<MeshType, TestSpaceType, ExpressionType, QRAdapterNeverAdapt> (request.volumeList(), request.indexList(), QRAdapterNeverAdapt (quadrature), testSpace, expression);
}

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateVectorVolumeID<MeshType, TestSpaceType, ExpressionType, QRAdapterType>
integrate ( const RequestLoopVolumeID<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapter,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression)
{
    return IntegrateVectorVolumeID<MeshType, TestSpaceType, ExpressionType, QRAdapterType>
           (request.volumeList(), request.indexList(), qrAdapter.implementation(), testSpace, expression);
}

/* Integration on the boundary of the domain */


template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorFaceID<MeshType, TestSpaceType, ExpressionType>
integrate ( const RequestLoopFaceID<MeshType>& request,
            const QuadratureBoundary& quadratureBoundary,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression)
{
    return IntegrateVectorFaceID<MeshType, TestSpaceType, ExpressionType>
           (request.mesh(), request.id(), quadratureBoundary, testSpace, expression);
}


template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixFaceID<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
integrate ( const RequestLoopFaceID<MeshType>& request,
            const QuadratureBoundary& quadratureBoundary,
            const boost::shared_ptr<TestSpaceType> testSpace,
            const boost::shared_ptr<SolutionSpaceType> solutionSpace,
            const ExpressionType& expression)
{
    return IntegrateMatrixFaceID<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
           (request.mesh(), request.id(), quadratureBoundary, testSpace, solutionSpace, expression);
}


template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>
integrate (const RequestLoopFaceID<MeshType>& request,
           const LevelSetBDQRAdapter<LSFESpaceType, LSVectorType>& quadratureAdapter,
           const boost::shared_ptr<TestSpaceType> testSpace,
           const boost::shared_ptr<SolutionSpaceType> solutionSpace,
           const ExpressionType& expression)
{
    return IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType> (request.mesh(), request.id(), quadratureAdapter, testSpace, solutionSpace, expression);
}

template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateVectorFaceIDLSAdapted < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>
integrate (const RequestLoopFaceID<MeshType>& request,
           const LevelSetBDQRAdapter<LSFESpaceType, LSVectorType>& quadratureAdapter,
           const boost::shared_ptr<TestSpaceType> testSpace,
           const ExpressionType& expression)
{
    return IntegrateVectorFaceIDLSAdapted < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType> (request.mesh(), request.id(), quadratureAdapter, testSpace, expression);
}



} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
