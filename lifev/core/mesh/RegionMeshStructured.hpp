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
    @brief Contains methods to generate nD regular meshes.

    @author Antonio Cervone <ant.cervone@gmail.com>
    @maintainer -

    @date 2013-03-15
 */

#ifndef REGIONMESHSTRUCTURED_HPP
#define REGIONMESHSTRUCTURED_HPP 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>

namespace LifeV
{

namespace MeshDetails
{

// 3D specialization
template <typename MeshType>
inline void regularMesh( MeshType& mesh, typename MeshType::threeD_Type, markerID_Type regionFlag,
                         const std::vector<UInt> numberOfElements, bool verbose, const Vector3D length, const Vector3D origin )
{
    ASSERT( numberOfElements.size() > 2, "the numberOfElements vector is too small" );

    regularMesh3D( mesh, regionFlag,
                   numberOfElements[0], numberOfElements[1], numberOfElements[2], verbose,
                   length[0], length[1], length[2],
                   origin[0], origin[1], origin[2] );
}

// 2D specialization
template <typename MeshType>
inline void regularMesh( MeshType& mesh, typename MeshType::twoD_Type, markerID_Type regionFlag,
                         const std::vector<UInt> numberOfElements, bool verbose, const Vector3D length, const Vector3D origin )
{
    ASSERT( numberOfElements.size() > 1, "the numberOfElements vector is too small" );

    regularMesh2D( mesh, regionFlag,
                   numberOfElements[0], numberOfElements[1], verbose,
                   length[0], length[1],
                   origin[0], origin[1] );
}

// 1D specialization
template <typename MeshType>
inline void regularMesh( MeshType& mesh, typename MeshType::oneD_Type, markerID_Type regionFlag,
                         const std::vector<UInt> numberOfElements, bool verbose, const Vector3D length, const Vector3D origin )
{
    ASSERT( numberOfElements.size() > 0, "the numberOfElements vector is too small" );

    regularMesh1D( mesh, regionFlag,
                   numberOfElements[0], verbose,
                   length[0], origin[0] );
}

} // namespace MeshDetails

template <typename MeshType>
inline void regularMesh ( MeshType& mesh,
                   markerID_Type regionFlag,
                   const std::vector<UInt> numberOfElements,
                   bool verbose = false,
                   const Vector3D length = Vector3D( 1., 1., 1. ),
                   const Vector3D origin = Vector3D( 0., 0., 0. ) )
{
    MeshDetails::regularMesh( mesh, mesh.geoDim(), regionFlag, numberOfElements, verbose, length, origin );
}

} // Namespace LifeV

#endif /* REGIONMESHSTRUCTURED_HPP */
