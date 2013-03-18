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
 *  @file
 *  @brief MeshTransformer class
 *
 *  @date 2013-03-15
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *
 *  @contributors Antonio Cervone <ant.cervone@gmail.com>
 *  @mantainer    Antonio Cervone <ant.cervone@gmail.com>
 */

#ifndef MESHTRANSFORMER_HPP
#define MESHTRANSFORMER_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MapEpetra.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/numeric/ublas/matrix.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


namespace LifeV
{

namespace MeshUtility
{

/** Class to transform a mesh.
 * A class that implements methods to transform a mesh without changing
 * mesh connectivities. It has a constructor that takes the mesh to be transformed
 * @note The Template MarkerCommonType is used to compile with IBM AIX compilers
 * @author Luca Formaggia
 * @date 2 August 2011
 */

// The Template MarkerCommonType is used to compile with IBM compilers
template <typename RegionMeshType, typename MarkerCommonType = typename RegionMeshType::markerCommon_Type >
class MeshTransformer
{
public:

    typedef RegionMeshType regionMesh_Type;
    typedef MarkerCommonType markerCommon_Type;

    /** the constructor may take a reference to the mesh to be manipulated */
    MeshTransformer (RegionMeshType& m);
    /** Move the mesh according to a given displacement.
    *
    *  It moves the mesh from the last position saved with savePoints()
    *  For backward compatibility, if it is called before without calling
    *  savePoints(), the first time it is called it will save the current mesh point and then
    *  apply the movement.
    *
    *  Displacement is a 3*numpoints() VectorType which stores the x-displacement first,
    *  then the y-displacements etc.
    *
    *  The VectorType object must comply with lifeV distributed vector concept EpetraVector
    *  in particular it must have the methods isGlobalIDPresent(Uint i).
    *
    *  @author Miguel Fernandez
    *  @date 11/2002
    *
    *  @param disp Displacement vector. In this version it must be an EpetraVector
    *  @param dim  Length of vector disp.
    */
    template <typename VectorType>
    void moveMesh ( const VectorType& disp, UInt dim);
    /** Transform the mesh. It uses  boost::numeric::ublas (3,3) matrices
     *  scale, rotate and translate to perform the mesh movement
     *  (operations performed in this order).
     *  @date   14/09/2009
     *  @author Cristiano Malossi
     *  @note - Rotation follows Paraview conventions: first rotate around z-axis,
     *          then around y-axis and finally around x-axis;
     *  @note - All the vectors must allow the command: operator[];
     *  @param scale        vector of three components for (x,y,z) scaling of the mesh
     *  @param rotate       vector of three components (radiants) for rotating the mesh
     *  @param translate    vector of three components for (x,y,z) translation the mesh
     *
     */
    template <typename VectorType>
    void transformMesh ( const VectorType& scale, const VectorType& rotate, const VectorType& translate );

    /** Transform the mesh according to a given mapping.
     *  Transform the mesh according to a given meshMapping(Real& x, Real& y, Real& z).
     *  @date   12/2010
     *  @author Mauro Perego
     *  @param meshMapping   function void meshMmapping(Real& x, Real& y, Real& z) which receive
     *                   x, y, z, and transform them according to a certain mapping
     */
    template <typename FunctionType>
    void transformMesh ( const FunctionType& meshMapping);

    //! Tells if we store old points
    /**
     * If true than we can interrogate the old point position
     */
    bool hasOldPoint() const
    {
        return ! (this->M_pointList.empty() );
    }
    //! Saves the mesh points
    /**
     * Useful for algorithms which require to remember the position of the mesh
     * before the movement
     */
    void savePoints();
    /** Resets movement. Next step is like the mesh has never moved
     */
    void resetMovement()
    {
        this->M_pointList.clear();
    }

    /** Returns the i-th mesh Point before the last movement.

     *
     *  If the mesh points have not been saved with a previous call to
     *  savePoints() the method returns the mesh point
     *
     *  @note Avoid extensive use: it is inefficient. use pointListInitial() instead
     *  @param i Id of the Point.
     *  @return i-th mesh Point before the last movement.
     */
    typename RegionMeshType::point_Type const& pointInitial ( ID const i ) const;
    /** Returns a constant reference to the list of Points before the last movement.
      *
      *  If the mesh points have not been saved with a previous call to
      *  savePoints() returns the current mesh Points
      *
      *  @return The list mesh Point before the last movement.
      */
    typename RegionMeshType::points_Type const& pointListInitial() const;

private:
    /** Appropriately sets internal switches
     *
     *  It must be called by any mesh transformation method
     *  to ensure that the handling of (possibly) stored points
     *  works;
     */
    RegionMeshType& M_mesh;
    typename RegionMeshType::points_Type M_pointList;
};

// *****   IMPLEMENTATIONS ****
template <typename RegionMeshType, typename MarkerCommonType >
MeshTransformer<RegionMeshType, MarkerCommonType >::MeshTransformer (RegionMeshType& m) : M_mesh (m), M_pointList() {}
/**
 * @todo this method should be changed to make sure not to generate invalid elements
 */
template <typename RegionMeshType, typename MarkerCommonType >
template <typename VectorType>
void MeshTransformer<RegionMeshType, MarkerCommonType >::moveMesh ( const VectorType& disp, UInt dim )
{
    // the method must be called with a Repeated vector
    if ( disp.mapType() == Unique )
    {
#ifdef HAVE_LIFEV_DEBUG
        std::cerr << "Info: moveMesh() requires a Repeated vector, a copy of the passed Unique vector will be created.\n"
                  << "To optimize your code, you should pass a repeated vector to avoid the conversion." << std::endl;
#endif
        this->moveMesh ( VectorType ( disp, Repeated ), dim );
        return;
    }

    if ( !this->hasOldPoint() )
    {
        this->savePoints();
    }

    typedef typename RegionMeshType::points_Type points_Type;
    points_Type& pointList ( M_mesh.pointList );
    for ( UInt i = 0; i < M_mesh.pointList.size(); ++i )
    {
        for ( UInt j = 0; j < nDimensions; ++j )
        {
            int globalId = pointList[i].id();
            ASSERT ( disp.isGlobalIDPresent ( globalId + dim * j ), "global ID missing" );
            pointList[ i ].coordinate ( j ) = M_pointList[ i ].coordinate ( j ) + disp[ j * dim + globalId ];
        }
    }
}

template<typename RegionMeshType, typename MarkerCommonType >
void MeshTransformer<RegionMeshType, MarkerCommonType >::savePoints()
{
    if (M_pointList.capacity() < M_mesh.pointList.size() )
    {
        // Create space and add
        M_pointList.clear();
        M_pointList.reserve ( M_mesh.numPoints() );
        std::copy (M_mesh.pointList.begin(), M_mesh.pointList.end(),
                   std::back_inserter (M_pointList) );
    }
    else
    {
        // Overwrite
        std::copy (M_mesh.pointList.begin(), M_mesh.pointList.end(), M_pointList.begin() );
    }
}

template <typename RegionMeshType, typename MarkerCommonType >
const typename RegionMeshType::point_Type&
MeshTransformer<RegionMeshType, MarkerCommonType >::pointInitial ( ID const i ) const
{
    ASSERT_BD ( i < M_mesh.pointList.size() );
    return M_pointList.empty() ? M_mesh.pointList[i] : this->M_pointList[i];
}

template <typename RegionMeshType, typename MarkerCommonType >
const typename RegionMeshType::points_Type&
MeshTransformer<RegionMeshType, MarkerCommonType >::pointListInitial() const
{
    return M_pointList.empty() ? M_mesh.points_Type : M_pointList;
}

//! @todo Change using homogeneous coordinates to make it more efficient.
template <typename RegionMeshType, typename MarkerCommonType >
template <typename VectorType>
void MeshTransformer<RegionMeshType, MarkerCommonType >::transformMesh ( const VectorType& scale, const VectorType& rotate, const VectorType& translate )
{
    // Make life easier
    typename RegionMeshType::points_Type& pointList (M_mesh.pointList);

    //Create the 3 planar rotation matrix and the scale matrix
    boost::numeric::ublas::matrix<Real> R (3, 3), R1 (3, 3), R2 (3, 3), R3 (3, 3), S (3, 3);

    R1 (0, 0) =  1.;
    R1 (0, 1) =  0.;
    R1 (0, 2) =  0.;
    R1 (1, 0) =  0.;
    R1 (1, 1) =  std::cos (rotate[0]);
    R1 (1, 2) = -std::sin (rotate[0]);
    R1 (2, 0) =  0.;
    R1 (2, 1) =  std::sin (rotate[0]);
    R1 (2, 2) =  std::cos (rotate[0]);

    R2 (0, 0) =  std::cos (rotate[1]);
    R2 (0, 1) =  0.;
    R2 (0, 2) =  std::sin (rotate[1]);
    R2 (1, 0) =  0.;
    R2 (1, 1) =  1.;
    R2 (1, 2) = 0.;
    R2 (2, 0) = -std::sin (rotate[1]);
    R2 (2, 1) =  0.;
    R2 (2, 2) =  std::cos (rotate[1]);

    R3 (0, 0) =  std::cos (rotate[2]);
    R3 (0, 1) = -std::sin (rotate[2]);
    R3 (0, 2) = 0.;
    R3 (1, 0) =  std::sin (rotate[2]);
    R3 (1, 1) =  std::cos (rotate[2]);
    R3 (1, 2) = 0.;
    R3 (2, 0) =  0;
    R3 (2, 1) =  0.;
    R3 (2, 2) = 1.;

    S (0, 0) = scale[0];
    S (0, 1) = 0.;
    S (0, 2) = 0.;
    S (1, 0) = 0.;
    S (1, 1) = scale[1];
    S (1, 2) = 0.;
    S (2, 0) = 0.;
    S (2, 1) = 0.;
    S (2, 2) = scale[2];

    //The total rotation is: R = R1*R2*R3 (as in Paraview we rotate first around z, then around y, and finally around x).
    //We also post-multiply by S to apply the scale before the rotation.
    R = prod ( R3, S );
    R = prod ( R2, R );
    R = prod ( R1, R );

    //Create the 3D translate vector
    boost::numeric::ublas::vector<Real> P (3), T (3);
    T (0) = translate[0];
    T (1) = translate[1];
    T (2) = translate[2];

    //Apply the transformation
    for ( UInt i (0); i < pointList.size(); ++i )
    {
        //P = pointList[ i ].coordinate(); // Try to avoid double copy if possible

        P ( 0 ) = pointList[ i ].coordinate ( 0 );
        P ( 1 ) = pointList[ i ].coordinate ( 1 );
        P ( 2 ) = pointList[ i ].coordinate ( 2 );

        P = T + prod ( R, P );

        pointList[ i ].coordinate ( 0 ) = P ( 0 );
        pointList[ i ].coordinate ( 1 ) = P ( 1 );
        pointList[ i ].coordinate ( 2 ) = P ( 2 );
    }
}

//  The Template MarkerCommonType is used to compile with IBM compilers
template <typename RegionMeshType, typename MarkerCommonType >
template <typename FunctionType>
void MeshTransformer<RegionMeshType, MarkerCommonType >::transformMesh ( const FunctionType& meshMapping)
{
    // Make life easier
    typename RegionMeshType::points_Type& pointList (M_mesh.pointList);

    for ( UInt i = 0; i < pointList.size(); ++i )
    {
        typename RegionMeshType::point_Type& p = pointList[ i ];
        meshMapping (p.coordinate (0), p.coordinate (1), p.coordinate (2) );
    }
}

} // namespace MeshUtility

} // namespace LifeV

#endif // MESHTRANSFORMER_HPP
