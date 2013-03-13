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
    @brief Base utilities operating on meshes

    @contributor Tiziano Passerini <tiziano@mathcs.emory.edu>
    @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>

    @sa MeshUtility.hpp
 */

#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

namespace LifeV
{
namespace MeshUtility
{
// ===================================================
// Constructors & Destructor
// ===================================================
GetCoordComponent::GetCoordComponent() : componentIndex ( -1 )
{}

GetCoordComponent::GetCoordComponent ( Int i ) : componentIndex ( i )
{}

// ===================================================
// Operators
// ===================================================
void GetCoordComponent::operator() ( Real const& x, Real const& y, Real const& z, Real ret[ 3 ] ) const
{
    switch ( componentIndex )
    {
        case ( 0 ) :
            ret[ 0 ] = x;
            ret[ 1 ] = 0.0;
            ret[ 2 ] = 0.0;
            break;
        case ( 1 ) :
            ret[ 0 ] = 0.0;
            ret[ 1 ] = y;
            ret[ 2 ] = 0.0;
            break;
        case ( 2 ) :
            ret[ 0 ] = 0.0;
            ret[ 1 ] = 0.0;
            ret[ 2 ] = z;
            break;
        default:
            ret[ 0 ] = x;
            ret[ 1 ] = y;
            ret[ 2 ] = z;
    }
}

void GetOnes::operator() ( Real const& /*x*/, Real const& /*y*/, Real const& /*z*/, Real ret[ 3 ] )
const
{
    ret[ 0 ] = 1.0;
    ret[ 1 ] = 1.0;
    ret[ 2 ] = 1.0;
}


//void
//printMeshInfos( boost::shared_ptr<RegionMesh<LinearTetra, defaultMarkerCommon_Type > > mesh )
//{
//#ifdef HAVE_MPI
//    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_MpiComm( MPI_COMM_WORLD ) );
//#else
//    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );
//#endif
//	Displayer displayer( Comm );
//    MeshUtility::MeshStatistics::meshSize meshSize = MeshUtility::MeshStatistics::computeSize( *mesh );
//    displayer.leaderPrint( "Mesh size (max): ", meshSize.maxH, "\n" );
//    displayer.leaderPrint( "Mesh size (min): ", meshSize.minH, "\n" );
//    displayer.leaderPrint( "Mesh size (av.): ", meshSize.meanH, "\n" );
//}

} // namespace MeshUtility

} // namespace LifeV
