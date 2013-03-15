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
 *  @brief mesh utilities to generate facets and ridges
 *
 *  @date 2013-03-15
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *
 *  @contributors Antonio Cervone <ant.cervone@gmail.com>
 *  @mantainer    Antonio Cervone <ant.cervone@gmail.com>
 */


#ifndef UPDATEMESHFACETSRIDGES_HPP
#define UPDATEMESHFACETSRIDGES_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>

namespace LifeV
{


//! Build localRidgeId table and optionally fills the list of Ridges
/**
 * @param mesh mesh object
 * @param createRidges is set true if we want also to create the actual list
 *  of edges. There is another utility (MeshChecks.hpp), which
 *  might be used for the same purpose if we want just to create the faces
 *  and not also the LocalRidgeID table.
 *  @param verbose If true, output is verbose.
 *  @param ridgeEstimate is a guess provided by the user of the total
 *  number of ridges. It is relevant only when createFacets=true. Setting it
 *  to a proper value helps in reducing time and memory.
 *  @param renumber Relevant only if createFacets=true.It makes sure that boundary edges are first
 *  if set to false possibly existing edges are never moved
 *
 *  @note This method does not assume that boundary edges are stores, since
 *  this condition is NOT a paradigm for a RegionMesh.
 */
template <typename RegionMeshType>
void updateMeshRidges ( RegionMeshType& mesh, bool createRidges = false, bool verbose = false, UInt ridgeEstimate = 0, bool renumber = true );


//! Build localFacetId table and optionally fills the list of Facets.
/**
 *  @param mesh mesh object
 *  @param createFacets is set true if we want also to create the actual list
 *  of internal facets. There is another utility (in mesh_util.hpp), which
 *  might be used for the same purpose if we want just to create the faces
 *  and not also the LocaFaceID table.
 *  @param verbose if true, output is verbose.
 *  @param estimateFacetNumber is a guess provided by the user of the total
 *  number of faces. It is relevant only when createFaces=true. Setting it
 *  to a proper value helps in reducing time and memory.
 *
 *  @note Facets are renumbered so that boundary facets are stored first
 *  @pre The routine assumes that the boundary facets are properly set, if not use the
 *  methods in MeshChecks.hpp
 *
 */
template <typename RegionMeshType>
void updateMeshFacets ( RegionMeshType& mesh, bool createFacets = false, bool verbose = false, UInt estimateFacets = 0 );

namespace MeshDetails
{

// 3D specialization
template <typename RegionMeshType>
void updateMeshRidges ( RegionMeshType& mesh, typename RegionMeshType::threeD_Type, bool createRidges, bool verbose, UInt ridgeEstimate, bool renumber )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );

    if ( RegionMeshType::S_geoDimensions != 3)
    {
        ERROR_MSG ("RegionMesh::updateElementRidges, It is not possible to use this method with 2D and 1D geometries.");
    }

    // If the counter is set we trust it! Otherwise we use Euler formula
    // this is ok for domains with at most 1 hole!

    if (verbose)
    {
        std::cout << "     Updating element ridges ... " << std::flush;
    }

    renumber = renumber && createRidges && !mesh.ridgeList().empty();
    if ( createRidges && ridgeEstimate == 0 )
    {
        ridgeEstimate = mesh.numEdges() > mesh.numBEdges() ? mesh.numEdges() : ( RegionMeshType::geoShape_Type::S_numFaces / 2 - 1 ) * mesh.numVolumes() + mesh.numBFaces() / 2 + mesh.numVertices();
    }


    if ( createRidges )
    {
        // We want to create the edges, we need to reserve space
        mesh.ridgeList().setMaxNumItems (ridgeEstimate);
    }
    MeshElementBareHandler<BareEdge> bareEdge;
    std::pair<UInt, bool> e;
    mesh.elemToRidge().reshape ( mesh.numLocalEdges(), mesh.numVolumes() ); // DIMENSION ARRAY

    UInt elemLocalID, i1, i2;
    std::pair<BareEdge, bool> _edge;
    typename RegionMeshType::geoShape_Type ele;
    typename RegionMeshType::facetShape_Type bele;
    // First We check if we have already Edges stored
    if ( ! mesh.ridgeList().empty() )
    {
        // dump first the existing edges, to maintain the correct numbering
        // if everything is correct the numbering in the bareedge
        // structure will reflect the actual edge numbering
        std::pair<UInt, bool> _check;
        for ( UInt j = 0; j < mesh.ridgeList().size(); ++j )
        {
            i1 = ( mesh.ridge ( j ).point ( 0 ) ).localId();
            i2 = ( mesh.ridge ( j ).point ( 1 ) ).localId();

            _edge  = makeBareEdge ( i1, i2 );
            _check = bareEdge.addIfNotThere ( _edge.first );
        }
    }

    typename RegionMeshType::ridge_Type edg;

    for ( typename RegionMeshType::faces_Type::iterator ifa = mesh.faceList.begin();
            ifa != mesh.faceList.begin() + mesh.numBFaces(); ++ifa )
    {
        for ( UInt j = 0; j < mesh.numLocalEdgesOfFace(); j++ )
        {
            i1 = bele.edgeToPoint ( j, 0 );
            i2 = bele.edgeToPoint ( j, 1 );
            // go to global
            i1 = ( ifa->point ( i1 ) ).localId();
            i2 = ( ifa->point ( i2 ) ).localId();

            _edge = makeBareEdge ( i1, i2 );

            e = bareEdge.addIfNotThere ( _edge.first );

            if ( createRidges && e.second )
            {
                //
                for ( UInt k = 0; k < 2 + RegionMeshType::facetShape_Type::S_numPointsPerEdge; k++ )
                {
                    UInt inode = bele.edgeToPoint (j, k);
                    edg.setPoint ( k, ifa->point ( inode ) );
                }
                MeshUtility::inheritPointsWeakerMarker ( edg );
                edg.setBoundary ( true );
                edg.setId ( mesh.ridgeList().size() );
                mesh.addRidge ( edg );
            }
        }

    }

    if ( createRidges )
    {
        mesh.setNumBEdges ( mesh.ridgeList().size() );
        mesh.setLinkSwitch ( "HAS_BOUNDARY_RIDGES" );
    }

    for ( typename RegionMeshType::elements_Type::iterator elemIt = mesh.elementList().begin();
            elemIt != mesh.elementList().end(); ++elemIt )
    {
        elemLocalID = elemIt->localId();

        for ( UInt j = 0; j < mesh.numLocalEdges(); j++ )
        {
            i1 = ele.edgeToPoint ( j, 0 );
            i2 = ele.edgeToPoint ( j, 1 );
            // go to global
            i1 = ( elemIt->point ( i1 ) ).localId();
            i2 = ( elemIt->point ( i2 ) ).localId();
            _edge = makeBareEdge ( i1, i2 );

            e = bareEdge.addIfNotThere ( _edge.first );
            mesh.elemToRidge() ( j, elemLocalID ) = e.first;
            if ( createRidges && e.second )
            {
                for ( UInt k = 0; k < 2 + RegionMeshType::geoShape_Type::S_numPointsPerEdge; k++ )
                {
                    UInt inode = ele.edgeToPoint (j, k);
                    edg.setPoint ( k, elemIt->point ( inode ) );
                }
                MeshUtility::inheritPointsWeakerMarker ( edg );
                edg.setBoundary ( true );
                edg.setId ( mesh.ridgeList().size() );
                mesh.addRidge ( edg );
            }
        }
    }

    if ( createRidges )
    {
        mesh.setNumEdges ( mesh.ridgeList().size() );
        mesh.setNumBEdges ( mesh.ridgeList().countElementsWithFlag (EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet) );
        mesh.setLinkSwitch ( "HAS_ALL_RIDGES" );
        if (mesh.numGlobalEdges() == 0)
        {
            mesh.setMaxNumGlobalEdges ( mesh.numEdges() );
        }
    }

    if (renumber && !mesh.ridgeList().empty() )
    {
        mesh.ridgeList().reorderAccordingToFlag (EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);
        std::vector<ID>newToOld = mesh.ridgeList().resetId(); //reset the ids so that they are in accord with position in the container.
        //Unfortunately I need oldToNew!
        std::vector<ID> oldToNew ( newToOld.size() );
        for (UInt j = 0; j < newToOld.size(); ++j)
        {
            oldToNew[ newToOld[j] ] = j;
        }
        // Save some memory annihilating newToOld
        std::vector<ID>().swap (newToOld);
        // Fix element to ridge array to reflect new ridge numbering
        // M_ElemToRidge is in fact a vector!
        std::vector<UInt> tmp ( mesh.elemToRidge().size() );
        std::vector<UInt>::iterator tmp_it = tmp.begin();
        for (std::vector<UInt>::iterator it = mesh.elemToRidge().begin(); it < mesh.elemToRidge().end(); ++it, ++tmp_it)
        {
            *tmp_it = oldToNew[*it];
        }
        std::copy (tmp.begin(), tmp.end(), mesh.elemToRidge().begin() );
    }

    UInt n = bareEdge.maxId();

    if (!createRidges)
    {
        if ( mesh.numEdges() == 0 || mesh.numEdges() == mesh.numBEdges() )
        {
            mesh.setNumEdges ( n );
        }
    }

    if (verbose)
    {
        std::cout << n << " edges found";
    }
    ASSERT_POS ( n == mesh.numEdges() , "#Edges found is not equal to that in RegionMesh" << n << " " << M_numEdges ) ;
    mesh.setLinkSwitch ( std::string ( "HAS_ELEMENT_TO_RIDGES" ) );

    if (verbose)
    {
        std::cout << " done." << std::endl;
    }
}

// 2D specialization
template <typename RegionMeshType>
inline void updateMeshRidges ( RegionMeshType&, typename RegionMeshType::twoD_Type, bool, bool, UInt, bool )
{}

// 1D specialization
template <typename RegionMeshType>
inline void updateMeshRidges ( RegionMeshType&, typename RegionMeshType::oneD_Type, bool, bool, UInt, bool )
{}

} // namespace MeshDetails

template <typename RegionMeshType>
inline void updateMeshRidges ( RegionMeshType& mesh, bool createRidges, bool verbose, UInt ridgeEstimate, bool renumber )
{
    MeshDetails::updateMeshRidges ( mesh, mesh.geoDim(), createRidges, verbose, ridgeEstimate, renumber );
}

template <typename RegionMeshType>
void updateMeshFacets ( RegionMeshType& mesh, bool createFacets, bool verbose, UInt facetEstimate )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );

    typedef BareEntitySelector<typename RegionMeshType::facetShape_Type::BasRefSha> bareEntitySelector_Type;
    typedef typename bareEntitySelector_Type::bareEntity_Type bareFacet_type;

    if (verbose)
    {
        std::cout << "     Updating element facets ... " << std::flush;
    }

    ASSERT0 ( ! createFacets || mesh.numBoundaryFacets() > 0, std::stringstream ( std::string ("Boundary Facets Must have been set") +
                                                                                  std::string ("in order to call updateElementFacets with createFacets=true") +
                                                                                  std::string ("\nUse buildBoundaryFacets(..) from mesh_util.h") ).str().c_str() );
    // If the counter is set we trust it! Otherwise we use Euler formula

    if ( createFacets && facetEstimate == 0 )
    {
        facetEstimate = mesh.numFacets() > mesh.numBoundaryFacets() ?
                        mesh.numFacets() : ( RegionMeshType::geoShape_Type::S_numFacets * mesh.numElements() + mesh.numBoundaryFacets() ) / 2;
    }

    ASSERT ( createFacets || mesh.numFacets() > 0 , "Mesh is not properly set!" );

    if ( createFacets )
    {
        mesh.facetList().setMaxNumItems ( facetEstimate );
    }



    typename RegionMeshType::facet_Type aFacet;

    MeshElementBareHandler<bareFacet_type> bareFacet;
    // Extra map for facets stored which are not boundary facets
    MeshElementBareHandler<bareFacet_type> extraBareFacet;
    std::pair<UInt, bool> e;
    mesh.elemToFacet().reshape ( RegionMeshType::element_Type::S_numLocalFacets, mesh.numElements() ); // DIMENSION ARRAY

    UInt elemLocalID;
    std::pair<bareFacet_type, bool>_facet;

    typename RegionMeshType::geoShape_Type ele;
    // If we have all facets and the facets store all adjacency info
    // everything is easier
    if ( (mesh.facetList().size() == mesh.numFacets() ) && mesh.getLinkSwitch ( "FACETS_HAVE_ADIACENCY" ) && mesh.getLinkSwitch ( "HAS_ALL_FACETS" ) )
    {
        for ( typename RegionMeshType::facets_Type::iterator itf = mesh.facetList().begin(); itf != mesh.facetList().end(); ++itf )
        {
            if ( itf->firstAdjacentElementPosition() != NotAnId && itf->firstAdjacentElementIdentity() != NotAnId)
            {
                mesh.elemToFacet() ( itf->firstAdjacentElementPosition() , itf->firstAdjacentElementIdentity() ) = itf->localId();
            }
            if ( itf->secondAdjacentElementPosition() != NotAnId && itf->secondAdjacentElementIdentity() != NotAnId)
            {
                mesh.elemToFacet() ( itf->secondAdjacentElementPosition(), itf->secondAdjacentElementIdentity() ) = itf->localId();
            }
        }
        // we finish here
        mesh.setLinkSwitch ( "HAS_ELEMENT_TO_FACETS" );
        if (verbose)
        {
            std::cout << " done." << std::endl;
        }

        return ;
    }

    // If I have only boundary facets I need to process them first to keep the correct numbering

    // First We check if we have already Facets stored
    UInt _numOriginalStoredFacets = mesh.facetList().size();
    ID points[RegionMeshType::facetShape_Type::S_numVertices];
    if ( ! mesh.facetList().empty() )
    {
        // dump all facets in the container, to maintain the correct numbering
        // if everything is correct the numbering in the bareFacet structure
        // will reflect the actual facet numbering. However, if I want to create
        // the internal facets I need to make sure that I am processing only the
        // boundary ones in a special way.
        std::pair<UInt, bool> _check;
        for ( UInt j = 0; j < mesh.facetList().size(); ++j )
        {
            for (UInt k = 0; k < RegionMeshType::facetShape_Type::S_numVertices; k++)
            {
                points[k] = ( mesh.facet ( j ).point ( k ) ).localId();
            }
            _facet = bareEntitySelector_Type::makeBareEntity ( points );
            _check = bareFacet.addIfNotThere ( _facet.first );
            if ( ! ( mesh.facet ( j ).boundary() ) )
            {
                extraBareFacet.addIfNotThere ( _facet.first, j);
            }
        }
    }
    UInt numFoundBoundaryFacets = bareFacet.size();
    UInt facetCount = numFoundBoundaryFacets;
    for ( typename RegionMeshType::elements_Type::iterator elemIt = mesh.elementList().begin();
            elemIt != mesh.elementList().end(); ++elemIt )
    {
        elemLocalID = elemIt->localId();
        for ( UInt j = 0; j < RegionMeshType::element_Type::S_numLocalFacets; j++ )
        {
            for (UInt k = 0; k < RegionMeshType::facetShape_Type::S_numVertices; k++)
            {
                UInt id = ele.facetToPoint ( j, k );
                points[k] = elemIt->point ( id ).localId();
            }
            _facet = bareEntitySelector_Type::makeBareEntity ( points );

            e = bareFacet.addIfNotThere ( _facet.first );
            mesh.elemToFacet() ( j, elemLocalID ) = e.first;
            bool _isBound = e.first < numFoundBoundaryFacets;
            // Is the facet an extra facet (not on the boundary but originally included in the list)?
            bool _isExtra = (e.first >= numFoundBoundaryFacets && e.first < _numOriginalStoredFacets);
            if ( _isBound )
            {
                typename RegionMeshType::facet_Type& _thisFacet ( mesh.facet (e.first) );
                _thisFacet.firstAdjacentElementIdentity()   = elemLocalID;
                _thisFacet.firstAdjacentElementPosition()   = j;
                _thisFacet.secondAdjacentElementIdentity()  = NotAnId;
                _thisFacet.secondAdjacentElementPosition()  = NotAnId;
            }
            else if (_isExtra)
            {
                // This is not a bfacets and I need to set up all info about adjacency properly
                typename RegionMeshType::facet_Type& _thisFacet ( mesh.facet (e.first) );
                // I need to check if it is the first time I meet it. Then I delete it from the
                // map: if it as there it means that it is the first time I am treating this face
                if (extraBareFacet.deleteIfThere (_facet.first) )
                {
                    // I need to be sure about orientation, the easiest thing is to rewrite the facet points
                    for ( UInt k = 0; k < RegionMeshType::facet_Type::S_numPoints; ++k )
                    {
                        _thisFacet.setPoint ( k, elemIt->point ( ele.facetToPoint ( j, k ) ) );
                    }
                    _thisFacet.firstAdjacentElementIdentity()  = elemLocalID;
                    _thisFacet.firstAdjacentElementPosition()  = j;

                }
                else
                {
                    _thisFacet.secondAdjacentElementIdentity()  = elemLocalID;
                    _thisFacet.secondAdjacentElementPosition()  = j;
                }
            }
            else if ( createFacets ) // A facet not contained in the original list.
                // I process it only if requested!
            {
                if ( e.second )
                {
                    // a new facet It must be internal.
                    for ( UInt k = 0; k < RegionMeshType::facet_Type::S_numPoints; ++k )
                    {
                        aFacet.setPoint ( k, elemIt->point ( ele.facetToPoint ( j, k ) ) );
                    }

                    aFacet.firstAdjacentElementIdentity()  = elemLocalID;
                    aFacet.firstAdjacentElementPosition() = j;

                    // gets the marker from the RegionMesh
                    aFacet.setMarkerID ( NotAnId );
                    aFacet.setBoundary (false);
                    aFacet.setId ( facetCount++ );
                    mesh.addFacet ( aFacet); //The id should be correct
                }
                else
                {
                    mesh.facet ( e.first ).secondAdjacentElementIdentity() = elemLocalID;
                    mesh.facet ( e.first ).secondAdjacentElementPosition() = j;
                }
            }
        }
    }

    UInt n = bareFacet.maxId();
    // LF Fix _numfacets. This part has to be checked. One may want to use
    // this method on a partitioned mesh, in which case the Global facets are there
    mesh.setNumFacets (n); // We have found the right total number of facets in the mesh
    if ( mesh.numGlobalFacets() == 0)
    {
        mesh.setMaxNumGlobalFacets (n);    // If not already set fix it.
    }

    if (verbose)
    {
        std::cout << n << " facets ";
    }
    //ASSERT_POS( n == M_numFacets , "#Facets found inconsistent with that stored in RegionMesh" ) ;
    mesh.setLinkSwitch ( "HAS_ELEMENT_TO_FACETS" );
    if ( createFacets )
    {
        mesh.setLinkSwitch ( "HAS_ALL_FACETS" );
    }
    //if ( cf ) Facets have adjacency in any case!
    mesh.setLinkSwitch ( "FACETS_HAVE_ADIACENCY" );
    if (verbose)
    {
        std::cout << " done." << std::endl;
    }
}

} // namespace LifeV

#endif //UPDATEMESHFACETSRIDGES_HPP
