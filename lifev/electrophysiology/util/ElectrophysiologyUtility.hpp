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
    @brief Utilities

    @contributor Simone Palamara <palamara.simone@gmail.com>
    @maintainer Simone Palamara <palamara.simone@gmail.com>

    This file contains a set of base utilities used to applied current on a specified point.
 */

#ifndef ELECTROPHYSIOLOGYUTILITY_H
#define ELECTROPHYSIOLOGYUTILITY_H 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>

namespace LifeV
{

// Predeclaration

namespace ElectrophysiologyUtility
{

//! HeartUtility - A string parser grammar based on \c boost::spirit::qi
/*!
 *  @author(s) Simone Palamara
 *
 *  \c ElectrophysiologyUtility contains methods for applied current on a specified point.
 *
 */

//! @name Methods
//@{

//! Find closest point within radius and applied a constant current
/*!
 * @param point Vector of real containing the coordinates of the point within the radius
 * @param radius Radius used to find point.
 * @param appliedCurrentVector    Vector epetra containing the applied current.
 * @param valueAppliedCurrent    Value of the current to apply at the specified point.
 * @param fullMesh   Pointer to the mesh.
 */
template<typename Mesh> inline void appliedCurrentClosestPointWithinRadius(std::vector<Real>& point, Real Radius,boost::shared_ptr<VectorEpetra> appliedCurrentVector, Real valueAppliedCurrent,  boost::shared_ptr< Mesh > fullMesh )
{
    int n = appliedCurrentVector -> epetraVector().MyLength();

    int ids;
    for( UInt i(0); i < n; i++)
    {
    	 int iGID = appliedCurrentVector -> blockMap().GID(i);
    	 Real px = fullMesh -> point ( iGID ).x();
    	 Real py = fullMesh -> point ( iGID ).y();
    	 Real pz = fullMesh -> point ( iGID ).z();

    	 Real distance = std::sqrt( ( point[0] - px) * (point[0] - px)
    			 	 	 	 	  + ( point[1] - py) * (point[1] - py)
    			 	 	 	 	  + ( point[2] - pz) * (point[2] - pz) );
    	 if(distance <= Radius)
	 {
		ids=iGID;
		Radius=distance;
         }
	
	
    }
    double localRadius=Radius;
    double globalRadius(0);
    MPI_Barrier(MPI_COMM_WORLD);
   
    MPI_Allreduce(&localRadius,&globalRadius,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

    if(globalRadius==localRadius)
    {
	appliedCurrentVector-> operator []( ids ) = valueAppliedCurrent;
    }

}

//! Find all the points within radius and applied a constant current
/*!
 * @param point Vector of real containing the coordinates of the point within the radius
 * @param radius Radius used to find point.
 * @param appliedCurrentVector    Vector epetra containing the applied current.
 * @param valueAppliedCurrent    Value of the current to apply at the specified point.
 * @param fullMesh   Pointer to the mesh.
 */
template<typename Mesh> inline void appliedCurrentPointsWithinRadius(std::vector<Real>& point, Real Radius, boost::shared_ptr<VectorEpetra> appliedCurrentVector, Real valueAppliedCurrent,  boost::shared_ptr< Mesh > fullMesh )
{
    int n = appliedCurrentVector -> epetraVector().MyLength();

    std::vector<UInt> ids;
    for( UInt i(0); i < n; i++)
    {
    	 int iGID = appliedCurrentVector -> blockMap().GID(i);
    	 Real px = fullMesh -> point ( iGID ).x();
    	 Real py = fullMesh -> point ( iGID ).y();
    	 Real pz = fullMesh -> point ( iGID ).z();

    	 Real distance = std::sqrt( ( point[0] - px) * (point[0] - px)
    			 	 	 	 	  + ( point[1] - py) * (point[1] - py)
    			 	 	 	 	  + ( point[2] - pz) * (point[2] - pz) );
    	 if(distance <= Radius)
	 {
		ids.push_back(iGID);
         }
	
	
    }
    
    for(int i(0); i< ids.size(); i++)
    {
    	appliedCurrentVector -> operator []( ids.at(i) ) = valueAppliedCurrent;
    }
}



//@}

} // namespace ElectrophysiologyUtility

} // namespace LifeV

#endif /* HEARTUTILITY_H */
