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
  @brief Cardiac Electrophysiology Test
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
  @date 11-2007
  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @last update 11-2010
 */

#ifndef __HEART_H
#define __HEART_H

#define MONODOMAIN

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#ifdef MONODOMAIN
#include <lifev/electrophysiology/solver/ElectroMonodomainSolver.hpp>
#else
#include <lifev/electrophysiology/solver/ElectroBidomainSolver.hpp>
#endif
#include <lifev/electrophysiology/solver/ElectroIonicSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroLuoRudy.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>



namespace LifeV
{
/*!
  \class Heart

  3D Action potential propagation class


 */
class Heart
{
public:

    //! @name Typedefs
    //@{

#ifdef MONODOMAIN
    typedef ElectroMonodomainSolver< RegionMesh<LinearTetra> >::vector_Type   vector_Type;
    typedef ElectroMonodomainSolver<RegionMesh<LinearTetra> >::matrix_Type        matrix_Type;
#else
    typedef ElectroBidomainSolver< RegionMesh<LinearTetra> >::vector_Type     vector_Type;
    typedef ElectroBidomainSolver<RegionMesh<LinearTetra> >::matrix_Type      matrix_Type;
#endif
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;
    typedef boost::shared_ptr<matrix_Type>                  matrixPtr_Type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Heart ( Int argc,
            char** argv );

    virtual ~Heart() {}

    //@}

    /** @name  Methods
     */
    //@{

    //! To build the system and iterate
    void run();

    //! To compute the righthand side of the system
#ifdef MONODOMAIN
    void computeRhs ( vector_Type& rhs,
                      ElectroMonodomainSolver< RegionMesh<LinearTetra> >& electricModel,
                      boost::shared_ptr< ElectroIonicSolver< RegionMesh<LinearTetra> > > ionicModel,
                      ElectroMonodomainData& dataMonodomain );
#else
    void computeRhs ( vector_Type& rhs,
                      ElectroBidomainSolver< RegionMesh<LinearTetra> >& electricModel,
                      boost::shared_ptr< ElectroIonicSolver< RegionMesh<LinearTetra> > > ionicModel,
                      ElectroBidomainData& dataBidomain );
#endif
    //@}


private:
    UInt ion_model;
    UInt nbeq;
    boost::shared_ptr<ElectroFunctors> M_heart_fct;
};
}
#endif /* __HEART_H */
