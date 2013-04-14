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

    This file contains a set of base utilities used to test mesh entities or
    operate on them
 */

#ifndef HEARTUTILITY_H
#define HEARTUTILITY_H 1

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

namespace HeartUtility
{

inline void importFibers( boost::shared_ptr<VectorEpetra> fiber, std::string& name, boost::shared_ptr< RegionMesh<LinearTetra> > mesh  )
{
    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef ExporterData<mesh_Type> 						   exporterData_Type;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >  hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                  hdf5FilterPtr_Type;


    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > fiberSpace( new FESpace< mesh_Type, MapEpetra > ( mesh, "P1", 3, comm ) );

    exporterData_Type impData (exporterData_Type::VectorField, "fibers.00000", fiberSpace,
                               fiber, UInt (0), exporterData_Type::UnsteadyRegime);

    //    filterPtr_Type importer( new hdf5Filter_Type(dataFile, name) );
    filterPtr_Type importer ( new hdf5Filter_Type() );
    importer -> setMeshProcId ( mesh, comm -> MyPID() );
    importer-> setPrefix (name);
    importer -> readVariable (impData);
    importer -> closeFile();

}

inline void importFibers( boost::shared_ptr<VectorEpetra> fiberVector, std::string filename, std::string filepath )
{

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef ExporterData<mesh_Type> 						   exporterData_Type;
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;

    std::ifstream fibers ( (filepath+filename).c_str() );

	UInt NumGlobalElements =  fiberVector -> size();
	std::vector<Real> fiber_global_vector (NumGlobalElements);

	//Importing fibers
	for ( UInt i = 0; i < NumGlobalElements; ++i)
	{
		fibers >> fiber_global_vector[i];
	}

	int n = (*fiberVector).epetraVector().MyLength();
	int d = n / 3;
	int i (0);
	int j (0);
	int k (0);

	//std::cout << "\nAssigning fibers to the vector epetra...";

    for (UInt l = 0; l < d; ++l)
	{
	i = (*fiberVector).blockMap().GID (l);
	j = (*fiberVector).blockMap().GID (l + d);
	k = (*fiberVector).blockMap().GID (l + 2 * d);
	(*fiberVector) [i] = fiber_global_vector[3 * i];
	(*fiberVector) [j] = fiber_global_vector[3 * i + 1];
	(*fiberVector) [k] = fiber_global_vector[3 * i + 2];

	//normalizing
    Real norm = std::sqrt( (*fiberVector) [i] * (*fiberVector) [i] + (*fiberVector) [j] * (*fiberVector) [j] + (*fiberVector) [k] * (*fiberVector) [k] );

    (*fiberVector) [i] = (*fiberVector) [i] / norm;
    (*fiberVector) [j] = (*fiberVector) [j] / norm;
    (*fiberVector) [k] = (*fiberVector) [k] / norm;
	}
	std::cout << std::endl;
	fiber_global_vector.clear();

}


inline void setupFibers ( VectorEpetra& vec, VectorSmall<3>& fibers)
{
    int n1 = vec.epetraVector().MyLength();
    int d1 = n1 / 3;
    vec *= 0;
    int i (0);
    int j (0);
    int k (0);

    for ( int l (0); l < d1; l++)
    {
        i = vec.blockMap().GID (l);
        j = vec.blockMap().GID (l + d1);
        k = vec.blockMap().GID (l + 2 * d1);
        vec [i] = fibers[0];
        vec [j] = fibers[1];
        vec [k] = fibers[2];
    }

}

inline void setupFibers ( VectorEpetra& vec, std::vector<Real>& fibers)
{
	VectorSmall<3> fiberVectorSmall;
	fiberVectorSmall[0]=fibers.at(0);
	fiberVectorSmall[1]=fibers.at(1);
	fiberVectorSmall[2]=fibers.at(2);
	setupFibers(vec, fiberVectorSmall);
}

inline void setupFibers ( VectorEpetra& vec, Real fx, Real fy, Real fz)
{
	VectorSmall<3> fiberVectorSmall;
	fiberVectorSmall[0] = fx;
	fiberVectorSmall[1] = fy;
	fiberVectorSmall[2] = fz;
	setupFibers(vec, fiberVectorSmall);
}

inline void setValueOnBoundary( VectorEpetra& vec, boost::shared_ptr<  RegionMesh<LinearTetra> > fullMesh, Real value, std::vector<UInt> flags)
{

	for( UInt j (0); j < vec.epetraVector().MyLength() ; ++j )
	{
		for ( UInt k(0); k < flags.size(); k++ )
		{
			if ( fullMesh -> point ( vec.blockMap().GID (j) ).markerID() == flags.at(k) )
			{
				if ( vec.blockMap().LID ( vec.blockMap().GID (j) ) != -1 )
				{
					(vec) ( vec.blockMap().GID (j) ) = value;
				}
			}
		}
    }
}

inline void setValueOnBoundary( VectorEpetra& vec, boost::shared_ptr<  RegionMesh<LinearTetra> > fullMesh, Real value, UInt flag)
{

	for( UInt j (0); j < vec.epetraVector().MyLength() ; ++j )
	{
		if ( fullMesh -> point ( vec.blockMap().GID (j) ).markerID() == flag )
		{
			if ( vec.blockMap().LID ( vec.blockMap().GID (j) ) != -1 )
			{
				(vec) ( vec.blockMap().GID (j) ) = value;
			}
		}
    }
}

inline void rescaleVector( VectorEpetra& vec, Real min, Real max, Real scaleFactor = 1.0 )
{
	vec -= min;
	cout << "\n\nmax - min:" << max - min << "\n\n";
	//assert((max-min)==0 && "\nRescaling is going to divide by zero!\n");
	vec *= ( scaleFactor / ( max - min ) );
}

inline void rescaleVector( VectorEpetra& vec, Real scaleFactor = 1.0 )
{
	Real max = vec.maxValue();
	Real min = vec.minValue();
	rescaleVector( vec, min, max, scaleFactor);
}

} // namespace HeartUtility

} // namespace LifeV

#endif /* HEARTUTILITY_H */

