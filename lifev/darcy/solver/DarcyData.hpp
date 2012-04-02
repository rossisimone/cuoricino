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
     @brief This file contains the data for all the Darcy solver

     @date 05/2010
     @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

     @contributor M. Kern <michel.kern@inria.fr>
     @maintainer M. Kern <michel.kern@inria.fr>
 */

#ifndef _DATADARCY_H_
#define _DATADARCY_H_ 1

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/Assembly.hpp>

// LifeV namespace
namespace LifeV
{
/*!
  @class DarcyData

  This class contain the basic data for the Darcy solver. In particoular it stores...
  @todo class not finished!
 */
template <typename Mesh>
class DarcyData
{
public:

    // Policies
    //! @name Policies
    //@{

    typedef GetPot                          Data_Type;
    typedef boost::shared_ptr< Data_Type >  Data_ptrType;

    typedef TimeData                        Time_Type;
    typedef boost::shared_ptr< Time_Type >  Time_ptrType;

    typedef MeshData                        Mesh_Type;
    typedef boost::shared_ptr< Mesh_Type >  Mesh_ptrType;

    //@}

    // Constructors.
    //! @name Constructors
    //@{

    //! Empty Constructor
    DarcyData();

    /*!
    Constructor using a data file.
      @param dataFile GetPot data file for setup the problem
      @param section the section for the Darcy data
    */
    DarcyData( const GetPot& dataFile, const std::string& section = "darcy" );

    /*!
    Copy constructor.
      @param darcyData object to take a copy
    */
    DarcyData( const DarcyData &darcyData );

    //@}

    // Methods.
    //! @name Methods
    //@{

    /*! Overloading of the operator =
        @param darcyData The DarcyData to be copied.
    */
    DarcyData& operator=( const DarcyData& darcyData );

    /*! External setup
        @param dataFile The data file with all the data.
        @param section The global section.
    */
    void setup( const Data_Type& dataFile, const std::string& section = "darcy"  );

    //@}

    // Set methods
    //! @name Set methods
    //@{

    /*! Set data time container
        @param TimeData Boost shared_ptr to TimeData container
    */
    inline void setTimeData( const Time_ptrType TimeData )
    {
        M_time = TimeData;
    }

    /*! Set mesh container
        @param MeshData Boost shared_ptr to meshData container
    */
    inline void setMeshData( const Mesh_ptrType MeshData )
    {
        M_mesh = MeshData;
    }

    // Get methods.
    //! @name Get methods
    //@{

    //! Get the level of verbosity of the problem.
    inline UInt verbose( void ) const
    {
        return M_verbose;
    }

    //! Get the main section of the data file.
    inline std::string section( void ) const
    {
        return M_section;
    }

    //! Get the data file of the problem.
    inline Data_ptrType dataFile( void ) const
    {
        return M_data;
    }

    /*! Get data time container.
        @return shared_ptr to TimeData container
    */
    inline Time_ptrType dataTime( void ) const
    {
        return M_time;
    }

    /*! Get mesh container
       @return shared_ptr to meshData container
    */
    inline Mesh_ptrType meshData( void ) const
    {
        return M_mesh;
    }

    //@}


private:

    //! Data containers for time and mesh
    Data_ptrType      M_data;
    Time_ptrType      M_time;
    Mesh_ptrType      M_mesh;

    //! Miscellaneous
    UInt              M_verbose;
    std::string       M_section;

};

// ===================================================
// Constructors
// ===================================================

template < typename Mesh >
DarcyData<Mesh>::DarcyData( ):
        // Data containers
        M_data          ( ),
        M_time          ( ),
        M_mesh          ( ),
        // Miscellaneous
        M_verbose       ( static_cast<UInt>(0) ),
        M_section       ( )
{

}

// Copy constructor
template < typename Mesh >
DarcyData<Mesh>::DarcyData( const DarcyData &darcyData ):
        // Data containers
        M_data                ( darcyData.M_data ),
        M_time                ( darcyData.M_time ),
        M_mesh                ( darcyData.M_mesh ),
        // Miscellaneous
        M_verbose             ( darcyData.M_verbose ),
        M_section             ( darcyData.M_section )
{

}

// Overloading of the operator =
template < typename Mesh >
DarcyData<Mesh>&
DarcyData<Mesh>::operator=( const DarcyData& darcyData )
{
    // Avoid auto-copy
    if ( this != &darcyData )
    {
        // Data containers
        M_data           = darcyData.M_data;
        M_time           = darcyData.M_time;
        M_mesh           = darcyData.M_mesh;
        // Mescellaneous
        M_verbose       = darcyData.M_verbose;
    }

    return *this;

}


// External set up method
template < typename Mesh >
void DarcyData<Mesh>::setup( const Data_Type& dataFile, const std::string& section )
{
    M_section = section;

    // If data has not been set
    if ( !M_data.get() )
        M_data.reset( new Data_Type( dataFile ) );

    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, M_section + "/time_discretization" ) );

    // If data mesh has not been set
    if ( !M_mesh.get() )
        M_mesh.reset( new Mesh_Type( dataFile, M_section + "/space_discretization" ) );

    // Miscellaneous
    M_verbose      = dataFile( ( M_section + "/miscellaneous/verbose" ).data(), 1 );
}


/////////////////////////////////////////////////////////////////////////////////

template< typename Mesh, typename SolverType = LifeV::SolverAztecOO >
class inversePermeability
{

public:

    // Policies.
    //! @name Policies
    //@{

    typedef boost::function<Matrix ( const Real&, const Real&,
                                     const Real&, const Real&,
                                     const std::vector<Real> & )> permeability_Type;

    typedef typename SolverType::vector_type      vector_Type;
    typedef boost::shared_ptr<vector_Type>        vectorPtr_Type;

    //@}

    // Constructors and destructor.
    //! @name Constructors and destructor
    //@{

    // Copy constructor
    inversePermeability ( const permeability_Type& invPerm, FESpace<Mesh, MapEpetra>& fESpace ):
	M_fields              ( std::vector< const vectorPtr_Type* >(0) ),
	M_inversePermeability ( invPerm ),
	M_fESpace             ( fESpace )
    {
    };

    //! Virtual destructor.
    // virtual ~inversePermeability ();

    //@}

    // Set methods
    //! @name Set methods
    //@{

    // Add one field
    inline void setField ( const vectorPtr_Type & field ) { M_fields.push_back( &field ); };

    inline void setFunction ( const permeability_Type & invPerm ) { M_inversePermeability = invPerm; };

    //@}

    // Get methods
    //! @name Get methods
    //@{

    Matrix operator() ( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& iElem );

    //@}

private:

    // Vector of pointers for the dependences of the permeability to an external field.
    std::vector< const vectorPtr_Type* > M_fields;

    // Inverse permeability function
    permeability_Type             M_inversePermeability;

    // Finite element space
    FESpace<Mesh, MapEpetra>&     M_fESpace;

};

template < typename Mesh, typename SolverType >
Matrix
inversePermeability < Mesh, SolverType >::
operator() ( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& iElem )
{
    std::vector<Real> values ( M_fields.size(), 0 );
    VectorElemental value ( M_fESpace.refFE().nbDof(), 1 );

    // Update the value of the current element
    M_fESpace.fe().update( M_fESpace.mesh()->element( iElem ),
                           UPDATE_QUAD_NODES | UPDATE_WDET );

    for ( UInt i(static_cast<UInt>(0)); i < values.size(); ++i )
    {
        extract_vec ( *( *(M_fields)[i] ),
                      value,
                      M_fESpace.refFE(),
                      M_fESpace.dof(),
                      M_fESpace.fe().currentLocalId(), 0 );

        values[i] = value[0];

    }

    return M_inversePermeability ( t, x, y, z, values );
}

}
#endif // _DATADARCY_H_

// -*- mode: c++ -*-
