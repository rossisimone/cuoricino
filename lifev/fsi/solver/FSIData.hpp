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
    @brief FSIData - File containing a data container for FSI problems

    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @author Gilles fourestey <gilles.fourestey@epfl.ch>
    @date 10-06-2010

    @maintainer Simone Deparis <simone.deparis@epfl.ch>
 */


#ifndef FSIDATA_H
#define FSIDATA_H

#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/structure/solver/VenantKirchhoffElasticData.hpp>

#include <boost/array.hpp>
#include <boost/scoped_ptr.hpp>

namespace LifeV
{

//! FSIData - Data container for FSI problems
/*!
 *  @author Cristiano Malossi
 */

class FSIData
{
public:

    //! @name Type definitions
    //@{

    typedef OseenData                               dataFluid_Type;
    typedef boost::shared_ptr< dataFluid_Type >     dataFluidPtr_Type;

    typedef VenantKirchhoffElasticData              dataSolid_Type;
    typedef boost::shared_ptr< dataSolid_Type >     dataSolidPtr_Type;

    typedef TimeData                                timeData_Type;
    typedef boost::shared_ptr< timeData_Type >      timeDataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit FSIData();

    //! Copy constructor
    /*!
     * @param FSIData - FSIData
     */
    explicit FSIData( const FSIData& FSIData );

    //! Destructor
    virtual ~FSIData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param FSIData - FSIData
     */
    FSIData& operator=( const FSIData& FSIData );

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile - data file
     */
    void setup( const GetPot& dataFile, const std::string& section = "problem" );

    //! Display the values
    void showMe( std::ostream& output = std::cout );

    bool isMonolithic();

    //@}


    //! @name Set methods
    //@{

    //! Set data fluid container
    /*!
     * @param dataFluid shared_ptr to dataFluid container
     */
    void setDataFluid( const dataFluidPtr_Type& dataFluid ) { M_dataFluid = dataFluid; }

    //! Set data solid container
    /*!
     * @param dataFluid shared_ptr to dataSolid container
     */
    void setDataSolid( const dataSolidPtr_Type& dataSolid ) { M_dataSolid = dataSolid; }


    //! Set data solid container
    /*!
     * @param dataFluid shared_ptr to dataSolid container
     */
    void setDataALE( const timeDataPtr_Type& dataALE ) { M_dataALE = dataALE; }

    const timeDataPtr_Type& dataALE ( ) const { return M_dataALE; }

    //@}


    //! @name Get methods
    //@{

    //! Get data fluid container
    /*!
     * @return shared_ptr to dataFluid container
     */
    const dataFluidPtr_Type& dataFluid() const { return M_dataFluid; }

    //! Get data solid container
    /*!
     * @return shared_ptr to dataSolid container
     */
    const dataSolidPtr_Type& dataSolid() const { return M_dataSolid; }

    //! Get maximum number of subiterations
    /*!
     * @return maximum number of subiterations
     */
    const UInt& maxSubIterationNumber() const { return M_maxSubIterationNumber; }

    //! Get absolute tolerance
    /*!
     * @return absolute tolerance
     */
    const Real& absoluteTolerance() const { return M_absoluteTolerance; }

    //! Get relative tolerance
    /*!
     * @return relative tolerance
     */
    const Real& relativeTolerance() const { return M_relativeTolerance; }

    //! Get error tolerance
    /*!
     * @return error tolerance
     */
    const Real& errorTolerance() const { return M_errorTolerance; }

    //! Get NonLinearLineSearch
    /*!
     * @return NonLinearLineSearch
     */
    const Int& NonLinearLineSearch() const { return M_NonLinearLineSearch; }

    //! Get method type
    /*!
     * @return method type
     */
    const std::string& method() const { return M_method; }

    //! Get algorithm type
    /*!
     * @return algorithm type
     */
    const std::string& algorithm() const { return M_algorithm; }

    //! Get default omega for Aitken iterations
    /*!
     * @return default omega for Aitken iterations
     */
    const Real& defaultOmega() const { return M_defaultOmega; }

    //! Get the range of omega for Aitken iterations
    /*!
     * @return range of omega for Aitken iterations
     */
    const boost::array< Real, 2 >& OmegaRange() const { return M_rangeOmega; }

    //! Get update every
    /*!
     * If M_updateEvery == 1, normal fixedPoint algorithm
     * If M_updateEvery  > 1, recompute computational domain every M_updateEvery iterations (transpiration)
     * If M_updateEvery <= 0, recompute computational domain and matrices only at first subiteration (semi-implicit)
     *
     * @return updateEvery value
     */
    const Int& updateEvery() const { return M_updateEvery; }

    //! Get the fluid Interface Flag
    /*!
     * @return Flag of the interface  on the fluid boundary side
     */
    const Int& fluidInterfaceFlag() const { return M_fluidInterfaceFlag; }

    //! Get the structure Interface Flag
    /*!
     * @return Flag of the interface  on the structure boundary side
     */
    const Int& structureInterfaceFlag() const { return M_structureInterfaceFlag; }

    //! Get the fluid Interface Flag (for Vertices)
    /*!
     * @return Flag of the vertex on the interface on the fluid boundary side
     */
    Int const* fluidInterfaceVertexFlag() const { return M_fluidInterfaceVertexFlag.get(); }

    //! Get the fluid Interface Flag (for Vertices)
    /*!
     * @return Flag of the vertex on the interface on the structure boundary side
     */
    Int const* structureInterfaceVertexFlag() const { return M_structureInterfaceVertexFlag.get(); }

    //! Get the tolerance for the Interface identification
    /*!
     * @return the tolerance for the Interface identification
     */
    const Real& interfaceTolerance() const { return M_interfaceTolerance; }

    //! Get the timestep to restart the simulation
    /*!
     * @return the timestep used in the previous simulation from which we want to restart, used for the initialization
     of the time discretization
     */
    inline Real restartTimeStep() const { return M_restartTimeStep; }


    //@}

private:

    dataFluidPtr_Type             M_dataFluid;
    dataSolidPtr_Type             M_dataSolid;
    timeDataPtr_Type              M_dataALE;

    // Problem - Non Linear Richardson parameters
    UInt                          M_maxSubIterationNumber;
    Real                          M_absoluteTolerance;
    Real                          M_relativeTolerance;
    Real                          M_errorTolerance;
    Int                           M_NonLinearLineSearch;

    // Problem - Methods
    std::string                   M_method;
    std::string                   M_algorithm;

    // Problem - FixedPoint / EJ
    Real                          M_defaultOmega;
    boost::array< Real, 2 >       M_rangeOmega;
    Int                           M_updateEvery;

    // Interface
    Int                           M_fluidInterfaceFlag;
    Int                           M_structureInterfaceFlag;

    boost::scoped_ptr<Int const>  M_fluidInterfaceVertexFlag;
    boost::scoped_ptr<Int const>  M_structureInterfaceVertexFlag;

    Real                          M_interfaceTolerance;

    Real                          M_restartTimeStep;
};

} // end namespace LifeV

#endif // end FSIDATA_H
