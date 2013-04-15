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
 *  @brief File containing the BCInterfaceFunctionParserFileSolver class
 *
 *  @date 26-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionParserFileSolver_H
#define BCInterfaceFunctionParserFileSolver_H 1

#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserFile.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserSolver.hpp>

namespace LifeV
{

//! BCInterfaceFunctionParserFileSolver - LifeV boundary condition function file wrapper for \c BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between the \c BCInterface, the \c Parser, and and a general
 *  LifeV physical solver (such as \c OseenSolver or \c FSISolver). It allows to construct LifeV
 *  functions type for boundary conditions, using a \c GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain physical solver parameters.
 *
 *  See \c BCInterfaceFunctionParser, \c BCInterfaceFunctionParserFile, and \c BCInterfaceFunctionParserSolver classes for more details.
 */
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFunctionParserFileSolver: public virtual BCInterfaceFunctionParserFile< BcHandlerType, PhysicalSolverType > ,
    public virtual BCInterfaceFunctionParserSolver< BcHandlerType, PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandlerType                                                          bcHandler_Type;
    typedef PhysicalSolverType                                                     physicalSolver_Type;

    typedef BCInterfaceFunction< bcHandler_Type, physicalSolver_Type >             function_Type;
    typedef BCInterfaceFunctionParser< bcHandler_Type, physicalSolver_Type >       functionParser_Type;
    typedef BCInterfaceFunctionParserFile< bcHandler_Type, physicalSolver_Type >   functionParserFile_Type;
    typedef BCInterfaceFunctionParserSolver< bcHandler_Type, physicalSolver_Type > functionParserSolver_Type;

    typedef typename function_Type::data_Type                                      data_Type;
    typedef typename function_Type::dataPtr_Type                                   dataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionParserFileSolver();

    //! Destructor
    virtual ~BCInterfaceFunctionParserFileSolver() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data for boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const dataPtr_Type& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionParserFileSolver ( const BCInterfaceFunctionParserFileSolver& function );

    BCInterfaceFunctionParserFileSolver& operator= ( const BCInterfaceFunctionParserFileSolver& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename BcHandlerType, typename PhysicalSolverType >
inline BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >* createBCInterfaceFunctionParserFileSolver()
{
    return new BCInterfaceFunctionParserFileSolver< BcHandlerType, PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
BCInterfaceFunctionParserFileSolver< BcHandlerType, PhysicalSolverType >::BCInterfaceFunctionParserFileSolver() :
    function_Type                (),
    functionParser_Type          (),
    functionParserFile_Type      (),
    functionParserSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5024 ) << "BCInterfaceFunctionFileSolver::BCInterfaceFunctionFileSolver()" << "\n";
#endif

}

template< typename BcHandlerType, typename PhysicalSolverType >
void
BCInterfaceFunctionParserFileSolver< BcHandlerType, PhysicalSolverType >::setData ( const dataPtr_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5024 ) << "BCInterfaceFunctionFileSolver::setData" << "\n";
#endif

    functionParserFile_Type::setData ( data );

    functionParserSolver_Type::M_boundaryID = data->boundaryID();

    functionParserSolver_Type::createAccessList ( data );
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionParserFileSolver_H */
