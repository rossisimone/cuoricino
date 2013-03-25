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
    @brief File for the definition of the expression used to interpolate FE functions.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_IF_CROSSED_HPP
#define EXPRESSION_IF_CROSSED_HPP

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionInterpolateValue Class representing an interpolation in an expression.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

*/
template<typename MeshType, typename MapType, UInt SpaceDim>
class ExpressionIfCrossed
    : public ExpressionBase<ExpressionIfCrossed < MeshType, MapType, SpaceDim > >
{
public:

    //! @name Public Types
    //@{

    // Base class, used only to make the code cleaner
    typedef ExpressionBase<ExpressionIfCrossed < MeshType, MapType, SpaceDim > > base_Type;

    //! Type of the finite element space
    typedef ETFESpace<MeshType, MapType, SpaceDim, 1> fespace_Type;

    //! Type for the pointer on the finite element space
    typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;

    //! Data vector type
    typedef VectorEpetra vector_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the finite element space and the data vector
    ExpressionIfCrossed (fespacePtr_Type fespace, const vector_Type& vector)
        : base_Type(), M_fespace (fespace), M_vector (vector) {}

    //! Copy constructor
    ExpressionIfCrossed (const ExpressionIfCrossed<MeshType, MapType, SpaceDim>& expr)
        : base_Type(), M_fespace (expr.M_fespace), M_vector (expr.M_vector) {}

    //! Destructor
    ~ExpressionIfCrossed() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "ifCrossed";
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the finite element space
    fespacePtr_Type fespace() const
    {
        return M_fespace;
    }

    //! Getter for the data vector
    const vector_Type vector() const
    {
        return M_vector;
    }

    // @}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    ExpressionIfCrossed();

    //@}

    // Storage for the finite element space
    fespacePtr_Type M_fespace;

    // Storage for the data vector
    vector_Type M_vector;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>MeshType</i>: The type of the mesh stored in the finite element space

  <i>MapType</i>: The type of map used in the finite element space

  <i>SpaceDim</i>: The dimension of the ambient space.

  <b> Template requirements </b>

  <i>MeshType</i>: Same as in LifeV::ETFESpace

  <i>MapType</i>: Same as in LifeV::ETFESpace

*/
template<typename MeshType, typename MapType, UInt SpaceDim>
inline ExpressionIfCrossed<MeshType, MapType, SpaceDim>
ifCrossed (
    boost::shared_ptr< ETFESpace<MeshType, MapType, SpaceDim, 1> > fespace,
    const VectorEpetra& vector
)
{
    return ExpressionIfCrossed<MeshType, MapType, SpaceDim> (fespace, vector);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
