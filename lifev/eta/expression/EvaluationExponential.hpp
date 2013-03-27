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
     @brief This file contains the definition of the EvaluationAddition class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUTATION_EXPONENTIAL_HPP
#define EVALUTATION_EXPONENTIAL_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/array/OperationSmallExponential.hpp>

#include <lifev/eta/expression/ExpressionExponential.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for the product of two other Evaluations
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing a product operation during the assembly

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <typename EvaluationBaseType>
class EvaluationExponential
{
public:

    //! @name Public Types
    //@{

    // //! Type of the value returned by the left operand
    typedef typename EvaluationBaseType::return_Type BaseReturn_Type;

    //! Type of the value returned by this class
    typedef typename OperationSmallExponential<BaseReturn_Type>::result_Type return_Type;
    //@}


    //! @name Static constants
    //@{

    //! Flag for the global current FE
    const static flag_Type S_globalUpdateFlag;

    //! Flag for the test current FE
    const static flag_Type S_testUpdateFlag;

    //! Flag for the solution current FE
    const static flag_Type S_solutionUpdateFlag;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Copy constructor
    EvaluationExponential (const EvaluationExponential& eval)
        : M_evaluationBase (eval.M_evaluationBase)
    {}

    //! Constructor from the corresponding expression
    template <typename BaseExpressionType>
    explicit EvaluationExponential (const ExpressionExponential<BaseExpressionType>& expression)
        : M_evaluationBase (expression.base() )
    {}

    //! Destructor
    ~EvaluationExponential() {}

    //@}


    //! @name Methods
    //@{

    //! Internal update method
    void update (const UInt& iElement)
    {
        M_evaluationBase.update (iElement);
    }

    //! Display method
    static void display (ostream& out = std::cout )
    {
        out << "Exponential";
        EvaluationBaseType::display (out);
    }

    //@}


    //! @name Set Methods
    //@{

    //! Setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* globalCFE)
    {
        M_evaluationBase.setGlobalCFE (globalCFE);
    }

    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        M_evaluationBase.setTestCFE (testCFE);
    }

    //! Setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_evaluationBase.setSolutionCFE (solutionCFE);
    }

    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        M_evaluationBase.setQuadrature (qr);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter a value
    return_Type value_q (const UInt& q) const
    {
        return std::exp (M_evaluationBase.value_q (q) );
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& i) const
    {
        return std::exp (M_evaluationBase.value_qi (q, i) );
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& i, const UInt& j) const
    {
        return std::exp (M_evaluationBase.value_qij (q, i, j) );
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationExponential();

    //@}

    //! Internal storage
    EvaluationBaseType M_evaluationBase;
};

template< typename EvaluationBaseType>
const flag_Type EvaluationExponential<EvaluationBaseType>::S_globalUpdateFlag
    = EvaluationBaseType::S_globalUpdateFlag;

template< typename EvaluationBaseType>
const flag_Type EvaluationExponential<EvaluationBaseType>::S_testUpdateFlag
    = EvaluationBaseType::S_testUpdateFlag;

template< typename EvaluationBaseType>
const flag_Type EvaluationExponential<EvaluationBaseType>::S_solutionUpdateFlag
    = EvaluationBaseType::S_solutionUpdateFlag;

} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
