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
     @brief This file contains the definition of the EvaluationRotateMatrix class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUATION_ROTATE_MATRIX_HPP
#define EVALUATION_ROTATE_MATRIX_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionRotateMatrix.hpp>

#include <boost/shared_ptr.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for the interpolation of a FE function
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing an interpolated FE value.

  This is the generic implementation, so representing a vectorial FE

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.

  !! NOT DEFINED YET !! (Reason: miss SimpleTensor)
 */
template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
class EvaluationRotateMatrix
{
public:

	//! @name Public Types
    //@{

    //@}


    //! @name Static constants
    //@{

    //@}


    //! @name Constructors, destructor
    //@{

    //@}


    //! @name Methods
    //@{

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //@}


private:

    //! @name Private Methods
    //@{

    //! No empty constructor
	EvaluationRotateMatrix();

	//! Copy constructor
    EvaluationRotateMatrix(const EvaluationRotateMatrix<MeshType,MapType,SpaceDim,FieldDim>& evaluation);

	//! Expression-based constructor
	explicit EvaluationRotateMatrix(const ExpressionRotateMatrix<MeshType,MapType,SpaceDim,FieldDim>& expression);

    //! Destructor
	~EvaluationRotateMatrix();

    //@}

};


//! Evaluation for the interpolation of the gradient of a FE function
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing an interpolated FE gradient.

  This is the specialized (partially) implementation representing a vector FE

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template<typename MeshType, typename MapType>
class EvaluationRotateMatrix<MeshType,MapType,3,3>
{

public:

	//! @name Public Types
    //@{

    //! Type of the value returned by this class
    typedef MatrixSmall<3,3> return_Type;

    //! Type of the FESpace to be used in this class
	typedef ETFESpace<MeshType,MapType,3,3> fespace_Type;

    //! Type of the pointer on the FESpace
	typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;

    //! Type of the vector to be used
	typedef VectorEpetra vector_Type;

	//!
	typedef VectorSmall<3> diagonalMatrix_Type;

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
    EvaluationRotateMatrix(const EvaluationRotateMatrix<MeshType,MapType,3,3>& evaluation)
	:
		M_fespace		( evaluation.M_fespace),
		M_vector		( evaluation.M_vector, Repeated),
		M_quadrature	( 0 ),
        M_currentFE		( evaluation.M_currentFE),
		M_rotatedMatrix	( evaluation.M_rotatedMatrix)
	{
		if (evaluation.M_quadrature !=0)
		{
			M_quadrature = new QuadratureRule(*(evaluation.M_quadrature));
		}
	}

	//! Expression-based constructor
        explicit EvaluationRotateMatrix(const ExpressionInterpolateGradient<MeshType,MapType,3,3>& expression)
	:
		M_fespace( expression.fespace()),
		M_vector( expression.vector(),Repeated ),
		M_quadrature(0),
		M_currentFE(M_fespace->refFE(),M_fespace->geoMap()),
		M_rotatedMatrix(0)
	{}

    //! Destructor
	~EvaluationRotateMatrix()
	{
		if (M_quadrature !=0) delete M_quadrature;
	}

    //@}


    //! @name Methods
    //@{

    //! Internal update: computes the interpolated gradients
	void update(const UInt& iElement)
	{
		zero();

		UInt totalDOF(M_vector.epetraVector().MyLength() / 3 );


//		int const nbq( M_quadrature->nbQuadPt() );
		std::vector<Real> u_x(M_quadrature->nbQuadPt(),0);
		std::vector<Real> u_y(M_quadrature->nbQuadPt(),0);
		std::vector<Real> u_z(M_quadrature->nbQuadPt(),0);
		VectorSmall< 3 > a_l;

		Real longitudinalComponent = M_matrix[0];
		Real transversalComponent  = M_matrix[1];

		M_currentFE.update(M_fespace->mesh()->element(iElement), ET_UPDATE_DPHI);
		Real nbFEDof(M_fespace->refFE().nbDof());


		for (UInt i(0); i< M_fespace->refFE().nbDof(); ++i)
		{
			for (UInt q(0); q< M_quadrature->nbQuadPt(); ++q)
			{
				UInt globalID(M_fespace->dof().localToGlobalMap(iElement,i));

               	u_x[q] += M_currentFE.phi(i,q) * M_vector[globalID];
               	u_y[q] += M_currentFE.phi(i,q) * M_vector[globalID + totalDOF];
               	u_z[q] += M_currentFE.phi(i,q) * M_vector[globalID + 2 * totalDOF];

			}

			for (UInt q(0); q< M_quadrature->nbQuadPt(); ++q)
			{
				a_l[0] = u_x[q];
				a_l[1] = u_y[q];
				a_l[2] = u_z[q];
				Real norm = std::sqrt( a_l[0] * a_l[0] + a_l[1] * a_l[1] + a_l[2] * a_l[2] );
				a_l[0] = a_l[0] / norm;
				a_l[1] = a_l[1] / norm;
				a_l[2] = a_l[2] / norm;
				for( UInt j(0); j < 3; j++){			///////  D = sigma_t * I + (sigma_l-sigma_t) * a_l * a_l^T
					M_rotatedMatrix[q][j][j] = transversalComponent;

					for( UInt k(0); k < 3; k++){
						M_rotatedMatrix[q][j][k] += (longitudinalComponent - transversalComponent ) * a_l[j] * a_l[k];
					}
				}
			}
		}
	}


    //! Erase the interpolated gradients stored internally
	void zero()
	{
	    for (UInt q(0); q<M_quadrature->nbQuadPt(); ++q)
	    {
	        for (UInt j(0); j<3; ++j) 
	        {
	            for (UInt i(0); i<3; ++i)
                    {
                        M_rotatedMatrix[q][j][i] = 0.0;
                    }
	        }
	    }
	}

    //! Show the values
	void showValues() const
	{
		std::cout << " Rotated Matrix : " << std::endl;

		for (UInt i(0); i<M_quadrature->nbQuadPt(); ++i)
		{
			std::cout << M_rotatedMatrix[i] << std::endl;
		}
	}

    //! Display method
	static void display(ostream& out=std::cout)
    {
        out << "rotated matrix [ " << 9 << " ]";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
	template< typename CFEType >
	void setGlobalCFE(const CFEType* /*globalCFE*/)
    {}

    //! Do nothing setter for the test current FE
	template< typename CFEType >
	void setTestCFE(const CFEType* /*testCFE*/)
    {}

    //! Do nothing setter for the solution current FE
	template< typename CFEType >
	void setSolutionCFE(const CFEType* /*solutionCFE*/)
    {}

    //! Setter for the quadrature rule
	void setQuadrature(const QuadratureRule& qr)
	{
		if (M_quadrature !=0) delete M_quadrature;
		M_quadrature = new QuadratureRule(qr);
		M_currentFE.setQuadratureRule(qr);
		M_rotatedMatrix.resize(qr.nbQuadPt());
	}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
	return_Type value_q(const UInt& q) const
    {
        return M_rotatedMatrix[q];
    }

    //! Getter for the value for a vector
	return_Type value_qi(const UInt& q, const UInt& /*i*/) const
    {
        return M_rotatedMatrix[q];
    }

    //! Getter for the value for a matrix
	return_Type value_qij(const UInt& q, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_rotatedMatrix[q];
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
	EvaluationRotateMatrix();

    //@}

    //! Data storage
	fespacePtr_Type M_fespace;
	vector_Type M_vector;
	diagonalMatrix_Type M_matrix;
	QuadratureRule* M_quadrature;

    //! Structure for the computations
	ETCurrentFE<3,3> M_currentFE;

    //! Storage for the temporary values
	std::vector<return_Type> M_rotatedMatrix;
};


template<typename MeshType, typename MapType >
const flag_Type
EvaluationRotateMatrix<MeshType,MapType,3,3>::
S_globalUpdateFlag=ET_UPDATE_NONE;

template<typename MeshType, typename MapType>
const flag_Type
EvaluationRotateMatrix<MeshType,MapType,3,3>::
S_testUpdateFlag=ET_UPDATE_NONE;

template<typename MeshType, typename MapType>
const flag_Type
EvaluationRotateMatrix<MeshType,MapType,3,3>::
S_solutionUpdateFlag=ET_UPDATE_NONE;


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
