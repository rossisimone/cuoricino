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
    @brief Rosenbrock method with variable change

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 01-04-2013
 */

#ifndef ROSENBROCKTRANSFORMED_HPP_
#define ROSENBROCKTRANSFORMED_HPP_

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <fstream>
#include <string>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
	typedef LinearSolver                      solver_Type;
    typedef MatrixEpetra<Real>                matrix_Type;
    typedef boost::shared_ptr<matrix_Type>    matrixPtr_Type;
    typedef VectorEpetra                      vector_Type;
    typedef boost::shared_ptr<VectorEpetra>   vectorPtr_Type;
    typedef MapEpetra						  map_Type;

template<UInt s> class RosenbrockTransformed
{
public:
	RosenbrockTransformed(){};

	RosenbrockTransformed(Real g, const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
						  const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat,
						  UInt order, boost::shared_ptr<Epetra_Comm>& Comm, const string& solvParam);

	virtual ~RosenbrockTransformed(){};

	template<typename RightHandSide>
	void solve(RightHandSide& Fun, vector_Type& y0, Real t0, Real TF, Real dt_init);

	template<typename RightHandSide>
	void solve(boost::shared_ptr<RightHandSide> Fun, vector_Type& y0, Real t0, Real TF, Real dt_init);

	template<typename RightHandSide>
	void solve(boost::shared_ptr<RightHandSide> Fun, vector<Real>& y0, Real t0, Real TF, Real dt_init);


public:

	matrix_Type getIdentity(const map_Type& map, Int n);
	vector_Type combLin(const vector<vector_Type>& U, const MatrixSmall<s,s>& M, const Int line);
	vector_Type combLin(const vector<vector_Type>& U, const VectorSmall<s>& v);
	void zero(vector<vector_Type>& U);
	void updateMatrix(matrixPtr_Type& B, const matrix_Type& J, const matrix_Type& aI);

	template<typename RightHandSide>
	void computeStages(vector<vector_Type>& U, const vector_Type& y,vector_Type& ytmp, vectorPtr_Type& rhs,
					   vectorPtr_Type& Utmp, boost::shared_ptr<RightHandSide> Fun, Real dt);
	bool computeError(const vector<vector_Type>& U, vectorPtr_Type& Utmp, Real& err_n, Real& err_n_1,
					  Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected);
	void initMembers();
	void setSolver(boost::shared_ptr<Epetra_Comm>& Comm, const string& solvParam);
	void setMethod(Real g, const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
				   const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat, UInt order);

	//LU stuff
	void LU(vector<vector<Real> >& A, vector<vector<Real> >& P, vector<vector<Real> >& Q,
			vector<vector<Real> >& L, vector<vector<Real> >& U, UInt n);
	vector<vector<Real> > getIdentity(UInt n);
	void pivot(const vector<vector<Real> >& A, UInt k, UInt n, vector<vector<Real> >& Pk, vector<vector<Real> >& Qk);
	void MGauss(const vector<vector<Real> >& A, UInt k, UInt n, vector<vector<Real> >& Mk);
	void mult(vector<vector<Real> >& B, Real a, UInt n);
	vector< vector<Real> > mult(const vector<vector<Real> >& A, const vector<vector<Real> >& B, UInt n);
	vector<Real> mult(const vector<vector<Real> >& A, const vector<Real>& b, UInt n);
	void minusequal(vector<vector<Real> >& A, const vector<vector<Real> >& B, UInt n);
	vector<vector<Real> > triu(const vector<vector<Real> >& B, UInt n);
	void zero(vector<vector<Real> >& U, UInt n);
	vector<Real> combLin(const vector<vector<Real> >& U, const MatrixSmall<s,s>& M, const Int line, UInt n);
	vector<Real> combLin(const vector<vector<Real> >& U, const VectorSmall<s>& v, UInt n);
	vector<Real> sum(const vector<Real>& a, const vector<Real>& b, UInt n);
	void multequal(vector<Real>& v, Real a, UInt n);
	vector<Real> solvel(const vector<vector<Real> >& L, const vector<Real>& b, UInt n);
	vector<Real> solveu(const vector<vector<Real> >& U, const vector<Real>& b, UInt n);
	Real norm2(const vector<Real>& v, UInt n);
	template<typename RightHandSide>
	void computeStages(vector< vector<Real> >& U, const vector<Real>& y, vector<Real>& ytmp,
			 vector<Real>& rhs, vector<Real>& Utmp, boost::shared_ptr<RightHandSide> Fun,
			 Real dt, vector< vector<Real> >& Lsys, vector< vector<Real> >& Usys,
			 vector< vector<Real> >& Psys, vector< vector<Real> >& Qsys, UInt n);
	bool computeError(const vector<vector<Real> >& U, vector<Real>& Utmp, Real& err_n, Real& err_n_1,
					  Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected, UInt n);
	void disp(vector<vector<Real> >& M, string name, UInt n);


	solver_Type 	 M_solver;
	Real 			 M_g;
	MatrixSmall<s,s> M_A;
	MatrixSmall<s,s> M_C;
	VectorSmall<s>   M_gammai;
	VectorSmall<s>   M_a;
	VectorSmall<s>   M_m;
	VectorSmall<s>   M_mhat;
	Real 			 M_S;
	Real			 M_D;
	Real			 M_p;
	Real			 M_p_1;
	Real			 M_absTol;
	Real			 M_relTol;

};


//KEEP LIKE THAT
//keep the following line and don't put RosenbrockTransformed.cpp in the CMakeLists.txt
//This is to avoid linking errors due to templates
#include <lifev/core/fem/RosenbrockTransformed.cpp>

} //namespace LifeV


#endif /* ROSENBROCKTRANSFORMED_HPP_ */
