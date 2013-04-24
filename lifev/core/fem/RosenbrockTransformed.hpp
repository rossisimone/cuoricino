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
#include <math.h>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixLU.hpp>
#include <lifev/core/array/VectorLU.hpp>
#include <lifev/core/array/MatrixStandard.hpp>
#include <lifev/core/array/VectorStandard.hpp>
#include <lifev/core/algorithm/MatrixStandardSolver.hpp>


namespace LifeV
{

template<UInt s> class RosenbrockTransformed
{
public:
	RosenbrockTransformed(){};

	RosenbrockTransformed(Real g, const MatrixStandard& A, const MatrixStandard& C, const VectorStandard& gammai,
						  const VectorStandard& a, const VectorStandard& m, const VectorStandard& mhat, UInt order);

	virtual ~RosenbrockTransformed(){};

	template<typename RightHandSide>
	void solve(RightHandSide& Fun, vector<Real>& y0, Real t0, Real TF, Real& dt_init);

	template<typename RightHandSide>
	void solve(boost::shared_ptr<RightHandSide> Fun, vector<Real>& y0, Real t0, Real TF, Real& dt_init);

	template<typename RightHandSide>
	void solve(RightHandSide& Fun, VectorStandard& y0, Real t0, Real TF, Real& dt_init);

	template<typename RightHandSide>
	void solve(boost::shared_ptr<RightHandSide> Fun, VectorStandard& y, Real t0, Real TF, Real& dt);


protected:

	void initMembers();
	void setMethod(Real g, const MatrixStandard& A, const MatrixStandard& C, const VectorStandard& gammai,
				   const VectorStandard& a, const VectorStandard& m, const VectorStandard& mhat, UInt order);

	template<typename RightHandSide>
	void computeStages(MatrixStandard& U, const VectorStandard& y, VectorStandard& ytmp, VectorStandard& rhs, VectorStandard& Utmp,
					   boost::shared_ptr<RightHandSide> Fun, Real dt, MatrixStandard& Lsys, MatrixStandard& Usys,
					   MatrixStandard& Psys, MatrixStandard& Qsys);
	bool computeError(const MatrixStandard& U, VectorStandard& Utmp, Real& err_n, Real& err_n_1,
					  Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected);


	MatrixStandardSolver M_solver;
	UInt 			 M_s;
	Real 			 M_g;
	MatrixStandard   M_A;
	MatrixStandard   M_C;
	VectorStandard   M_gammai;
	VectorStandard   M_a;
	VectorStandard   M_m;
	VectorStandard   M_mdiff;
	Real 			 M_S;
	Real			 M_D;
	Real			 M_p;
	Real			 M_p_1;
	Real			 M_absTol;
	Real			 M_relTol;

};


//keep the following line and don't put RosenbrockTransformed.cpp in the CMakeLists.txt
//This is to avoid linking errors due to templates
#include <lifev/core/fem/RosenbrockTransformed.cpp>

} //namespace LifeV


#endif /* ROSENBROCKTRANSFORMED_HPP_ */
