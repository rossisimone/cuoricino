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

namespace LifeV
{

template<UInt s> class RosenbrockTransformed
{
public:
	RosenbrockTransformed(){};

	RosenbrockTransformed(Real g, const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
						  const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat, UInt order);

	virtual ~RosenbrockTransformed(){};

	template<typename RightHandSide>
	void solve(RightHandSide& Fun, vector<Real>& y0, Real t0, Real TF, Real& dt_init);

	template<typename RightHandSide>
	void solve(boost::shared_ptr<RightHandSide> Fun, vector<Real>& y0, Real t0, Real TF, Real& dt_init);

	template<typename RightHandSide>
	void solve(RightHandSide& Fun, VectorLU& y0, Real t0, Real TF, Real& dt_init);

	template<typename RightHandSide>
	void solve(boost::shared_ptr<RightHandSide> Fun, VectorLU& y0, Real t0, Real TF, Real& dt_init);


public:

	void initMembers();
	void setMethod(Real g, const MatrixSmall<s,s>& A, const MatrixSmall<s,s>& C, const VectorSmall<s>& gammai,
				   const VectorSmall<s>& a, const VectorSmall<s>& m, const VectorSmall<s>& mhat, UInt order);

	template<typename RightHandSide>
	void computeStages(MatrixLU& U, const VectorLU& y, VectorLU& ytmp, VectorLU& rhs, VectorLU& Utmp,
					   boost::shared_ptr<RightHandSide> Fun, Real dt, MatrixLU& Lsys, MatrixLU& Usys,
					   MatrixLU& Psys, MatrixLU& Qsys);
	bool computeError(const MatrixLU& U, VectorLU& Utmp, Real& err_n, Real& err_n_1,
					  Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected);


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


//keep the following line and don't put RosenbrockTransformed.cpp in the CMakeLists.txt
//This is to avoid linking errors due to templates
#include <lifev/core/fem/RosenbrockTransformed.cpp>

} //namespace LifeV


#endif /* ROSENBROCKTRANSFORMED_HPP_ */
