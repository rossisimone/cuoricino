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
    @brief Implementation of the Rosenbrock method ROS3P

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 01-04-2013
 */


#ifndef ROS3P_HPP_
#define ROS3P_HPP_

#include <lifev/core/fem/RosenbrockTransformed.hpp>

#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>


namespace LifeV
{
	typedef LinearSolver                      solver_Type;
    typedef MatrixEpetra<Real>                matrix_Type;
    typedef boost::shared_ptr<matrix_Type>    matrixPtr_Type;
    typedef VectorEpetra                      vector_Type;
    typedef boost::shared_ptr<VectorEpetra>   vectorPtr_Type;
    typedef MapEpetra						  map_Type;

class ROS3P : public RosenbrockTransformed<3>
{
public:
	ROS3P(boost::shared_ptr<Epetra_Comm>& Comm, const string& solvParam);
/*
	template<typename RightHandSide>
	void solve(RightHandSide Fun, vector_Type& y0, Real t0, Real TF, Real dt_init);

	template<typename RightHandSide>
	void solve(boost::shared_ptr<RightHandSide>& Fun, vector_Type& y0, Real t0, Real TF, Real dt_init);
*/
};

ROS3P::ROS3P(boost::shared_ptr<Epetra_Comm>& Comm, const string& solvParam)
: RosenbrockTransformed<3>()
{
	MatrixSmall<3,3> A;
	MatrixSmall<3,3> C;
	VectorSmall<3> gammai;
	VectorSmall<3> a;
	VectorSmall<3> m;
	VectorSmall<3> mhat;
	Real g(0.7886751345948129);

	A(1,0) = 1.267949192431123;			A(2,0) = 1.267949192431123;
	C(1,0) = -1.607695154586736;		C(2,0) = -3.464101615137755;		C(2,1) = -1.732050807568877;
	gammai(0) = 0.7886751345948129;		gammai(1) = -0.2113248654051871;	gammai(2) = -1.077350269189626;
	a(0) = 0.0;							a(1) = 1.0;							a(2) = 1.0;
	m(0) = 2.0;							m(1) = 0.5773502691896258;			m(2) = 0.4226497308103742;
	mhat(0) = 2.113248654051871;		mhat(1) = 1.0;						mhat(2) = 0.4226497308103742;

	this->setMethod(g, A, C, gammai, a, m, mhat, 3);
	this->initMembers();
	this->setSolver(Comm, solvParam);
}
/*
template<typename RightHandSide>
void ROS3P::solve(RightHandSide Fun, vector_Type& y0, Real t0, Real TF, Real dt_init)
{
	cout<<"solve1...\n";
	RosenbrockTransformed<3>::solve(Fun, y0, t0, TF, dt_init);
}

template<typename RightHandSide>
void ROS3P::solve(boost::shared_ptr<RightHandSide>& Fun, vector_Type& y0, Real t0, Real TF, Real dt_init)
{
	cout<<"solve2...\n";
	RosenbrockTransformed<3>::solve(Fun, y0, t0, TF, dt_init);
}

*/





} //namesapce LifeV



#endif /* ROS3P_HPP_ */
