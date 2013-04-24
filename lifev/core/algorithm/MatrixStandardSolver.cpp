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
    @brief Implementation of LU for MatrixStandard

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 22-04-2013
 */

#include <lifev/core/algorithm/MatrixStandardSolver.hpp>


namespace LifeV
{

void MatrixStandardSolver::LU(MatrixStandard& A, MatrixStandard& P, MatrixStandard& Q, MatrixStandard& L, MatrixStandard& U) const
{
	MatrixStandard I(A.size());
	LU(A, P, Q, L, U, I);
}

void MatrixStandardSolver::LU(MatrixStandard& A, MatrixStandard& P, MatrixStandard& Q,
							  MatrixStandard& L, MatrixStandard& U, const MatrixStandard& I) const
{
	UInt n = I.size();
	MatrixStandard Minv(I), Pk(I), Qk(I), Mk(I), Mkinv(I), I2(I*2.0), TMP(I);
	P = I;
	Q = I;

	for(UInt k = 1; k<n; k++)
	{
		Pk = I;
		Qk = I;
		Pivot(A, Pk, Qk, k);
		Pk.times(A, TMP);
		TMP.times(Qk, A);
		Mk = I;

		for(UInt i = k; i<n; i++)
			Mk[i][k-1] = -A[i][k-1]/A[k-1][k-1];

		Mkinv = I2;
		Mkinv -= Mk;
		A = Mk*A;
		P = Pk*P;
		Q = Q*Qk;
		Minv.times(Pk, TMP);
		TMP.times(Mkinv, Minv);
	}

	A.triU(U);
	L = P*Minv;
}

void MatrixStandardSolver::Pivot(const MatrixStandard& A, MatrixStandard& Pk, MatrixStandard& Qk, UInt k) const
{
	UInt i = k-1;
	UInt j = k-1;
	UInt n = A.size();
	Real max = std::abs(A[i][j]);

	for(UInt l = k-1; l<n; l++)
	{
		for(UInt c = k-1; c<n; c++)
		{
			if( std::abs(A[l][c]) > max)
			{
				max = std::abs(A[l][c]);
				i = l;
				j = c;
			}
		}
	}

	Pk[i][i] = 0.;	Pk[k-1][k-1] = 0.;	Pk[k-1][i] = 1.;	Pk[i][k-1] = 1.;
	Qk[j][j] = 0.;	Qk[k-1][k-1] = 0.;	Qk[k-1][j] = 1.;	Qk[j][k-1] = 1.;
}

void MatrixStandardSolver::solveU(const MatrixStandard& A, VectorStandard& b) const
{
	UInt n = b.size();

	b[n-1] = b[n-1]/A[n-1][n-1];

	for(UInt i=2; i<=n; i++)
	{
		for(UInt k = n-i+1; k<=n-1; k++)
			b[n-i] -= A[n-i][k]*b[k];

		b[n-i] /= A[n-i][n-i];
	}
}

void MatrixStandardSolver::solveL(const MatrixStandard& A, VectorStandard& b) const
{
	UInt n = b.size();

	b[0] = b[0]/A[0][0];

	for(UInt i=1; i<n; i++)
	{
		for(UInt k = 0; k<=i-1; k++)
			b[i] -= A[i][k]*b[k];

		b[i] /= A[i][i];
	}
}

} //namespace LifeV
