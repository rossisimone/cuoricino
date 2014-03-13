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
    @brief Rosenbrock method with change of variables

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

class RosenbrockTransformed
{
public:
    RosenbrockTransformed() {};

    RosenbrockTransformed (Real g, const MatrixStandard& A, const MatrixStandard& C, const VectorStandard& gammai,
                           const VectorStandard& a, const VectorStandard& m, const VectorStandard& mhat, UInt order);

    virtual ~RosenbrockTransformed() {};

    template<typename RightHandSide>
    void solve (RightHandSide& Fun, std::vector<Real>& y0, Real t0, Real TF, Real& dt_init);

    template<typename RightHandSide>
    void solve (boost::shared_ptr<RightHandSide> Fun, std::vector<Real>& y0, Real t0, Real TF, Real& dt_init);

    template<typename RightHandSide>
    void solve (RightHandSide& Fun, VectorStandard& y0, Real t0, Real TF, Real& dt_init);

    template<typename RightHandSide>
    void solve (boost::shared_ptr<RightHandSide> Fun, VectorStandard& y, Real t0, Real TF, Real& dt);


protected:

    void initMembers();
    void setMethod (Real g, const MatrixStandard& A, const MatrixStandard& C, const VectorStandard& gammai,
                    const VectorStandard& a, const VectorStandard& m, const VectorStandard& mhat, UInt order);

    template<typename RightHandSide>
    void computeStages (MatrixStandard& U, const VectorStandard& y, VectorStandard& ytmp, VectorStandard& rhs, VectorStandard& Utmp,
                        boost::shared_ptr<RightHandSide> Fun, Real dt, MatrixStandard& Lsys, MatrixStandard& Usys,
                        MatrixStandard& Psys, MatrixStandard& Qsys);
    bool computeError (const MatrixStandard& U, VectorStandard& Utmp, Real& err_n, Real& err_n_1,
                       Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected);


    MatrixStandardSolver M_solver;
    UInt             M_s;
    Real             M_g;
    MatrixStandard   M_A;
    MatrixStandard   M_C;
    VectorStandard   M_gammai;
    VectorStandard   M_a;
    VectorStandard   M_m;
    VectorStandard   M_mdiff;
    Real             M_S;
    Real             M_D;
    Real             M_p;
    Real             M_p_1;
    Real             M_absTol;
    Real             M_relTol;

};


template<typename RightHandSide>
void RosenbrockTransformed::solve (RightHandSide& Fun, VectorStandard& y0, Real t0, Real TF, Real& dt_init)
{
    boost::shared_ptr<RightHandSide> FunPtr (new RightHandSide (Fun) );

    solve (FunPtr, y0, t0, TF, dt_init);
}

template<typename RightHandSide>
void RosenbrockTransformed::solve (RightHandSide& Fun, std::vector<Real>& y0, Real t0, Real TF, Real& dt_init)
{
    boost::shared_ptr<RightHandSide> FunPtr (new RightHandSide (Fun) );

    VectorStandard y0LU (y0);
    solve (FunPtr, y0LU, t0, TF, dt_init);
    y0 = y0LU;
}

template<typename RightHandSide>
void RosenbrockTransformed::solve ( boost::shared_ptr<RightHandSide> Fun, std::vector<Real>& y0, Real t0, Real TF, Real& dt_init)
{
    VectorStandard y0LU (y0);
    solve (Fun, y0LU, t0, TF, dt_init);
    y0 = y0LU;
}

template<typename RightHandSide>
void RosenbrockTransformed::solve ( boost::shared_ptr<RightHandSide> Fun, VectorStandard& y, Real t0, Real TF, Real& dt)
{

    Real t (t0);                        // time t_k
    Real dt_old (dt);

    UInt n = y.size();
    MatrixStandard U (n, M_s);
    MatrixStandard I (n);
    MatrixStandard Psys (n, n), Qsys (n, n), Usys (n, n), Lsys (n, n);
    MatrixStandard B (n);               // Linear system matrix
    VectorStandard ytmp (y);            // temporary variable
    VectorStandard Utmp (y);
    VectorStandard rhs (y);             // rhs will be the right hand side
    Real err_n;                         // error at step n
    Real err_n_1;                       // error at step n-1
    Real fac_max = 5.0;                 // maximal value for this factor, dt(k+1) < dt(k)*fac_max
    Int k = 1;                          // iteration counter
    bool rejected = false;              // used to know if a step is rejected two times consecutively

    // First step, to set err_n_1

    B = I / (dt * M_g);
    B -= Fun->getJac (y);
    M_solver.LU (B, Psys, Qsys, Lsys, Usys, I);
    computeStages<RightHandSide> (U, y, ytmp, rhs, Utmp, Fun, dt, Lsys, Usys, Psys, Qsys);

    U.times (M_m, ytmp);
    y += ytmp ;
    U.times (M_mdiff, Utmp);
    err_n_1 = Utmp.norm2();
    t += dt;

    while (t < TF)
    {

        U *= 0.0;
        B = I / (dt * M_g);
        B -=  Fun->getJac (y);
        M_solver.LU (B, Psys, Qsys, Lsys, Usys, I);// Computing the inverse, which will be used s times
        computeStages<RightHandSide> (U, y, ytmp, rhs, Utmp, Fun, dt, Lsys, Usys, Psys, Qsys);

        if ( computeError (U, Utmp, err_n, err_n_1, fac_max, dt, dt_old, TF - t, y.norm2(), rejected) )
        {
            rejected = true;
            continue;
        }
        else
        {
            U.times (M_m, ytmp);
            y += ytmp;
            t += dt;                    // updating the time
            k++;

            rejected = false;
        }

    }

}

template<typename RightHandSide>
void RosenbrockTransformed::computeStages (MatrixStandard& U, const VectorStandard& y, VectorStandard& ytmp, VectorStandard& rhs, VectorStandard& Utmp,
                                           boost::shared_ptr<RightHandSide> Fun, Real dt, MatrixStandard& Lsys, MatrixStandard& Usys,
                                           MatrixStandard& Psys, MatrixStandard& Qsys)
{
    for (UInt i = 0; i < M_s; i++)
    {
        U.times (M_A.getLine (i), Utmp);
        ytmp = y + Utmp;                // ytmp = y0 + sum_{j=1}^{i-1} A(i,j) * U(:,j)
        U.times (M_C.getLine (i), Utmp);// Utmp = sum_{j=1}^{i-1} C(i,j) * U(:,j) / dt
        Utmp /= dt;
        Fun->computeRhs ( ytmp, rhs);
        Utmp += rhs;
        Psys.times (Utmp, rhs);
        M_solver.solveL (Lsys, rhs);
        M_solver.solveU (Usys, rhs);
        Qsys.times (rhs, Utmp);
        U.setCol (i, Utmp);
    }
}


} //namespace LifeV


#endif /* ROSENBROCKTRANSFORMED_HPP_ */
