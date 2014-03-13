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

#include <lifev/core/fem/RosenbrockTransformed.hpp>


namespace LifeV
{

RosenbrockTransformed::RosenbrockTransformed (Real g, const MatrixStandard& A, const MatrixStandard& C, const VectorStandard& gammai,
                                              const VectorStandard& a, const VectorStandard& m, const VectorStandard& mhat,
                                              UInt order)
    : M_g (g), M_A (A), M_C (C), M_gammai (gammai), M_a (a), M_m (m), M_mdiff (m - mhat), M_p (order)
{
    initMembers();
}


void RosenbrockTransformed::initMembers()
{
    M_s = M_m.size();
    M_D = 1.5;
    M_S = 0.95;
    M_absTol = 0.0000001;
    M_relTol = 0.0000001;
    M_p_1 = 1.0 / M_p;
}

void RosenbrockTransformed::setMethod (Real g, const MatrixStandard& A, const MatrixStandard& C, const VectorStandard& gammai,
                                       const VectorStandard& a, const VectorStandard& m, const VectorStandard& mhat, UInt order)
{
    M_g = g;
    M_A = A;
    M_C = C;
    M_gammai = gammai;
    M_a = a;
    M_m = m;
    M_mdiff = m - mhat;
    M_p = (double) (order);
}


bool RosenbrockTransformed::computeError (const MatrixStandard& U, VectorStandard& Utmp, Real& err_n, Real& err_n_1,
                                          Real fac_max, Real& dt, Real& dt_old, Real Trem, Real ynorm, bool& rejected)
{
    Real Tol = M_absTol + M_relTol * ynorm; // Tol = atol + rtol*|y_k|
    Real fac;
    U.times (M_mdiff, Utmp);                // difference with the embedded method
    err_n = Utmp.norm2();                   // norm of the error
    if (err_n == 0.0)                       // here we set fac, where dt(k+1) = fac*dt(k)
    {
        fac = fac_max;                      // if the actual error is zero we set fac to its maximal value
    }
    else if (err_n_1 == 0.0)                // if the previous error was zero and the actual is not then fac~1
    {
        fac = M_S;
    }
    else                                    // formula to compute fac, takes in account Tol, errors and time steps
    {
        fac = M_S * std::pow ( (Tol * err_n_1) / (err_n * err_n) , M_p_1 ) * ( dt / dt_old ) ;
    }

    if (err_n > Tol)                        // the step is rejected
    {
        if (rejected)                       // if it is the second time that it is rejected we
        {
            dt /= M_D;                      // divide the time step by M_D
        }
        else                                // else dt = dt*fac, if the previous step has not grown too much then
        {
            dt *= fac;                      // fac<1, if fac>1 it will be rejected one more time and dt will be divided by 10.
        }

        return true;                        // this timestep has been rejected
    }
    else
    {
        err_n_1 = err_n;
        dt_old = dt;
        dt = std::min<Real> ( Trem, std::min<Real> ( fac, fac_max ) * dt ); // dt(k+1) = min( TF-t, fac*dt(k), fac_max*dt(k) )
        return false;
    }
}

} //namespace LifeV
