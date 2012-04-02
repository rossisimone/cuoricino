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
 *  @file
 *  @brief File containing the boundary conditions for the Monolithic Test
 *
 *  @date 2009-04-09
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @maintainer Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 */

#ifndef __FLOWCONDITIONS_HPP
#define __FLOWCONDITIONS_HPP

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_SerialDenseVector.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV includes
#include <life/lifecore/LifeV.hpp>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSIOperator.hpp>

namespace LifeV
{

class FlowConditions
{
public:

    FlowConditions();

    static void setParamsFromGetPot ( const GetPot& dataFile );

    static void initParameters      ( FSIOperator& oper, const int& outflowFlag);

    static void renewParameters     ( FSISolver&   oper, const int& outflowFlag);

    static void renewLumpedParameters  (const int&    Flag, const Real& flux);
    static Real fextvessel             (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real force_cardium          (const Real& t);

    static Real fZero               (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    static Real inPressure          (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    static Real outPressure0         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure1         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure2         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure3         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure4         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure5         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure6         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    Real inDeltaRadius          (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    Real outDeltaRadius          (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    static Real valuePout           ( );
    static Real valuePin            ( );

private:
    static Real pi;

    static Real mu;
    static Real nu;
    static Real rho;
    static Real dt;
    static Real Pin;
    static Real Pout;
    static Real Rext_d; //external distal resistance
    static Real Rext_p; //external proximal resistance
    static Real Cp; //capacity for the pressure term in Windkessel
    static int BDType; //boundary condition type (explicit Resistance, Windkessel RC, RCR)
    static Real Flux;
    static Real Flux_old;
    static Real fExt;
    static Real tAppl;
    static Real periode;

    static bool bcOnFluid;

    static Real M_outflux;
    static Real M_influx;
    static Real M_outP;
    static Real M_inP;

    static Real M_area0;
    static Real M_inRadius0;
    static Real M_outRadius0;
    static Real M_inDeltaRadius;
    static Real M_outDeltaRadius;

    static Real M_beta;
    static Real M_rhos;

    static std::vector<Real> outputVector;
    static UInt conditionNumber;
};
}

#endif /* __FLOWCONDITIONS_HPP */

