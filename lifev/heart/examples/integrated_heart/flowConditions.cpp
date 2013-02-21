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
 *  @brief File containing the lumped parameter models for the Integrated Heart example
 *
 *  @date 2012-09-25
 *  @author Toni Lassila <toni.lassila@epfl.ch>
 *          Paolo Crosetto <crosetto@iacspc70.epfl.ch>

 *  @maintainer Toni Lassila <toni.lassila@epfl.ch>
 *
*/

#include "flowConditions.hpp"

namespace LifeV
{
FlowConditions::FlowConditions()
{
    outputVector.push_back (0);
    conditionNumber = FlowConditions::outputVector.size() - 1;
}

void FlowConditions::setParamsFromGetPot ( const GetPot& dataFile )
{
    mu      = dataFile ( "fluid/physics/viscosity", 1. );
    rho     = dataFile ( "fluid/physics/density", 1. );
    dt      = dataFile ( "fluid/time_discretization/timestep", 1. );
    nu      = mu / rho;
    Pin     = dataFile ( "parameters/Pin", 1. );
    Pout    = dataFile ( "parameters/Pout", 1. );
    M_outflux = 0;
    Rext_d  = dataFile ( "parameters/RExtD", 0. );
    Rext_p  = dataFile ( "parameters/RExtP", 0. );
    Cp      = dataFile ( "parameters/Cp", 0. );
    BDType  = dataFile ( "parameters/BDType", 0.);
    fExt    = dataFile ( "parameters/fExt", 1. );
    tAppl   = dataFile ( "parameters/tAppl", 1. );
    periode = dataFile ( "parameters/periode", 1. );

    // if true impose absorbing boundary conditions on fluid,
    // otherwise on structure.
    bcOnFluid = dataFile ( "parameters/bcOnFluid", true );;
}


void FlowConditions::initParameters ( FSIOperator&  Oper,
                                      const int&    outflowFlag)
{

    Epetra_SerialDenseVector fluidQuantities (1); // M_area0
    Epetra_SerialDenseVector solidQuantities (2); // M_beta and M_rhos

    if (Oper.isFluid() )
    {
        fluidQuantities (0) = Oper.fluid().area (outflowFlag);
    }

    M_area0      = fluidQuantities (0);
    M_outRadius0 = std::sqrt (M_area0 / pi);
    M_inRadius0 = M_outRadius0;

    Oper.displayer().leaderPrint ( "  Outflow BC : area0     = ", M_area0 );
    Oper.displayer().leaderPrint ( "  Outflow BC : radius    = ", M_outRadius0 );


    if (Oper.isSolid() )
    {
        solidQuantities (0) =  ( ( Oper.solid().thickness() * Oper.solid().young (1)     ) / ( 1 - Oper.solid().poisson (1) * Oper.solid().poisson (1) ) * pi / M_area0 );

        solidQuantities (1) = Oper.solid().rho();

        Oper.displayer().leaderPrint ( "  Outflow BC : thickness = " , Oper.solid().thickness() );
        Oper.displayer().leaderPrint ( "  Outflow BC : young     = " , Oper.solid().young (1) );
        Oper.displayer().leaderPrint ( "  Outflow BC : poisson   = " , Oper.solid().poisson (1) );

    }

    M_beta  = solidQuantities (0);
    M_rhos  = solidQuantities (1);
    Oper.displayer().leaderPrint ( "  Outflow BC : beta      = " , M_beta );
    Oper.displayer().leaderPrint ( "  Outflow BC : rho       = " , M_rhos );


}

void FlowConditions::renewParameters ( FSISolver&  oper_,
                                       const int&    outflowFlag)
{

    Epetra_SerialDenseVector fluidQuantities (2); // Flux and Area
    FSIOperator* Oper (oper_.FSIOper().get() );

    if (Oper->isFluid() )
    {
        fluidQuantities (0) = Oper->fluid().flux (outflowFlag, oper_.displacement() );
        fluidQuantities (1) = Oper->fluid().area (outflowFlag);
    }

    Oper->worldComm()->Broadcast ( fluidQuantities.Values(), fluidQuantities.Length(),
                                   Oper->getFluidLeaderId() );


    Real qn;
    Real area;

    qn   = fluidQuantities (0);
    area = fluidQuantities (1);

    // Setting parameters for our simulation:
    // If imposing the absorbing boundary condition through the pressure:

    if (bcOnFluid)
    {
        M_outP =  std::pow ( (M_rhos / (2.*std::sqrt (2.) ) * qn / area + std::sqrt (M_beta * std::sqrt (M_area0) ) ), 2.)
                  - M_beta * std::sqrt (M_area0);
        FlowConditions::outputVector[conditionNumber] = M_outP;

        Oper->displayer().leaderPrint ( " Flow rate = " , qn );
        Oper->displayer().leaderPrint ( " outflow pressure   = " , M_outP );

        M_outDeltaRadius = 0;


    }
    else
    {
        // If imposing the absorbing boundary condition through change in radius: --> Not ready
#ifdef  TESTING
        M_outP = Pout;

        area = qn * std::sqrt (M_rhos) / ( (2.*std::sqrt (2) ) *
                                           std::sqrt ( M_outP + M_beta * std::sqrt (M_area0) ) - std::sqrt ( M_beta * std::sqrt (M_area0) ) );

        assert (area >= 0 );
        if (area < 1e-8 * M_area0)
        {
            area = M_area0;
        }

        M_outDeltaRadius = std::sqrt ( area / pi  ) - M_outRadius0;

        Oper->displayer().leaderPrint ( " outflow A = " , area );
        Oper->displayer().leaderPrint ( " outflow dr = " , M_outDeltaRadius );
        Oper->displayer().leaderPrint ( " Flow rate = " , qn );
        Oper->displayer().leaderPrint ( " outflow pressure   = " , M_outP );
#endif

    }

    // For now applying absBC only at outflow
    M_inDeltaRadius = 0;

}

void FlowConditions::renewLumpedParameters ( const int&    Flag , const Real& Pvalve )
{

    // Setting parameters for our simulation:

    if (Pvalve > Pout)
    {
        M_outflux = (Pvalve - Pout) / Rext_d; // Aortic outflow according to Q_av = (P_av - P_ao) / R_av
    }
    else
    {
        M_outflux = 0;
    }

}


Real FlowConditions::fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}


Real FlowConditions::inPressure (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -Pin;
}


Real FlowConditions::outPressure (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    return Pout;

}

Real FlowConditions::outFlux (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    /* Imposes flow rate at the aortic valve */
    return (M_outflux > 0 ? M_outflux : 0); // On-off valve

}

Real FlowConditions::outProfile (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    /* Imposes flat velocity profile at the aortic valve */
    Real outVelocity = (M_outflux > 0 ? M_outflux : 0) / 5.29414; // Scale velocity to aortic valve area
    Real normalVector[3] = {0, 0, 1};

    switch (i)
    {
        case 0:
            return normalVector[0] * outVelocity;
            break;
        case 1:
            return normalVector[1] * outVelocity;
            break;
        case 2:
            return normalVector[2] * outVelocity;
            break;
        default:
            ERROR_MSG ("This entry is not allowed: flowConditions.hpp");
            break;
    }

}


Real FlowConditions::force_cardium (const Real& t)
{
    int n = 5;

    Real coeffA[5] = { 0.7120, -0.154736633422335,  -0.152048143771071,  -0.037433333333333,  -0.011781889473260  };

    Real coeffB[5] = {  0 ,  0.464877195821670,  -0.053906336561117,  -0.007409328454600 , -0.024744907037225};

    Real force (0);
    force = coeffA[0] / 2;
    for (int k = 1; k < n; ++k)
    {
        force += coeffA[k] * cos (k * t) + coeffB[k] * sin (k * t);
    }

    force *= (force > 0);
    return force;

}

Real FlowConditions::fextvessel (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{

    Real t_loc = t * 2 * acos (-1) / periode;

    return -fExt * force_cardium (t_loc) * 2 * (8 - z) / 9.5;

}

Real FlowConditions::outPressure0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[0];
}

Real FlowConditions::outPressure1 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[1];
}

Real FlowConditions::outPressure2 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[2];
}
Real FlowConditions::outPressure3 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[4];
}

Real FlowConditions::outPressure5 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[5];
}

Real FlowConditions::outPressure6 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[6];
}


Real FlowConditions::inDeltaRadius (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
    if (i == 2)
    {
        return 0;
    }

    Real r ( std::sqrt (x * x + y * y) );

    switch (i)
    {
        case 0:
            return M_inDeltaRadius * x / r;
        case 1:
            return M_inDeltaRadius * y / r;
        default:
            ERROR_MSG ("This entry is not allowed: flowConditions.hpp");
    };
    return 0.;
}

Real FlowConditions::outDeltaRadius (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
    if (i == 2)
    {
        return 0;
    }

    Real r ( std::sqrt (x * x + y * y) );

    switch (i)
    {
        case 0:
            return M_outDeltaRadius * x / r;
        case 1:
            return M_outDeltaRadius * y / r;
        default:
            ERROR_MSG ("This entry is not allowed: flowConditions.hpp");
    };
    return 0.;
}

Real FlowConditions::nu (-1);
Real FlowConditions::mu (-1);
Real FlowConditions::dt (1);
Real FlowConditions::rho (-1);
Real FlowConditions::Pin (-1);
Real FlowConditions::Pout (-1);

//For Explicit Resistance or Windkessel

Real FlowConditions::Rext_d (0);
Real FlowConditions::Rext_p (0);
Real FlowConditions::Cp (0);
int FlowConditions::BDType (0);
Real FlowConditions::Flux (0);
Real FlowConditions::Flux_old (0);

Real FlowConditions::fExt (-1);
Real FlowConditions::tAppl (-1);
Real FlowConditions::periode (-1);
Real FlowConditions::pi (3.141592635);

bool FlowConditions::bcOnFluid (true);

Real FlowConditions::M_outflux (0);
Real FlowConditions::M_influx (0);
Real FlowConditions::M_outP (0);
Real FlowConditions::M_inP (0);

Real FlowConditions::M_area0 (0);
Real FlowConditions::M_inRadius0 (0);
Real FlowConditions::M_outRadius0 (0);
Real FlowConditions::M_inDeltaRadius (0);
Real FlowConditions::M_outDeltaRadius (0);

Real FlowConditions::M_beta (0);
Real FlowConditions::M_rhos (0);

UInt FlowConditions::conditionNumber (0);

std::vector<Real> FlowConditions::outputVector;
}
