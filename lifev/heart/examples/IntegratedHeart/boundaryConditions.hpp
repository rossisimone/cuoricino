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
 *  @brief File containing the boundary conditions for the Integrated Heart example
 *
 *  @date 2012-09-25
 *  @author Toni Lassila <toni.lassila@epfl.ch>
 *          Paolo Crosetto <crosetto@iacspc70.epfl.ch>

 *  @maintainer Toni Lassila <toni.lassila@epfl.ch>
 *
*/

#ifndef BC_HPP
#define BC_HPP

// LifeV includes
#include "lifev/core/LifeV.hpp"
#include "lifev/core/fem/BCHandler.hpp"

// Mathcard includes
#include "lifev/fsi/solver/FSIMonolithicGE.hpp"
#include "lifev/fsi/solver/FSIMonolithicGI.hpp"

#include "flowConditions.hpp"
#include "ud_functions.hpp"

#define OUTLET 3
#define INLET  2
#define OUTLETRING 30
#define INLETRING  20
#define FLUIDINTERFACE 5
#define TOPCARDIUM  4
#define TOPCARDIUMRING  40
#define AORTICROOT 50

#define SOLIDINTERFACE 5
#define OUTERRING 20
#define INNERRING 30
#define OUTERWALL 10
#define TOP 40

namespace LifeV
{

typedef FSIOperator::fluid_Type fluid;
typedef FSIOperator::solid_Type solid;

FSIOperator::fluidBchandlerPtr_Type BCh_harmonicExtension (FSIOperator& _oper)
{

    BCFunctionBase bcf (fZero);

    FSISolver::fluidBchandlerPtr_Type BCh_he (new FSIOperator::fluidBchandler_Type );

    BCh_he->addBC ("HE1", INLET,          Essential,         Full, bcf, 3);
    BCh_he->addBC ("HE2", OUTLET,         Essential,         Full, bcf, 3);
    BCh_he->addBC ("HE3", TOPCARDIUM,     Essential,         Full, bcf, 3);
    BCh_he->addBC ("HE4", INLETRING,      EssentialVertices, Full, bcf, 3);
    BCh_he->addBC ("HE5", OUTLETRING,     EssentialVertices, Full, bcf, 3);
    BCh_he->addBC ("HE6", TOPCARDIUMRING, EssentialVertices, Full, bcf, 3);
    BCh_he->addBC ("HE7", AORTICROOT,     EssentialVertices, Full, bcf, 3);

    if (_oper.data().method() == "monolithicGE")
    {
        FSIMonolithicGE* MOper = dynamic_cast<FSIMonolithicGE*> (&_oper);
        MOper->setStructureDispToHarmonicExtension (_oper.lambdaFluidRepeated() );
        BCh_he->addBC ("Interface", SOLIDINTERFACE, Essential, Full,
                       *MOper->bcvStructureDispToHarmonicExtension(), 3);
    }
    else if (_oper.data().method() == "monolithicGI")
    {

    }

    return BCh_he;
}


FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFlux ( bool /*AorticValveisOpen*/, bool /*MitralValveisOpen*/ )
{

    FSIOperator::fluidBchandlerPtr_Type BCh_fluid ( new FSIOperator::fluidBchandler_Type );

    /* Defective version of outflow b.c. (unstable) */
    // BCFunctionBase out_flux (LifeV::FlowConditions::outFlux);
    // BCh_fluid->addBC("FL2", OUTLET, Flux, Normal, out_flux);

    return BCh_fluid;
}

FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFluid (FSIOperator& _oper )
{
    if (! _oper.isFluid() )
    {
        return FSIOperator::fluidBchandlerPtr_Type();
    }

    FSIOperator::fluidBchandlerPtr_Type BCh_fluid ( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase bcNoSlip (fZero);
    BCFunctionBase mitral_flow (u_mitral );

    /*  Note: cannot fix both fluid and structure on a coupled interface ring */
    BCh_fluid->addBC ("FL1", INLET,      Essential, Full,   mitral_flow, 3);
    BCh_fluid->addBC ("FL3", TOPCARDIUM, EssentialVertices, Full,   bcNoSlip,    3);
    BCh_fluid->addBC ("FL4", INLETRING,  Essential, Full,   bcNoSlip,    3);
    BCh_fluid->addBC ("FL5", OUTLETRING, Essential, Full,   bcNoSlip,    3);
    BCh_fluid->addBC ("FL6", AORTICROOT, EssentialVertices, Full,   bcNoSlip,    3);

    /* Dirichlet version of outflow b.c. (stable?) */
    BCFunctionBase out_profile (LifeV::FlowConditions::outProfile);
    BCh_fluid->addBC ("FL2", OUTLET, Essential, Full, out_profile, 3);

    return BCh_fluid;
}

FSIOperator::solidBchandlerPtr_Type BCh_monolithicSolid (FSIOperator& _oper)
{

    if (! _oper.isSolid() )
    {
        return FSIOperator::solidBchandlerPtr_Type();
    }

    // Boundary conditions for the solid displacement
    FSIOperator::solidBchandlerPtr_Type BCh_solid ( new FSIOperator::solidBchandler_Type );

    BCFunctionBase bcf (fZero);
    BCFunctionBase aroundheart (FlowConditions::fextvessel);

    std::vector<LifeV::ID> zComp (1);
    zComp[0] = 3;

    BCh_solid->addBC ("OuterWall", OUTERWALL, Natural,           Normal, aroundheart);
    BCh_solid->addBC ("Top",       TOP,       EssentialVertices, Full,   bcf,         3);
    BCh_solid->addBC ("InnerRing", INNERRING, EssentialVertices, Full,   bcf,         3);

    return BCh_solid;
}

}

#endif
