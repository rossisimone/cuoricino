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
    @brief

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @date 00-00-0000
 */

#ifndef BC_HPP
#define BC_HPP

#include "life/lifecore/LifeV.hpp"
#include "ud_functions.hpp"
#include "life/lifefem/BCHandler.hpp"
#include "life/lifefem/BCFunction.hpp"

#include "life/lifesolver/FSIExactJacobian.hpp"
#include "life/lifesolver/FSIFixedPoint.hpp"


//#define FLUX = 

namespace LifeV
{

typedef FSIOperator::fluid_Type fluid;
typedef FSIOperator::solid_Type solid;

FSIOperator::fluidBchandlerPtr_Type BCh_harmonicExtension(FSIOperator &_oper)
{

//     debugStream(10000) << "SP harmonic extension\n";
//     fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

//     FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluid());

    if (! _oper.isFluid() )
        return FSIOperator::fluidBchandlerPtr_Type();

//    FPOper->bcvStructureDispToHarmonicExtension()->showMe(true,std::cout);

    // Boundary condition for the mesh
    debugStream( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);

    FSISolver::fluidBchandlerPtr_Type BCh_he(new FSIOperator::fluidBchandler_Type );

    BCh_he->addBC("Top",         3, Essential, Full, bcf,   3);
    BCh_he->addBC("Base",        2, Essential, Full, bcf,   3);
    BCh_he->addBC("Base",        4, Essential, Full, bcf,   3);
    BCh_he->addBC("Top",         30, Essential, Full, bcf,   3);
    BCh_he->addBC("Top",         20, Essential, Full, bcf,   3);

    if (_oper.data().method() == "steklovPoincare")
    {
//         debugStream(10000) << "SP harmonic extension\n";
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setFluidInterfaceDisp((LifeV::Vector&) _oper.lambdaFluidRepeated());
//         BCh_he->addBC("Interface", 1, Essential, Full,
//                       *SPOper->bcvFluidInterfaceDisp(), 3);
    }
    else if (_oper.data().method() == "exactJacobian")
    {
        debugStream(10000) << "EJ harmonic extension\n";
        FSIExactJacobian *EJOper = dynamic_cast<FSIExactJacobian *>(&_oper);
        EJOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", 1, Essential, Full,
                      *EJOper->bcvStructureDispToHarmonicExtension(), 3);
    }
    else if (_oper.data().method() == "fixedPoint")
    {
        debugStream(10000) << "FP harmonic extension\n";
        FSIFixedPoint *FPOper = dynamic_cast<FSIFixedPoint *>(&_oper);

        FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", 1, Essential, Full,
                      *FPOper->bcvStructureDispToHarmonicExtension(), 3);
    }



    return BCh_he;
}


FSIOperator::fluidBchandlerPtr_Type BCh_fluid(FSIOperator &_oper)
{
    // Boundary conditions for the fluid velocity
    debugStream( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSIOperator::fluidBchandlerPtr_Type();

    FSIOperator::fluidBchandlerPtr_Type BCh_fluid( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase bcf           (fZero);
    BCFunctionBase in_flow       (u2);
    BCFunctionBase in_vel        (u2vel);
    BCFunctionBase in_flow_pr       (pressure);
    BCFunctionBase in_flow_flux  (PhysFlux);
    BCFunctionBase out_flow      (fZero);


    //  #ifdef FLUX
    // BCh_fluid->addBC("InFlow" ,   2,  Flux,   Full, in_flow_flux, 3);
    // #else
     BCh_fluid->addBC("InFlow" , 2,  Natural,   Full, in_flow, 3);
     // #endif

    BCh_fluid->addBC("EdgesIn",  20, Essential, Full, bcf,  3);

    BCh_fluid->addBC("OutFlowBrain",   3,  Natural,   Full, out_flow, 3);
    BCh_fluid->addBC("EdgesFace",  30, Essential, Full, bcf,  3);
    BCh_fluid->addBC("OutFlowBrain",   4,  Natural,   Full, out_flow, 3);
    BCh_fluid->addBC("EdgesBrain",  40, Essential, Full, bcf,  3);

    _oper.setStructureToFluid(_oper.veloFluidMesh());
    // _oper.setHarmonicExtensionVelToFluid(_oper.veloFluidMesh());

    if (_oper.data().algorithm()=="RobinNeumann")
    {
        // _oper.setAlphafbcf(alpha); // if alpha is bcFunction define in ud_function.cpp

        _oper.setSolidLoadToStructure( _oper.minusSigmaFluidRepeated());
        _oper.setStructureToFluidParameters();

        BCh_fluid->addBC("Interface",   1,  Robin, Full,
                         *_oper.bcvStructureToFluid(),  3);
        BCh_fluid->addBC("Interface",   1,  Natural, Full,
                         *_oper.bcvSolidLoadToStructure(), 3);
    }
    else
    {
        BCh_fluid->addBC("Interface",   1,  Essential, Full,
                         *_oper.bcvStructureToFluid(),  3);
    }
    return BCh_fluid;
}


FSIOperator::fluidBchandlerPtr_Type BCh_fluidInv(FSIOperator &_oper)
{

    if (! _oper.isFluid() )
        return FSIOperator::fluidBchandlerPtr_Type();

    // Boundary conditions for the fluid velocity
    debugStream( 10000 ) << "Boundary condition for the inverse fluid\n";
    FSIOperator::fluidBchandlerPtr_Type BCh_fluidInv( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase bcf(fZero);
    BCFunctionBase in_flow(u2);
    BCFunctionBase in_vel(u2vel);
    BCFunctionBase in_flow_pr(pressure);


    BCh_fluidInv->addBC("InFlow", 2,  Essential,  Normal, in_vel);
    //BCh_fluidInv->addBC("InFlow", 5,  Natural,   Normal, in_flow_pr);
    BCh_fluidInv->addBC("EdgesIn",  20, Essential, Full, bcf,     3);

    return BCh_fluidInv;
}




FSIOperator::fluidBchandlerPtr_Type BCh_fluidLin(FSIOperator &_oper)
{
    if (! _oper.isFluid() )
        return FSIOperator::fluidBchandlerPtr_Type();

    // Boundary conditions for the fluid velocity
    debugStream( 10000 ) << "Boundary condition for the linearized fluid\n";
    FSIOperator::fluidBchandlerPtr_Type BCh_fluidLin( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase bcf(fZero);
    BCFunctionBase in_flow(u2);

    #ifdef FLUX
     BCh_fluidLin->addBC("InFlow",   2,       Flux, Full, bcf,     3);
    #else
    BCh_fluidLin->addBC("InFlow",  2,  Essential,  Full, bcf, 3);
    // BCh_fluidLin->addBC("InFlow",  5,  Natural, Normal, bcf);
    #endif

    BCh_fluidLin->addBC("InletFace",  20,  Essential, Full, bcf,     3);

    BCh_fluidLin->addBC("outFlowFace",  3,    Natural, Full, bcf,     3);
    BCh_fluidLin->addBC("edgesFace",  30,  Essential, Full, bcf,     3);
    BCh_fluidLin->addBC("outFlowBrain",  4,    Natural, Full, bcf,     3);
    BCh_fluidLin->addBC("InletFace",  40,  Essential, Full, bcf,     3);


    //BCh_fluidLin->addBC("ainterface",  1,  Essential,   Full, bcf,     3);


    if (_oper.data().method() == "steklovPoincare")
    {
//             steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//             SPOper->setDerHarmonicExtensionVelToFluid(_oper.dw());
//             BCh_fluidLin->addBC("Wall"     , 1, Essential  , Full,
//                                 *SPOper->bcvDerHarmonicExtensionVelToFluid(), 3);
    }
    if (_oper.data().method() == "exactJacobian")
    {
        FSIExactJacobian* EJOper = dynamic_cast<FSIExactJacobian *>(&_oper);
        EJOper->setDerHarmonicExtensionVelToFluid(_oper.derVeloFluidMesh());
        BCh_fluidLin->addBC("Interface", 1, Essential  , Full,
                            *_oper.bcvDerHarmonicExtensionVelToFluid(), 3);
    }

    return BCh_fluidLin;
}



FSIOperator::solidBchandlerPtr_Type BCh_solid(FSIOperator &_oper)
{

    if (! _oper.isSolid() )
        return FSIOperator::solidBchandlerPtr_Type();

    // Boundary conditions for the solid displacement
    debugStream( 10000 ) << "Boundary condition for the solid\n";
    FSIOperator::solidBchandlerPtr_Type BCh_solid( new FSIOperator::solidBchandler_Type );

    BCFunctionBase bcf(fZero);


    BCh_solid->addBC("BaseRingF5",      3, Essential, Full, bcf,  3);
    BCh_solid->addBC("BaseRingF3",      2, Essential, Full, bcf,  3);
    BCh_solid->addBC("Ring6",      20, Essential, Full, bcf,  3);
    BCh_solid->addBC("Ring6",      30, Essential, Full, bcf,  3);

    std::vector<ID> zComp(1);
    zComp[0] = 3;
//     debugStream(10000) << "SP harmonic extension\n";

    if (_oper.data().method() == "steklovPoincare")
    {
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setSolidInterfaceDisp((LifeV::Vector&) _oper.displacement());

//         BCh_solid->addBC("Interface", 1, Essential, Full,
//                          *SPOper->bcvSolidInterfaceDisp(), 3);
    }
    else if (_oper.data().method() == "exactJacobian")
    {
        FSIExactJacobian  *EJOper = dynamic_cast<FSIExactJacobian *>(&_oper);
        EJOper->setFluidLoadToStructure(_oper.sigmaSolidRepeated());

        BCh_solid->addBC("Interface", 1, Natural,   Full,
                         *EJOper->bcvFluidLoadToStructure(), 3);
    }
    else if (_oper.data().method() == "fixedPoint")
    {
        FSIFixedPoint *FPOper = dynamic_cast<FSIFixedPoint *>(&_oper);

        FPOper->setFluidLoadToStructure(_oper.sigmaSolidRepeated());

        BCh_solid->addBC("Interface", 1, Natural, Full,
                         *FPOper->bcvFluidLoadToStructure(), 3);
    }

    return BCh_solid;
}


FSIOperator::solidBchandlerPtr_Type BCh_solidLin(FSIOperator &_oper)
{
    if (! _oper.isSolid() )
        return FSIOperator::solidBchandlerPtr_Type();

    // Boundary conditions for the solid displacement
    debugStream( 10000 ) << "Boundary condition for the linear solid\n";
    FSIOperator::solidBchandlerPtr_Type BCh_solidLin( new FSIOperator::solidBchandler_Type );

    BCFunctionBase bcf(fZero);

    BCh_solidLin->addBC("BaseRingF5",      3, Essential, Full, bcf,  3);
    BCh_solidLin->addBC("BaseRingF3",      2, Essential, Full, bcf,  3);
    BCh_solidLin->addBC("Ring6",      20, Essential, Full, bcf,  3);


    std::vector<ID> zComp(1);
    zComp[0] = 3;

    if (_oper.data().method() == "steklovPoincare")
    {
//             steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//             SPOper->setSolidLinInterfaceDisp((LifeV::Vector&) _oper.displacement());
//             BCh_solidLin->addBC("Interface", 1, Essential, Full,
//                                 *SPOper->bcvSolidLinInterfaceDisp(), 3);
    }
    else if (_oper.data().method() == "exactJacobian")
    {
        FSIExactJacobian  *EJOper = dynamic_cast<FSIExactJacobian *>(&_oper);
        EJOper->setDerFluidLoadToStructure(_oper.sigmaSolidRepeated());
        BCh_solidLin->addBC("Interface", 1, Natural,   Full,
                            *EJOper->bcvDerFluidLoadToStructure(), 3);
    }

    return BCh_solidLin;
}

FSIOperator::solidBchandlerPtr_Type BCh_solidInvLin(FSIOperator &_oper)
{

    if (! _oper.isSolid() )
        return FSIOperator::solidBchandlerPtr_Type();

    // Boundary conditions for the solid displacement
    debugStream( 10000 ) << "Boundary condition for the inverse linear solid\n";
    FSIOperator::solidBchandlerPtr_Type BCh_solidLinInv( new FSIOperator::solidBchandler_Type );

    BCFunctionBase bcf(fZero);


    BCh_solidLinInv->addBC("BasRingF5",     3, Essential, Full, bcf,  3);
    BCh_solidLinInv->addBC("BaseRingF3",   2, Essential, Full, bcf,  3);
    BCh_solidLinInv->addBC("Ring6",      20, Essential, Full, bcf,  3);

//     if (_oper.method() == "steklovPoincare")
//     {
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setSolidInvLinInterfaceStress((LifeV::Vector&) _oper.residual());

//         BCh_solidLinInv->addBC("Interface", 100, Natural, Full,
//                                *SPOper->bcvSolidInvLinInterfaceStress(), 3);
//     }
//     else
//     {
//         exactJacobian  *EJOper = dynamic_cast<exactJacobian *>(&_oper);
//         EJOper->setDerFluidLoadToStructure(_oper.fluid().residual());
//         BCh_solidLinInv->addBC("Interface", 1, Natural,   Full,
//                             *EJOper->bcvDerFluidLoadToStructure(), 3);
//     }

    return BCh_solidLinInv;
}









}

#endif
