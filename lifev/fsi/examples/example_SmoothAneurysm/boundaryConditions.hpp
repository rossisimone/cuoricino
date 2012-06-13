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
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#ifndef BC_HPP
#define BC_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/BCHandler.hpp>

#include <lifev/fsi/solver/FSIMonolithicGE.hpp>
#include <lifev/fsi/solver/FSIMonolithicGI.hpp>

#include "flowConditions.hpp"
#include "ud_functions.hpp"

/*
//Fluid mesh 100170_cm_N3H0.008_fluid.mesh
#define OUTLET_2 2
#define OUTLET_3 3
#define OUTLET_4 4
#define OUTLET_5 5
#define INLET 7
#define FLUIDINTERFACE 1

//Solid mesh 100170_cm_N3H0.008_solid.mesh
#define OUTERWALL 10
#define SOLIDINTERFACE 1
#define BORDERSURFACE 20
#define BORDERPROBLEM 7
*/

//Regular mesh
#define INLET 3
#define INLETRING 30
#define OUTLET 2
#define OUTLETRING 20
#define FLUIDINTERFACE 200

#define SOLIDINTERFACE 200
#define INLETWALL 3
#define INLETWALL_INTRING 30
#define INLETWALL_OUTRING 31

#define OUTLETWALL 2
#define OUTLETWALL_INTRING 20
#define OUTLETWALL_OUTRING 21

#define OUTERWALL 210



namespace LifeV
{

typedef FSIOperator::fluid_Type fluid;
typedef FSIOperator::solid_Type solid;

FSIOperator::fluidBchandlerPtr_Type BCh_harmonicExtension(FSIOperator &_oper)
{

    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);

    FSISolver::fluidBchandlerPtr_Type BCh_he(new FSIOperator::fluidBchandler_Type );

    BCh_he->addBC("in", INLET, Essential, Full, bcf,   3);
    BCh_he->addBC("in", INLETRING, Essential, Full, bcf,   3);
    BCh_he->addBC("in", OUTLET, Essential, Full, bcf,   3);
    BCh_he->addBC("out3", OUTLETRING, Essential, Full, bcf,   3);


    if (_oper.data().method() == "monolithicGE")
    {
        Debug(10000) << "FSIMonolithic GCE harmonic extension\n";
        FSIMonolithicGE *MOper = dynamic_cast<FSIMonolithicGE *>(&_oper);
        MOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", SOLIDINTERFACE, Essential, Full,
                      *MOper->bcvStructureDispToHarmonicExtension(), 3);
    }
    else if (_oper.data().method() == "monolithicGI")
    {

    }

    return BCh_he;
}


FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFlux(bool /*isOpen=true*/)
{
    FSIOperator::fluidBchandlerPtr_Type BCh_fluid( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase flowAneurysm (fluxFunctionAneurysm);
    BCFunctionBase bcf      (fZero);
    //uncomment  to use fluxes

    //BCh_fluid->addBC("InFlow" , INLET,  Flux, Normal, flowAneurysm);
//   if(!isOpen)
//       BCh_fluid->addBC("InFlow" , INLET,  Flux,   Normal, bcf);

    //uncomment  to use fluxes
    //BCh_fluid->addBC("InFlow" , INLET,  Flux, Normal, flowAneurysm);

    return BCh_fluid;
}

FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFluid(FSIOperator &_oper, bool const & /*isOpen=true*/)
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSIOperator::fluidBchandlerPtr_Type();

    FSIOperator::fluidBchandlerPtr_Type BCh_fluid( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase bcf      (fZero);
    BCFunctionBase in_flow  (uInterpolated);
    //    BCFunctionBase out_flow (fZero);

    BCFunctionBase out_press3 (FlowConditions::outPressure0);

    BCFunctionBase InletVect (aneurismFluxInVectorial);
    //BCFunctionBase bcfw0 (w0);

    //Inlets
    BCh_fluid->addBC("InFlow" , INLET,  EssentialVertices, Full, InletVect, 3);

    //Outlets

    //Absorbing BC seemed not to work
    //Absorbing BC on outlet 2and3 caused instabilities
    BCh_fluid->addBC("out3", OUTLET, Natural,  Normal, out_press3);
    //BCh_fluid->addBC("out3", OUTLET, Natural,  Normal, bcf);

    return BCh_fluid;
}

FSIOperator::solidBchandlerPtr_Type BCh_monolithicSolid(FSIOperator &_oper)
{

    if (! _oper.isSolid() )
        return FSIOperator::solidBchandlerPtr_Type();

    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the solid\n";
    FSIOperator::solidBchandlerPtr_Type BCh_solid( new FSIOperator::solidBchandler_Type );

    BCFunctionBase bcf(fZero);

    //Inlets & Outlets
    BCh_solid->addBC("BORDERS",   INLETWALL, Essential, Full, bcf,  3);
    BCh_solid->addBC("BORDERS-RIN",   INLETWALL_INTRING, Essential, Full, bcf,  3);
    BCh_solid->addBC("BORDERS-ROUT",   INLETWALL_OUTRING, Essential, Full, bcf,  3);
    BCh_solid->addBC("BORDERS",   OUTLETWALL, Essential, Full, bcf,  3);
    BCh_solid->addBC("BORDERS-rin",   OUTLETWALL_INTRING, Essential, Full, bcf,  3);
    BCh_solid->addBC("BORDERS-rout",   OUTLETWALL_OUTRING, Essential, Full, bcf,  3);

    //aortaVelIn::S_timestep = _oper.dataFluid()->dataTime()->timeStep();


    //Robin BC
    BCFunctionBase hyd(fZero);
    BCFunctionBase young (E);
    //robin condition on the outer wall
    _oper.setRobinOuterWall(hyd, young);
    //BCh_solid->addBC("OuterWall", OUTERWALL, Robin, Normal, _oper.bcfRobinOuterWall());


    //First try: Homogeneous Neumann
    BCh_solid->addBC("OuterWall", OUTERWALL, Natural, Normal, bcf);

    return BCh_solid;
}

}

#endif
