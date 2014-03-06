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
 *  @date 2014-02-18
 *  @author Davide Forti <crosetto@iacspc70.epfl.ch>
 *
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

#include "ud_functions.hpp"

#define OUTLET 3
#define INLET 2
#define FLUIDINTERFACE 1
#define SOLIDINTERFACE 1
#define OUTERWALL 10
#define RING  2
#define RING2 3
#define INOUTEDGE 20
#define INEDGE 30

namespace LifeV
{

    typedef BCHandler                  bc_Type;
    typedef boost::shared_ptr<bc_Type> bcPtr_Type;
    
class boundaryConditions
{
public:
    
    boundaryConditions();
    
    ~boundaryConditions();
    
    bcPtr_Type fluidBC();
    
    bcPtr_Type structureBC();
    
    bcPtr_Type aleBC();
    
    void setup();
    
private:
    
    bcPtr_Type M_bcFluid;
    bcPtr_Type M_bcStructure;
    bcPtr_Type M_bcAle;
    
};
    
boundaryConditions::boundaryConditions()
{}

boundaryConditions::~boundaryConditions()
{}

bcPtr_Type boundaryConditions::fluidBC()
{
    return M_bcFluid;
}

bcPtr_Type boundaryConditions::structureBC()
{
    return M_bcStructure;
}

bcPtr_Type boundaryConditions::aleBC()
{
    return M_bcAle;
}
    
void boundaryConditions::setup()
{
    BCFunctionBase bcf (fZero);
    BCFunctionBase vel (w0);
    BCFunctionBase disp (vinit);
    
    vector <ID> compy (1);
    compy[0] = 1;

    // Fluid BC
    M_bcFluid.reset(new bc_Type());
    M_bcFluid->addBC ("InFlow" , INLET,   Essential, Full,   vel, 3);
    M_bcFluid->addBC ("OutFlow", OUTLET,  Natural,   Normal, bcf);
    
    // Structure BC
    M_bcStructure.reset(new bc_Type());
    M_bcStructure->addBC ("Top",   RING,  Essential, Full, bcf, 3);
    M_bcStructure->addBC ("Base",  RING2, Essential, Full, bcf, 3);
    M_bcStructure->addBC ("Interface",  SOLIDINTERFACE, Essential, Component, disp, compy);

    // ALE BC (displacement on gamma to be prescibed later)
    M_bcAle.reset(new bc_Type());
    M_bcAle->addBC ("Edges", INOUTEDGE, EssentialVertices, Full, bcf, 3);
    M_bcAle->addBC ("Edges", INEDGE,    EssentialVertices, Full, bcf, 3);
    M_bcAle->addBC ("Base",  INLET,     EssentialVertices, Full, bcf, 3);
}

} // end Namespace LifeV

#endif
