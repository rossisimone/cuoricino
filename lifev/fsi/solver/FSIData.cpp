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
    @brief DataFSI - File containing a data container for FSI problems

    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @author Gilles fourestey <gilles.fourestey@epfl.ch>
    @date 10-06-2010

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>
 */


#include <lifev/fsi/solver/FSIData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
FSIData::FSIData( ) :
    M_dataFluid                     ( new dataFluid_Type() ),
    M_dataSolid                     ( new dataSolid_Type() ),
    M_maxSubIterationNumber         (),
    M_absoluteTolerance             (),
    M_relativeTolerance             (),
    M_errorTolerance                (),
    M_NonLinearLineSearch           (),
    M_method                        (),
    M_algorithm                     (),
    M_defaultOmega                  (),
    M_rangeOmega                    (),
    M_updateEvery                   (),
    M_fluidInterfaceFlag            (),
    M_structureInterfaceFlag        (),
    M_fluidInterfaceVertexFlag      (),
    M_structureInterfaceVertexFlag  (),
    M_interfaceTolerance            ()
{
}

FSIData::FSIData ( const FSIData& FSIData ) :
    M_dataFluid                     ( FSIData.M_dataFluid ),
    M_dataSolid                     ( FSIData.M_dataSolid ),
    M_maxSubIterationNumber         ( FSIData.M_maxSubIterationNumber ),
    M_absoluteTolerance             ( FSIData.M_absoluteTolerance ),
    M_relativeTolerance             ( FSIData.M_relativeTolerance ),
    M_errorTolerance                ( FSIData.M_errorTolerance ),
    M_NonLinearLineSearch           ( FSIData.M_NonLinearLineSearch ),
    M_method                        ( FSIData.M_method ),
    M_algorithm                     ( FSIData.M_algorithm ),
    M_defaultOmega                  ( FSIData.M_defaultOmega ),
    M_rangeOmega                    ( FSIData.M_rangeOmega ),
    M_updateEvery                   ( FSIData.M_updateEvery ),
    M_fluidInterfaceFlag            ( FSIData.M_fluidInterfaceFlag ),
    M_structureInterfaceFlag        ( FSIData.M_structureInterfaceFlag ),
    M_fluidInterfaceVertexFlag      ( new Int const ( *FSIData.M_fluidInterfaceVertexFlag ) ),
    M_structureInterfaceVertexFlag  ( new Int const ( *FSIData.M_structureInterfaceVertexFlag ) ),
    M_interfaceTolerance            ( FSIData.M_interfaceTolerance ),
    M_restartTimeStep               ( 0. )
{
}


// ===================================================
// Methods
// ===================================================
FSIData&
FSIData::operator= ( const FSIData& FSIData )
{
    if ( this != &FSIData )
    {
        M_dataFluid                     = FSIData.M_dataFluid;
        M_dataSolid                     = FSIData.M_dataSolid;
        M_maxSubIterationNumber         = FSIData.M_maxSubIterationNumber;
        M_absoluteTolerance             = FSIData.M_absoluteTolerance;
        M_relativeTolerance             = FSIData.M_relativeTolerance;
        M_errorTolerance                = FSIData.M_errorTolerance;
        M_NonLinearLineSearch           = FSIData.M_NonLinearLineSearch;
        M_method                        = FSIData.M_method;
        M_algorithm                     = FSIData.M_algorithm;
        M_defaultOmega                  = FSIData.M_defaultOmega;
        M_rangeOmega                    = FSIData.M_rangeOmega;
        M_updateEvery                   = FSIData.M_updateEvery;
        M_fluidInterfaceFlag            = FSIData.M_fluidInterfaceFlag;
        M_structureInterfaceFlag        = FSIData.M_structureInterfaceFlag;

        M_fluidInterfaceVertexFlag.reset    ( new Int const ( *FSIData.M_fluidInterfaceVertexFlag ) );
        M_structureInterfaceVertexFlag.reset ( new Int const ( *FSIData.M_structureInterfaceVertexFlag ) );

        M_interfaceTolerance            = FSIData.M_interfaceTolerance;
    }

    return *this;
}

void
FSIData::setup ( const GetPot& dataFile, const std::string& section )
{
    if ( !M_timeALE.get() )
    {
        M_timeALE.reset ( new time_Type ( dataFile, "mesh_motion/time_discretization" ) );
    }

    if ( !M_timeAdvanceALE.get() )
    {
        M_timeAdvanceALE.reset ( new timeAdvance_Type ( dataFile, "mesh_motion/time_discretization" ) );
    }

    M_dataFluid->setup ( dataFile );
    M_dataSolid->setup ( dataFile );

    // Problem - Non Linear Richardson Parameters
    M_maxSubIterationNumber = dataFile ( ( section + "/maxSubIter" ).data(), 300 );
    M_absoluteTolerance = dataFile ( ( section + "/abstol" ).data(), 1.e-07 );
    M_relativeTolerance = dataFile ( ( section + "/reltol" ).data(), 1.e-04 );
    M_errorTolerance = dataFile ( ( section + "/etamax" ).data(), 1.e-03 );
    M_NonLinearLineSearch = static_cast<Int> ( dataFile ( ( section + "/NonLinearLineSearch" ).data(), 0 ) );

    // Problem - Methods
    M_method = dataFile ( ( section + "/method" ).data(), "steklovPoincare" );
    M_algorithm = dataFile ( ( section + "/algorithm" ).data(), "DirichletNeumann" );

    // Problem - FixPoint / EJ
    M_defaultOmega = dataFile ( ( section + "/defOmega" ).data(), 0.001);
    M_rangeOmega[0] = dataFile ( ( section + "/defOmega" ).data(), std::fabs ( M_defaultOmega ) /**1024.*/, 0);
    M_rangeOmega[1] = dataFile ( ( section + "/defOmega" ).data(), std::fabs ( M_defaultOmega ) /*/1024.*/, 1);
    M_updateEvery = dataFile ( ( section + "/updateEvery" ).data(), 1);

    // Interface
    M_fluidInterfaceFlag     = dataFile ( "interface/fluid_flag", 1 );
    M_structureInterfaceFlag = dataFile ( "interface/solid_flag", M_fluidInterfaceFlag );

    Int vertexFlag;
    vertexFlag               = dataFile ( "interface/edgeFlag",      -1 );
    vertexFlag               = dataFile ( "interface/fluid_vertex_flag", vertexFlag );

    if (vertexFlag >= 0)
    {
        M_fluidInterfaceVertexFlag.reset    ( new Int const ( vertexFlag ) );
    }

    vertexFlag               = dataFile ( "interface/structure_vertex_flag", -1 );
    if (vertexFlag >= 0)
    {
        M_structureInterfaceVertexFlag.reset ( new Int const ( vertexFlag ) );
    }

    M_interfaceTolerance = dataFile ( "interface/tolerance",      0. );

    M_restartTimeStep  = dataFile ( "importer/restart_timestep",      0. );
}

bool
FSIData::isMonolithic()
{
    return ! ( M_method.compare ( "monolithicGE" ) && M_method.compare ( "monolithicGI" ) );
}

void
FSIData::showMe ( std::ostream& output )
{
    output << "\n*** Values for data fluid\n\n";
    M_dataFluid->showMe();

    output << "\n*** Values for data solid\n\n";
    M_dataSolid->showMe();

    output << "\n*** Values for problem\n\n";
    output << "Max subiteration number          = " << M_maxSubIterationNumber << std::endl;
    output << "Absolute tolerance               = " << M_absoluteTolerance << std::endl;
    output << "Relative tolerance               = " << M_relativeTolerance << std::endl;
    output << "Max error tolerance              = " << M_errorTolerance << std::endl;
    output << "NonLinearLineSearch                       = " << M_NonLinearLineSearch << std::endl;

    output << "Method                           = " << M_method << std::endl;
    output << "Algorithm                        = " << M_algorithm << std::endl;

    output << "Default Omega                    = " << M_defaultOmega << std::endl;
    output << "Omega range                      = " << "(" << M_rangeOmega[0] << " " << M_rangeOmega[1] << ")" << std::endl;
    output << "Update every                     = " << M_updateEvery << std::endl;

    output << "\n*** Values for interface\n\n";
    output << "Interface fluid                  = " << M_fluidInterfaceFlag << std::endl;
    output << "Interface structure              = " << M_structureInterfaceFlag << std::endl;
    if (M_fluidInterfaceVertexFlag.get() != 0)
    {
        output << "Interface fluid vertices         = " << *M_fluidInterfaceVertexFlag << std::endl;
    }
    if (M_structureInterfaceVertexFlag.get() != 0)
    {
        output << "Interface structure vertices     = " << *M_structureInterfaceVertexFlag << std::endl;
    }
    output << "Interface tolerance              = " << M_interfaceTolerance << std::endl;
}

}
