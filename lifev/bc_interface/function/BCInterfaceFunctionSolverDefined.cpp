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
 *  @brief File containing the BCInterfaceFunctionSolverDefined class and specializations
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/bc_interface/function/BCInterfaceFunctionSolverDefined.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceFunctionSolverDefined< FSIOperator >::BCInterfaceFunctionSolverDefined() :
        M_FSIFunction           (),
        M_physicalSolver        (),
        M_name                  (),
        M_flag                  (),
        M_type                  (),
        M_mode                  (),
        M_componentsVector      (),
        M_vectorFunctionRobin   (),
        M_robinRHS              (),
        M_robinAlphaCoefficient (),
        M_robinBetaCoefficient  ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::BCInterfaceFunctionSolverDefined()" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< FSIOperator >::exportData( BCInterfaceData3D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::exportData" << "\n";
#endif

    data.setName( M_name );
    data.setFlag( M_flag );
    data.setType( M_type );
    data.setMode( M_mode );
    data.setComponentsVector( M_componentsVector );
}

void
BCInterfaceFunctionSolverDefined< FSIOperator >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::updatePhysicalSolverVariables" << "\n";
#endif

    switch ( M_FSIFunction )
    {
    case RobinWall:
    {
        if ( !M_physicalSolver->isSolid() )
            return;

        // Update the physical solver variables
        for ( UInt i( 0 ); i < M_vectorFunctionRobin.size(); ++i )
        {
            boost::shared_ptr< BCInterfaceFunctionParserSolver< physicalSolver_Type > > castedFunctionSolver =
                boost::dynamic_pointer_cast< BCInterfaceFunctionParserSolver< physicalSolver_Type > > ( M_vectorFunctionRobin[i] );

            if ( castedFunctionSolver != 0 )
                castedFunctionSolver->updatePhysicalSolverVariables();
        }

        // Set coefficients
        Int gid;
        Real x, y, z;
        Real alpha, beta;
        Real t( M_physicalSolver->dataSolid()->dataTime()->time() );
        Real timeStep( M_physicalSolver->dataSolid()->dataTime()->timeStep() );

        // Update Time advance
        M_physicalSolver->solidTimeAdvance()->updateRHSFirstDerivative( timeStep );

        Int verticesGlobalNumber( M_physicalSolver->solidMeshPart().meshPartition()->numGlobalVertices() );
        for ( UInt i(0) ; i < M_physicalSolver->solidMeshPart().meshPartition()->numVertices() ; ++i )
        {
            gid = M_physicalSolver->solidMeshPart().meshPartition()->meshTransformer().pointInitial( i ).id();

            x   = M_physicalSolver->solidMeshPart().meshPartition()->meshTransformer().pointInitial( i ).x();
            y   = M_physicalSolver->solidMeshPart().meshPartition()->meshTransformer().pointInitial( i ).y();
            z   = M_physicalSolver->solidMeshPart().meshPartition()->meshTransformer().pointInitial( i ).z();

            alpha = M_vectorFunctionRobin[0]->functionTimeSpace( t, x, y, z, 0 );
            beta  = M_vectorFunctionRobin[1]->functionTimeSpace( t, x, y, z, 0 );

            alpha += M_physicalSolver->solidTimeAdvance()->coefficientFirstDerivative( 0 ) / timeStep * beta;

            (*M_robinAlphaCoefficient)[gid] = alpha;
            (*M_robinBetaCoefficient)[gid]  = beta;

            (*M_robinAlphaCoefficient)[gid + verticesGlobalNumber] = alpha;
            (*M_robinBetaCoefficient)[gid + verticesGlobalNumber]  = beta;

            (*M_robinAlphaCoefficient)[gid + verticesGlobalNumber * 2] = alpha;
            (*M_robinBetaCoefficient)[gid + verticesGlobalNumber * 2]  = beta;
        }

        // Set displacement and velocity at time tn (mid-point scheme for the solid)
        FSIOperator::vector_Type displacementTn( M_physicalSolver->dFESpace().map(), Repeated, Zero );
        FSIOperator::vector_Type velocityTn( M_physicalSolver->dFESpace().map(), Repeated, Zero );

        M_physicalSolver->exportSolidDisplacement( displacementTn );
        M_physicalSolver->exportSolidVelocity( velocityTn );

        *M_robinRHS = M_physicalSolver->solidTimeAdvance()->rhsContributionFirstDerivative();

        break;
    }
    default:
        break;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< FSIOperator >::setData( const BCInterfaceData3D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::setData" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, FSIFunction > mapFunction;
    mapFunction["DerFluidLoadToFluid"]              = DerFluidLoadToFluid;
    mapFunction["DerFluidLoadToStructure"]          = DerFluidLoadToStructure;
    mapFunction["DerHarmonicExtensionVelToFluid"]   = DerHarmonicExtensionVelToFluid;
    mapFunction["DerStructureDispToSolid"]          = DerStructureDispToSolid;
    mapFunction["FluidInterfaceDisp"]               = FluidInterfaceDisp;
    mapFunction["FluidLoadToStructure"]             = FluidLoadToStructure;
    mapFunction["HarmonicExtensionVelToFluid"]      = HarmonicExtensionVelToFluid;
    mapFunction["SolidLoadToStructure"]             = SolidLoadToStructure;
    mapFunction["StructureDispToHarmonicExtension"] = StructureDispToHarmonicExtension;
    mapFunction["StructureDispToSolid"]             = StructureDispToSolid;
    mapFunction["StructureToFluid"]                 = StructureToFluid;
    mapFunction["RobinWall"]                        = RobinWall;

    // Retrieving the strings
    M_FSIFunction = mapFunction[ data.baseString() ];

    M_name = data.name();
    M_flag = data.flag();
    M_type = data.type();
    M_mode = data.mode();
    M_componentsVector = data.componentsVector();

    if ( M_FSIFunction == RobinWall )
    {
        factory_Type factory;
        M_vectorFunctionRobin.reserve(2);
        BCInterfaceData3D temporaryData ( data );

        // Create the mass term function
        temporaryData.setRobinBaseAlpha();
        M_vectorFunctionRobin.push_back( factory.createFunctionParser( temporaryData ) );

        // Create the RHS
        temporaryData.setRobinBaseBeta();
        M_vectorFunctionRobin.push_back( factory.createFunctionParser( temporaryData ) );
    }
}





// ===================================================
// Constructors
// ===================================================
BCInterfaceFunctionSolverDefined< OneDFSISolver >::BCInterfaceFunctionSolverDefined() :
        M_defaultFunction (),
        M_function        ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::BCInterfaceFunctionSolverDefined()" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< OneDFSISolver >::assignFunction( OneDFSIFunction& base )
{
    switch ( M_defaultFunction )
    {
    case Riemann:

        base.setFunction( boost::bind( &OneDFSIFunctionSolverDefinedRiemann::operator(),
                                       dynamic_cast<OneDFSIFunctionSolverDefinedRiemann *> ( &( *M_function ) ), _1, _2 ) );

        break;

    case Compatibility:

        base.setFunction( boost::bind( &OneDFSIFunctionSolverDefinedCompatibility::operator(),
                                       dynamic_cast<OneDFSIFunctionSolverDefinedCompatibility *> ( &( *M_function ) ), _1, _2 ) );

        break;

    case Absorbing:

        base.setFunction( boost::bind( &OneDFSIFunctionSolverDefinedAbsorbing::operator(),
                                       dynamic_cast<OneDFSIFunctionSolverDefinedAbsorbing *> ( &( *M_function ) ), _1, _2 ) );

        break;

    case Resistance:

        base.setFunction( boost::bind( &OneDFSIFunctionSolverDefinedResistance::operator(),
                                       dynamic_cast<OneDFSIFunctionSolverDefinedResistance *> ( &( *M_function ) ), _1, _2 ) );

        break;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< OneDFSISolver >::setData( const BCInterfaceData1D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::setData( data )" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, solverDefinedFunctions > mapFunction;
    mapFunction["Riemann"]       = Riemann;
    mapFunction["Compatibility"] = Compatibility;
    mapFunction["Absorbing"]     = Absorbing;
    mapFunction["Resistance"]    = Resistance;

    M_defaultFunction = mapFunction[data.baseString()];

    switch ( M_defaultFunction )
    {
    case Riemann:

        M_function.reset( new OneDFSIFunctionSolverDefinedRiemann( data.side(), data.quantity() ) );

        break;

    case Compatibility:

        M_function.reset( new OneDFSIFunctionSolverDefinedCompatibility( data.side(), data.quantity() ) );

        break;

    case Absorbing:

        M_function.reset( new OneDFSIFunctionSolverDefinedAbsorbing( data.side(), data.quantity() ) );

        break;

    case Resistance:

        M_function.reset( new OneDFSIFunctionSolverDefinedResistance( data.side(), data.quantity(), data.resistance()[0] ) );

        break;
    }
}

} // Namespace LifeV
