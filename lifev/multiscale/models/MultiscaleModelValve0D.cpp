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
 *  @brief File containing the Multiscale Valve 0D
 *
 *  @date 17-01-2014
 *  @author Toni Lassila <toni.lassila@epfl.ch>
 *
 *  @mantainer Toni Lassila <toni.lassila@epfl.ch>
 */

#include <lifev/multiscale/models/MultiscaleModelValve0D.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelValve0D::MultiscaleModelValve0D() :
                    multiscaleModel_Type           (),
                    MultiscaleInterface            (),
                    M_outputFile                   (),
                    M_bc                           ( new bcInterface_Type() ),
                    M_data                         ( new data_Type() ),
                    M_pressureLeft_tn              (),
                    M_flowRateLeft_tn              (),
                    M_pressureLeft                 (),
                    M_flowRateLeft                 (),
                    M_pressureRight                (),
                    M_tangentPressureLeft          (),
                    M_tangentFlowRateLeft          (),
                    M_openingAngle                 (),
                    M_openingAngle_tn              (),
                    M_thetaVel                     (),
                    M_thetaVel_tn                  (),
                    M_minimumOpeningAngle          (),
                    M_maximumOpeningAngle          (),
                    M_flowDischargeCoefficient     (),
                    M_frictionalMomentCoefficient  (),
                    M_resistiveMomentCoefficient   ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::MultiscaleModelValve0D() \n";
#endif

    M_type = Valve0D;
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelValve0D::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData ( fileName );

    GetPot dataFile ( fileName );

    // All parameters should be in cm - s - dyn/cm^2 - g units
    // Reference values are for human mitral valves, see Korakianitis-Shi 2006

    M_minimumOpeningAngle         = dataFile ( "Coefficients/MinimumOpeningAngle",         0.0 );   // theta_min
    M_maximumOpeningAngle         = dataFile ( "Coefficients/MaximumOpeningAngle",         75.0 );  // theta_max
    M_flowDischargeCoefficient    = dataFile ( "Coefficients/FlowDischargeCoefficient",    10.95 ); // C_D
    M_frictionalMomentCoefficient = dataFile ( "Coefficients/FrictionalMomentCoefficient", 4.125 ); // K_P
    M_resistiveMomentCoefficient  = dataFile ( "Coefficients/ResistiveMomentCoefficient",  50.0 );  // K_f

    if ( M_globalData.get() )
    {
        setupGlobalData ( fileName );
    }

    // We need to create the BCHandler before using it
    M_bc->createHandler();

    // Exporter/Importer
    setupExporterImporter();
}

void
MultiscaleModelValve0D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::setupModel() \n";
#endif

    initializeSolution();

    M_bc->setPhysicalSolver ( M_data );

    // Safety check
    if ( M_bc->handler()->bc ( 1 ).bcType() != Voltage )
    {
        std::cout << "!!! Error: the Windkessel model support only stress boundary conditions on the right at the present time !!!" << std::endl;
    }
}

void
MultiscaleModelValve0D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::buildModel() \n";
#endif

    updateModel();
}

void
MultiscaleModelValve0D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::updateModel() \n";
#endif

    M_pressureLeft_tn = M_pressureLeft;
    M_flowRateLeft_tn = M_flowRateLeft;
    M_openingAngle_tn = M_openingAngle;
    M_thetaVel_tn     = M_thetaVel;

    // Update BCInterface solver variables
    M_bc->updatePhysicalSolverVariables();

    M_pressureRight   = -M_bc->handler()->bc ( 1 ).evaluate ( M_globalData->dataTime()->time() );
}

void
MultiscaleModelValve0D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::solveModel() \n";
#endif

    displayModelStatus ( "Solve" );

    solveForOpeningAngle();

    switch ( M_bc->handler()->bc ( 0 ).bcType() )
    {
        case Current:

            M_flowRateLeft = M_bc->handler()->bc ( 0 ).evaluate ( M_globalData->dataTime()->time() );
            M_pressureLeft = solveForPressure();

            break;

        case Voltage:

            M_pressureLeft = -M_bc->handler()->bc ( 0 ).evaluate ( M_globalData->dataTime()->time() );
            M_flowRateLeft = solveForFlowRate();

            break;

        default:

            std::cout << "Warning: bcType \"" << M_bc->handler()->bc ( 0 ).bcType() << "\"not available!" << std::endl;

            break;
    }
}

void
MultiscaleModelValve0D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::updateSolution() \n";
#endif

}

void
MultiscaleModelValve0D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::saveSolution() \n";
#endif

    M_outputFile << "    " << M_globalData->dataTime()->time()
                                 << "    " << M_flowRateLeft
                                 << "    " << M_pressureLeft
                                 << "    " << M_openingAngle << std::endl;

    if ( M_globalData->dataTime()->isLastTimeStep() )
    {
        M_outputFile.close();
    }
}

void
MultiscaleModelValve0D::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();

        std::cout << "Minimum opening angle         = " << M_minimumOpeningAngle << std::endl
                        << "Maximum opening angle         = " << M_maximumOpeningAngle << std::endl
                        << "Frictional moment coefficient = " << M_frictionalMomentCoefficient << std::endl
                        << "Resistive moment coefficient  = " << M_resistiveMomentCoefficient << std::endl
                        << "Flow discharge coefficient    = " << M_flowDischargeCoefficient << std::endl << std::endl;

    }
}

Real
MultiscaleModelValve0D::checkSolution() const
{
    return M_pressureLeft + M_flowRateLeft;
}

// ===================================================
// MultiscaleInterface Methods
// ===================================================
void
MultiscaleModelValve0D::imposeBoundaryFlowRate ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    ZeroDimensionalFunction base;
    base.setFunction ( boost::bind ( function, _1, _1, _1, _1, _1 ) );

    M_bc->handler()->setBC ( boundaryFlag ( boundaryID ), Current, base );
}

void
MultiscaleModelValve0D::imposeBoundaryMeanNormalStress ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    ZeroDimensionalFunction base;
    base.setFunction ( boost::bind ( function, _1, _1, _1, _1, _1 ) );

    M_bc->handler()->setBC ( boundaryFlag ( boundaryID ), Voltage, base );
}

Real
MultiscaleModelValve0D::boundaryDeltaFlowRate ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    if ( boundaryFlag ( boundaryID ) == 1 )
    {
        return 0;
    }

    solveLinearModel ( solveLinearSystem );

    return M_tangentFlowRateLeft;
}

Real
MultiscaleModelValve0D::boundaryDeltaMeanNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    if ( boundaryFlag ( boundaryID ) == 1 )
    {
        return 0;
    }

    solveLinearModel ( solveLinearSystem );

    return -M_tangentPressureLeft;
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModelValve0D::setupGlobalData ( const std::string& fileName )
{
    GetPot dataFile ( fileName );

    //Global data time
    M_data->setTimeData ( M_globalData->dataTime() );

    //if ( !dataFile.checkVariable ( "Coefficients/VenousPressure" ) )
    //{
    //    M_pressureRight = M_globalData->fluidVenousPressure();
    //}
    //M_data->setVenousPressure ( M_pressureRight );
}

void
MultiscaleModelValve0D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::initializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        std::string fileName = multiscaleProblemFolder + multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) + "_" + number2string ( multiscaleProblemStep - 1 ) + ".mfile";

        std::ifstream inputFile;
        inputFile.open ( fileName.c_str(), std::ios::in );

        if ( inputFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;
            Real deltaT (1e15);

            // Read the first line with comments
            std::getline ( inputFile, line, '\n' );

            // Read one-by-one all the others lines of the fileName
            while ( std::getline ( inputFile, line, '\n' ) )
            {
                // Split the four entries
                boost::split ( stringsVector, line, boost::is_any_of ( " " ), boost::token_compress_on );

                // Import values
                if ( std::abs ( string2number ( stringsVector[1] ) - M_globalData->dataTime()->initialTime() ) <= deltaT )
                {
                    deltaT = std::abs ( string2number ( stringsVector[1] ) - M_globalData->dataTime()->initialTime() );

                    M_flowRateLeft = string2number ( stringsVector[2] );
                    M_pressureLeft = string2number ( stringsVector[3] );
                    M_openingAngle = string2number ( stringsVector[4] );
                }
            }

            // Close file
            inputFile.close();
        }
        else
        {
            std::cerr << " !!! Error: cannot open fileName: " << fileName.c_str() << " !!!" << std::endl;
        }
    }
    else
    {
        M_flowRateLeft = 0;
        M_pressureLeft = M_globalData->solidExternalPressure();
        M_openingAngle = 0; // Valve assumed to be closed in the start

        switch ( M_bc->handler()->bc ( 0 ).bcType() )
        {
            case Current:

                M_flowRateLeft = M_bc->handler()->bc ( 0 ).evaluate ( M_globalData->dataTime()->time() );

                break;

            case Voltage:

                //            if ( std::abs( M_globalData->solidExternalPressure() + M_bc->handler()->bc( 0 ).evaluate( M_globalData->dataTime()->time() ) ) > 1e-14 )
                //                std::cout << "!!! Warning: external pressure should be equal to the initial pressure !!! " << std::endl;

                break;

            default:

                std::cout << "Warning: bcType \"" << M_bc->handler()->bc ( 0 ).bcType() << "\"not available!" << std::endl;

                break;
        }
    }
}

void
MultiscaleModelValve0D::setupExporterImporter()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::setupExporterImporter() \n";
#endif

    std::string file = multiscaleProblemFolder + multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) + "_" + number2string ( multiscaleProblemStep ) + ".mfile";
    M_outputFile.open ( file.c_str(), std::ios::trunc );
    M_outputFile << std::scientific << std::setprecision ( 15 )
    << "%   MODEL: " << M_modelName << std::endl
    << "%   TIME                     FLOW RATE                PRESSURE                 OPENING ANGLE" << std::endl;
}

Real
MultiscaleModelValve0D::solveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::solveForFlowRate() \n";
#endif

    // Orifice model of Korakianitis-Shi '06
    Real A ( (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) * (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) / ((1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.)) * (1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.))) );

    if (M_pressureLeft > M_pressureRight)
    {
        // Outflow
        return M_flowDischargeCoefficient * A * sqrt(M_pressureLeft - M_pressureRight);
    }
    else
    {
        // Inflow
        return M_flowDischargeCoefficient * A * sqrt(M_pressureRight - M_pressureLeft);
    }

}

Real
MultiscaleModelValve0D::solveForPressure()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::solveForPressure() \n";
#endif

    // Orifice model of Korakianitis-Shi '06
    Real A ( (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) * (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) / ((1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.)) * (1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.))) );

    if ( M_flowRateLeft > 0)
    {
        // Outflow
        return M_pressureRight + M_flowRateLeft * M_flowRateLeft / ( (M_flowDischargeCoefficient * A) * (M_flowDischargeCoefficient * A) ) ;
    }
    else
    {
        // Inflow
        return M_pressureRight - M_flowRateLeft * M_flowRateLeft / ( (M_flowDischargeCoefficient * A) * (M_flowDischargeCoefficient * A) );
    }

}

void
MultiscaleModelValve0D::solveLinearModel ( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::solveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
    {
        return;
    }

    //Solve the linear problem
    displayModelStatus ( "Solve linear" );

    //Solve the opening angle first
    solveForOpeningAngle();

    //Solve the flow rate/pressure next
    switch ( M_bc->handler()->bc ( 0 ).bcType() )
    {
        case Current: // dP/dQ

            M_tangentFlowRateLeft = 1.;
            M_tangentPressureLeft = tangentSolveForPressure();

            break;

        case Voltage: // dQ/dS

            M_tangentPressureLeft = 1.;
            M_tangentFlowRateLeft = -tangentSolveForFlowRate();

            break;

        default:

            std::cout << "Warning: bcType \"" << M_bc->handler()->bc ( 0 ).bcType() << "\"not available!" << std::endl;

            break;
    }

    // This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

Real
MultiscaleModelValve0D::tangentSolveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::tangentSolveForFlowRate() \n";
#endif

    // Orifice model of Korakianitis-Shi '06
    Real A ( (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) * (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) / ((1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.)) * (1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.))) );

    if ( M_flowRateLeft > 0)
    {
        // Outflow
        return 2 * M_flowRateLeft / ( M_flowDischargeCoefficient * A * M_flowDischargeCoefficient * A ) ;
    }
    else
    {
        // Inflow
        return - 2 * M_flowRateLeft / ( M_flowDischargeCoefficient * A * M_flowDischargeCoefficient * A );
    }
}

Real
MultiscaleModelValve0D::tangentSolveForPressure()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::tangentSolveForPressure() \n";
#endif

    // Orifice model of Korakianitis-Shi '06
    Real A ( (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) * (1 - cos(M_openingAngle_tn * 2. * M_PI / 360.)) / ((1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.)) * (1 - cos(M_maximumOpeningAngle * 2. * M_PI / 360.))) );

    if (M_pressureLeft > M_pressureRight)
    {
        // Outflow
        return M_flowDischargeCoefficient * A * 0.5 / sqrt(M_pressureLeft - M_pressureRight);
    }
    else
    {
        // Inflow
        return - M_flowDischargeCoefficient * A * 0.5 * sqrt(M_pressureRight - M_pressureLeft);
    }

}

void
MultiscaleModelValve0D::solveForOpeningAngle()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelValve0D::solveForOpeningAngle() \n";
#endif

    Real dt = M_globalData->dataTime()->timeStep();

    // Leaflet moment model of Korakianitis-Shi '06
    Real F_tn ( M_frictionalMomentCoefficient * (M_pressureLeft - M_pressureRight) * std::cos(M_openingAngle_tn * 2. * M_PI / 360.) - M_resistiveMomentCoefficient * M_thetaVel_tn );
    M_thetaVel = M_thetaVel_tn + dt * F_tn;
    M_openingAngle = M_openingAngle_tn + dt * M_thetaVel;

    if (M_openingAngle > M_maximumOpeningAngle)
    {
        M_openingAngle = M_maximumOpeningAngle;
        M_thetaVel     = 0;
    }
    else if (M_openingAngle < M_minimumOpeningAngle)
    {
        M_openingAngle = M_minimumOpeningAngle;
        M_thetaVel     = 0;
    }

}

} // Namespace multiscale
} // Namespace LifeV
