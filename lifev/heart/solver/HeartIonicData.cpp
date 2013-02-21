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
    @brief File containing a class for handling Ionic model data with GetPot

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>
    @author Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>, J.Castelneau (INRIA)
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#include <lifev/heart/solver/HeartIonicData.hpp>


namespace LifeV
{


// ===================================================
// Constructors & Destructor
// ===================================================
//! Constructors
HeartIonicData::HeartIonicData ( const GetPot& dataFile ) :
    MeshData ( dataFile, "electric/space_discretization" ),
    TimeData ( dataFile, "electric/time_discretization" )
{
    setup (dataFile);
}

HeartIonicData::HeartIonicData() :
    MeshData                        ( ),
    TimeData                        ( ),
    M_verbose                       ( ),
    M_RMCParameterA                 ( ),
    M_RMCParameterB                 ( ),
    M_RMCParameterC1                ( ),
    M_RMCParameterC2                ( ),
    M_RMCParameterD                 ( ),
    M_RMCTimeUnit                   ( ),
    M_RMCPotentialAmplitude         ( ),
    M_RMCRestPotential              ( ),
    M_RMCInitialRepolarization      ( ),
    // Mitchell & Schaeffer
    M_MSTauIn                       ( ),
    M_MSTauOut                       ( ),
    M_MSTauOpen                     ( ),
    M_MSTauClose                    ( ),
    M_MSCriticalPotential           ( ),
    M_MSPotentialMinimum            ( ),
    M_MSPotentialMaximum            ( ),
    M_MSReactionAmplitude           ( ),
    M_MSInitialTime                 ( ),
    M_MSTend                        ( ),
    M_MSBDForder                    ( ),
    M_MSHasHeterogeneousTauClose    ( ),
    //Minimal model:
    M_MinimalEpitheta0 ( ),
    M_MinimalEpitheta1 ( ),
    M_MinimalEpitheta2 ( ),
    M_MinimalEpitheta1minus ( ),
    M_MinimalEpitau301 ( ),
    M_MinimalEpitau302 ( ),
    M_MinimalEpitau31 ( ),
    M_MinimalEpitau32 ( ),
    M_MinimalEpitau01 ( ),
    M_MinimalEpitau02 ( ),
    M_MinimalEpiwStar2inf ( ),
    M_MinimalEpivv ( ),
    M_MinimalEpitau11minus ( ),
    M_MinimalEpitau12minus ( ),
    M_MinimalEpitau1plus ( ),
    M_MinimalEpitau21minus ( ),
    M_MinimalEpitau22minus ( ),
    M_MinimalEpitau2plus ( ),
    M_MinimalEpik2minus ( ),
    M_MinimalEpiv2minus ( ),
    M_MinimalEpitauFi ( ),
    M_MinimalEpik30 ( ),
    M_MinimalEpiv30 ( ),
    M_MinimalEpik3 ( ),
    M_MinimalEpiv3 ( ),
    M_MinimalEpitauSi ( ),
    M_MinimalEpitau2inf ( )
{
}

HeartIonicData::HeartIonicData ( const HeartIonicData& dataIonic ) :
    MeshData                        ( dataIonic ),
    TimeData                        ( dataIonic ),
    M_verbose                       ( dataIonic.M_verbose ),
    M_RMCParameterA                 ( dataIonic.M_RMCParameterA ),
    M_RMCParameterB                 ( dataIonic.M_RMCParameterB ),
    M_RMCParameterC1                ( dataIonic.M_RMCParameterC1 ),
    M_RMCParameterC2                ( dataIonic.M_RMCParameterC2 ),
    M_RMCParameterD                 ( dataIonic.M_RMCParameterD ),
    M_RMCTimeUnit                   ( dataIonic.M_RMCTimeUnit ),
    M_RMCPotentialAmplitude         ( dataIonic.M_RMCPotentialAmplitude ),
    M_RMCRestPotential              ( dataIonic.M_RMCRestPotential ),
    M_RMCInitialRepolarization      ( dataIonic.M_RMCInitialRepolarization ),
    // Mitchell & Schaeffer
    M_MSTauIn                       ( dataIonic.M_MSTauIn ),
    M_MSTauOut                       ( dataIonic.M_MSTauOut ),
    M_MSTauOpen                     ( dataIonic.M_MSTauOpen ),
    M_MSTauClose                    ( dataIonic.M_MSTauClose ),
    M_MSCriticalPotential           ( dataIonic.M_MSCriticalPotential ),
    M_MSPotentialMinimum            ( dataIonic.M_MSPotentialMinimum ),
    M_MSPotentialMaximum            ( dataIonic.M_MSPotentialMaximum ),
    M_MSReactionAmplitude           ( dataIonic.M_MSReactionAmplitude ),
    M_MSInitialTime                 ( dataIonic.M_MSInitialTime ),
    M_MSTend                        ( dataIonic.M_MSTend ),
    M_MSBDForder                    ( dataIonic.M_MSBDForder ),
    M_MSHasHeterogeneousTauClose    ( dataIonic.M_MSHasHeterogeneousTauClose ),
    M_MinimalEpitheta0 (dataIonic.M_MinimalEpitheta0 ),
    M_MinimalEpitheta1 (dataIonic.M_MinimalEpitheta1 ),
    M_MinimalEpitheta2 (dataIonic.M_MinimalEpitheta2 ),
    M_MinimalEpitheta1minus (dataIonic.M_MinimalEpitheta1minus ),
    M_MinimalEpitau301 (dataIonic.M_MinimalEpitau301 ),
    M_MinimalEpitau302 (dataIonic.M_MinimalEpitau302 ),
    M_MinimalEpitau31 (dataIonic.M_MinimalEpitau31 ),
    M_MinimalEpitau32 (dataIonic.M_MinimalEpitau32 ),
    M_MinimalEpitau01 (dataIonic.M_MinimalEpitau01 ),
    M_MinimalEpitau02 (dataIonic.M_MinimalEpitau02 ),
    M_MinimalEpiwStar2inf (dataIonic.M_MinimalEpiwStar2inf ),
    M_MinimalEpivv (dataIonic.M_MinimalEpivv ),
    M_MinimalEpitau11minus (dataIonic.M_MinimalEpitau11minus ),
    M_MinimalEpitau12minus (dataIonic.M_MinimalEpitau12minus ),
    M_MinimalEpitau1plus (dataIonic.M_MinimalEpitau1plus ),
    M_MinimalEpitau21minus (dataIonic.M_MinimalEpitau21minus ),
    M_MinimalEpitau22minus (dataIonic.M_MinimalEpitau22minus ),
    M_MinimalEpitau2plus (dataIonic.M_MinimalEpitau2plus ),
    M_MinimalEpik2minus (dataIonic.M_MinimalEpik2minus ),
    M_MinimalEpiv2minus (dataIonic.M_MinimalEpiv2minus ),
    M_MinimalEpitauFi (dataIonic.M_MinimalEpitauFi ),
    M_MinimalEpik30 (dataIonic.M_MinimalEpik30 ),
    M_MinimalEpiv30 (dataIonic.M_MinimalEpiv30 ),
    M_MinimalEpik3 (dataIonic.M_MinimalEpik3 ),
    M_MinimalEpiv3 (dataIonic.M_MinimalEpiv3 ),
    M_MinimalEpitauSi (dataIonic.M_MinimalEpitauSi ),
    M_MinimalEpitau2inf (dataIonic.M_MinimalEpitau2inf )

{
}


// ===================================================
// Methods
// ===================================================
HeartIonicData&
HeartIonicData::operator= ( const HeartIonicData& dataIonic )
{
    if ( this != &dataIonic )
    {
        M_MSHasHeterogeneousTauClose    = dataIonic.M_MSHasHeterogeneousTauClose;
        M_RMCParameterA                           = dataIonic.M_RMCParameterA;
        M_RMCParameterB                           = dataIonic.M_RMCParameterB;
        M_RMCParameterC1                          = dataIonic.M_RMCParameterC1;
        M_RMCParameterC2                          = dataIonic.M_RMCParameterC2;
        M_RMCParameterD                           = dataIonic.M_RMCParameterD;
        M_RMCTimeUnit                    = dataIonic.M_RMCTimeUnit;
        M_RMCPotentialAmplitude          = dataIonic.M_RMCPotentialAmplitude;
        M_RMCRestPotential               = dataIonic.M_RMCRestPotential;
        M_RMCInitialRepolarization       = dataIonic.M_RMCInitialRepolarization;
        // Mitchell & Schaeffer
        M_MSTauIn                      = dataIonic.M_MSTauIn;
        M_MSTauOut                     = dataIonic.M_MSTauOut;
        M_MSTauOpen                    = dataIonic.M_MSTauOpen;
        M_MSTauClose                   = dataIonic.M_MSTauClose;
        M_MSCriticalPotential           = dataIonic.M_MSCriticalPotential;
        M_MSPotentialMinimum            = dataIonic.M_MSPotentialMinimum;
        M_MSPotentialMaximum            = dataIonic.M_MSPotentialMaximum;
        M_MSReactionAmplitude           = dataIonic.M_MSReactionAmplitude;
        M_MSInitialTime                 = dataIonic.M_MSInitialTime;
        M_MSTend                        = dataIonic.M_MSTend;
        M_MSBDForder                    = dataIonic.M_MSBDForder;

        /*/Minimal model
        M_MinimalEpitheta0 = dataIonic.M_MinimalEpitheta0;
        M_MinimalEpitheta1 = dataIonic.M_MinimalEpitheta1;
        M_MinimalEpitheta2 = dataIonic.M_MinimalEpitheta2;
        M_MinimalEpitheta1minus = dataIonic.M_MinimalEpitheta1minus;
        M_MinimalEpitau301  = dataIonic.M_MinimalEpitau301;
        M_MinimalEpitau302 = dataIonic.M_MinimalEpitau302;
        M_MinimalEpitau31 = dataIonic.M_MinimalEpitau31;
        M_MinimalEpitau32 = dataIonic.M_MinimalEpitau32;
        M_MinimalEpitau01 = dataIonic.M_MinimalEpitau01;
        M_MinimalEpitau02 = dataIonic.M_MinimalEpitau02;
        M_MinimalEpiwStar2inf = dataIonic.M_MinimalEpiwStar2inf;
        M_MinimalEpivv = dataIonic.M_MinimalEpivv;
        M_MinimalEpitau11minus = dataIonic.M_MinimalEpitau11minus;
        M_MinimalEpitau12minus = dataIonic.M_MinimalEpitau12minus;
        M_MinimalEpitau1plus = dataIonic.M_MinimalEpitau1plus;
        M_MinimalEpitau21minus = dataIonic.M_MinimalEpitau21minus;
        M_MinimalEpitau22minus = dataIonic.M_MinimalEpitau22minus;
        M_MinimalEpitau2plus = dataIonic.M_MinimalEpitau2plus;
        M_MinimalEpik2minus = dataIonic.M_MinimalEpik2minus;
        M_MinimalEpiv2minus = dataIonic.M_MinimalEpiv2minus;
        M_MinimalEpitauFi = dataIonic.M_MinimalEpitauFi;
        M_MinimalEpik30 = dataIonic.M_MinimalEpik30;
        M_MinimalEpiv30 = dataIonic.M_MinimalEpiv30;
        M_MinimalEpik3 = dataIonic.M_MinimalEpik3;
        M_MinimalEpiv3 = dataIonic.M_MinimalEpiv3;
        M_MinimalEpitauSi = dataIonic.M_MinimalEpitauSi;
        M_MinimalEpitau2inf = dataIonic.M_MinimalEpitau2inf;
        */
    }
    return *this;
}

void
HeartIonicData::setup (  const GetPot& dataFile )
{

    M_MSHasHeterogeneousTauClose = dataFile ( "electric/physics/hasHeteroTauClose", 1 );
    M_RMCParameterA                        = dataFile ( "electric/physics/a", 0.13 ); // 0.13  adim  //RogersMcCulloch1994
    M_RMCParameterB                        = dataFile ( "electric/physics/b", 0.013 ); // 0.013 adim //RogersMcCulloch1994
    M_RMCParameterC1                       = dataFile ( "electric/physics/c1", 0.26 ); // 0.26  adim //RogersMcCulloch1994
    M_RMCParameterC2                       = dataFile ( "electric/physics/c2", 0.1 ); //0.1    adim //RogersMcCulloch1994
    M_RMCParameterD                        = dataFile ( "electric/physics/d", 1 );    //1      adim //RogersMcCulloch1994
    M_RMCTimeUnit                  = dataFile ( "electric/physics/T", 0.63 );  //0.63ms    //RogersMcCulloch1994
    M_RMCPotentialAmplitude        = dataFile ( "electric/physics/A", 110 );  //130mV    //RogersMcCulloch1994
    M_RMCRestPotential            = dataFile ( "electric/physics/u0", -84.0 );    //-84mV    //RogersMcCulloch1994
    M_RMCInitialRepolarization    = dataFile ( "electric/physics/winit", 0 );
    // Mitchell & Schaeffer
    M_MSTauIn                   = dataFile ( "electric/physics/tau_in", 0.8 );
    M_MSPotentialMinimum         = dataFile ( "electric/physics/v_min", -80.0 );
    M_MSPotentialMaximum         = dataFile ( "electric/physics/v_max", 20.0 );
    M_MSReactionAmplitude        = dataFile ( "electric/physics/reac_amp", 0.2 );
    M_MSTauOut                  = dataFile ( "electric/physics/tau_out", 18.0 );
    M_MSTauOpen                 = dataFile ( "electric/physics/tau_open", 100.0 );
    M_MSTauClose                = dataFile ( "electric/physics/tau_close", 100.0 );
    M_MSCriticalPotential        = dataFile ( "electric/physics/vcrit", -67.0 );
    M_MSInitialTime              = dataFile ( "electric/physics/init_time", 0.0 );
    M_MSTend                     = dataFile ( "electric/physics/end_time", 1000.0 );
    M_MSBDForder                 = dataFile ( "electric/time_discretization/BDF_order", 1 );

    //Minimal model
    M_MinimalEpitheta0 = dataFile ("electric/physics/Epitheta0", 0.005);
    M_MinimalEpitheta1 = dataFile ("electric/physics/Epitheta1", 0.3 );
    M_MinimalEpitheta2 = dataFile ("electric/physics/Epitheta2", 0.13 );
    M_MinimalEpitheta1minus = dataFile ("electric/physics/Epitheta1minus", 0.1 );
    M_MinimalEpitau301  = dataFile ("electric/physics/Epitau301", 91. );
    M_MinimalEpitau302 = dataFile ("electric/physics/Epitau302", 0.8 );
    M_MinimalEpitau31 = dataFile ("electric/physics/Epitau31", 2.7342 );
    M_MinimalEpitau32 = dataFile ("electric/physics/Epitau32", 4. );
    M_MinimalEpitau01 = dataFile ("electric/physics/Epitau01", 410. );
    M_MinimalEpitau02 = dataFile ("electric/physics/Epitau02", 7. );
    M_MinimalEpiwStar2inf = dataFile ("electric/physics/EpiwStar2inf", 0.5 );
    M_MinimalEpivv = dataFile ("electric/physics/Epivv", 1.61);
    M_MinimalEpitau11minus = dataFile ("electric/physics/Epitau11minus", 80. );
    M_MinimalEpitau12minus = dataFile ("electric/physics/Epitau12minus", 1.4506);
    M_MinimalEpitau1plus = dataFile ("electric/physics/Epitau1plus", 1.4506 );
    M_MinimalEpitau21minus = dataFile ("electric/physics/Epitau21minus", 70.);
    M_MinimalEpitau22minus = dataFile ("electric/physics/Epitau22minus", 8. );
    M_MinimalEpitau2plus = dataFile ("electric/physics/Epitau2plus", 280. );
    M_MinimalEpik2minus = dataFile ("electric/physics/Epik2minus", 200. );
    M_MinimalEpiv2minus = dataFile ("electric/physics/Epiv2minus", 0.016 );
    M_MinimalEpitauFi = dataFile ("electric/physics/EpitauFi", 0.078);
    M_MinimalEpik30 = dataFile ("electric/physics/Epik30", 2.1 );
    M_MinimalEpiv30 = dataFile ("electric/physics/Epiv30", 0.6 );
    M_MinimalEpik3 = dataFile ("electric/physics/Epik3", 2.0994 );
    M_MinimalEpiv3 = dataFile ("electric/physics/Epiv3", 0.9087 );
    M_MinimalEpitauSi = dataFile ("electric/physics/EpitauSi", 3.3849 );
    M_MinimalEpitau2inf = dataFile ("electric/physics/Epitau2inf", 0.01 );

    M_subiter = dataFile ( "electric/time_discretization/subiter", 100);

    M_SLinitial = dataFile ("electric/physics/SLinitial", 1.9 );
    M_SLrest = dataFile ("electric/physics/SLinitial", 1.9 );
    M_SLset = dataFile ("electric/physics/SLinitial", 1.9 );
}



void HeartIonicData::showMe ( std::ostream& /*output*/ )
{
}

}

