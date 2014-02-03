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
 @brief Class for applying cardiac stimulus at a single point according to a pacing protocol

 @date 02-2014
 @author Simone Palamara <palamara.simone@gmail.com>

 @last update 02-2014
 */

#include <lifev/electrophysiology/util/CardiacStimulusSingleSource.hpp>

namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
CardiacStimulusPacingProtocol::CardiacStimulusPacingProtocol() :
    M_radius ( 0 ),
    M_totalCurrent ( 0 ),
    M_pacingSite_X ( 0 ),
    M_pacingSite_Y ( 0 ),
    M_pacingSite_Z ( 0 ),
    M_stimulusValue ( 0 ),
    M_startingActivationTime ( 0 ),
    M_endingActivationTime ( 0 ),
{

}


// ===================================================
//! Methods
// ===================================================


void CardiacStimulusPacingProtocol::pacingProtocolChoice ( const Real& t)
{

    if ( M_pacingProtocol == "FCL" )
    {
        fixedCycleLength ( t );
    }

    else if ( M_pacingProtocol == "FCL-ExtraSt" )
    {
        fixedCycleLengthwExtraStim ( t );
    }

    else if ( M_pacingProtocol == "S1S2Pro" )
    {
        standardS1S2Protocol ( t );
    }

    else if ( M_pacingProtocol == "DynPro" )
    {
        dynamicProtocol ( t );
    }

    else
    {
        fixedCycleLength ( t );
    }


}


void CardiacStimulusPacingProtocol::fixedCycleLength ( const Real& t )
{
    if ( M_numberStimulus < M_nbStimMax )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            M_totalCurrent = M_Istim;
            if ( t >= M_timeSt + M_StimDuration - M_dt && t <= M_timeSt + M_StimDuration )
            {
                M_numberStimulus++;
                M_timeSt = M_timeSt + M_stInt;
            }
        }
        else
        {
            M_totalCurrent = 0;
        }
    }
    else
    {
        M_totalCurrent = 0;
    }
}

void CardiacStimulusPacingProtocol::fixedCycleLengthwExtraStim ( const Real& t, Real& Iapp )
{
    if ( M_numberStimulus < M_nbStimMax )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            Iapp = M_Istim;

            if ( t >= M_timeSt + M_StimDuration - M_dt && t <= M_timeSt + M_StimDuration )
            {
                if ( M_numberStimulus < M_repeatSt )
                {
                    M_timeSt = M_timeSt + M_stInt;
                    M_numberStimulus++;
                }
                else
                {
                    if ( M_pacingProtocolType == "S1-S2" )
                    {
                        M_timeSt = M_timeSt + M_stIntS1S2;
                        M_numberStimulus = 0;
                    }
                    else if ( M_pacingProtocolType == "S1-S2-S3" )
                    {
                        if ( M_numberStimulus == M_repeatSt )
                        {
                            M_timeSt = M_timeSt + M_stIntS1S2;
                            M_numberStimulus++;
                        }
                        else
                        {
                            M_timeSt = M_timeSt + M_stIntS2S3;
                            M_numberStimulus = 0;
                        }
                    }
                    else if ( M_pacingProtocolType == "S1-S2-S3-S4" )
                    {
                        if ( M_numberStimulus == M_repeatSt )
                        {
                            M_timeSt = M_timeSt + M_stIntS1S2;
                            M_numberStimulus++;
                        }
                        else if ( M_numberStimulus == M_repeatSt + 1 )
                        {
                            M_timeSt = M_timeSt + M_stIntS2S3;
                            M_numberStimulus++;
                        }
                        else
                        {
                            M_timeSt = M_timeSt + M_stIntS3S4;
                            M_numberStimulus = 0;
                        }
                    }
                    else
                    {
                        M_timeSt = M_timeSt + M_stInt;
                    }
                }
            }
        }
        else
        {
            Iapp = 0;
        }
    }
    else
    {
        Iapp = 0;
    }
}

void CardiacStimulusPacingProtocol::standardS1S2Protocol ( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
{
    if ( t < M_nbStimMax * M_stInt )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            Iapp = M_Istim;

            if ( t >= M_timeSt + M_StimDuration - dt && t <= M_timeSt + M_StimDuration )
            {
                M_timeSt = M_timeSt + M_stInt;

                if ( t > ( M_nbStimMax - 1 ) * M_stInt && t < M_nbStimMax * M_stInt )
                {
                    NbStimulus = 0;
                }

                else
                {
                    NbStimulus++;
                }
            }
        }
        else
        {
            Iapp = 0;
        }
    }
    else
    {
        if ( M_stIntS1S2 >= M_stIntS1S2Min )
        {
            if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration)
            {
                Iapp = M_Istim;

                if ( t >= M_timeSt + M_StimDuration - dt && t <= M_timeSt + M_StimDuration )
                {
                    NbStimulus++;

                    if ( NbStimulus < M_repeatSt )
                    {
                        M_timeSt = M_timeSt + M_stInt;
                    }

                    else if ( NbStimulus == M_repeatSt )
                    {
                        M_timeSt = M_timeSt + M_stIntS1S2;
                    }

                    else if ( NbStimulus == M_repeatSt + 1 )
                    {
                        M_timeSt = M_timeSt + M_stInt;
                        NbStimulus = 0;

                        if ( M_stIntS1S2 > 1000 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 1000;
                        }

                        else if ( M_stIntS1S2 <= 1000 && M_stIntS1S2 > 300 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 50;
                        }

                        else if (  M_stIntS1S2 <= 300 &&  M_stIntS1S2 > 200 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 10;
                        }

                        else if ( M_stIntS1S2 <= 200 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 5;
                        }
                    }
                }
            }
            else
            {
                Iapp = 0;
            }
        }
        else
        {
            Iapp = 0;
        }
    }
}

void CardiacStimulusPacingProtocol::dynamicProtocol ( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
{
    if ( M_stInt >= M_stIntMin )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            Iapp = M_Istim;

            if ( t >= M_timeSt + M_StimDuration - dt && t <= M_timeSt + M_StimDuration )
            {
                NbStimulus++;
                M_timeSt = M_timeSt + M_stInt;
            }
        }
        else
        {
            Iapp = 0;
        }

        if ( t > M_tShortS1S1 )
        {
            if ( M_stInt > 1000 )
            {
                M_stInt      = M_stInt - 1000;
                M_tShortS1S1 = M_tShortS1S1 + M_stInt * 20;
            }
            else if ( M_stInt <= 1000 && M_stInt > 300 )
            {
                M_stInt      = M_stInt - 50;
                M_tShortS1S1 = M_tShortS1S1 + M_stInt * 20;
            }
            else if (  M_stInt <= 300 &&  M_stInt > 200 )
            {
                M_stInt      = M_stInt - 10;
                M_tShortS1S1 = M_tShortS1S1 + M_stInt * 20;
            }
            else if ( M_stInt <= 200 )
            {
                M_stInt      = M_stInt - 5;
                M_tShortS1S1 = M_tShortS1S1 + M_stInt * 20;
            }
        }
    }
    else
    {
        Iapp = 0;
    }
}

Real CardiacStimulusPacingProtocol::appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/ )
{

    Real current = 0.0;
    const Real volumeOfBall = (4. / 3.) * M_PI * M_radius * M_radius * M_radius;
    Real distance = std::sqrt ( (x - M_pacingSite_X) * (x - M_pacingSite_X) + (y - M_pacingSite_Y) * (y - M_pacingSite_Y) + (z - M_pacingSite_Z) * (z - M_pacingSite_Z) );
    if (distance <= M_radius && t >= (M_startingActivationTime) && t <= (M_endingActivationTime) )
    {
       pacingProtocolChoice( t );
       current += M_totalCurrent / volumeOfBall;
    }
    return current;

}

