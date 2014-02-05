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

#include <lifev/electrophysiology/util/CardiacStimulusPacingProtocol.hpp>

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
    M_stimulusValue ( 0 )
{

}


// ===================================================
//! Methods
// ===================================================


Real CardiacStimulusPacingProtocol::pacingProtocolChoice ( const Real& t)
{

    if ( M_pacingProtocol == "FCL" )
    {
        return fixedCycleLength ( t );
    }

    else if ( M_pacingProtocol == "FCL-ExtraSt" )
    {
        return fixedCycleLengthwExtraStim ( t );
    }

    else if ( M_pacingProtocol == "S1S2Pro" )
    {
        return standardS1S2Protocol ( t );
    }

    else if ( M_pacingProtocol == "DynPro" )
    {
        return dynamicProtocol ( t );
    }

    else
    {
        return fixedCycleLength ( t );
    }


}


Real CardiacStimulusPacingProtocol::fixedCycleLength ( const Real& t )
{
    Real current = 0;
    if ( M_numberStimulus < M_nbStimMax )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            current = M_totalCurrent;
            if ( t >= M_timeSt + M_StimDuration - M_dt && t <= M_timeSt + M_StimDuration )
            {
                M_numberStimulus++;
                M_timeSt = M_timeSt + M_stInt;
            }
        }
        else
        {
            current = 0;
        }
    }
    else
    {
        current = 0;
    }
    return current;
}

Real CardiacStimulusPacingProtocol::fixedCycleLengthwExtraStim ( const Real& t )
{
    Real current = 0;
    if ( M_numberStimulus < M_nbStimMax )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            current = M_totalCurrent;

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
            current = 0;
        }
    }
    else
    {
        current = 0;
    }
    return current;
}

Real CardiacStimulusPacingProtocol::standardS1S2Protocol ( const Real& t)
{
    Real current = 0;
    if ( t < M_nbStimMax * M_stInt )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            current = M_totalCurrent;

            if ( t >= M_timeSt + M_StimDuration - M_dt && t <= M_timeSt + M_StimDuration )
            {
                M_timeSt = M_timeSt + M_stInt;

                if ( t > ( M_nbStimMax - 1 ) * M_stInt && t < M_nbStimMax * M_stInt )
                {
                    M_numberStimulus = 0;
                }

                else
                {
                    M_numberStimulus++;
                }
            }
        }
        else
        {
            current = 0;
        }
    }
    else
    {
        if ( M_stIntS1S2 >= M_stIntS1S2Min )
        {
            if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration)
            {
                current = M_totalCurrent;

                if ( t >= M_timeSt + M_StimDuration - M_dt && t <= M_timeSt + M_StimDuration )
                {
                    M_numberStimulus++;

                    if ( M_numberStimulus < M_repeatSt )
                    {
                        M_timeSt = M_timeSt + M_stInt;
                    }

                    else if ( M_numberStimulus == M_repeatSt )
                    {
                        M_timeSt = M_timeSt + M_stIntS1S2;
                    }

                    else if ( M_numberStimulus == M_repeatSt + 1 )
                    {
                        M_timeSt = M_timeSt + M_stInt;
                        M_numberStimulus = 0;

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
                current = 0;
            }
        }
        else
        {
            current = 0;
        }
    }
    return current;
}

Real CardiacStimulusPacingProtocol::dynamicProtocol ( const Real& t )
{
    Real current = 0;
    if ( M_stInt >= M_stIntMin )
    {
        if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
        {
            current = M_totalCurrent;

            if ( t >= M_timeSt + M_StimDuration - M_dt && t <= M_timeSt + M_StimDuration )
            {
                M_numberStimulus++;
                M_timeSt = M_timeSt + M_stInt;
            }
        }
        else
        {
            current = 0;
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
        current = 0;
    }
    return current;
}

Real CardiacStimulusPacingProtocol::appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/ )
{

    Real current = 0.0;
    const Real volumeOfBall = (4. / 3.) * M_PI * M_radius * M_radius * M_radius;
    Real distance = std::sqrt ( (x - M_pacingSite_X) * (x - M_pacingSite_X) + (y - M_pacingSite_Y) * (y - M_pacingSite_Y) + (z - M_pacingSite_Z) * (z - M_pacingSite_Z) );
    if (distance <= M_radius )
    {
       current += pacingProtocolChoice( t ) / volumeOfBall;
    }
    return current;

}

void CardiacStimulusPacingProtocol::showMe()
{
    std::cout << "\n\n\t\tStimulation protocol Informations\n\n";

    std::cout << "\n\t\tList of parameters:\n\n";

    std::cout << "Istim: " << M_totalCurrent << std::endl;
    std::cout << "StimDuration: " << M_StimDuration << std::endl;
    std::cout << "1st stimuli time: " << M_timeSt << std::endl;
    std::cout << "Pacing protocol: " << M_pacingProtocol << std::endl;
    std::cout << "NumberStimMax: " << M_nbStimMax << std::endl;

    if ( M_pacingProtocol == "FCL" )
    {
        std::cout << "S1-S1 interval: " << M_tShortS1S1 << std::endl;
        std::cout << "NbStimuliMax: " << M_numberStimulus << std::endl;
    }
//    else if ( M_pacingProtocol == "FCL-ExtraSt" )
//    {
//        std::cout << "Pacing protocol type: " << this->pacProTyp() << std::endl;
//        std::cout << "S1-S1 interval: " << this->stInt() << std::endl;
//        std::cout << "NbStimuliMax: " << this ->nbStimMax() << std::endl;
//        std::cout << "Repeat S1 stimuli: " << this->repeatSt() << std::endl;
//
//        if ( M_pacingProtocolType == "S1-S2" )
//        {
//            std::cout << "S1-S2 interval: " << this->stIntS1S2() << std::endl;
//        }
//        else if ( M_pacingProtocolType == "S1-S2-S3" )
//        {
//            std::cout << "S1-S2 interval: " << this->stIntS1S2() << std::endl;
//            std::cout << "S2-S3 interval: " << this->stIntS2S3() << std::endl;
//        }
//        else if ( M_pacingProtocolType == "S1-S2-S3-S4" )
//        {
//            std::cout << "S1-S2 interval: " << this->stIntS1S2() << std::endl;
//            std::cout << "S2-S3 interval: " << this->stIntS2S3() << std::endl;
//            std::cout << "S3-S4 interval: " << this->stIntS3S4() << std::endl;
//        }
//
//    }
//    else if ( M_pacingProtocol == "S1S2Pro" )
//    {
//        std::cout << "S1-S1 interval: " << this->stInt() << std::endl;
//        std::cout << "NbStimuliMax for stabilisation: " << this ->nbStimMax() << std::endl;
//        std::cout << "First S1-S2 interval: " << this->stIntS1S2() << std::endl;
//        std::cout << "Minimum S1-S2 interval: " << this->stIntS1S2Min() << std::endl;
//        std::cout << "Repeat S1 stimuli: " << this->repeatSt() << std::endl;
//    }
//    else if ( M_pacingProtocol == "DynPro" )
//    {
//        std::cout << "First S1-S1 interval: " << this->stInt() << std::endl;
//        std::cout << "Minimum S1-S1 interval: " << this->stIntMin() << std::endl;
//        std::cout << "First time S1-S1 interval decrease: " << this->timeShortS1S1() << std::endl;
//    }
//    else
//    {
//        std::cout << "S1-S1 interval: " << this->stInt() << std::endl;
//        std::cout << "NbStimuliMax: " << this ->nbStimMax() << std::endl;
//    }



    std::cout << "\n\t\t End of Stimulation protocol Informations\n\n\n";
}

}
