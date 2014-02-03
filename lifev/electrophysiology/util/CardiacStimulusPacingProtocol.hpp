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

#ifndef CARDIACSTIMULUSPACINGPROTOCOL_HPP_
#define CARDIACSTIMULUSPACINGPROTOCOL_HPP_

#include <lifev/electrophysiology/util/CardiacStimulus.hpp>

namespace LifeV
{

class CardiacStimulusPacingProtocol : public CardiacStimulus
{

public:


    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    CardiacStimulusPacingProtocol();

    //! Destructor
    virtual ~CardiacStimulusPacingProtocol() {}

    //@}

    //! @name Set Methods
    //@{

    inline void setRadius ( Real r )
    {
        ASSERT (r > 0, "Invalid radius value.");
        M_radius = r;
    }

    inline void setTotalCurrent ( Real I )
    {
        ASSERT (I >= 0, "Invalid current value.");
        M_totalCurrent = I;
    }

    inline void setPacingSite ( Real x , Real y , Real z )
    {
        M_pacingSite_X = x;
        M_pacingSite_Y = y;
        M_pacingSite_Z = z;
    }

    inline void setStimulusValue ( Real stimulusValue )
    {
        M_stimulusValue = stimulusValue;
    }


    inline void setPacingProtocol ( std::string PacingProtocol )
    {
        M_pacingProtocol = PacingProtocol;
    }


    inline const Real& timeSt() const
    {
        return M_timeSt;
    }
    inline void setTimeSt (const Real& timeSt)
    {
        this->M_timeSt = timeSt;
    }

    inline const Real& timeShortS1S1() const
    {
        return M_tShortS1S1;
    }
    inline void setTimeShortS1S1 (const Real& tShortS1S1)
    {
        this->M_tShortS1S1 = tShortS1S1;
    }

    inline const Real& stInt() const
    {
        return M_stInt;
    }
    inline void setStInt (const Real& stInt)
    {
        this->M_stInt = stInt;
    }

    inline const Real& stIntMin() const
    {
        return M_stIntMin;
    }
    inline void setStIntMin (const Real& stIntMin)
    {
        this->M_stInt = stIntMin;
    }

    inline const Real& stIntS1S2() const
    {
        return M_stIntS1S2;
    }
    inline void setStIntS1S2 (const Real& stIntS1S2 )
    {
        this->M_stIntS1S2 = stIntS1S2;
    }

    inline const Real& stIntS1S2Min() const
    {
        return M_stIntS1S2Min;
    }
    inline void setStIntS1S2Min (const Real& stIntS1S2Min )
    {
        this->M_stIntS1S2Min = stIntS1S2Min;
    }

    inline const Real& stIntS2S3() const
    {
        return M_stIntS2S3;
    }
    inline void setStIntS2S3 (const Real& stIntS2S3)
    {
        this->M_stIntS2S3 = stIntS2S3;
    }

    inline const Real& stIntS3S4() const
    {
        return M_stIntS3S4;
    }
    inline void setStIntS3S4 (const Real& stIntS3S4)
    {
        this->M_stIntS3S4 = stIntS3S4;
    }

    inline const int& nbStimMax() const
    {
        return M_nbStimMax;
    }
    inline void setNbStimMax (const int& nbStimMax)
    {
        this->M_nbStimMax = nbStimMax;
    }

    inline const int& repeatSt() const
    {
        return M_repeatSt;
    }
    inline void setRepeatSt (const int& repeatSt)
    {
        this->M_repeatSt = repeatSt;
    }

    inline const Real& currStim() const
    {
        return M_Istim;
    }
    inline void setIstim (const Real& currStim)
    {
        this->M_Istim = currStim;
    }

    inline const Real& stimDuration() const
    {
        return M_StimDuration;
    }
    inline void setStimDuration (const Real& stimDur)
    {
        this->M_StimDuration = stimDur;
    }

    inline void setPacingProtocolType (const std::string& pacProTyp)
    {
        this->M_pacingProtocolType = pacProTyp;
    }

    //@}

    //! @name Copy Methods
    //@{

    //@}

    //! @name Methods
    //@{
    Real appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i );


    //@}

private:

    Real                 M_radius;
    Real                 M_totalCurrent;
    Real 		 M_pacingSite_X;
    Real 		 M_pacingSite_Y;
    Real  		 M_pacingSite_Z;
    Real 		 M_stimulusValue;

    // Values of the stimulation interval used in the protocols.
    Real M_timeSt;
    Real M_tShortS1S1;
    Real M_stInt;
    Real M_stIntMin;
    Real M_stIntS1S2;
    Real M_stIntS1S2Min;
    Real M_stIntS2S3;
    Real M_stIntS3S4;

    // Number of stimulation information
    int M_nbStimMax;
    int M_repeatSt;
    int M_numberStimulus;

    // Choice of the protocol
    std::string M_pacingProtocol;
    std::string M_pacingProtocolType;

    // Value of the current stimulation
    Real M_Istim;
    Real M_StimDuration;

    //Discretization time step used to understand when a stimulus in
    //the pacing train
    Real M_dt;

    // Stimulation at constant frequency
    void fixedCycleLength ( const Real& t );

    // Stimulation at constant frequency with extra stimulation ( respective Si-Sj interval = constant )
    void fixedCycleLengthwExtraStim ( const Real& t );

    // Standard stimulation protocol
    void standardS1S2Protocol ( const Real& t );

    // Dynamic stimulation protocol
    void dynamicProtocol ( const Real& t );

    void pacingProtocolChoice ( const Real& t );

};

} // namespace LifeV

#endif /* CARDIACSTIMULUSPACINGPROTOCOL_HPP_ */
