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
  @brief Stimulation protocol class with different protocol methods
  @ List of protocols parameters used in the different models:
  @ PACPRO
  *********
  FCL -> Fixed Cycle Length
  FCL-ExtraSt -> Fixed Cycle Length with premature or long-short
                 extrastimulation.
                 The choice of the type depends only on the parameters in the xml.
  S1S2Pro -> Standard S1-S2 protocol
  DynPro  -> Dynamic protocol

  PACPROTYPE
  *************
  S1-S2       -> Single Extrastimuli interval
  S1-S2-S3    -> Double Extrastimuli interval
  S1-S2-S3-S4 -> Triple Extrastimuli interval

  @
  @date 04-2013
  @author Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

  @contributors
  @mantainer Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
  @last update 03-2013
 */

#ifndef STIMULATIONPROTOCOL_HPP_
#define STIMULATIONPROTOCOL_HPP_

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/core/LifeV.hpp>

#include <string>

namespace LifeV
{
//! StimulatioProtocol - This class implements the different types of stimulation protocol

class StimulationProtocol
{

public:
    //! @name Constructors & Destructor
    //@{

    //! Constructor

    StimulationProtocol();

    StimulationProtocol ( Teuchos::ParameterList& parameterList );

    /*!
    * @param StimulationProtocol object
    */
    StimulationProtocol ( const StimulationProtocol& protocol );

    //! Destructor
    virtual ~StimulationProtocol() {}

    //@}

    //! @name Overloads
    //@{

    StimulationProtocol& operator= ( const StimulationProtocol& protocol );

    //@}


    //! @name Setters and getters
    //@{

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

    inline const std::string& pacPro() const
    {
        return M_pacingProtocol;
    }
    inline void setPacingProtocol (const std::string& pacPro)
    {
        this->M_pacingProtocol = pacPro;
    }

    inline const std::string& pacProTyp() const
    {
        return M_pacingProtocolType;
    }
    inline void setPacingProtocolType (const std::string& pacProTyp)
    {
        this->M_pacingProtocolType = pacProTyp;
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

    //@}

    //! @name Methods
    //@{

    // Choice of the type or protocol use. It will choose one of the protocol defined
    // in the methods below

    void pacingProtocolChoice ( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

    void showMe();

    //@}

private:

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

    // Choice of the protocol
    std::string M_pacingProtocol;
    std::string M_pacingProtocolType;

    // Value of the current stimulation
    Real M_Istim;
    Real M_StimDuration;

    // Protocol methods

    // Stimulation at constant frequency
    void fixedCycleLength ( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

    // Stimulation at constant frequency with extra stimulation ( respective Si-Sj interval = constant )
    void fixedCycleLengthwExtraStim ( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

    // Standard stimulation protocol
    void standardS1S2Protocol ( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

    // Dynamic stimulation protocol
    void dynamicProtocol ( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

};

}



#endif /* STIMULATIONPROTOCOL_HPP_ */
