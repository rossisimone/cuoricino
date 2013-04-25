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
	  @date 04-2013
	  @author Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

	  @contributors
	  @mantainer Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
	  @last update 03-2013
	 */

#ifndef STIMULATION_HPP_
#define STIMULATION_HPP_

	#include <Teuchos_RCP.hpp>
	#include <Teuchos_ParameterList.hpp>
	#include "Teuchos_XMLParameterListHelpers.hpp"

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
		inline void setTimeSt(const Real& timeSt)
		{
			this->M_timeSt = timeSt;
		}

		inline const Real& stInt() const
		{
			return M_stInt;
		}
		inline void setStInt(const Real& stInt)
		{
			this->M_stInt = stInt;
		}

		inline const Real& stIntS1S2() const
		{
			return M_stIntS1S2;
		}
		inline void setStIntS1S2(const Real& stIntS1S2 )
		{
			this->M_stIntS1S2 = stIntS1S2;
		}

		inline const Real& stIntS2S3() const
		{
			return M_stIntS2S3;
		}
		inline void setStIntS2S3(const Real& stIntS2S3)
		{
			this->M_stIntS2S3 = stIntS2S3;
		}

		inline const Real& stIntS3S4() const
		{
			return M_stIntS3S4;
		}
		inline void setStIntS3S4(const Real& stIntS3S4)
		{
			this->M_stIntS3S4 = stIntS3S4;
		}

		inline const std::string& pacPro() const
		{
			return M_pacingProtocol;
		}
		inline void setPacingProtocol(const std::string& pacPro)
		{
			this->M_pacingProtocol = pacPro;
		}

		inline const std::string& pacProTyp() const
		{
			return M_pacingProtocolType;
		}
		inline void setPacingProtocolType(const std::string& pacProTyp)
		{
			this->M_pacingProtocolType = pacProTyp;
		}

		inline const Real& currStim() const
		{
			return M_Istim;
		}
		inline void setIstim(const Real& currStim)
		{
			this->M_Istim = currStim;
		}


		//@}

		//! @name Methods
		//@{

		// Choice of the type or protocol use. It will choose one of the protocol defined
		// in the methods below

		void pacingProtocolChoice( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

		// Protocol methods

		void fixedCycleLength( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

		void fixedCycleLengthwExtraStim( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp );

		void showMe();

		//@}

	private:

		Real M_timeSt;
		Real M_stInt;
		Real M_stIntS1S2;
		Real M_stIntS2S3;
		Real M_stIntS3S4;

		// Choice of the protocol
		std::string M_pacingProtocol;
		std::string M_pacingProtocolType;

		// Value of the current stimulation
		Real M_Istim;

	};

	// ===================================================
	//! Constructors
	// ===================================================

	StimulationProtocol::StimulationProtocol() :
		M_timeSt             ( 20.0 ),
		M_stInt              ( 600.0 ),
		M_stIntS1S2          ( 350.0 ),
		M_stIntS2S3          ( 250.0 ),
		M_stIntS3S4          ( 150.0 ),
		M_pacingProtocol     ( "FCL" ) ,
		M_pacingProtocolType ( "" ),
		M_Istim              ( 0 )
	{}

	StimulationProtocol::StimulationProtocol( Teuchos::ParameterList& parameterList )
	{
		M_timeSt             = parameterList.get ( "stimuliTime", 20.0 );
		M_stInt     	     = parameterList.get ( "stInt", 600.0 );
		M_stIntS1S2          = parameterList.get ( "stIntS1S2", 350.0 );
		M_stIntS2S3          = parameterList.get ( "stIntS1S2S3", 250.0 );
		M_stIntS3S4          = parameterList.get ( "stIntS2S3S4", 150.0 );
		M_pacingProtocol     = parameterList.get ( "pacPro", 150.0 );
		M_pacingProtocolType = parameterList.get ( "pacProType", 150.0 );
		M_Istim              = parameterList.get ( "iStim", 0 );
	}

	StimulationProtocol::StimulationProtocol ( const StimulationProtocol& protocol )
	{
		M_timeSt             = protocol.M_timeSt;
		M_stInt              = protocol.M_stInt;
		M_stIntS1S2          = protocol.M_stIntS1S2;
		M_stIntS2S3          = protocol.M_stIntS2S3;
		M_stIntS3S4          = protocol.M_stIntS3S4;
		M_pacingProtocol     = protocol.M_pacingProtocol;
		M_pacingProtocolType = protocol.M_pacingProtocolType;
		M_Istim              = protocol.M_Istim;
	}

	// ===================================================
	//! Operator
	// ===================================================

	StimulationProtocol& StimulationProtocol::operator= ( const StimulationProtocol& protocol )
	{
		M_timeSt             = protocol.M_timeSt;
		M_stInt              = protocol.M_stInt;
		M_stIntS1S2          = protocol.M_stIntS1S2;
		M_stIntS2S3          = protocol.M_stIntS2S3;
		M_stIntS3S4          = protocol.M_stIntS3S4;
		M_pacingProtocol     = protocol.M_pacingProtocol;
		M_pacingProtocolType = protocol.M_pacingProtocolType;
		M_Istim              = protocol.M_Istim;

		return *this;
	}

	// ===================================================
	//! Methods
	// ===================================================


	void StimulationProtocol::pacingProtocolChoice( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
	{
		// Choice the protocol according to the parameter passed.
		// Is this method that must be used in the main.cpp

		if ( M_pacingProtocol == "FCL" )
			fixedCycleLength ( t, dt, NbStimulus, Iapp );
		if ( M_pacingProtocol == "FCL-ExtraSt" )
			fixedCycleLengthwExtraStim ( t, dt, NbStimulus, Iapp );
	}

	void StimulationProtocol::fixedCycleLength( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
	{
		if ( t >= M_timeSt && t <= M_timeSt + 1.0 )
		{
			Iapp = M_Istim;

			if ( t >= M_timeSt + 1.0 - dt && t <= M_timeSt + 1.0 )
			{
				NbStimulus++;
				cout << "\n NbStimulus: " << NbStimulus << "\n";
				M_timeSt = M_timeSt + M_stInt;
			}
			else
				Iapp = 0;
		}
	}

	void StimulationProtocol::fixedCycleLengthwExtraStim( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
	{
		if ( t >= M_timeSt && t <= M_timeSt + 1.0 )
		{
			Iapp = M_Istim;

			if ( t >= M_timeSt + 1.0 - dt && t <= M_timeSt + 1.0 )
			{
				NbStimulus++;

				if ( NbStimulus < 8 )
					M_timeSt = M_timeSt + M_stInt;

				else
				{
					if ( M_pacingProtocolType == "S1-S2" )
						M_timeSt = M_timeSt + M_stIntS1S2;

					else if ( M_pacingProtocolType == "S1-S2-S3" )
					{
						if ( NbStimulus < 9 )
							M_timeSt = M_timeSt + M_stIntS1S2;
						else
							M_timeSt = M_timeSt + M_stIntS2S3;
					}

					else if ( M_pacingProtocolType == "S1-S2-S3-S4" )
					{
						if ( NbStimulus < 9 )
							M_timeSt = M_timeSt + M_stIntS1S2;
						else if ( NbStimulus < 10 )
							M_timeSt = M_timeSt + M_stIntS2S3;
						else
							M_timeSt = M_timeSt + M_stIntS3S4;
					}

					else
						M_timeSt = M_timeSt + M_stInt;
				}

			}
		}
		else
			Iapp = 0;
	}



	void StimulationProtocol::showMe()
		{
			std::cout << "\n\n\t\tStimulation protocol Informations\n\n";

			std::cout << "\n\t\tList of parameters:\n\n";
			std::cout << "1st stimuli time: " << this->timeSt() << std::endl;
			std::cout << "Stimulation interval: " << this->stInt() << std::endl;
			std::cout << "S1-S2 interval: " << this->stIntS1S2() << std::endl;
			std::cout << "S1-S2-S3 interval: " << this->stIntS2S3() << std::endl;
			std::cout << "S1-S2-S3-S4 interval: " << this->stIntS3S4() << std::endl;
			std::cout << "Pacing protocol: " << this->pacPro() << std::endl;
			std::cout << "Pacing protocol type: " << this->pacProTyp() << std::endl;
			std::cout << "Istim: " << this->currStim() << std::endl;


			std::cout << "\n\t\t End of Stimulation protocol Informations\n\n\n";
		}


	}



#endif /* STIMULATION_HPP_ */
