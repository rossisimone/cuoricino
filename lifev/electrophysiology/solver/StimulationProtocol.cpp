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


	#include <lifev/electrophysiology/solver/StimulationProtocol.hpp>

	namespace LifeV
	{


		// ===================================================
		//! Constructors
		// ===================================================

		StimulationProtocol::StimulationProtocol() :
			M_timeSt             ( 20.0 ),
			M_tShortS1S1         ( 1000.0 ),
			M_stInt              ( 600.0 ),
			M_stIntMin           ( 600.0 ),
			M_stIntS1S2          ( 350.0 ),
			M_stIntS1S2Min       ( 200.0 ),
			M_stIntS2S3          ( 250.0 ),
			M_stIntS3S4          ( 150.0 ),
			M_nbStimMax          ( 20 ),
			M_repeatSt           ( 10 ),
			M_pacingProtocol     ( "FCL" ) ,
			M_pacingProtocolType ( "" ),
			M_Istim              ( 0.0 ),
			M_StimDuration       ( 1.0 )
		{}

		StimulationProtocol::StimulationProtocol( Teuchos::ParameterList& parameterList )
		{
			M_timeSt             = parameterList.get ( "stimuliTime", 20.0 );
			M_tShortS1S1         = parameterList.get ( "tShortS1S1", 1000.0);
			M_stInt     	     = parameterList.get ( "stInt", 600.0 );
			M_stIntMin     	     = parameterList.get ( "stIntMin", 600.0 );
			M_stIntS1S2          = parameterList.get ( "stIntS1S2", 350.0 );
			M_stIntS1S2Min       = parameterList.get ( "stIntS1S2Min", 200.0 );
			M_stIntS2S3          = parameterList.get ( "stIntS2S3", 250.0 );
			M_stIntS3S4          = parameterList.get ( "stIntS3S4", 150.0 );
			M_nbStimMax          = parameterList.get ( "nbStimMax", 20 );
			M_repeatSt           = parameterList.get ( "repeatSt", 10 );
			M_pacingProtocol     = parameterList.get ( "pacPro", "FCL" );
			M_pacingProtocolType = parameterList.get ( "pacProType", "" );
			M_Istim              = parameterList.get ( "iStim", 0.0 );
			M_StimDuration       = parameterList.get ( "stimDur", 1.0);
		}

		StimulationProtocol::StimulationProtocol ( const StimulationProtocol& protocol )
		{
			M_timeSt             = protocol.M_timeSt;
			M_tShortS1S1         = protocol.M_tShortS1S1;
			M_stInt              = protocol.M_stInt;
			M_stIntMin           = protocol.M_stIntMin;
			M_stIntS1S2          = protocol.M_stIntS1S2;
			M_stIntS1S2Min       = protocol.M_stIntS1S2Min;
			M_stIntS2S3          = protocol.M_stIntS2S3;
			M_stIntS3S4          = protocol.M_stIntS3S4;
			M_nbStimMax          = protocol.M_nbStimMax;
			M_repeatSt           = protocol.M_repeatSt;
			M_pacingProtocol     = protocol.M_pacingProtocol;
			M_pacingProtocolType = protocol.M_pacingProtocolType;
			M_Istim              = protocol.M_Istim;
			M_StimDuration       = protocol.M_StimDuration;
		}

		// ===================================================
		//! Operator
		// ===================================================

		StimulationProtocol& StimulationProtocol::operator= ( const StimulationProtocol& protocol )
		{
			M_timeSt             = protocol.M_timeSt;
			M_tShortS1S1         = protocol.M_tShortS1S1;
			M_stInt              = protocol.M_stInt;
			M_stIntMin           = protocol.M_stIntMin;
			M_stIntS1S2          = protocol.M_stIntS1S2;
			M_stIntS1S2Min       = protocol.M_stIntS1S2Min;
			M_stIntS2S3          = protocol.M_stIntS2S3;
			M_stIntS3S4          = protocol.M_stIntS3S4;
			M_nbStimMax          = protocol.M_nbStimMax;
			M_repeatSt           = protocol.M_repeatSt;
			M_pacingProtocol     = protocol.M_pacingProtocol;
			M_pacingProtocolType = protocol.M_pacingProtocolType;
			M_Istim              = protocol.M_Istim;
			M_StimDuration       = protocol.M_StimDuration;

			return *this;
		}

		// ===================================================
		//! Methods
		// ===================================================

		// Change the number of stimulus already applied at a given time and the appropriate value of Iapp

		void StimulationProtocol::pacingProtocolChoice( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
		{
			// Choice the protocol according to the parameter passed.
			// Is this method that must be used in the main.cpp

			if ( M_pacingProtocol == "FCL" )
				fixedCycleLength ( t, dt, NbStimulus, Iapp );

			else if ( M_pacingProtocol == "FCL-ExtraSt" )
				fixedCycleLengthwExtraStim ( t, dt, NbStimulus, Iapp );

			else if ( M_pacingProtocol == "S1S2Pro" )
				standardS1S2Protocol( t, dt, NbStimulus, Iapp );

			else if ( M_pacingProtocol == "DynPro" )
				dynamicProtocol( t, dt, NbStimulus, Iapp );

			else
				fixedCycleLength ( t, dt, NbStimulus, Iapp );


		}

		void StimulationProtocol::fixedCycleLength( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
		{
			if ( NbStimulus < M_nbStimMax )
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
					Iapp = 0;
			}
			else
				Iapp = 0;
		}

		void StimulationProtocol::fixedCycleLengthwExtraStim( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
		{
			if ( NbStimulus < M_nbStimMax )
			{
				if ( t >= M_timeSt && t <= M_timeSt + M_StimDuration )
				{
					Iapp = M_Istim;

					if ( t >= M_timeSt + M_StimDuration - dt && t <= M_timeSt + M_StimDuration )
					{
						if ( NbStimulus < M_repeatSt )
						{
							M_timeSt = M_timeSt + M_stInt;
							NbStimulus++;
						}
						else
						{
							if ( M_pacingProtocolType == "S1-S2" )
							{
								M_timeSt = M_timeSt + M_stIntS1S2;
								NbStimulus = 0;
							}
							else if ( M_pacingProtocolType == "S1-S2-S3" )
							{
								if ( NbStimulus == M_repeatSt )
								{
									M_timeSt = M_timeSt + M_stIntS1S2;
									NbStimulus++;
								}
								else
								{
									M_timeSt = M_timeSt + M_stIntS2S3;
									NbStimulus = 0;
								}
							}
							else if ( M_pacingProtocolType == "S1-S2-S3-S4" )
							{
								if ( NbStimulus == M_repeatSt )
								{
									M_timeSt = M_timeSt + M_stIntS1S2;
									NbStimulus++;
								}
								else if ( NbStimulus == M_repeatSt + 1 )
								{
									M_timeSt = M_timeSt + M_stIntS2S3;
									NbStimulus++;
								}
								else
								{
									M_timeSt = M_timeSt + M_stIntS3S4;
									NbStimulus = 0;
								}
							}
							else
								M_timeSt = M_timeSt + M_stInt;
						}
					}
				}
				else
					Iapp = 0;
			}
			else
				Iapp = 0;
		}

		void StimulationProtocol::standardS1S2Protocol( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
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
							NbStimulus = 0;
						
						else
							NbStimulus++;
					}
				}
				else
					Iapp = 0;
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
								M_timeSt = M_timeSt + M_stInt;

							else if ( NbStimulus == M_repeatSt )
								M_timeSt = M_timeSt + M_stIntS1S2;

							else if ( NbStimulus == M_repeatSt + 1 )
							{
								M_timeSt = M_timeSt + M_stInt;
								NbStimulus = 0;

								if ( M_stIntS1S2 > 1000 )
									M_stIntS1S2 = M_stIntS1S2 - 1000;

								else if ( M_stIntS1S2 <= 1000 && M_stIntS1S2 > 300 )
									M_stIntS1S2 = M_stIntS1S2 - 50;

								else if (  M_stIntS1S2 <= 300 &&  M_stIntS1S2 > 200 )
									M_stIntS1S2 = M_stIntS1S2 - 10;

								else if ( M_stIntS1S2 <= 200 )
									M_stIntS1S2 = M_stIntS1S2 - 5;
							}
						}
					}
					else 
						Iapp = 0;
				}
				else
					Iapp = 0;
			}
		}

		void StimulationProtocol::dynamicProtocol( const Real& t, const Real& dt, int& NbStimulus, Real& Iapp )
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
					Iapp = 0;

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
				Iapp = 0;
		}

		void StimulationProtocol::showMe()
		{
			std::cout << "\n\n\t\tStimulation protocol Informations\n\n";

			std::cout << "\n\t\tList of parameters:\n\n";

			std::cout << "Istim: " << this->currStim() << std::endl;
			std::cout << "StimDuration: " << this->stimDuration() << std::endl;
			std::cout << "1st stimuli time: " << this->timeSt() << std::endl;
			std::cout << "Pacing protocol: " << this->pacPro() << std::endl;

			if ( M_pacingProtocol == "FCL" )
			{
				std::cout << "S1-S1 interval: " << this->stInt() << std::endl;
				std::cout << "NbStimuliMax: " << this ->nbStimMax() << std::endl;
			}
			else if ( M_pacingProtocol == "FCL-ExtraSt" )
			{
				std::cout << "Pacing protocol type: " << this->pacProTyp() << std::endl;
				std::cout << "S1-S1 interval: " << this->stInt() << std::endl;
				std::cout << "NbStimuliMax: " << this ->nbStimMax() << std::endl;
				std::cout << "Repeat S1 stimuli: " << this->repeatSt() << std::endl;

				if ( M_pacingProtocolType == "S1-S2" )
					std::cout << "S1-S2 interval: " << this->stIntS1S2() << std::endl;
				else if ( M_pacingProtocolType == "S1-S2-S3" )
				{
					std::cout << "S1-S2 interval: " << this->stIntS1S2() << std::endl;
					std::cout << "S2-S3 interval: " << this->stIntS2S3() << std::endl;
				}
				else if ( M_pacingProtocolType == "S1-S2-S3-S4" )
				{
					std::cout << "S1-S2 interval: " << this->stIntS1S2() << std::endl;
					std::cout << "S2-S3 interval: " << this->stIntS2S3() << std::endl;
					std::cout << "S3-S4 interval: " << this->stIntS3S4() << std::endl;
				}

			}
			else if ( M_pacingProtocol == "S1S2Pro" )
			{
				std::cout << "S1-S1 interval: " << this->stInt() << std::endl;
				std::cout << "NbStimuliMax for stabilisation: " << this ->nbStimMax() << std::endl;
				std::cout << "First S1-S2 interval: " << this->stIntS1S2() << std::endl;
				std::cout << "Minimum S1-S2 interval: " << this->stIntS1S2Min() << std::endl;
				std::cout << "Repeat S1 stimuli: " << this->repeatSt() << std::endl;
			}
			else if ( M_pacingProtocol == "DynPro" )
			{
				std::cout << "First S1-S1 interval: " << this->stInt() << std::endl;
				std::cout << "Minimum S1-S1 interval: " << this->stIntMin() << std::endl;
				std::cout << "First time S1-S1 interval decrease: " << this->timeShortS1S1() << std::endl;
			}
			else
			{
				std::cout << "S1-S1 interval: " << this->stInt() << std::endl;
				std::cout << "NbStimuliMax: " << this ->nbStimMax() << std::endl;
			}



			std::cout << "\n\t\t End of Stimulation protocol Informations\n\n\n";
		}

	}
