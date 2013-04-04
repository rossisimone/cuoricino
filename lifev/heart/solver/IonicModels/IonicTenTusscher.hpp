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
	  @brief Ionic model based on ten Tusscher model.
	  @date 03-2013
	  @author Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

	  @contributors
	  @mantainer Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
	  @last update 03-2013
	 */


	#ifndef TENTUSSCHER_HPP_INCLUDED
	#define TENTUSSCHER_HPP_INCLUDED

	#include <lifev/heart/solver/IonicModels/HeartIonicModel.hpp>
	#include <Teuchos_RCP.hpp>
	#include <Teuchos_ParameterList.hpp>
	#include "Teuchos_XMLParameterListHelpers.hpp"

	#include <cmath>
	#include <string>


	namespace LifeV
	{
	//! XbModel - This class implements a mean field model.

	class IonicTenTusscher : public virtual HeartIonicModel
	{

	public:
		//! @name Type definitions
		//@{
		typedef HeartIonicModel super;
		typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
		typedef RegionMesh<LinearTetra> mesh_Type;
		//@}



		//! @name Constructors & Destructor
		//@{

		//! Constructor
		IonicTenTusscher();

		/*!
		 * @param Epetra communicator
		 * @param list of parameters in an xml file
		 */
		IonicTenTusscher ( Teuchos::ParameterList& parameterList );

		/*!
		 * @param IonicTenTusscher object
		 */
		IonicTenTusscher ( const IonicTenTusscher& model );

		//! Destructor
		virtual ~IonicTenTusscher() {}

		//@}

		//! @name Overloads
		//@{

		IonicTenTusscher& operator= ( const IonicTenTusscher& model );

		//@}

		//! @name Setters and getters
		//@{

		inline const Real& gasConst() const
		{
			return M_R;
		}
		inline void setR(const Real& gasConst)
		{
			this->M_R = gasConst;
		}

		inline const Real& temp() const
		{
			return M_T;
		}
		inline void setTemp(const Real& temp)
		{
			this->M_T = temp;
		}

		inline const Real& farad() const
		{
			return M_F;
		}
		inline void setFarad(const Real& farad)
		{
			this->M_F = farad;
		}

		inline const Real& capMem() const
		{
			return M_Cm;
		}
		inline void setCapMem(const Real& capMem)
		{
			this->M_Cm = capMem;
		}

		inline const Real& svRatio() const
		{
			return M_S;
		}
		inline void setSVRatio(const Real& svRatio)
		{
			this->M_S = svRatio;
		}

		inline const Real& resCell() const
		{
			return M_rho;
		}
		inline void setResCell(const Real& resCell)
		{
			this->M_rho = resCell;
		}

		inline const Real& volCyt() const
		{
			return M_VCyt;
		}
		inline void setVJsr(const Real& volCyt)
		{
			this->M_VCyt = volCyt;
		}

		inline const Real& volSR() const
		{
			return M_VSr;
		}
		inline void setVSr(const Real& volSR)
		{
			this->M_VSr = volSR;
		}

		inline const Real& concKO() const
		{
			return M_KO;
		}
		inline void setKO(const Real& concKO)
		{
			this->M_KO = concKO;
		}

		inline const Real& concNaO() const
		{
			return M_NaO;
		}
		inline void setNaO(const Real& concNaO)
		{
			this->M_NaO = concNaO;
		}

		inline const Real& concCaO() const
		{
			return M_CaO;
		}
		inline void setCaO(const Real& concCaO)
		{
			this->M_CaO = concCaO;
		}

		inline const Real& maxCondNa() const
		{
			return M_GNa;
		}
		inline void setGNa(const Real& maxCondNa)
		{
			this->M_GNa = maxCondNa;
		}

		inline const Real& maxCondK1() const
		{
			return M_GK1;
		}
		inline void setGK1(const Real& maxCondK1)
		{
			this->M_GK1 = maxCondK1;
		}

		inline const Real& maxCondToEpiM() const
		{
			return M_GToEpiM;
		}
		inline void setGToEpiM(const Real& maxCondToEpiM)
		{
			this->M_GToEpiM = maxCondToEpiM;
		}

		inline const Real& maxCondToEndo() const
		{
			return M_GToEndo;
		}
		inline void setGToEndo(const Real& maxCondToEndo)
		{
			this->M_GToEndo = maxCondToEndo;
		}

		inline const Real& maxCondKr() const
		{
			return M_GKr;
		}
		inline void setGKr(const Real& maxCondKr)
		{
			this->M_GKr = maxCondKr;
		}

		inline const Real& maxCondKs() const
		{
			return M_GKs;
		}
		inline void setGKs(const Real& maxCondKs)
		{
			this->M_GKs = maxCondKs;
		}

		inline const Real& maxCondKsM() const
		{
			return M_GKsM;
		}
		inline void setGKsM(const Real& maxCondKsM)
		{
			this->M_GKsM = maxCondKsM;
		}

		inline const Real& relPermKNa() const
		{
			return M_pKNa;
		}
		inline void setRelPermKNa(const Real& relPermKNa)
		{
			this->M_pKNa = relPermKNa;
		}

		inline const Real& maxCondCaL() const
		{
			return M_GCaL;
		}
		inline void setGCaL(const Real& maxCondCaL)
		{
			this->M_GCaL = maxCondCaL;
		}

		inline const Real& maxCourNaCa() const
		{
			return M_kNaCa;
		}
		inline void setMaxCourNaCa(const Real& maxCourNaCa)
		{
			this->M_kNaCa = maxCourNaCa;
		}

		inline const Real& gamma() const
		{
			return M_gamma;
		}
		inline void setGamma(const Real& gamma)
		{
			this->M_gamma = gamma;
		}

		inline const Real& constmCa() const
		{
			return M_KmCa;
		}
		inline void setKmCa(const Real& constmCa)
		{
			this->M_KmCa = constmCa;
		}

		inline const Real& constmNai() const
		{
			return M_KmNai;
		}
		inline void setKmNai(const Real& constmNai)
		{
			this->M_KmNai = constmNai;
		}

		inline const Real& kSat() const
		{
			return M_kSat;
		}
		inline void setKSat(const Real& kSat)
		{
			this->M_kSat = kSat;
		}

		inline const Real& alpha() const
		{
			return M_alpha;
		}
		inline void setAlpha(const Real& alpha)
		{
			this->M_alpha = alpha;
		}

		inline const Real& permNaK() const
		{
			return M_PNaK;
		}
		inline void setPNaK(const Real& permNaK)
		{
			this->M_PNaK = permNaK;
		}

		inline const Real& constmK() const
		{
			return M_KmK;
		}
		inline void setKmK(const Real& constmK)
		{
			this->M_KmK = constmK;
		}

		inline const Real& constmNa() const
		{
			return M_KmNa;
		}
		inline void setKmNa(const Real& constmNa)
		{
			this->M_KmNa = constmNa;
		}

		inline const Real& maxCondKp() const
		{
			return M_GKp;
		}
		inline void setGKp(const Real& maxCondKp)
		{
			this->M_GKp = maxCondKp;
		}

		inline const Real& maxCondCap() const
		{
			return M_GCap;
		}
		inline void setGCap(const Real& maxCondCap)
		{
			this->M_GCap = maxCondCap;
		}

		inline const Real& constpCa() const
		{
			return M_KpCa;
		}
		inline void setKmPCa(const Real& constpCa)
		{
			this->M_KpCa = constpCa;
		}

		inline const Real& maxCondNab() const
		{
			return M_GNab;
		}
		inline void setGNab(const Real& maxCondNab)
		{
			this->M_GNab = maxCondNab;
		}

		inline const Real& maxCondCab() const
		{
			return M_GCab;
		}
		inline void setGCab(const Real& maxCondCab)
		{
			this->M_GCab = maxCondCab;
		}

		inline const Real& maxCourUp() const
		{
			return M_VMaxUp;
		}
		inline void setVmaxUp(const Real& maxCourUp)
		{
			this->M_VMaxUp = maxCourUp;
		}

		inline const Real& constUp() const
		{
			return M_Kup;
		}
		inline void setKup(const Real& constUp)
		{
			this->M_Kup = constUp;
		}

		inline const Real& aRel() const
		{
			return M_aRel;
		}
		inline void setARel(const Real& aRel)
		{
			this->M_aRel = aRel;
		}

		inline const Real& bRel() const
		{
			return M_bRel;
		}
		inline void setBRel(const Real& bRel)
		{
			this->M_bRel = bRel;
		}

		inline const Real& cRel() const
		{
			return M_cRel;
		}
		inline void setCRel(const Real& cRel)
		{
			this->M_cRel = cRel;
		}

		inline const Real& maxCourLeak() const
		{
			return M_VLeak;
		}
		inline void setVLeak(const Real& maxCourLeak)
		{
			this->M_VLeak = maxCourLeak;
		}

		inline const Real& buffCyt() const
		{
			return M_Buffc;
		}
		inline void setBuffc(const Real& buffCyt)
		{
			this->M_Buffc = buffCyt;
		}

		inline const Real& constBuffc() const
		{
			return M_KBuffc;
		}
		inline void setKBuffc(const Real& constBuffc)
		{
			this->M_KBuffc = constBuffc;
		}

		inline const Real& buffSR() const
		{
			return M_BuffSR;
		}
		inline void setBuffSR(const Real& buffSR)
		{
			this->M_BuffSR = buffSR;
		}

		inline const Real& constBuffSR() const
		{
			return M_KBuffSR;
		}
		inline void setKBuffSR(const Real& constBuffSR)
		{
			this->M_KBuffSR = constBuffSR;
		}

		inline const std::string& typeCell() const
		{
			return M_typeCell;
		}
		inline void setKBuffSR(const std::string& typeCell)
		{
			this->M_typeCell = typeCell;
		}

		//@}

		//! @name Methods
		//@{

		//Compute the rhs on a single node or for the 0D case
		void computeRhs ( const std::vector<Real>&  v, std::vector<Real>& rhs );

		void computeRhs ( const std::vector<Real>&  v, const Real& Istim, std::vector<Real>& rhs );


		//Compute the rhs with state variable interpolation
		Real computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Istim );

		std::vector<Real> computeLocalGatingRhs ( const std::vector<Real>& v );

		std::vector<Real> computeLocalConcRhs ( const std::vector<Real>& v );

		std::vector<Real> computeLocalConcRhs ( const std::vector<Real>& v, const Real& Istim );

		std::vector<Real> computeLocalSubSysCaRhs( const std::vector<Real>& v );

		std::vector<Real> fastINa( const std::vector<Real>& v );

		std::vector<Real> transientIto( const std::vector<Real>& v );

		std::vector<Real> slowIKs( const std::vector<Real>& v );

		std::vector<Real> rapDelIKr( const std::vector<Real>& v );

		Real inwardIK1( const std::vector<Real>& v );

		Real exINaCa( const std::vector<Real>& v );

		Real pumpINaK( const std::vector<Real>& v );

		Real pumpIpCa( const std::vector<Real>& v );

		Real pumpIpK( const std::vector<Real>& v );

		Real backICab( const std::vector<Real>& v );

		Real backINab( const std::vector<Real>& v );

		//! Display information about the model
		void showMe();

		//@}

	private:

        Real M_R;
        Real M_T;
		Real M_F;
		Real M_Cm;
		Real M_S;
		Real M_rho;
		Real M_VCyt;
		Real M_VSr;
		Real M_KO;
		Real M_NaO;
		Real M_CaO;
		Real M_GNa;
		Real M_GK1;
		Real M_GToEpiM;
		Real M_GToEndo;
		Real M_GKr;
		Real M_GKs;
		Real M_GKsM;
		Real M_pKNa;
		Real M_GCaL;
		Real M_kNaCa;
		Real M_gamma;
		Real M_KmCa;
		Real M_KmNai;
		Real M_kSat;
		Real M_alpha;
		Real M_PNaK;
		Real M_KmK;
		Real M_KmNa;
		Real M_GKp;
		Real M_GCap;
		Real M_KpCa;
		Real M_GNab;
		Real M_GCab;
		Real M_VMaxUp;
		Real M_Kup;
		Real M_aRel;
		Real M_bRel;
		Real M_cRel;
		Real M_VLeak;
		Real M_Buffc;
		Real M_KBuffc;
		Real M_BuffSR;
		Real M_KBuffSR;
		std::string M_typeCell;

		//@}

	}; // class IonicTenTusscher

	// ===================================================
	//! Constructors
	// ===================================================
	IonicTenTusscher::IonicTenTusscher()    :
		super      ( 17 ),
		M_R        ( 8314.472 ),
		M_T        ( 310.0 ),
		M_F        ( 96485.3415 ),
		M_Cm       ( 0.185 ),
		M_S        ( 0.2 ),
		M_rho      ( 162.0 ),
		M_VCyt     ( 0.016404 ),
		M_VSr      ( 0.001094 ),
		M_KO       ( 5.4 ),
		M_NaO      ( 140.0 ),
		M_CaO      ( 2.0 ),
		M_GNa      ( 14.838 ),
		M_GK1      ( 5.405 ),
		M_GToEpiM  ( 0.294 ),
		M_GToEndo  ( 0.073 ),
		M_GKr      ( 0.096 ),
		M_GKs      ( 0.245 ),
		M_GKsM     ( 0.062 ),
		M_pKNa     ( 0.03 ),
		M_GCaL     ( 0.000175 ),
		M_kNaCa    ( 1000 ),
		M_gamma    ( 0.35 ),
		M_KmCa     ( 1.38 ),
		M_KmNai    ( 87.5 ),
		M_kSat     ( 0.1 ),
		M_alpha    ( 2.5 ),
		M_PNaK     ( 1.362 ),
		M_KmK      ( 1.0 ),
		M_KmNa     ( 40.0 ),
		M_GKp      ( 0.0146 ),
		M_GCap     ( 0.825 ),
		M_KpCa     ( 0.0005 ),
		M_GNab     ( 0.00029 ),
		M_GCab     ( 0.000592 ),
		M_VMaxUp   ( 0.000425 ),
		M_Kup      ( 0.00025 ),
		M_aRel     ( 0.016464 ),
		M_bRel     ( 0.25 ),
		M_cRel     ( 0.008232 ),
		M_VLeak    ( 0.00008 ),
		M_Buffc    ( 0.15 ),
		M_KBuffc   ( 0.001 ),
		M_BuffSR   ( 10.0 ),
		M_KBuffSR  ( 0.3 ),
		M_typeCell ("epicardial")
	{}

	IonicTenTusscher::IonicTenTusscher ( Teuchos::ParameterList& parameterList ) :
		super       ( 17 )
	{
		M_R        = parameterList.get ( "gasConst", 8314.472 );
		M_T        = parameterList.get ( "temp", 310.0 );
		M_F        = parameterList.get ( "farad", 96485.3415 );
		M_Cm       = parameterList.get ( "capMem", 0.185 );
		M_S        = parameterList.get ( "svRatio", 0.2 );
		M_rho      = parameterList.get ( "resCell", 162.0 );
		M_VCyt     = parameterList.get ( "volCyt", 0.016404 );
		M_VSr      = parameterList.get ( "volSR", 0.001094  );
		M_KO       = parameterList.get ( "concKO", 5.4  );
		M_NaO      = parameterList.get ( "concNaO", 140.0 );
		M_CaO      = parameterList.get ( "concCaO", 2.0 );
		M_GNa      = parameterList.get ( "maxCondNa", 14.838 );
		M_GK1      = parameterList.get ( "maxCondK1", 5.405 );
		M_GToEpiM  = parameterList.get ( "maxCondToEpiM", 0.294 );
		M_GToEndo  = parameterList.get ( "maxCondToEndo", 0.073 );
		M_GKr      = parameterList.get ( "maxCondKr", 0.096 );
		M_GKs      = parameterList.get ( "maxCondKs", 0.245 );
		M_GKsM      = parameterList.get ( "maxCondKsM", 0.062 );
		M_pKNa     = parameterList.get ( "relPermKNa", 0.03 );
		M_GCaL     = parameterList.get ( "maxCondCaL", 0.000175 );
		M_kNaCa    = parameterList.get ( "maxCourNaCa", 1000.0 );
		M_gamma    = parameterList.get ( "gamma", 0.35 );
		M_KmCa     = parameterList.get ( "constmCa", 1.38 );
		M_KmNai    = parameterList.get ( "constmNai", 87.5 );
		M_kSat     = parameterList.get ( "kSat", 0.1 );
		M_alpha    = parameterList.get ( "alpha", 2.5 );
		M_PNaK     = parameterList.get ( "permNaK", 1.362 );
		M_KmK      = parameterList.get ( "constmK", 1.0 );
		M_KmNa     = parameterList.get ( "constmNa", 40.0 );
		M_GKp      = parameterList.get ( "maxCondKp", 0.0146 );
		M_GCap     = parameterList.get ( "maxCondCap", 0.825 );
		M_KpCa     = parameterList.get ( "constpCa", 0.0005 );
		M_GNab     = parameterList.get ( "maxCondNab", 0.00029 );
		M_GCab     = parameterList.get ( "maxCondCab", 0.000592 );
		M_VMaxUp   = parameterList.get ( "maxCourUp", 0.000425 );
		M_Kup      = parameterList.get ( "constUp", 0.00025 );
		M_aRel     = parameterList.get ( "aRel", 0.016464 );
		M_bRel     = parameterList.get ( "bRel", 0.25 );
		M_cRel     = parameterList.get ( "cRel", 0.008232 );
		M_VLeak    = parameterList.get ( "maxCourLeak", 0.00008 );
		M_Buffc    = parameterList.get ( "buffCyt", 0.15 );
		M_KBuffc   = parameterList.get ( "constBuffc", 0.001 );
		M_BuffSR   = parameterList.get ( "buffSR", 10.0 );
		M_KBuffSR  = parameterList.get ( "constBuffSR", 0.3 );
		M_typeCell = parameterList.get ( "typeCell", "epicardial");
	}
	
	IonicTenTusscher::IonicTenTusscher ( const IonicTenTusscher& model )
	{
		M_R        = model.M_R;
		M_T        = model.M_T;
		M_F        = model.M_F;
		M_Cm       = model.M_Cm;
		M_S        = model.M_S;
		M_rho      = model.M_rho;
		M_VCyt     = model.M_VCyt;
		M_VSr      = model.M_VSr;
		M_KO       = model.M_KO;
		M_NaO      = model.M_NaO;
		M_CaO      = model.M_CaO;
		M_GNa      = model.M_GNa;
		M_GK1      = model.M_GK1;
		M_GToEpiM  = model.M_GToEpiM;
		M_GToEndo  = model.M_GToEndo;
		M_GKr      = model.M_GKr;
		M_GKs      = model.M_GKs;
		M_GKsM     = model.M_GKsM;
		M_pKNa     = model.M_pKNa;
		M_GCaL     = model.M_GCaL;
		M_kNaCa    = model.M_kNaCa;
		M_gamma    = model.M_gamma;
		M_KmCa     = model.M_KmCa;
		M_KmNai    = model.M_KmNai;
		M_kSat     = model.M_kSat;
		M_alpha    = model.M_alpha;
		M_PNaK     = model.M_PNaK;
		M_KmK      = model.M_KmK;
		M_KmNa     = model.M_KmNa;
		M_GKp      = model.M_GKp;
		M_GCap     = model.M_GCap;
		M_KpCa     = model.M_KpCa;
		M_GNab     = model.M_GNab;
		M_GCab     = model.M_GCab;
		M_VMaxUp   = model.M_VMaxUp;
		M_Kup      = model.M_Kup;
		M_aRel     = model.M_aRel;
		M_bRel     = model.M_bRel;
		M_cRel     = model.M_cRel;
		M_VLeak    = model.M_VLeak;
		M_Buffc    = model.M_Buffc;
		M_KBuffc   = model.M_KBuffc;
		M_BuffSR   = model.M_BuffSR;
		M_KBuffSR  = model.M_KBuffSR;
		M_typeCell = model.M_typeCell;
	
		M_numberOfEquations = model.M_numberOfEquations;
	}

	// ===================================================
	//! Operator
	// ===================================================
	IonicTenTusscher& IonicTenTusscher::operator= ( const IonicTenTusscher& model )
	{
		M_R        = model.M_R;
		M_T        = model.M_T;
		M_F        = model.M_F;
		M_Cm       = model.M_Cm;
		M_S        = model.M_S;
		M_rho      = model.M_rho;
		M_VCyt     = model.M_VCyt;
		M_VSr      = model.M_VSr;
		M_KO       = model.M_KO;
		M_NaO      = model.M_NaO;
		M_CaO      = model.M_CaO;
		M_GNa      = model.M_GNa;
		M_GK1      = model.M_GK1;
		M_GToEpiM  = model.M_GToEpiM;
		M_GToEndo  = model.M_GToEndo;
		M_GKr      = model.M_GKr;
		M_GKs      = model.M_GKs;
		M_GKsM     = model.M_GKsM;
		M_pKNa     = model.M_pKNa;
		M_GCaL     = model.M_GCaL;
		M_kNaCa    = model.M_kNaCa;
		M_gamma    = model.M_gamma;
		M_KmCa     = model.M_KmCa;
		M_KmNai    = model.M_KmNai;
		M_kSat     = model.M_kSat;
		M_alpha    = model.M_alpha;
		M_PNaK     = model.M_PNaK;
		M_KmK      = model.M_KmK;
		M_KmNa     = model.M_KmNa;
		M_GKp      = model.M_GKp;
		M_GCap     = model.M_GCap;
		M_KpCa     = model.M_KpCa;
		M_GNab     = model.M_GNab;
		M_GCab     = model.M_GCab;
		M_VMaxUp   = model.M_VMaxUp;
		M_Kup      = model.M_Kup;
		M_aRel     = model.M_aRel;
		M_bRel     = model.M_bRel;
		M_cRel     = model.M_cRel;
		M_VLeak    = model.M_VLeak;
		M_Buffc    = model.M_Buffc;
		M_KBuffc   = model.M_KBuffc;
		M_BuffSR   = model.M_BuffSR;
		M_KBuffSR  = model.M_KBuffSR;
		M_typeCell = model.M_typeCell;

		M_numberOfEquations = model.M_numberOfEquations;

		return *this;
	}


	// ===================================================
	//! Methods
	// ===================================================
	//Only gating variables
	void IonicTenTusscher::computeRhs ( const std::vector<Real>&  v,
											 std::vector<Real>& rhs )
	{
		std::vector<Real> gatingRhs     ( computeLocalGatingRhs(v) );
		std::vector<Real> concRhs       ( computeLocalConcRhs(v) );
		
		
		std::copy( gatingRhs.begin(), gatingRhs.end(), rhs.begin() );

		std::copy( concRhs.begin(), concRhs.end(), rhs.begin() + 12 );

	}

	//Potential and gating variables
	void IonicTenTusscher::computeRhs (const   std::vector<Real>&  v,
											 const   Real& Istim,
											 std::vector<Real>& rhs )
	{
		std::vector<Real> gatingRhs     ( computeLocalGatingRhs(v) );
		std::vector<Real> concRhs       ( computeLocalConcRhs(v, Istim) );
		
		
		rhs[0] = computeLocalPotentialRhs(v, Istim);

		std::copy( gatingRhs.begin(), gatingRhs.end(), rhs.begin() + 1 );

		std::copy( concRhs.begin(), concRhs.end(), rhs.begin() + 13 );
		
	}

	Real IonicTenTusscher::computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Istim )
	{
		std::vector<Real> courSubSysCa ( computeLocalSubSysCaRhs(v) );
		std::vector<Real> courINa	   ( fastINa(v) );
		std::vector<Real> courIto	   ( transientIto(v) );
		std::vector<Real> courIKs      ( slowIKs(v) );
		std::vector<Real> courIKr      ( rapDelIKr(v) );
	
		Real iIon = courINa[0] + inwardIK1(v) + courIto[0] + courIKr[0] 
					+ courIKs[0] + courSubSysCa[0] + exINaCa(v) + pumpINaK(v) 
					+ pumpIpCa(v) + backICab(v) + backINab(v);
		
		return  - ( iIon + Istim ) ;
	}

	std::vector<Real> IonicTenTusscher::computeLocalGatingRhs ( const std::vector<Real>& v )
	{
		std::vector<Real> gatingINa (fastINa(v));
		std::vector<Real> gatingIKr (rapDelIKr(v));
		std::vector<Real> gatingIKs (slowIKs(v));
		std::vector<Real> gatingICa (computeLocalSubSysCaRhs(v));
		std::vector<Real> gatingIto (transientIto(v));
	
		std::vector<Real> gatingRhs (12);

		std::copy( gatingINa.begin() + 1, gatingINa.end(), gatingRhs.begin() );
		std::copy( gatingIKr.begin() + 1, gatingIKr.end(), gatingRhs.begin() + 3 );
		gatingRhs[5] = gatingIKs[1];
		std::copy( gatingICa.begin() + 6, gatingICa.end() - 1, gatingRhs.begin() + 6 );
		std::copy( gatingIto.begin() + 1, gatingIto.end(), gatingRhs.begin() + 9 );
		gatingRhs[12] = gatingICa[9];

		return gatingRhs;
	}

	std::vector<Real> IonicTenTusscher::computeLocalConcRhs ( const std::vector<Real>& v )
	{

		std::vector<Real> courSubSysCa ( computeLocalSubSysCaRhs(v) );
		std::vector<Real> courINa	   ( fastINa(v) );
		std::vector<Real> courIto	   ( transientIto(v) );
		std::vector<Real> courIKs      ( slowIKs(v) );
		std::vector<Real> courIKr      ( rapDelIKr(v) );
		
		std::vector<Real> concRhs(4);

		concRhs[0] = courSubSysCa[4] *( courSubSysCa[1] - courSubSysCa[2] + courSubSysCa[3] -
						M_Cm * ( courSubSysCa[0] + backICab(v) + pumpIpCa(v) - 2 * exINaCa(v) ) / ( 2 * M_VCyt * M_F ) );
		concRhs[1] = courSubSysCa[5] * ( M_VCyt / M_VSr ) * ( - courSubSysCa[1] + courSubSysCa[2] - courSubSysCa[3] );
		concRhs[2] = - ( courINa[0] + backINab(v) + 3 * exINaCa(v) + 3 * pumpINaK(v) ) * M_Cm / ( M_VCyt * M_F );
		concRhs[3] = - ( inwardIK1(v) + courIto[0] + courIKr[0] + courIKs[0] - 2 * pumpINaK(v) + pumpIpK(v) ) * M_Cm / ( M_VCyt * M_F );
		
		return concRhs;
	}
	
	std::vector<Real> IonicTenTusscher::computeLocalConcRhs ( const std::vector<Real>& v, const Real& Istim )
	{

		std::vector<Real> courSubSysCa ( computeLocalSubSysCaRhs(v) );
		std::vector<Real> courINa	   ( fastINa(v) );
		std::vector<Real> courIto	   ( transientIto(v) );
		std::vector<Real> courIKs      ( slowIKs(v) );
		std::vector<Real> courIKr      ( rapDelIKr(v) );
		
		std::vector<Real> concRhs(4);

		concRhs[0] = courSubSysCa[4] *( courSubSysCa[1] - courSubSysCa[2] + courSubSysCa[3] -
						M_Cm * ( courSubSysCa[0] + backICab(v) + pumpIpCa(v) - 2 * exINaCa(v) ) / ( 2 * M_VCyt * M_F ) );
		concRhs[1] = courSubSysCa[5] * ( M_VCyt / M_VSr ) * ( - courSubSysCa[1] + courSubSysCa[2] - courSubSysCa[3] );
		concRhs[2] = - ( courINa[0] + backINab(v) + 3 * exINaCa(v) + 3 * pumpINaK(v) ) * M_Cm / ( M_VCyt * M_F );
		concRhs[3] = - ( inwardIK1(v) + courIto[0] + courIKr[0] + courIKs[0] - 2 * pumpINaK(v) + pumpIpK(v) + Istim ) * M_Cm / ( M_VCyt * M_F );

		return concRhs;
	}
	
	//! Ca2+ Subsystem

	std::vector<Real> IonicTenTusscher::computeLocalSubSysCaRhs( const std::vector<Real>& v )
	{
		std::vector<Real> subSysCaRHS(10);

		Real V 	   ( v[0] );
		Real d     ( v[7] );
		Real f     ( v[8] );
		Real fCa   ( v[9] );
		Real g     ( v[12] );
		Real cCa   ( v[13] );
		Real cCaSR ( v[14] );

		// Internal Parameters

		Real jLeak = M_VLeak * ( cCaSR - cCa );
		Real jUp   = M_VMaxUp * ( pow(cCa, 2) / ( pow(M_Kup, 2) + pow(cCa, 2) ) );
		Real jRel  = ( M_aRel * pow(cCaSR, 2) / ( pow(M_bRel, 2) + pow(cCaSR, 2) ) + M_cRel ) * d * g ;
		
		Real g_inf (0);

		if ( cCa < 0.00035 )
			g_inf = 1.0 / ( 1.0 + pow(cCa / 0.00035, 6) );
		else
			g_inf = 1.0 / ( 1.0 + pow(cCa / 0.00035, 16) );
		
		Real tau_g (2.0);
		
		int k_g (0);
		if ( ( g_inf > g ) && ( V > -37 ) )
			k_g = 0;
		else 
			k_g = 1;
		
		//Real cCai_bufc   = 1.0 / ( 1.0 + M_Buffc * M_KBuffc / pow(cCa + M_KBuffc, 2) );
		//Real cCaSR_bufsr = 1.0 / ( 1.0 + M_BuffSR * M_KBuffSR / pow(cCaSR + M_KBuffSR, 2) );
		Real cCai_bufc   = cCa * M_Buffc / ( cCa + M_KBuffc );
		Real cCaSR_bufsr = cCaSR * M_BuffSR / ( cCaSR + M_KBuffSR );

		Real d_inf   = 1.0 / ( 1.0 + exp( ( -5.0 - V ) / 7.5 ) );
		Real alpha_d = 1.4 / ( 1.0 + exp( ( -35.0 - V ) / 13.0 ) ) + 0.25;
		Real beta_d  = 1.4 / ( 1.0 + exp( ( 5.0 + V ) / 5.0 ) );
		Real gamma_d = 1.0 / ( 1.0 + exp( ( 50.0 - V ) / 20.0 ) );	
		Real tau_d   = alpha_d * beta_d + gamma_d;
		
		Real f_inf = 1.0 / ( 1.0 + exp( ( 20.0 + V ) / 7.0 ) );
		Real tau_f = 1125.0 * exp( -pow(27.0 + V, 2 ) / 300.0 ) + 80.0 + 165.0 / ( 1.0 + exp( ( 25.0 - V ) / 10.0 ) );
		
		Real alpha_fCa = 1.0 / ( 1.0 + pow(cCa / 0.000325, 8) );
		Real beta_fCa  = 0.1 / ( 1.0 + exp( ( cCa - 0.0005 ) / 0.0001 ) );
		Real gamma_fCa = 0.2 / ( 1.0 + exp( ( cCa - 0.00075 ) / 0.0008 ) );
		Real fCa_inf   = ( alpha_fCa + beta_fCa + gamma_fCa + 0.23 ) / 1.46;
		Real tau_fCa   ( 2 );

		int k_f (0);
		if ( ( f_inf > f ) && ( V > -37 ) )
			k_f = 0;
		else 
			k_f = 1;
		
		// RHS of the Ca2+ subsystem

		subSysCaRHS[0] = M_GCaL * d * f * fCa * 4.0 * ( ( V * pow(M_F,2) ) / ( M_R * M_T ) ) * (cCa * exp( 2 * V * M_F / ( M_R * M_T ) ) - 0.341 * M_CaO ) / ( exp( 2 * V * M_F / ( M_R * M_T ) ) - 1.0 );
		subSysCaRHS[1] = jLeak;
		subSysCaRHS[2] = jUp;
		subSysCaRHS[3] = jRel;
		subSysCaRHS[4] = cCai_bufc;
		subSysCaRHS[5] = cCaSR_bufsr;
		subSysCaRHS[6] = ( d_inf - d ) / tau_d;
		subSysCaRHS[7] = ( f_inf - f ) / tau_f;
		subSysCaRHS[8] = k_f * ( fCa_inf - fCa ) / tau_fCa;
		subSysCaRHS[9] = k_g * ( g_inf - g ) / tau_g;

		return subSysCaRHS;
	}

	//! Ionic Currents

	// Fast Na+ Current INa
	std::vector<Real> IonicTenTusscher::fastINa( const std::vector<Real>& v )
	{
		std::vector<Real> fastNa(4);

		Real V   ( v[0] );
		Real m   ( v[1] );
		Real j   ( v[2] );
		Real h   ( v[3] );
		Real cNa ( v[15] );

		Real potNa = ( M_R * M_T / M_F ) * log( M_NaO / cNa );

		fastNa[0] = M_GNa * pow(m, 3) * h * j * ( V - potNa );
		
		// Fast Na+ current m gate
		Real m_inf   = 1.0 / pow( 1.0 + exp( ( -58.6 - V ) / 9.03), 2);
		Real alpha_m = 1.0 / ( 1.0 + exp( ( -60.0 - V ) / 5.0 ) );
		Real beta_m  = 0.1 / ( 1.0 + exp( ( V + 35.0 ) / 5.0 ) ) +  0.1 / ( 1.0 + exp( ( V - 50.0 ) / 200.0 ) );
		Real tau_m   = alpha_m * beta_m;
		
		fastNa[1] = ( m_inf - m ) / tau_m;

		// Fast Na+ current h and j gate
		Real h_inf   = 1.0 / pow( 1 + exp( ( V + 71.55 ) / 7.43), 2);
		Real j_inf   = 1.0 / pow( 1 + exp( ( V + 71.55 ) / 7.43), 2);
		

		Real alpha_h (0);
		Real alpha_j (0);
		Real beta_h  (0);
		Real beta_j  (0);

		if (V >= -40)
		{
			alpha_h = 0.0;
			alpha_j = 0.0;
			beta_h  = 0.77 / ( 0.13 * ( 1.0 + exp( -( V + 10.66 ) / 11.1 ) ) );
			beta_j  = 0.6 * exp( 0.057 * V) / ( 1.0 + exp( -0.1 * ( V + 32.0 ) ) );
		}
		else
		{
			alpha_h = 0.057 * exp( -( 80 + V ) / 6.8 );
			alpha_j = ( -25428 * exp( 0.2444 * V) - 6.948e-6 * exp( -0.04391 * V ) ) * ( V + 37.78 ) / ( 1.0 + exp ( 0.311 * ( V + 79.23 ) ) );
			beta_h  = 2.7 * exp( 0.079 * V ) + 3.1e5 * exp( 0.3485 * V ) ;
			beta_j  = 0.02424 * exp( -0.01052 * V ) / ( 1.0 + exp( -0.1378 * ( V + 40.14 ) ) );
		}
		
		Real tau_h 	= 1.0 / ( alpha_h + beta_h );
		Real tau_j  = 1.0 / ( alpha_j + beta_j );
	
		fastNa[2] = ( j_inf - j ) / tau_j;
		fastNa[3] = ( h_inf - h ) / tau_h;

		return fastNa;
	}

	// Transient outward current Ito
	std::vector<Real> IonicTenTusscher::transientIto( const std::vector<Real>& v )
	{
		std::vector<Real> transIto(3);

		Real V   ( v[0] );
		Real r   ( v[10] );
		Real s   ( v[11] );
		Real cKi ( v[16] );
		
		Real potK     = ( M_R * M_T / M_F ) * log( M_KO / cKi );
		
		if ( ( M_typeCell == "epicardial" ) || ( M_typeCell == "M cell" ) )
			transIto[0] = M_GToEpiM * r * s * ( V - potK );
		else if ( M_typeCell == "endocardial" )
			transIto[0] = M_GToEndo * r * s * ( V - potK );
		else
			transIto[0] = M_GToEpiM * r * s * ( V - potK );


		Real r_inf = 1.0 / ( 1.0 + exp( ( 20.0 - V ) / 6.0 ) );
		Real tau_r = 9.5 * exp( - pow( V + 40, 2) / 1800.0 ) + 0.8;
		

		Real s_inf (0);
		Real tau_s (0);

		if ( ( M_typeCell == "epicardial" ) || ( M_typeCell == "M cell" ) )
		{
			s_inf = 1.0 / ( 1.0 + exp( ( 20.0 + V ) / 5.0 ) );
			tau_s = 85.0 * exp( - pow( V + 45.0, 2) / 320.0 ) + 5.0 / ( 1.0 + exp( ( -20.0 + V ) / 5.0 ) ) + 3.0;
		}
		else if ( M_typeCell == "endocardial" )
		{
			s_inf = 1.0 / ( 1.0 + exp( ( 28.0 + V ) / 5.0 ) );
			tau_s = 1000.0 * exp( - pow( V + 67.0, 2) / 1000.0 ) + 8.0;
		}
		else
		{
			s_inf = 1.0 / ( 1.0 + exp( ( 20.0 + V ) / 5.0 ) );
			tau_s = 85.0 * exp( - pow( V + 45.0, 2) / 320.0 ) + 5.0 / ( 1.0 + exp( ( -20.0 + V ) / 5.0 ) ) + 3.0;
		}


		transIto[1] = ( r_inf -r ) / tau_r;
		transIto[2] = ( s_inf -r ) / tau_s;

		return transIto;
	}

	// Slow Delayed Rectifier Current IKs
	std::vector<Real> IonicTenTusscher::slowIKs( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real Xs  ( v[6] );
		Real cNa ( v[15] );
		Real cKi ( v[16] );
		
		std::vector<Real> slowIKs(2);
		
		Real potKs    = ( M_R * M_T / M_F ) * log( ( M_KO + M_pKNa * M_NaO ) / ( cKi + M_pKNa * cNa ) );
		Real alpha_xs = 1100.0 / sqrt( 1.0 + exp( ( -10.0 - V ) / 6.0 ) );
		Real beta_xs  = 1.0 / ( 1.0 + exp( ( V - 60.0 ) / 20.0 ) );
		

		if ( ( M_typeCell == "epicardial" ) || ( M_typeCell == "endocardial" ) )
			slowIKs[0] = M_GKs * pow(Xs, 2) * ( V - potKs );
		else if ( M_typeCell == "M cell" )
			slowIKs[0] = M_GKsM * pow(Xs, 2) * ( V - potKs );
		else
			slowIKs[0] = M_GKs * pow(Xs, 2) * ( V - potKs );

		Real Xs_inf = 1.0 / ( 1.0 + exp( ( -5.0 - V ) / 14.0 ) );
		Real tau_xs = alpha_xs * beta_xs;
		
		slowIKs[1] = ( Xs_inf - Xs ) / tau_xs;
		
		return  slowIKs;
	}

	// Rapid Delayed Rectifier Current
	std::vector<Real> IonicTenTusscher::rapDelIKr( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real Xr1 ( v[4] );
		Real Xr2 ( v[5] );
		Real cKi ( v[16] );

		std::vector<Real> rapidIKr(3);
		
		Real potK    = ( M_R * M_T / M_F ) * log( M_KO / cKi );
		
		Real alpha_xr1 = 450.0 / ( 1.0 + exp( ( -45.0 - V ) / 10.0 ) );
		Real beta_xr1  = 6.0 / ( 1.0 + exp( ( 30.0 + V ) / 11.5 ) );
		Real alpha_xr2 = 3.0 / ( 1.0 + exp( ( -60.0 - V ) / 20.0 ) );
		Real beta_xr2  = 1.12 / ( 1.0 + exp( ( -60.0 + V ) / 20.0 ) );
		
		rapidIKr[0] = M_GKr * sqrt( M_KO / 5.4 ) * Xr1 * Xr2 * ( V - potK );
		
		Real Xr1_inf = 1.0 / ( 1.0 + exp( ( -26.0 - V ) / 7.0 ) );
		Real tau_xr1 = alpha_xr1 * beta_xr1;
		Real Xr2_inf = 1.0 / ( 1.0 + exp( ( 88.0 + V ) / 24.0 ) );
		Real tau_xr2 = alpha_xr2 * beta_xr2;
		
		rapidIKr[1] = ( Xr1_inf - Xr1 ) / tau_xr1;
		rapidIKr[2] = ( Xr2_inf - Xr2 ) / tau_xr2;

		return rapidIKr;
	}

	// Inward Rectifier K+ Current
	Real IonicTenTusscher::inwardIK1( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cKi ( v[16] );
		
		Real potK    = ( M_R * M_T / M_F ) * log( M_KO / cKi );
		
		Real alpha_K1 = 0.1 / ( 1.0 + exp( 0.06 * ( V - potK - 200.0 ) ) );
		Real beta_K1  = ( 3.0 * exp( 0.0002 * ( V - potK + 100.0 ) ) + exp( 0.1 * ( V - potK - 10.0 ) ) )
						/ ( 1.0 + exp( -0.5 * ( V - potK ) ) );
		Real XK1_inf   = alpha_K1 / ( alpha_K1 + beta_K1 );
		
		return  M_GK1 * XK1_inf * sqrt( M_KO / 5.4 ) * ( V - potK );
	}
	
	
	// Na+/Ca2+ exchanger current INaCa
	Real IonicTenTusscher::exINaCa( const std::vector<Real>& v )
	{

		Real V   ( v[0] );
		Real cNa ( v[15] );
		Real cCa ( v[13] );


		return M_kNaCa * ( 1.0 / ( pow(M_KmNai, 3) + pow(M_NaO, 3) ) ) * ( 1.0 / ( M_KmCa + M_CaO ) )
			* (1.0 / ( 1.0 + M_kSat * exp( ( M_gamma - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) ) )
			* ( exp( M_gamma * ( V * M_F ) / ( M_R * M_T ) ) * pow(cNa, 3) * M_CaO -
						exp( ( M_gamma - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) * pow(M_NaO, 3) * cCa * 2.5 );
	}

	// Na+/K+ pump INaK
	Real IonicTenTusscher::pumpINaK( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cNa ( v[15] );

		Real fNaK  = 1.0 / ( 1.0 + 0.1245 * exp( -0.1 * ( V * M_F ) / ( M_R * M_T ) ) + 0.0353 * exp( -( V * M_F ) / ( M_R * M_T ) ) );

		return M_PNaK * fNaK * ( cNa / ( M_KmNa + cNa ) ) * ( M_KO / ( M_KO + M_KmK ) );
	}

	// Sarcolemmal Ca2+ pump current IpCa
	Real IonicTenTusscher::pumpIpCa( const std::vector<Real>& v )
	{
		return M_GCap * v[13] / ( M_KpCa + v[13] );
	}
	
	// K+ pump current IpK
	Real IonicTenTusscher::pumpIpK( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cKi ( v[16] );
		
		Real potK    = ( M_R * M_T / M_F ) * log( M_KO / cKi );
		
		return M_GKp * ( V - potK ) / ( 1 + exp( ( 25.0 - V ) / 5.98 ) );
	}

	// Ca2+ background current ICab
	Real IonicTenTusscher::backICab( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cCa ( v[13] );

		Real potCaN = ( ( M_R * M_T ) / ( 2 * M_F ) ) * log( M_CaO / cCa );

		return M_GCab * ( V - potCaN );
	}

	// Na+ background current INab
	Real IonicTenTusscher::backINab( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cNa ( v[15] );

		Real potNaN = M_R * M_T / M_F * log( M_NaO / cNa );

		return M_GNab * ( V - potNaN );
	}

	void IonicTenTusscher::showMe()
	{
		std::cout << "\n\n\t\tIonicJafriRiceWinslow Informations\n\n";
		std::cout << "number of unkowns: "  << this->Size() << std::endl;

		std::cout << "\n\t\tList of model parameters:\n\n";
		std::cout << "gasConst: " << this->gasConst() << std::endl;
		std::cout << "temp: " << this->temp() << std::endl;
		std::cout << "farad: " << this->farad() << std::endl;
		std::cout << "capMem: " << this->capMem() << std::endl;
		std::cout << "svRatio: " << this->svRatio() << std::endl;
		std::cout << "resCell: " << this->resCell() << std::endl;
		std::cout << "volCyt: " << this->volCyt() << std::endl;
		std::cout << "volSR: " << this->volSR() << std::endl;
		std::cout << "concKO: " << this->concKO() << std::endl;
		std::cout << "concNaO: " << this->concNaO() << std::endl;
		std::cout << "concCaO: " << this->concCaO() << std::endl;
		std::cout << "maxCondNa: " << this->maxCondNa() << std::endl;
		std::cout << "maxCondK1: " << this->maxCondK1() << std::endl;
		std::cout << "maxCondToEpiM: " << this->maxCondToEpiM() << std::endl;
		std::cout << "maxCondToEndo: " << this->maxCondToEndo() << std::endl;
		std::cout << "maxCondKr: " << this->maxCondKr() << std::endl;
		std::cout << "maxCondKs: " << this->maxCondKs() << std::endl;
		std::cout << "maxCondKsM: " << this->maxCondKsM() << std::endl;
		std::cout << "relPermKNa: " << this->relPermKNa() << std::endl;
		std::cout << "maxCondCaL: " << this->maxCondCaL() << std::endl;
		std::cout << "maxCourNaCa: " << this->maxCourNaCa() << std::endl;
		std::cout << "gamma: " << this->gamma() << std::endl;
		std::cout << "constmCa: " << this->constmCa() << std::endl;
		std::cout << "constmNai: " << this->constmNai() << std::endl;
		std::cout << "kSat: " << this->kSat() << std::endl;
		std::cout << "alpha: " << this->alpha() << std::endl;
		std::cout << "permNaK: " << this->permNaK() << std::endl;
		std::cout << "constmK: " << this->constmK() << std::endl;
		std::cout << "constmNa: " << this->constmNa() << std::endl;
		std::cout << "maxCondKp: " << this->maxCondKp() << std::endl;
		std::cout << "maxCondCap: " << this->maxCondCap() << std::endl;
		std::cout << "constpCa: " << this->constpCa() << std::endl;
		std::cout << "maxCondNab: " << this->maxCondNab() << std::endl;
		std::cout << "maxCondCab: " << this->maxCondCab() << std::endl;
		std::cout << "constUp: " << this->constUp() << std::endl;
		std::cout << "aRel: " << this->aRel() << std::endl;
		std::cout << "bRel: " << this->bRel() << std::endl;
		std::cout << "cRel: " << this->cRel() << std::endl;
		std::cout << "maxCourLeak: " << this->maxCourLeak() << std::endl;
		std::cout << "buffCyt: " << this->buffCyt() << std::endl;
		std::cout << "constBuffc: " << this->constBuffc() << std::endl;
		std::cout << "buffSR: " << this->buffSR() << std::endl;
		std::cout << "constBuffSR: " << this->constBuffSR() << std::endl;
		std::cout << "typeCell: " << this->typeCell() << std::endl;

		std::cout << "\n\t\t End of IonicTenTusscher Informations\n\n\n";
	}


	}

	#endif // TENTUSSCHER_HPP_INCLUDED
