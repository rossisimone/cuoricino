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

	#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

	#include <Teuchos_RCP.hpp>
	#include <Teuchos_ParameterList.hpp>
	#include "Teuchos_XMLParameterListHelpers.hpp"

	#include <cmath>
	#include <string>

	namespace LifeV
	{
	//! XbModel - This class implements a mean field model.

	class IonicTenTusscher : public virtual ElectroIonicModel
	{

	public:
		//! @name Type definitions
		//@{
		typedef ElectroIonicModel super;
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

		std::vector<Real> gateInf( const std::vector<Real>& v );

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


	}

	#endif // TENTUSSCHER_HPP_INCLUDED
