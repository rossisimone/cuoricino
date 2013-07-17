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
	  @brief Ionic model based on Jafri, Rice And Winslow model.
	  @date 03-2013
	  @author Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

	  @contributors
	  @mantainer Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
	  @last update 03-2013
	 */


	#ifndef JAFRIRICEWINSLOW_HPP_INCLUDED
	#define JAFRIRICEWINSLOW_HPP_INCLUDED

	#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>
	#include <Teuchos_RCP.hpp>
	#include <Teuchos_ParameterList.hpp>
	#include "Teuchos_XMLParameterListHelpers.hpp"

	#include <cmath>


	namespace LifeV
	{
	//! XbModel - This class implements a mean field model.

	class IonicJafriRiceWinslow : public virtual ElectroIonicModel
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
		IonicJafriRiceWinslow();

		/*!
		 * @param Epetra communicator
		 * @param list of parameters in an xml file
		 */
		IonicJafriRiceWinslow ( Teuchos::ParameterList& parameterList );

		/*!
		 * @param IonicJafriRiceWinslow object
		 */
		IonicJafriRiceWinslow ( const IonicJafriRiceWinslow& model );

		//! Destructor
		virtual ~IonicJafriRiceWinslow() {}

		//@}

		//! @name Overloads
		//@{

		IonicJafriRiceWinslow& operator= ( const IonicJafriRiceWinslow& model );

		//@}

		//! @name Setters and getters
		//@{

		inline const Real& areaCap() const
		{
			return M_ACap;
		}
		inline void setACap(const Real& areaCap)
		{
			this->M_ACap = areaCap;
		}

		inline const Real& volMyo() const
		{
			return M_VMyo;
		}
		inline void setVMyo(const Real& volMyo)
		{
			this->M_VMyo = volMyo;
		}

		inline const Real& volJSR() const
		{
			return M_VJsr;
		}
		inline void setVJsr(const Real& volJSR)
		{
			this->M_VJsr = volJSR;
		}

		inline const Real& volNSR() const
		{
			return M_VNsr;
		}

		inline void setVNsr(const Real& volNSR)
		{
			this->M_VNsr = volNSR;
		}

		inline const Real& volSS() const
		{
			return M_VSs;
		}
		inline void setVSs(const Real& volSS)
		{
			this->M_VSs = volSS;
		}

		inline const Real& concNa0() const
		{
			return M_NaO;
		}
		inline void setNa0(const Real& concNa0)
		{
			this->M_NaO = concNa0;
		}

//		inline const Real& concK0() const
//		{
//			return M_KO;
//		}
//		inline void setK0(const Real& concK0)
//		{
//			this->M_KO = concK0;
//		}

		inline const Real& concCa0() const
		{
			return M_CaO;
		}
		inline void setCa0(const Real& concCa0)
		{
			this->M_CaO = concCa0;
		}

		inline const Real& lTrpnTot() const
		{
			return M_LtrpnTot;
		}
		inline void setLtrpnTot(const Real& lTrpnTot)
		{
			this->M_LtrpnTot = lTrpnTot;
		}

		inline const Real& hTrpnTot() const
		{
			return M_HtrpnTot;
		}
		inline void setHtrpnTot(const Real& hTrpnTot)
		{
			this->M_HtrpnTot = hTrpnTot;
		}

		inline const Real& kpHtrpn() const
		{
			return M_kPHtrpn;
		}
		inline void setkPHtrpn(const Real& kpHtrpn)
		{
			this->M_kPHtrpn = kpHtrpn;
		}

		inline const Real& knHtrpn() const
		{
			return M_kNHtrpn;
		}
		inline void setKNHtrpn(const Real& knHtrpn)
		{
			this->M_kNHtrpn = knHtrpn;
		}

		inline const Real& kpLtrpn() const
		{
			return M_kPLtrpn;
		}
		inline void setKPLtrpn(const Real& kpLtrpn)
		{
			this->M_kPLtrpn = kpLtrpn;
		}

		inline const Real& knLtrpn() const
		{
			return M_kNLtrpn;
		}
		inline void setKNLtrpn(const Real& knLtrpn)
		{
			this->M_kNLtrpn = knLtrpn;
		}

		inline const Real& cmdnTot() const
		{
			return M_CmdnTot;
		}
		inline void setCmdnTot(const Real& cmdnTot)
		{
			this->M_CmdnTot = cmdnTot;
		}

		inline const Real& csqnTot() const
		{
			return M_CsqnTot;
		}
		inline void setCsqnTot(const Real& csqnTot)
		{
			this->M_CsqnTot = csqnTot;
		}

		inline const Real& constmCmdn() const
		{
			return M_KmCmdn;
		}
		inline void setKmCmdn(const Real& constmCmdn)
		{
			this->M_KmCmdn = constmCmdn;
		}

		inline const Real& constmCsqn() const
		{
			return M_KmCsqn;
		}
		inline void setKmCsqn(const Real& constmCsqn)
		{
			this->M_KmCsqn = constmCsqn;
		}

		inline const Real& capMem() const
		{
			return M_Cm;
		}
		inline void setCapMem(const Real& capMem)
		{
			this->M_Cm = capMem;
		}

		inline const Real& farad() const
		{
			return M_F;
		}
		inline void setFarad(const Real& farad)
		{
			this->M_F = farad;
		}

		inline const Real& temp() const
		{
			return M_T;
		}
		inline void setTemp(const Real& temp)
		{
			this->M_T = temp;
		}

		inline const Real& gasConst() const
		{
			return M_R;
		}
		inline void setR(const Real& gasConst)
		{
			this->M_R = gasConst;
		}

		inline const Real& maxCondNa() const
		{
			return M_GNa;
		}
		inline void setGNa(const Real& maxCondNa)
		{
			this->M_GNa = maxCondNa;
		}

		inline const Real& maxCondKp() const
		{
			return M_GKp;
		}
		inline void setGkp(const Real& maxCondKp)
		{
			this->M_GKp = maxCondKp;
		}

		inline const Real& permNaK() const
		{
			return M_PNaK;
		}
		inline void setPNaK(const Real& permNaK)
		{
			this->M_PNaK = permNaK;
		}

		inline const Real& kNaCa() const
		{
			return M_kNaCa;
		}
		inline void setKNaCa(const Real& kNaCa)
		{
			this->M_kNaCa = kNaCa;
		}

		inline const Real& constmNa() const
		{
			return M_KmNa;
		}
		inline void setKmNa(const Real& constmNa)
		{
			this->M_KmNa = constmNa;
		}

		inline const Real& constmCa() const
		{
			return M_KmCa;
		}
		inline void setKmCa(const Real& constmCa)
		{
			this->M_KmCa = constmCa;
		}

		inline const Real& kSat() const
		{
			return M_kSat;
		}
		inline void setKSat(const Real& kSat)
		{
			this->M_kSat = kSat;
		}

		inline const Real& eta() const
		{
			return M_eta;
		}
		inline void setEta(const Real& eta)
		{
			this->M_eta = eta;
		}

		inline const Real& courNaK() const
		{
			return M_INaK;
		}
		inline void setINaK(const Real& courNaK)
		{
			this->M_INaK = courNaK;
		}

		inline const Real& constmNai() const
		{
			return M_KmNai;
		}
		inline void setKmNai(const Real& constmNai)
		{
			this->M_KmNai = constmNai;
		}

		inline const Real& constmK0() const
		{
			return M_KmK0;
		}
		inline void setKmK0(const Real& constmK0)
		{
			this->M_KmK0 = constmK0;
		}

		inline const Real& permNsCa() const
		{
			return M_PnsCa;
		}
		inline void setPnsCa(const Real& permNsCa)
		{
			this->M_PnsCa = permNsCa;
		}

		inline const Real& constmNsCa() const
		{
			return M_KmNsCa;
		}
		inline void setKmNsCa(const Real& constmNsCa)
		{
			this->M_KmNsCa = constmNsCa;
		}

		inline const Real& courpCa() const
		{
			return M_IpCa;
		}
		inline void setIpCa(const Real& courpCa)
		{
			this->M_IpCa = courpCa;
		}

		inline const Real& constmpCa() const
		{
			return M_KmPCa;
		}
		inline void setKmPCa(const Real& constmpCa)
		{
			this->M_KmPCa = constmpCa;
		}

		inline const Real& maxCondCab() const
		{
			return M_GCab;
		}
		inline void setGCab(const Real& maxCondCab)
		{
			this->M_GCab = maxCondCab;
		}

		inline const Real& maxCondNab() const
		{
			return M_GNab;
		}
		inline void setGNab(const Real& maxCondNab)
		{
			this->M_GNab = maxCondNab;
		}

		inline const Real& maxRyRPerm() const
		{
			return M_v1;
		}
		inline void setV1(const Real& maxRyRPerm)
		{
			this->M_v1 = maxRyRPerm;
		}

		inline const Real& leakRateConst() const
		{
			return M_v2;
		}
		inline void setV2(const Real& leakRateConst)
		{
			this->M_v2 = leakRateConst;
		}

		inline const Real& pumpRateATPase() const
		{
			return M_v3;
		}
		inline void setV3(const Real& pumpRateATPase)
		{
			this->M_v3 = pumpRateATPase;
		}

		inline const Real& constmUp() const
		{
			return M_KmUp;
		}
		inline void setKmUp(const Real& constmUp)
		{
			this->M_KmUp = constmUp;
		}

		inline const Real& timeConstNsrJsr() const
		{
			return M_tauTr;
		}
		inline void setTauTr(const Real& timeConstNsrJsr)
		{
			this->M_tauTr = timeConstNsrJsr;
		}

		inline const Real& timeConstSubMyo() const
		{
			return M_tauXFer;
		}
		inline void setTauXFer(const Real& timeConstSubMyo)
		{
			this->M_tauXFer = timeConstSubMyo;
		}

		inline const Real& kAPlus() const
		{
			return M_kap;
		}
		inline void setKap(const Real& kAPlus)
		{
			this->M_kap = kAPlus;
		}

		inline const Real& kANeg() const
		{
			return M_kan;
		}
		inline void setKan(const Real& kANeg)
		{
			this->M_kan = kANeg;
		}

		inline const Real& kBPlus() const
		{
			return M_kbp;
		}
		inline void setKbp(const Real& kBPlus)
		{
			this->M_kbp = kBPlus;
		}

		inline const Real& kBNeg() const
		{
			return M_kbn;
		}
		inline void setKbn(const Real& kBNeg)
		{
			this->M_kbn = kBNeg;
		}

		inline const Real& kCPlus() const
		{
			return M_kcp;
		}
		inline void setKcp(const Real& kCPlus)
		{
			this->M_kcp = kCPlus;
		}

		inline const Real& kCNeg() const
		{
			return M_kcn;
		}
		inline void setKcn(const Real& kCNeg)
		{
			this->M_kcn = kCNeg;
		}

		inline const Real& coopParamN() const
		{
			return M_n;
		}
		inline void setN(const Real& coopParamN)
		{
			this->M_n = coopParamN;
		}

		inline const Real& coopParamM() const
		{
			return M_m;
		}
		inline void setM(const Real& coopParamM)
		{
			this->M_m = coopParamM;
		}

		inline const Real& intoOpenSt() const
		{
			return M_f;
		}
		inline void setF(const Real& intoOpenSt)
		{
			this->M_f = intoOpenSt;
		}

		inline const Real& outOpenSt() const
		{
			return M_g;
		}
		inline void setG(const Real& outOpenSt)
		{
			this->M_g = outOpenSt;
		}

		inline const Real& intoOpenStCa() const
		{
			return M_fprime;
		}
		inline void setFprime(const Real& intoOpenStCa)
		{
			this->M_fprime = intoOpenStCa;
		}

		inline const Real& outOpenSt2() const
		{
			return M_gprime;
		}
		inline void setGprime(const Real& outOpenSt2)
		{
				this->M_gprime = outOpenSt2;
		}

		inline const Real& modeTParamA() const
		{
			return M_a;
		}
		inline void setA(const Real& modeTParamA)
		{
			this->M_a = modeTParamA;
		}

		inline const Real& modeTParamB() const
		{
			return M_b;
		}
		inline void setB(const Real& modeTParamB)
		{
			this->M_b = modeTParamB;
		}

		inline const Real& modeTParamO() const
		{
			return M_omega;
		}
		inline void setOmega(const Real& modeTParamO)
		{
			this->M_kcn = modeTParamO;
		}

		inline const Real& permCa() const
		{
			return M_PCa;
		}
		inline void setPCa(const Real& permCa)
		{
			this->M_PCa = permCa;
		}

		inline const Real& permK() const
		{
			return M_PK;
		}
		inline void setPK(const Real& permK)
		{
			this->M_PK = permK;
		}

		inline const Real& courCaHalf() const
		{
			return M_ICaHalf;
		}
		inline void setICaHalf(const Real& courCaHalf)
		{
			this->M_ICaHalf = courCaHalf;
		}

		//@}

		//! @name Methods
		//@{

		//Compute the rhs on a single node or for the 0D case
		void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

		void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );


		//Compute the rhs with state variable interpolation
		Real computeLocalPotentialRhs ( const std::vector<Real>& v );

		std::vector<Real> computeLocalGatingRhs ( const std::vector<Real>& v );

		std::vector<Real> computeLocalSubSysCaRhs( const std::vector<Real>& v );

		std::vector<Real> computeLocalChannelRyrRhs( const std::vector<Real>& v );

		std::vector<Real> computeLocalConcRhs ( const std::vector<Real>& v );

		//Compute the ionic currents (Luo and Rudy)
		std::vector<Real> fastINa( const std::vector<Real>& v );

		std::vector<Real> timeDIK( const std::vector<Real>& v );

		Real timeIIK1( const std::vector<Real>& v );

		Real plaIKp( const std::vector<Real>& v );

		Real exINaCa( const std::vector<Real>& v );

		Real pumpINaK( const std::vector<Real>& v );

		std::vector<Real> noSpecInsCa( const std::vector<Real>& v );

		Real pumpIpCa( const std::vector<Real>& v );

		Real backICab( const std::vector<Real>& v );

		Real backINab( const std::vector<Real>& v);

		//! For resolution with implicit numerical integrator

//		std::vector<Real> computeYParameters( const std::vector<Real>& v);
//		std::vector<Real> channelLCaMatrix( const std::vector<Real>& v );
        Real computeNewtonNa   ( const std::vector<Real>& v, const Real& dt, const int& nitermax );
        Real computeNewtonKi   ( const std::vector<Real>& v, const Real& dt, const int& nitermax );
        Real computeNewtonKo   ( const std::vector<Real>& v, const Real& dt, const int& nitermax );
        Real computeNewtonCai  ( const std::vector<Real>& v, const Real& dt, const int& nitermax );
        Real computeNewtonCaNSR( const std::vector<Real>& v, const Real& dt, const int& nitermax );
		Real computeNewtonCaSS ( const std::vector<Real>& v, const Real& dt, const int& nitermax );
		Real computeNewtonCaJSR( const std::vector<Real>& v, const Real& dt, const int& nitermax );

		//! For resolution with Rush and Larsen numerical integrator
		std::vector<Real> gateInf( const std::vector<Real>& v );
		std::vector<Real> otherVarInf( const std::vector<Real>& v );

		//! Display information about the model
		void showMe();

		//@}

	private:
		//! Cell Geometry Parameters (5)

		Real M_ACap;
		Real M_VMyo;
		Real M_VJsr;
		Real M_VNsr;
		Real M_VSs;

		//! Standard Ionic Concentrations (3)

		Real M_NaO;
//		Real M_KO;
		Real M_CaO;

		//! Buffering Parameters (10)

		Real M_LtrpnTot;
		Real M_HtrpnTot;
		Real M_kPHtrpn;
		Real M_kNHtrpn;
		Real M_kPLtrpn;
		Real M_kNLtrpn;
		Real M_CmdnTot;
		Real M_CsqnTot;
		Real M_KmCmdn;
		Real M_KmCsqn;

		//! Membrane Current Parameters (20)

		Real M_Cm;
		Real M_F;
		Real M_T;
		Real M_R;
		Real M_GNa;
		Real M_GKp;
		Real M_PNaK;
		Real M_kNaCa;
		Real M_KmNa;
		Real M_KmCa;
		Real M_kSat;
		Real M_eta;
		Real M_INaK;
		Real M_KmNai;
		Real M_KmK0;
		Real M_PnsCa;
		Real M_KmNsCa;
		Real M_IpCa;
		Real M_KmPCa;
		Real M_GCab;
		Real M_GNab;

		//! SR Parameters (14)

		Real M_v1;
		Real M_v2;
		Real M_v3;
		Real M_KmUp;
		Real M_tauTr;
		Real M_tauXFer;
		Real M_kap;
		Real M_kan;
		Real M_kbp;
		Real M_kbn;
		Real M_kcp;
		Real M_kcn;
		Real M_n;
		Real M_m;

		//! L-type Ca2+ Channel Parameters (10)

		Real M_f;
		Real M_g;
		Real M_fprime;
		Real M_gprime;
		Real M_a;
		Real M_b;
		Real M_omega;
		Real M_PCa;
		Real M_PK;
		Real M_ICaHalf;

		//@}



	}; // class IonicJafriRiceWinslow


	}

	#endif // IonicJafriRiceWinslow_HPP_INCLUDED
