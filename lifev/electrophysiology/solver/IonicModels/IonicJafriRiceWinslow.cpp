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
	#include <lifev/electrophysiology/solver/IonicModels/IonicJafriRiceWinslow.hpp>


	#include <cmath>


	namespace LifeV
	{
	// ===================================================
	//! Constructors
	// ===================================================
	IonicJafriRiceWinslow::IonicJafriRiceWinslow()    :
		super       ( 31 ),
		M_ACap      ( 1.5343e-2 ),
		M_VMyo      ( 25.84e-6 ),
		M_VJsr      ( 0.12e-6 ),
		M_VNsr      ( 2.098e-6 ),
		M_VSs       ( 1.485e-9 ),
		M_NaO       ( 140.0 ),
//		M_KO        ( 5.4 ),
		M_CaO       ( 1.8 ),
		M_LtrpnTot  ( 0.07 ),
		M_HtrpnTot  ( 0.14 ),
		M_kPHtrpn   ( 20.0 ),
		M_kNHtrpn   ( 0.066e-3 ),
		M_kPLtrpn   ( 40.0 ),
		M_kNLtrpn   ( 0.04 ),
		M_CmdnTot   ( 0.05 ),
		M_CsqnTot   ( 15.0 ),
		M_KmCmdn    ( 2.38e-3 ),
		M_KmCsqn    ( 0.8 ),
		M_Cm   		( 0.01 ),
		M_F         ( 9.6485e4 ),
		M_T         ( 310.0 ),
		M_R         ( 8.314e3 ),
		M_GNa       ( 0.128 ),
		M_GKp       ( 8.28e-5 ),
		M_PNaK      ( 0.01833 ),
		M_kNaCa     ( 10.0 ),
		M_KmNa      ( 87.5 ),
		M_KmCa      ( 1.38 ),
		M_kSat      ( 0.1 ),
		M_eta       ( 0.35 ),
		M_INaK      ( 2.724e-2 ),
		M_KmNai     ( 10.0 ),
		M_KmK0      ( 1.5 ),
		M_PnsCa     ( 0.0 ),
		M_KmNsCa    ( 1.2e-3 ),
		M_IpCa      ( 1.15e-2 ),
		M_KmPCa     ( 0.5e-3 ),
		M_GCab      ( 6.032e-5 ),
		M_GNab      ( 1.41e-5 ),
		M_v1       	( 1.8 ),
		M_v2        ( 3.6e-4 ),
		M_v3        ( 6.375e-3 ),
		M_KmUp      ( 0.5e-3 ),
		M_tauTr     ( 34.48 ),
		M_tauXFer   ( 3.125 ),
		M_kap       ( 1.215e10 ),
		M_kan       ( 0.1425 ),
		M_kbp       ( 4.05e6 ),
		M_kbn       ( 1.93 ),
		M_kcp       ( 0.018 ),
		M_kcn       ( 0.0008 ),
		M_n 		( 4.0 ),
		M_m         ( 3.0 ),
		M_f         ( 0.3 ),
		M_g         ( 2.0 ),
		M_fprime    ( 0.0 ),
		M_gprime    ( 0.0 ),
		M_a         ( 2.0 ),
		M_b         ( 2.0 ),
		M_omega     ( 0.01 ),
		M_PCa       ( 33.75e-6 ),
		M_PK        ( 1.0e-9 ),
		M_ICaHalf   ( -4.58e-3 )
	{}

	IonicJafriRiceWinslow::IonicJafriRiceWinslow ( Teuchos::ParameterList& parameterList ) :
		super       ( 31 )
	{
		M_ACap      = parameterList.get ( "areaCap", 1.5343e-2);
		M_VMyo      = parameterList.get ( "volMyo", 25.84e-6 );
		M_VJsr      = parameterList.get ( "volJSR", 0.12e-6 );
		M_VNsr      = parameterList.get ( "volNSR", 2.098e-6 );
		M_VSs       = parameterList.get ( "volSS", 1.485e-9 );
		M_NaO       = parameterList.get ( "concNa0", 140.0 );
//		M_KO        = parameterList.get ( "concK0", 5.4 );
		M_CaO       = parameterList.get ( "concCa0", 1.8 );
		M_LtrpnTot  = parameterList.get ( "lTrpnTot", 0.07 );
		M_HtrpnTot  = parameterList.get ( "hTrpnTot", 0.14 );
		M_kPHtrpn   = parameterList.get ( "kpHtrpn", 20.0 );
		M_kNHtrpn   = parameterList.get ( "knHtrpn", 0.066e-3 );
		M_kPLtrpn   = parameterList.get ( "kpLtrpn", 40.0 );
		M_kNLtrpn   = parameterList.get ( "knLtrpn", 0.04 );
		M_CmdnTot   = parameterList.get ( "cmdnTot", 0.05 );
		M_CsqnTot   = parameterList.get ( "csqnTot", 15.0 );
		M_KmCmdn    = parameterList.get ( "constmCmdn", 2.38e-3 );
		M_KmCsqn    = parameterList.get ( "constmCsqn", 0.8 );
		M_Cm		= parameterList.get ( "capMem", 0.01 );
		M_F         = parameterList.get ( "farad", 9.6485e4 );
		M_T         = parameterList.get ( "temp", 310.0 );
		M_R         = parameterList.get ( "gasConst", 8.3145e3 );
		M_GNa       = parameterList.get ( "maxCondNa", 0.128 );
		M_GKp       = parameterList.get ( "maxCondKp", 8.28e-5 );
		M_PNaK      = parameterList.get ( "permNaK", 0.01833 );
		M_kNaCa     = parameterList.get ( "kNaCa", 10.0 );
		M_KmNa      = parameterList.get ( "constmNa", 87.5 );
		M_KmCa      = parameterList.get ( "constmCa", 1.38 );
		M_kSat      = parameterList.get ( "kSat", 0.1 );
		M_eta       = parameterList.get ( "eta", 0.35 );
		M_INaK      = parameterList.get ( "courNaK", 2.742e-2 );
		M_KmNai     = parameterList.get ( "constmNai", 10.0 );
		M_KmK0      = parameterList.get ( "constmK0", 1.5 );
		M_PnsCa     = parameterList.get ( "permNsCa", 0.0 );
		M_KmNsCa    = parameterList.get ( "constmNsCa", 1.2e-3 );
		M_IpCa      = parameterList.get ( "courpCa", 1.15e-2 );
		M_KmPCa     = parameterList.get ( "constmpCa", 0.5e-3 );
		M_GCab      = parameterList.get ( "maxCondCab", 6.032e-5 );
		M_GNab      = parameterList.get ( "maxCondNab", 1.41e-5 );
		M_v1       	= parameterList.get ( "maxRyRPerm", 1.8 );
		M_v2        = parameterList.get ( "leakRateConst", 3.6e-4 );
		M_v3        = parameterList.get ( "pumpRateATPase", 6.375e-3 );
		M_KmUp      = parameterList.get ( "constmUp", 0.5e-3 );
		M_tauTr     = parameterList.get ( "timeConstNsrJsr", 34.48 );
		M_tauXFer   = parameterList.get ( "timeConstSubMyo", 3.125 );
		M_kap       = parameterList.get ( "kAPlus", 1.215e10 );
		M_kan       = parameterList.get ( "kANeg", 0.1425 );
		M_kbp       = parameterList.get ( "kBPlus", 4.05e6 );
		M_kbn       = parameterList.get ( "kBNeg", 1.93 );
		M_kcp       = parameterList.get ( "kCPlus", 0.018 );
		M_kcn       = parameterList.get ( "kCNeg", 0.0008 );
		M_n 		= parameterList.get ( "coopParamN", 4.0 );
		M_m         = parameterList.get ( "coopParamM", 3.0 );
		M_f         = parameterList.get ( "intoOpenSt", 0.3 );
		M_g         = parameterList.get ( "outOpenSt", 2.0 );
		M_fprime    = parameterList.get ( "intoOpenStCa", 0.0 );
		M_gprime    = parameterList.get ( "outOpenSt2", 0.0 );
		M_a         = parameterList.get ( "modeTParamA", 2.0 );
		M_b         = parameterList.get ( "modeTParamB", 2.0 );
		M_omega     = parameterList.get ( "modeTParamO", 0.01 );
		M_PCa       = parameterList.get ( "permCa", 33.75e-6 );
		M_PK        = parameterList.get ( "permK", 1.0e-9 );
		M_ICaHalf   = parameterList.get ( "courCaHalf", -4.58e-3 );
	}

	IonicJafriRiceWinslow::IonicJafriRiceWinslow ( const IonicJafriRiceWinslow& model )
	{
		M_ACap		= model.M_ACap;
		M_VMyo		= model.M_VMyo;
		M_VJsr		= model.M_VJsr;
		M_VNsr		= model.M_VNsr;
		M_VSs		= model.M_VSs;
		M_NaO		= model.M_NaO;
//		M_KO		= model.M_KO;
		M_CaO 		= model.M_CaO;
		M_LtrpnTot	= model.M_LtrpnTot;
		M_HtrpnTot	= model.M_HtrpnTot;
		M_kPHtrpn	= model.M_kPHtrpn;
		M_kNHtrpn	= model.M_kNHtrpn;
		M_kPLtrpn	= model.M_kPLtrpn;
		M_kNLtrpn	= model.M_kNLtrpn;
		M_CmdnTot	= model.M_CmdnTot;
		M_CsqnTot	= model.M_CsqnTot;
		M_KmCmdn	= model.M_KmCmdn;
		M_KmCsqn	= model.M_KmCsqn;
		M_Cm		= model.M_Cm;
		M_F		    = model.M_F;
		M_T		    = model.M_T;
		M_R		    = model.M_R;
		M_GNa		= model.M_GNa;
		M_GKp		= model.M_GKp;
		M_PNaK		= model.M_PNaK;
		M_kNaCa	    = model.M_kNaCa;
		M_KmNa		= model.M_KmNa;
		M_KmCa		= model.M_KmCa;
		M_kSat		= model.M_kSat;
		M_eta		= model.M_eta;
		M_INaK		= model.M_INaK;
		M_KmNai	    = model.M_KmNai;
		M_KmK0		= model.M_KmK0;
		M_PnsCa	    = model.M_PnsCa;
		M_KmNsCa	= model.M_KmNsCa;
		M_IpCa		= model.M_IpCa;
		M_KmPCa	    = model.M_KmPCa;
		M_GCab		= model.M_GCab;
		M_GNab		= model.M_GNab;
		M_v1       	= model.M_v1;
		M_v2        = model.M_v2;
		M_v3        = model.M_v3;
		M_KmUp      = model.M_KmUp;
		M_tauTr     = model.M_tauTr;
		M_tauXFer   = model.M_tauXFer;
		M_kap       = model.M_kap;
		M_kan       = model.M_kan;
		M_kbp       = model.M_kbp;
		M_kbn       = model.M_kbn;
		M_kcp       = model.M_kcp;
		M_kcn       = model.M_kcn;
		M_n 		= model.M_n;
		M_m         = model.M_m;
		M_f         = model.M_f;
		M_g         = model.M_g;
		M_fprime    = model.M_fprime;
		M_gprime    = model.M_gprime;
		M_a         = model.M_a;
		M_b         = model.M_b;
		M_omega     = model.M_omega;
		M_PCa       = model.M_PCa;
		M_PK        = model.M_PK;
		M_ICaHalf   = model.M_ICaHalf;

		M_numberOfEquations = model.M_numberOfEquations;
	}

	// ===================================================
	//! Operator
	// ===================================================
	IonicJafriRiceWinslow& IonicJafriRiceWinslow::operator= ( const IonicJafriRiceWinslow& model )
	{
		M_ACap      = model.M_ACap;
		M_VMyo      = model.M_VMyo;
		M_VJsr      = model.M_VJsr;
		M_VNsr      = model.M_VNsr;
		M_VSs       = model.M_VSs;
		M_NaO       = model.M_NaO;
//		M_KO		= model.M_KO;
		M_CaO       = model.M_CaO;
		M_LtrpnTot  = model.M_LtrpnTot;
		M_HtrpnTot  = model.M_HtrpnTot;
		M_kPHtrpn   = model.M_kPHtrpn;
		M_kNHtrpn   = model.M_kNHtrpn;
		M_kPLtrpn   = model.M_kPLtrpn;
		M_kNLtrpn   = model.M_kNLtrpn;
		M_CmdnTot   = model.M_CmdnTot;
		M_CsqnTot   = model.M_CsqnTot;
		M_KmCmdn    = model.M_KmCmdn;
		M_KmCsqn    = model.M_KmCsqn;
		M_F         = model.M_F;
		M_T         = model.M_T;
		M_R         = model.M_R;
		M_GNa       = model.M_GNa;
		M_GKp       = model.M_GKp;
		M_PNaK      = model.M_PNaK;
		M_kNaCa     = model.M_kNaCa;
		M_KmNa      = model.M_KmNa;
		M_KmCa      = model.M_KmCa;
		M_kSat      = model.M_kSat;
		M_eta       = model.M_eta;
		M_INaK      = model.M_INaK;
		M_KmNai     = model.M_KmNai;
		M_KmK0      = model.M_KmK0;
		M_PnsCa     = model.M_PnsCa;
		M_KmNsCa    = model.M_KmNsCa;
		M_IpCa      = model.M_IpCa;
		M_KmPCa     = model.M_KmPCa;
		M_GCab      = model.M_GCab;
		M_GNab      = model.M_GNab;
		M_v1       	= model.M_v1;
		M_v2        = model.M_v2;
		M_v3        = model.M_v3;
		M_KmUp      = model.M_KmUp;
		M_tauTr     = model.M_tauTr;
		M_tauXFer   = model.M_tauXFer;
		M_kap       = model.M_kap;
		M_kan       = model.M_kan;
		M_kbp       = model.M_kbp;
		M_kbn       = model.M_kbn;
		M_kcp       = model.M_kcp;
		M_kcn       = model.M_kcn;
		M_n 		= model.M_n;
		M_m         = model.M_m;
		M_f         = model.M_f;
		M_g         = model.M_g;
		M_fprime    = model.M_fprime;
		M_gprime    = model.M_gprime;
		M_a         = model.M_a;
		M_b         = model.M_b;
		M_omega     = model.M_omega;
		M_PCa       = model.M_PCa;
		M_PK        = model.M_PK;
		M_ICaHalf   = model.M_ICaHalf;

		M_numberOfEquations = model.M_numberOfEquations;

		return *this;
	}


	// ===================================================
	//! Methods
	// ===================================================
	//Only gating variables
	void IonicJafriRiceWinslow::computeRhs ( const std::vector<Real>&  v,
											 std::vector<Real>& rhs )
	{
		std::vector<Real> gatingRhs     ( computeLocalGatingRhs(v) );
		std::vector<Real> concRhs       ( computeLocalConcRhs(v) );
		std::vector<Real> subSysCaRhs   ( computeLocalSubSysCaRhs(v) );
		std::vector<Real> channelRyrRhs ( computeLocalChannelRyrRhs(v) );

		std::copy( gatingRhs.begin(), gatingRhs.end(), rhs.begin() );

		std::copy( concRhs.begin(), concRhs.end() - 2, rhs.begin() + gatingRhs.size() );

		std::copy( subSysCaRhs.begin() + 2, subSysCaRhs.begin() + 6 , rhs.begin() + gatingRhs.size() + ( concRhs.size() - 2 ) );

		std::copy( channelRyrRhs.begin(), channelRyrRhs.end(), rhs.begin() + gatingRhs.size() + ( concRhs.size() - 2 ) + 4 );

		std::copy( subSysCaRhs.begin() + 6, subSysCaRhs.end(), rhs.begin() + gatingRhs.size() + ( concRhs.size() - 2 ) + 4 + channelRyrRhs.size() );

		std::copy( concRhs.begin() + 3, concRhs.end(), rhs.end() - 2);

	}

	//Potential and gating variables
	void IonicJafriRiceWinslow::computeRhs (const   std::vector<Real>&  v,
											 const   Real&           Iapp,
											 std::vector<Real>& rhs )
	{
		std::vector<Real> gatingRhs     ( computeLocalGatingRhs(v) );
		std::vector<Real> concRhs       ( computeLocalConcRhs(v) );
		std::vector<Real> subSysCaRhs   ( computeLocalSubSysCaRhs(v) );
		std::vector<Real> channelRyrRhs ( computeLocalChannelRyrRhs(v) );


		rhs[0] = computeLocalPotentialRhs(v, Iapp);

		std::copy( gatingRhs.begin(), gatingRhs.end(), rhs.begin() + 1 );

		std::copy( concRhs.begin(), concRhs.end() - 2, rhs.begin() + 1 + gatingRhs.size() );

		std::copy( subSysCaRhs.begin() + 2, subSysCaRhs.begin() + 6 , rhs.begin() + 1 + gatingRhs.size() + ( concRhs.size() - 2 ) );

		std::copy( channelRyrRhs.begin(), channelRyrRhs.end(), rhs.begin() + 1 + gatingRhs.size() + ( concRhs.size() - 2 ) + 4 );

		std::copy( subSysCaRhs.begin() + 6, subSysCaRhs.end(), rhs.begin() + 1 + gatingRhs.size() + ( concRhs.size() - 2 ) + 4 + channelRyrRhs.size() );

		std::copy( concRhs.begin() + 3, concRhs.end(), rhs.end() - 2);

	}


	Real IonicJafriRiceWinslow::computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp )
	{
		std::vector<Real> courSubSysCa ( computeLocalSubSysCaRhs(v) );
		std::vector<Real> courINa	   ( fastINa(v) );
		std::vector<Real> courDIK 	   ( timeDIK(v) );
		std::vector<Real> courInsCa    ( noSpecInsCa(v) );

		return  - ( 1 / M_Cm ) * ( courINa[0] + courSubSysCa[0] + courDIK[0] + timeIIK1(v) + plaIKp(v) + exINaCa(v)
				+ pumpINaK(v) + courInsCa[3] + pumpIpCa(v) + backICab(v) + courSubSysCa[1] + backINab(v) + Iapp );
	}

	std::vector<Real> IonicJafriRiceWinslow::computeLocalGatingRhs ( const std::vector<Real>& v )
	{
		Real m ( v[1] );
		Real h ( v[2] );
		Real j ( v[3] );
		Real x ( v[4] );

		std::vector<Real> gatingRhs (4);
		std::vector<Real> param1    ( fastINa(v) );
		std::vector<Real> param2    ( timeDIK(v) );

		gatingRhs[0] = param1[1] * ( 1 - m ) - param1[2] * m;
		gatingRhs[1] = param1[3] * ( 1 - h ) - param1[4] * h;
		gatingRhs[2] = param1[5] * ( 1 - j ) - param1[6] * j;
		gatingRhs[3] = param2[1] * ( 1 - x ) - param2[2] * x;
//		gatingRhs[0] = param1[1];
//		gatingRhs[1] = param1[2];
//		gatingRhs[2] = param1[3];
//		gatingRhs[3] = param2[1];

		return gatingRhs;
	}

	std::vector<Real> IonicJafriRiceWinslow::computeLocalConcRhs ( const std::vector<Real>& v )
	{

		std::vector<Real> courSubSysCa ( computeLocalSubSysCaRhs(v) );
		std::vector<Real> courINa	   ( fastINa(v) );
		std::vector<Real> courDIK 	   ( timeDIK(v) );
		std::vector<Real> courInsCa    ( noSpecInsCa(v) );

		std::vector<Real> concRhs(5);

		concRhs[0] = - ( courINa[0] + backINab(v) + courInsCa[0] + 3 * exINaCa(v) + 3 * pumpINaK(v) ) * M_ACap / ( M_VMyo * M_F );
		concRhs[1] = - ( courDIK[0] + timeIIK1(v) + plaIKp(v) + courInsCa[1] - 2 * pumpINaK(v) + courSubSysCa[1] ) * M_ACap / ( M_VMyo * M_F );
		concRhs[2] =   ( courDIK[0] + timeIIK1(v) + plaIKp(v) + courInsCa[1] - 2 * pumpINaK(v) + courSubSysCa[1] ) * M_ACap / ( M_VMyo * M_F );
		concRhs[3] = M_kPLtrpn * v[8] * ( M_LtrpnTot - v[29] ) - M_kNLtrpn * v[29];
		concRhs[4] = M_kPHtrpn * v[8] * ( M_HtrpnTot - v[30] ) - M_kNHtrpn * v[30];
//		concRhs[3] = - ( M_kNLtrpn + M_kPLtrpn * std::pow(v[8], M_n) );
//		concRhs[4] = - ( M_kNHtrpn + M_kPHtrpn * std::pow(v[8], M_n) );

		return concRhs;
	}



	//! Ca2+ Subsystem

	std::vector<Real> IonicJafriRiceWinslow::computeLocalSubSysCaRhs( const std::vector<Real>& v )
	{
		std::vector<Real> subSysCaRHS(19);

		Real V 	      ( v[0] );
		Real cKi 	  ( v[6] );
		Real cKo	  ( v[7] );
		Real cCa      ( v[8] );
		Real cCaNSR   ( v[9] );
		Real cCaSS    ( v[10] );
		Real cCaJSR   ( v[11] );

		Real fracPO1  ( v[13] );
		Real fracPO2  ( v[14] );

		Real c0       ( v[16] );
		Real c1 	  ( v[17] );
		Real c2       ( v[18] );
		Real c3       ( v[19] );
		Real c4       ( v[20] );
		Real o        ( v[21] );
		Real cCa0     ( v[22] );
		Real cCa1     ( v[23] );
		Real cCa2     ( v[24] );
		Real cCa3     ( v[25] );
		Real cCa4     ( v[26] );
		Real oCa      ( v[27] );
		Real y        ( v[28] );
		Real cLTrpnCa ( v[29] );
		Real cHTrpnCa ( v[30] );


		// Internal Parameters


		Real jRel  = M_v1 * ( fracPO1 + fracPO2 ) * ( cCaJSR - cCaSS );
		Real jLeak = M_v2 * ( cCaNSR - cCa );
		Real jUp   = M_v3 * ( cCa * cCa ) / ( M_KmUp * M_KmUp + cCa * cCa );
		Real jTr   = ( cCaNSR - cCaJSR ) / M_tauTr;
		Real jXFer = ( cCaSS - cCa ) / M_tauXFer;
		Real jTRPN = M_kPHtrpn * cCa * ( M_HtrpnTot - cHTrpnCa ) - M_kNHtrpn * cHTrpnCa;
//						+ M_kPLtrpn * cCa * ( M_LtrpnTot - cLTrpnCa ) - M_kNLtrpn * cLTrpnCa;

		Real bI   = 1 / ( 1 + M_CmdnTot * M_KmCmdn / ( ( M_KmCmdn + cCa ) * ( M_KmCmdn + cCa ) ) );
		Real bSS  = 1 / ( 1 + M_CmdnTot * M_KmCmdn / ( ( M_KmCmdn + cCaSS ) * ( M_KmCmdn + cCaSS ) ) );
		Real bJSR = 1 / ( 1 + M_CsqnTot * M_KmCsqn / ( ( M_KmCsqn + cCaJSR ) * ( M_KmCsqn + cCaJSR ) ) );

		Real alpha      = 0.4 * std::exp( ( V + 12.0 ) / 10.0 );
		Real beta       = 0.05 * std::exp( - ( V + 12 ) / 13.0 );
		Real alphaprime = M_a * alpha;
		Real betaprime  = beta / M_b;
		Real gamma      = 0.1875 * cCaSS;

		Real yinf   = 1 / ( 1 + std::exp( ( V + 55.0 ) / 7.5 ) ) + 0.1 / ( 1 + std::exp( ( -V + 21.0 ) / 6.0 ) );
		Real tauY   = 20 + 600 / ( 1 + std::exp( ( V + 30 ) / 9.5 ) );
		Real iCaMax = M_PCa * 4 * ( V * M_F * M_F ) / ( M_R * M_T ) * ( 0.001 * std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 0.341 * M_CaO )
						/ ( std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );
		Real iCa    = y * ( o + oCa ) * iCaMax;
		Real pK     = M_PK / ( 1.0 + iCaMax / M_ICaHalf );
		Real iCaK   = pK * y * ( o + oCa ) * ( V * M_F * M_F ) / ( M_R * M_T ) * ( cKi * std::exp( ( V * M_F ) / ( M_R * M_T ) ) - cKo )
						/ ( std::exp( ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );


		// RHS of the Ca2+ subsystem

		subSysCaRHS[0] 	= iCa;
		subSysCaRHS[1] 	= iCaK;

		subSysCaRHS[2] 	= bI * ( jLeak + jXFer - ( jUp + jTRPN
							+ ( backICab(v) - 2 * exINaCa(v) + pumpIpCa(v) ) * M_ACap / ( 2 * M_VMyo * M_F ) ) );

		subSysCaRHS[3] 	= ( jUp - jLeak ) * M_VMyo / M_VNsr - jTr * M_VJsr / M_VNsr;
		subSysCaRHS[4] 	= bSS * ( jRel * M_VJsr / M_VSs - jXFer * M_VMyo / M_VSs - iCa * M_ACap / ( 2 * M_VSs * M_F ) );
		subSysCaRHS[5]	= bJSR * ( jTr - jRel );

		subSysCaRHS[6] 	= beta * c1 + M_omega * cCa0 - ( 4 * alpha + gamma ) * c0;
		subSysCaRHS[7] 	= 4 * alpha * c0 + 2 * beta * c2 + cCa1 * M_omega / M_b - ( beta + 3 * alpha + gamma * M_a ) * c1;
		subSysCaRHS[8] 	= 3 * alpha * c1  + 3 * beta * c3 + cCa2 * M_omega * ( M_b * M_b ) - ( 2 * beta + 2 * alpha + gamma * ( M_a * M_a ) ) * c2;
		subSysCaRHS[9] 	= 2 * alpha * c2 + 4 * beta * c4 + cCa3 * M_omega / std::pow(M_b, 3) - ( 3 * beta + alpha + gamma * std::pow(M_a, 3) ) * c3;
		subSysCaRHS[10] = alpha *c3 + M_g * o + cCa4 * M_omega / std::pow(M_b, 4) - ( 4 * beta + M_f + gamma * std::pow(M_a, 4) ) * c4;
		subSysCaRHS[11] = M_f * c4 - M_g * o;
		subSysCaRHS[12] = betaprime * cCa1 + gamma * c0 - ( 4 * alphaprime + M_omega ) * cCa0;
		subSysCaRHS[13] = 4 * alphaprime * cCa0 + 2 * betaprime * cCa2 + gamma * M_a * c1 - ( betaprime + 3 * alphaprime + M_omega / M_b ) * cCa1;
		subSysCaRHS[14] = 3 * alphaprime * cCa1 +  3 * betaprime * cCa3 + gamma * ( M_a * M_a ) * c2 - ( 2* betaprime + 2 * alphaprime + M_omega / ( M_b * M_b ) ) * cCa2;
		subSysCaRHS[15] = 2 * alphaprime * cCa2 + 4 * betaprime * cCa4 + gamma * std::pow(M_a, 3) * c3  - ( 3 * betaprime + alphaprime + M_omega / std::pow(M_b, 3) ) * cCa3;
		subSysCaRHS[16] = alphaprime * cCa3 + M_gprime * oCa + gamma * std::pow(M_a, 4) * c4 - ( 4 * betaprime + M_fprime + M_omega / std::pow(M_b, 4) ) * cCa4;
		subSysCaRHS[17] = M_fprime * cCa4 - M_gprime * oCa;
		subSysCaRHS[18] = ( yinf - y ) / tauY ;

//		subSysCaRHS[6] 	= - ( 4 * alpha + gamma );
//		subSysCaRHS[7] 	= - ( beta + 3 * alpha + gamma * M_a );
//		subSysCaRHS[8] 	= - ( 2 * beta + 2 * alpha + gamma * ( M_a * M_a ) );
//		subSysCaRHS[9] 	= - ( 3 * beta + alpha + gamma * std::pow(M_a, 3) );
//		subSysCaRHS[10] = - ( 4 * beta + M_f + gamma * std::pow(M_a, 4) );
//		subSysCaRHS[11] = - M_g ;
//		subSysCaRHS[12] = - ( 4 * alphaprime + M_omega );
//		subSysCaRHS[13] = - ( betaprime + 3 * alphaprime + M_omega / M_b );
//		subSysCaRHS[14] = - ( 2* betaprime + 2 * alphaprime + M_omega / ( M_b * M_b ) );
//		subSysCaRHS[15] = - ( 3 * betaprime + alphaprime + M_omega / std::pow(M_b, 3) );
//		subSysCaRHS[16] = - ( 4 * betaprime + M_fprime + M_omega / std::pow(M_b, 4) );
//		subSysCaRHS[17] = - M_gprime;
//		subSysCaRHS[18] = - 1.0 / tauY ;


		return subSysCaRHS;
	}


	//! Ionic Currents (Luo and Rudy)

	// Fast Na+ Current INa
	std::vector<Real> IonicJafriRiceWinslow::fastINa( const std::vector<Real>& v )
	{
		std::vector<Real> fastNa(7);
//		std::vector<Real> fastNa(4);

		Real V   ( v[0] );
		Real m   ( v[1] );
		Real h   ( v[2] );
		Real j   ( v[3] );
		Real cNa ( v[5] );

		Real potNa = ( M_R * M_T / M_F ) * std::log( M_NaO / cNa );

		fastNa[0] = M_GNa * std::pow(m, 3) * h * j * ( V - potNa );

		Real alpha_m = 0.32 * ( V + 47.13 ) / ( 1.0 - std::exp( - 0.1 * ( V + 47.13 ) ) );
		Real beta_m  = 0.08 * std::exp(- V / 11.0 );

		Real alpha_h ( 0 );
		Real alpha_j ( 0 );
		Real beta_h  ( 0 );
		Real beta_j  ( 0 );

		if (V >= -40)
		{
			alpha_h = 0.0;
			alpha_j = 0.0;
			beta_h  = 1.0 / ( 0.13 * ( 1.0 + std::exp( ( V + 10.66 ) / -11.1 ) ) );
			beta_j  = 0.3 * std::exp( -2.535e-7 * V ) / ( 1.0 + std::exp( -0.1 * ( V + 32.0 ) ) );
		}
		else
		{
			alpha_h = 0.135 * std::exp( ( 80 + V ) / -6.8 );
			alpha_j = ( -127140 * std::exp( 0.2444 * V ) - 3.474e-5 * std::exp( -0.04391 * V ) ) * ( V + 37.78 ) / ( 1.0 + std::exp ( 0.311 * ( V + 79.23 ) ) );
			beta_h  = 3.56 * std::exp( 0.079 * V ) + 3.1e5 * std::exp( 0.35 * V ) ;
			beta_j  = 0.1212 * std::exp( -0.01052 * V ) / ( 1.0 + std::exp( -0.1378 * ( V + 40.14 ) ) );
		}

		fastNa[1] = alpha_m;
		fastNa[2] = beta_m;
		fastNa[3] = alpha_h;
		fastNa[4] = beta_h;
		fastNa[5] = alpha_j;
		fastNa[6] = beta_j;

//		fastNa[1] = - ( alpha_m + beta_m );
//		fastNa[2] = - ( alpha_h + beta_h );
//		fastNa[3] = - ( alpha_j + beta_j );


		return fastNa;
	}

	// Time dependent K+ currrent IK
	std::vector<Real> IonicJafriRiceWinslow::timeDIK( const std::vector<Real>& v )
	{
		std::vector<Real> timeDK(3);
//		std::vector<Real> timeDK(2);

		Real V   ( v[0] );
		Real x   ( v[4] );
		Real cNa ( v[5] );
		Real cKi ( v[6] );
		Real cKo ( v[7] );


		Real potK     = ( M_R * M_T / M_F ) * std::log( ( cKo + M_PNaK * M_NaO ) / ( cKi + M_PNaK * cNa ) );
		Real maxCondK = 0.1128e-2 * sqrt( cKo / 5.4 );
		Real xi       = 1.0 / ( 1.0 + std::exp( ( V - 56.26 ) / 32.1 ) );

		timeDK[0] = maxCondK * xi * x * x * ( V - potK );
		timeDK[1] = 7.19e-5 * ( V + 30.0 ) / ( 1.0 - std::exp( -0.148 * ( V + 30.0 ) ) );
		timeDK[2] = 1.31e-4 * ( V + 30.0 ) / ( -1.0 + std::exp( 0.0687 * ( V + 30.0 ) ) );
//		Real alpha_x   = 7.19e-5 * ( V + 30.0 ) / ( 1.0 - std::exp( -0.148 * ( V + 30.0 ) ) );
//		Real beta_x    = 1.31e-4 * ( V + 30.0 ) / ( -1.0 + std::exp( 0.0687 * ( V + 30.0 ) ) );
//		timeDK[1] = - ( alpha_x + beta_x );

		return timeDK;
	}

	// Time-independent K+ Current IK1
	Real IonicJafriRiceWinslow::timeIIK1( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cKi ( v[6] );
		Real cKo ( v[7] );

		Real potK1      = ( M_R * M_T / M_F ) * std::log( cKo / cKi );
		Real maxCondK1  = 0.75e-2 * sqrt( cKo / 5.4 );
		Real alphaK1    = 1.02 / ( 1.0 + std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) );
		Real betaK1     = ( 0.4912 * std::exp( 0.08032 * ( V - potK1 + 5.476 ) ) + std::exp( 0.06175 * ( V - potK1 - 594.31 ) ) ) / ( 1.0 + std::exp( -0.5143 * ( V - potK1 + 4.753 ) ) );
		Real constK1inf = alphaK1 / ( alphaK1 + betaK1 );

		return  maxCondK1 * constK1inf * ( V - potK1 );
	}

	// Plateau K+ current IKp
	Real IonicJafriRiceWinslow::plaIKp( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cKi ( v[6] );
		Real cKo ( v[7] );

		Real potKp   = M_R * M_T / M_F * std::log( cKo / cKi );
		Real constKp = 1 / ( 1 + std::exp( ( 7.488 - V ) / 5.98 ) );

		return M_GKp * constKp *( V - potKp );
	}

	// Na+/Ca2+ exchanger current INaCa
	Real IonicJafriRiceWinslow::exINaCa( const std::vector<Real>& v )
	{

		Real V   ( v[0] );
		Real cNa ( v[5] );
		Real cCa ( v[8] );

		 return M_kNaCa * ( 1.0 / ( std::pow(M_KmNa, 3) + std::pow(M_NaO, 3) ) ) * ( 1.0 / ( M_KmCa + M_CaO ) )
				* ( 1.0 / ( 1.0 + M_kSat * std::exp( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) ) )
				* ( std::exp( M_eta * ( V * M_F ) / ( M_R * M_T ) ) * std::pow(cNa, 3) * M_CaO -
						std::exp( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) * std::pow(M_NaO, 3) * cCa );
	}

	// Na+/K+ pump INaK
	Real IonicJafriRiceWinslow::pumpINaK( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cNa ( v[5] );
		Real cKo ( v[7] );

		Real sigma = ( 1.0 / 7.0 ) * ( std::exp( M_NaO / 67.3 ) - 1.0 );
		Real fNak  = 1.0 / ( 1.0 + 0.1245 * std::exp( -0.1 * ( V * M_F ) / ( M_R * M_T ) ) + 0.0365 * sigma * std::exp( -( V * M_F ) / ( M_R * M_T ) ) );

		return M_INaK * fNak * ( 1.0 / ( 1.0 + std::pow( M_KmNai / cNa, 1.5) ) ) * ( cKo / ( cKo + M_KmK0 ) );

	}

	// Nonspecific Ca2+ activated current InsCa
	std::vector<Real> IonicJafriRiceWinslow::noSpecInsCa( const std::vector<Real>& v )
	{
		std::vector<Real> iNsCa(3);

		Real V   ( v[0] );
		Real cNa ( v[5] );
		Real cKi ( v[6] );
		Real cKo ( v[7] );
		Real cCa ( v[8] );

		Real Vns      = V - ( M_R * M_T / M_F ) * std::log( ( cKo + M_NaO ) / ( cKi + cNa ) );

		Real maxInsNa = M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) *
							( ( 0.75 * cNa * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * M_NaO ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) );

		Real maxInsK  = M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) *
							( ( 0.75 * cKi * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * cKo ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1 ) );

		iNsCa[0] = maxInsNa * ( 1.0 / ( 1.0 + std::pow(M_KmNsCa / cCa, 3) ) );
		iNsCa[1] = maxInsK * ( 1.0 / ( 1.0 + std::pow(M_KmNsCa / cCa, 3) ) );
		iNsCa[2] = iNsCa[0] + iNsCa[1];

		return iNsCa;
	}

	// Sarcolemmal Ca2+ pump current IpCa
	Real IonicJafriRiceWinslow::pumpIpCa( const std::vector<Real>& v )
	{
		return M_IpCa * v[8] / ( M_KmPCa + v[8] );
	}

	// Ca2+ background current ICab
	Real IonicJafriRiceWinslow::backICab( const std::vector<Real>& v )
	{
		Real V   ( v[0] );
		Real cCa ( v[8] );

		Real potCaN = ( ( M_R * M_T ) / ( 2 * M_F ) ) * std::log( M_CaO / cCa );

		return M_GCab * ( V - potCaN );
	}

	// Na+ background current INab
	Real IonicJafriRiceWinslow::backINab( const std::vector<Real>& v)
	{
		Real V   ( v[0] );
		Real cNa ( v[5] );

		Real potNaN = M_R * M_T / M_F * std::log( M_NaO / cNa );

		return M_GNab * ( V - potNaN );
	}


	//! RyR Channel States (Keizer and Levine)

	std::vector<Real> IonicJafriRiceWinslow::computeLocalChannelRyrRhs( const std::vector<Real>& v )
	{
		std::vector<Real> channelRyrRHS(4);

		Real cCaSS   ( v[10] );
		Real fracPC1 ( v[12] );
		Real fracPO1 ( v[13] );
		Real fracPO2 ( v[14] );
		Real fracPC2 ( v[15] );


		channelRyrRHS[0] = -M_kap * std::pow(cCaSS, M_n) * fracPC1 + M_kan * fracPO1;
		channelRyrRHS[1] = M_kap * std::pow(cCaSS, M_n) * fracPC1 - M_kan * fracPO1 -
                			M_kbp * std::pow(cCaSS, M_m) * fracPO1 + M_kbn * fracPO2 - M_kcp * fracPO1 + M_kcn * fracPC2; 
		channelRyrRHS[2] = M_kbp * std::pow(cCaSS, M_m) * fracPO1 - M_kbn * fracPO2; 
		channelRyrRHS[3] = M_kcp * fracPO1 - M_kcn * fracPC2;
//		channelRyrRHS[0]  = - ( M_kap * std::pow(cCaSS, M_n) );
//		channelRyrRHS[1]  = - ( M_kbp * std::pow(cCaSS, M_m) + M_kcp + M_kan );
//		channelRyrRHS[2]  = - M_kbn;
//		channelRyrRHS[3]  = - M_kcn;

		return channelRyrRHS;
	}

	Real IonicJafriRiceWinslow::computeNewtonNa( const std::vector<Real>& v, const Real& dt, const int& nitermax )
	{
		Real V   ( v[0] );
		Real m   ( v[1] );
		Real h   ( v[2] );
		Real j   ( v[3] );
		Real cNa ( v[5] );
		Real cKi ( v[6] );
		Real cKo ( v[7] );
		Real cCa ( v[8] );


		Real constINa   = M_GNa * std::pow(m, 3) * h * j;
		Real constInsCa = 1.0 / ( 1.0 + std::pow(M_KmNsCa / cCa, 3) );
		Real constINaCa = M_kNaCa * ( 1.0 / ( std::pow(M_KmNa, 3) + std::pow(M_NaO, 3) ) ) * ( 1.0 / ( M_KmCa + M_CaO ) )
								* ( 1.0 / ( 1.0 + M_kSat * std::exp( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) ) );
		Real sigma      = ( 1.0 / 7.0 ) * ( std::exp( M_NaO / 67.3 ) - 1.0 );
		Real fNak       = 1.0 / ( 1.0 + 0.1245 * std::exp( -0.1 * ( V * M_F ) / ( M_R * M_T ) ) +
								0.0365 * sigma * std::exp( -( V * M_F ) / ( M_R * M_T ) ) );
		Real constINaK  = M_INaK * fNak * ( cKo / ( cKo + M_KmK0 ) );

		Real cNa0      (cNa);

		Real potNa     ( 0 );
		Real INa       ( 0 );
		Real dINa      ( 0 );
		Real INab      ( 0 );
		Real dINab     ( 0 );
		Real Vns       ( 0 );
		Real maxInsNa  ( 0 );
		Real InsNa     ( 0 );
		Real dInsNa    ( 0 );
		Real dmaxInsNa ( 0 );
		Real INaCa     ( 0 );
		Real dINaCa    ( 0 );
		Real INaK      ( 0 );
		Real dINaK     ( 0 );

		Real newtF  ( 0 );
		Real newtdF ( 0 );

		for ( int i(0); i < nitermax; ++i )
		{
			potNa     = ( M_R * M_T / M_F ) * std::log( M_NaO / cNa );
			INa       = constINa * ( V - potNa );
			dINa      = constINa * M_R * M_T / ( cNa * M_F );
			INab      = M_GNab * ( V - potNa );
			dINab     = M_GNab * M_R * M_T / ( cNa * M_F );

			Vns       = V - ( M_R * M_T / M_F ) * std::log( ( cKo + M_NaO ) / ( cKi + cNa ) );
			maxInsNa  = M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) *
							( ( 0.75 * cNa * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * M_NaO ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) );
			InsNa     = maxInsNa * constInsCa;
			dmaxInsNa = ( M_PnsCa * M_F / ( cNa +cKi ) ) * ( ( 0.75 * cNa * std::exp( ( Vns * M_F ) / ( M_R * M_T ) )
							- 0.75 * M_NaO ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) )
							+ M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) * (
							0.75 * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) * ( 1.0 + cNa / ( cNa + cKi ) ) * ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 )
							- std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) * ( 1.0 / ( cNa + cKi ) ) *
							( ( 0.75 * cNa * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * M_NaO ) ) )
							/ ( ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) * ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) );
			dInsNa    = dmaxInsNa * constInsCa;

			INaCa     = constINaCa * ( std::exp( M_eta * ( V * M_F ) / ( M_R * M_T ) ) * std::pow(cNa, 3) * M_CaO -
								std::exp( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) * std::pow(M_NaO, 3) * cCa );
			dINaCa    = constINaCa * 3 * cNa * cNa * M_CaO * std::exp( M_eta * ( V * M_F ) / ( M_R * M_T ) );

			INaK      = constINaK * ( 1.0 / ( 1.0 + std::pow( M_KmNai / cNa, 1.5) ) );
			dINaK     = constINaK * ( 3 / 2 * ( std::pow( M_KmNai, 1.5 ) / ( std::pow(cNa, 2.5) ) ) ) / ( ( 1.0 + std::pow( M_KmNai / cNa, 1.5) ) * ( 1.0 + std::pow( M_KmNai / cNa, 1.5) ) );

			newtF     = cNa - cNa0 + dt * ( INa + INab + InsNa + 3 * INaCa + 3 * INaK ) * M_ACap / ( M_VMyo * M_F );
			newtdF    = 1.0 + dt * ( dINa + dINab + dInsNa + 3 * dINaCa + 3 * dINaK ) * M_ACap / ( M_VMyo * M_F );

			cNa = cNa - newtF / newtdF;
		}

		return cNa;
	}

	Real IonicJafriRiceWinslow::computeNewtonKi( const std::vector<Real>& v, const Real& dt, const int& nitermax )
	{
		Real V   ( v[0] );
		Real x   ( v[4] );
		Real cNa ( v[5] );
		Real cKi ( v[6] );
		Real cKo ( v[7] );
		Real cCa ( v[8] );
		Real o   ( v[21] );
		Real oCa ( v[27] );
		Real y   ( v[28] );


		Real maxCondK  = 0.001128 * sqrt( cKo / 5.4 );
		Real xi        = 1.0 / ( 1.0 + std::exp( ( V - 56.26 ) / 32.1 ) );
		Real constIK   =  maxCondK * xi * x * x;
		Real maxCondK1 = 7.5e-3 * sqrt( cKo / 5.4 );
        Real constKp   = 1 / ( 1 + std::exp( ( 7.488 - V ) / 5.98 ) );
        Real sigma     = ( 1.0 / 7.0 ) * ( std::exp( M_NaO / 67.3 ) - 1.0 );
		Real fNak      = 1.0 / ( 1.0 + 0.1245 * std::exp( -0.1 * ( V * M_F ) / ( M_R * M_T ) ) + 0.0365 * sigma * std::exp( -( V * M_F ) / ( M_R * M_T ) ) );
        Real iCaMax    = M_PCa * 4 * ( V * M_F * M_F ) / ( M_R * M_T ) * ( 0.001 * std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 0.341 * M_CaO )
                            / ( std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );
		Real constInsCa = 1.0 / ( 1.0 + std::pow(M_KmNsCa / cCa, 3) );
		Real pK        = M_PK / ( 1.0 + iCaMax / M_ICaHalf );
		Real constICaK = pK * y * ( o + oCa ) * (V * M_F * M_F ) / ( M_R * M_T );

        Real cKi0   (cKi);

		Real potK        ( 0 );
		Real IK          ( 0 );
		Real dIK         ( 0 );
		Real potK1       ( 0 );
		Real alphaK1     ( 0 );
		Real betaK1      ( 0 );
		Real constK1inf  ( 0 );
		Real IK1         ( 0 );
		Real dalphaK1    ( 0 );
		Real numebetaK1  ( 0 );
		Real denobetaK1  ( 0 );
		Real dbetaK1     ( 0 );
		Real dconstK1inf ( 0 );
		Real dIK1        ( 0 );
		Real Vns		 ( 0 );
		Real maxInsK     ( 0 );
        Real InsK        ( 0 );
        Real dInsK       ( 0 );
        Real dmaxInsK    ( 0 );
        Real IKp         ( 0 );
        Real dIKp        ( 0 );
        Real INaK        = M_INaK * fNak * ( cKo / ( cKo + M_KmK0 ) ) * ( 1.0 / ( 1.0 + std::pow( M_KmNai / cNa, 1.5) ) );
        Real dINaK       ( 0 );
        Real ICaK        ( 0 );
        Real dICaK       ( 0 );

		Real newtF  ( 0 );
		Real newtdF ( 0 );

		for ( int i(0); i < nitermax; ++i )
		{
			potK        = ( M_R * M_T / M_F ) * std::log( ( cKo + M_PNaK * M_NaO ) / ( cKi + M_PNaK * cNa ) );
			IK          = constIK * ( V - potK );
			dIK         = constIK * ( M_R * M_T / M_F ) * 1.0 / ( cKi + M_PNaK * cNa );

			potK1       = ( M_R * M_T / M_F ) * std::log( cKo / cKi );
			alphaK1     = 1.02 / ( 1.0 + std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) );
			betaK1      = ( 0.4912 * std::exp( 0.08032 * ( V - potK1 + 5.476 ) ) + std::exp( 0.06175 * ( V - potK1 - 594.31 ) ) ) / ( 1.0 + std::exp( -0.5143 * (V - potK1 + 4.753 ) ) );
			constK1inf  = alphaK1 / ( alphaK1 + betaK1 );
			IK1         = maxCondK1 * constK1inf * ( V - potK1 );
            dalphaK1    = - 1.02 / ( ( 1.0 + std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) ) * ( 1.0 + std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) ) )
                            * std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) * 0.2385 * ( M_R * M_T / M_F ) * 1.0 / cKi;
            numebetaK1  = 0.4912 * std::exp( 0.08032 * ( V - potK1 + 5.476 ) ) + std::exp( 0.06175 * ( V - potK1 - 594.31 ) );
            denobetaK1  =  1.0 + std::exp( -0.5143 * ( V - potK1 + 4.753 ) );
            dbetaK1     = ( denobetaK1 * ( 0.4912 * std::exp( 0.08032 * ( V - potK1 + 5.476 ) ) * 0.08032 * ( M_R * M_T / M_F ) * 1.0 / cKi
                            + std::exp( 0.06175 * ( V - potK1 - 594.31 ) ) * 0.06175 * ( M_R * M_T / M_F ) * 1.0 / cKi )
                            + numebetaK1 * std::exp( -0.5143 * ( V - potK1 + 4.753 ) ) * 0.5413 * 1.0 / cKi )
                            / ( denobetaK1 * denobetaK1 );
            dconstK1inf = ( dalphaK1 * ( alphaK1 + betaK1 ) - ( dalphaK1 + dbetaK1 ) * alphaK1 ) / ( ( alphaK1 + betaK1 ) * ( alphaK1 + betaK1 ) );
            dIK1        = maxCondK1 * dconstK1inf * ( V - potK1 ) + maxCondK1 * constK1inf * ( M_R * M_T / M_F ) * 1.0 / cKi;

            Vns         = V - ( M_R * M_T / M_F ) * std::log( ( cKo + M_NaO ) / ( cKi + cNa ) );
			maxInsK     = M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) *
									( ( 0.75 * cKi * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * cKo ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) );
            InsK        = maxInsK * constInsCa;
            dmaxInsK    = ( M_PnsCa * M_F / ( cNa +cKi ) ) * ( ( 0.75 * cKi * std::exp( ( Vns * M_F ) / ( M_R * M_T ) )
									- 0.75 * cKo ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) )
									+ M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) * (
									0.75 * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) * ( 1.0 + cKi / ( cNa + cKi ) ) * ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 )
									- std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) * ( 1.0 / ( cNa + cKi ) ) *
									( ( 0.75 * cKi * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * cKo ) ) )
									/ ( ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) * ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) );
            dInsK       = dmaxInsK * constInsCa;

            IKp         = M_GKp * constKp *( V - potK1 );
            dIKp        = M_GKp * constKp * ( M_R * M_T / M_F ) * 1.0 / cKi;

            ICaK        = constICaK * ( cKi * std::exp( ( V * M_F ) / ( M_R * M_T ) ) - cKo )
                            / ( std::exp( ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );
            dICaK       = constICaK * std::exp( ( V * M_F ) / ( M_R * M_T ) ) / ( std::exp( ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );

			newtF = cKi - cKi0 + dt * ( IK + IK1 + InsK + IKp - 2 * INaK + ICaK ) * M_ACap / ( M_VMyo * M_F );
			newtdF = 1.0 + dt * ( dIK + dIK1 + dInsK + dIKp - 2 * dINaK + dICaK ) * M_ACap / ( M_VMyo * M_F );

            cKi = cKi - newtF / newtdF;
        }

		return cKi;
    }

    Real IonicJafriRiceWinslow::computeNewtonKo( const std::vector<Real>& v, const Real& dt, const int& nitermax )
	{
		Real V   ( v[0] );
		Real x   ( v[4] );
		Real cNa ( v[5] );
		Real cKi ( v[6] );
		Real cKo ( v[7] );
		Real cCa ( v[8] );
		Real o   ( v[21] );
		Real oCa ( v[27] );
		Real y   ( v[28] );


        Real xi        = 1.0 / ( 1.0 + std::exp( ( V - 56.26 ) / 32.1 ) );
		Real constIK   = xi * x * x;
        Real constKp   = 1 / ( 1 + std::exp( ( 7.488 - V ) / 5.98 ) );
       	Real sigma     = ( 1.0 / 7.0 ) * ( std::exp( M_NaO / 67.3 ) - 1.0 );
		Real fNak      = 1.0 / ( 1.0 + 0.1245 * std::exp( -0.1 * ( V * M_F ) / ( M_R * M_T ) ) + 0.0365 * sigma * std::exp( -( V * M_F ) / ( M_R * M_T ) ) );
        Real constINaK = M_INaK * fNak * ( 1.0 / ( 1.0 + std::pow( M_KmNai / cNa, 1.5) ) );
        Real iCaMax    = M_PCa * 4 * ( V * M_F * M_F ) / ( M_R * M_T ) * ( 0.001 * std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 0.341 * M_CaO )
                            / ( std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );
		Real iCa       = y * ( o + oCa ) * iCaMax;
		Real constInsCa = 1.0 / ( 1.0 + std::pow(M_KmNsCa / cCa, 3) );
		Real pK        = M_PK / ( 1.0 + iCaMax / M_ICaHalf );
		Real constICaK = pK * y * ( o + oCa ) * ( V * M_F * M_F ) / ( M_R * M_T );

        Real cKo0        (cKo);

        Real maxCondK    ( 0 );
        Real dmaxCondK   ( 0 );
		Real potK        ( 0 );
		Real IK          ( 0 );
		Real dIK         ( 0 );
		Real maxCondK1   ( 0 );
		Real dmaxCondK1  ( 0 );
		Real potK1       ( 0 );
		Real alphaK1     ( 0 );
		Real betaK1      ( 0 );
		Real constK1inf  ( 0 );
		Real IK1         ( 0 );
		Real dalphaK1    ( 0 );
		Real numebetaK1  ( 0 );
		Real denobetaK1  ( 0 );
		Real dbetaK1     ( 0 );
		Real dconstK1inf ( 0 );
		Real dIK1        ( 0 );
		Real Vns         ( 0 );
		Real maxInsK     ( 0 );
        Real InsK        ( 0 );
        Real dmaxInsK    ( 0 );
        Real dInsK       ( 0 );
        Real IKp         ( 0 );
        Real dIKp        ( 0 );
        Real INaK        ( 0 );
        Real dINaK       ( 0 );
        Real ICaK        ( 0 );
        Real dICaK       ( 0 );

		Real newtF  ( 0 );
		Real newtdF ( 0 );

		for ( int i(0); i < nitermax; ++i )
		{
		    maxCondK    = 0.001128 * sqrt( cKo / 5.4 );
            dmaxCondK   = 0.001128 * sqrt( 1.0 / ( 5.4 * cKo ) ) / 2;
			potK        = ( M_R * M_T / M_F ) * std::log( ( cKo + M_PNaK * M_NaO ) / ( cKi + M_PNaK * cNa ) );
			IK          = maxCondK * constIK * ( V - potK );
			dIK         = dmaxCondK * constIK * ( V - potK ) - maxCondK * constIK * ( M_R * M_T / M_F ) * 1.0 / ( cKo + M_PNaK * M_NaO );

            maxCondK1   = 7.5e-3 * sqrt( cKo / 5.4 );
            dmaxCondK1  = 7.5e-3 * sqrt( 1.0 / ( 5.4 * cKo ) ) / 2;
			potK1       = ( M_R * M_T / M_F ) * std::log( cKo / cKi );
			alphaK1     = 1.02 / ( 1.0 + std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) );
			betaK1      = ( 0.4912 * std::exp( 0.08032 * ( V - potK1 + 5.476 ) ) + std::exp( 0.06175 * ( V - potK1 - 594.31 ) ) ) / ( 1.0 + std::exp( -0.5143 * (V - potK1 + 4.753 ) ) );
			constK1inf  = alphaK1 / ( alphaK1 + betaK1 );
			IK1         = maxCondK1 * constK1inf * ( V - potK1 );
            dalphaK1    = 1.02 / ( ( 1.0 + std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) ) * ( 1.0 + std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) ) )
                            * std::exp( 0.2385 * ( V - potK1 - 59.215 ) ) * 0.2385 * ( M_R * M_T / M_F ) * 1.0 / cKo;
            numebetaK1  = 0.4912 * std::exp( 0.08032 * ( V - potK1 + 5.476 ) ) + std::exp( 0.06175 * ( V - potK1 - 594.31 ) );
            denobetaK1  =  1.0 + std::exp( -0.5143 * ( V - potK1 + 4.753 ) );
            dbetaK1     = - ( denobetaK1 * ( 0.4912 * std::exp( 0.08032 * ( V - potK1 + 5.476 ) ) * 0.08032 * ( M_R * M_T / M_F ) * 1.0 / cKo
                            + std::exp( 0.06175 * ( V - potK1 - 594.31 ) ) * 0.06175 * ( M_R * M_T / M_F ) * 1.0 / cKo )
                            + numebetaK1 * std::exp( -0.5143 * ( V - potK1 + 4.753 ) ) * 0.5413 * 1.0 / cKo )
                            / ( denobetaK1 * denobetaK1 );
            dconstK1inf = ( dalphaK1 * ( alphaK1 + betaK1 ) - ( dalphaK1 + dbetaK1 ) * alphaK1 ) / ( ( alphaK1 + betaK1 ) * ( alphaK1 + betaK1 ) );
            dIK1        = maxCondK1 * dconstK1inf * ( V - potK1 )
                            + dmaxCondK1 * constK1inf * ( V - potK1 )
                            - maxCondK1 * constK1inf * ( M_R * M_T / M_F ) * 1.0 / cKo;

            Vns         = V - ( M_R * M_T / M_F ) * std::log( ( cKo + M_NaO ) / ( cKi + cNa ) );
			maxInsK     = M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) *
									( ( 0.75 * cKi * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * cKo ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) );
            InsK        = maxInsK * constInsCa;
            dmaxInsK    = - ( M_PnsCa * M_F / ( M_NaO +cKo ) ) * ( ( 0.75 * cKi * std::exp( ( Vns * M_F ) / ( M_R * M_T ) )
									- 0.75 * cKo ) / ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) )
									+ M_PnsCa * ( ( Vns * M_F * M_F ) / ( M_R * M_T ) ) * (
									- 0.75 * ( 1.0 + std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) * ( cKi / ( M_NaO + cKo ) ) )
                                    * ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) +
                                    ( 0.75 * cKi * std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 0.75 * cKo ) * std::exp( ( Vns * M_F ) / ( M_R * M_T ) )
                                    * 1.0 / ( M_NaO + cKo ) ) / ( ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) *
                                    ( std::exp( ( Vns * M_F ) / ( M_R * M_T ) ) - 1.0 ) );
            dInsK       = dmaxInsK * constInsCa;

            IKp         = M_GKp * constKp *( V - potK1 );
            dIKp        = - M_GKp * constKp * ( M_R * M_T / M_F ) * 1.0 / cKo;

            INaK        = constINaK * cKo / ( cKo + M_KmK0 );
            dINaK       = constINaK * M_KmK0 / ( ( cKo + M_KmK0 ) * ( cKo + M_KmK0 ) );

            ICaK        =  constICaK * ( cKi * std::exp( ( V * M_F ) / ( M_R * M_T ) ) - cKo )
                            / ( std::exp( ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );
            dICaK       = - constICaK * 1.0 / ( std::exp( ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );

			newtF = cKo - cKo0 - dt * ( IK + IK1 + InsK + IKp - 2 * INaK + ICaK ) * M_ACap / ( M_VMyo * M_F );
			newtdF = 1.0 - dt * ( dIK + dIK1 + dInsK + dIKp - 2 * dINaK + dICaK ) * M_ACap / ( M_VMyo * M_F );

            cKo = cKo - newtF / newtdF;
        }

		return cKo;
    }

    Real IonicJafriRiceWinslow::computeNewtonCai( const std::vector<Real>& v, const Real& dt, const int& nitermax )
	{
		Real V 	      ( v[0] );
		Real cNa      ( v[5] );
		Real cCa      ( v[8] );
		Real cCaNSR   ( v[9] );
		Real cCaSS    ( v[10] );
		Real cLTrpnCa ( v[29] );
		Real cHTrpnCa ( v[30] );

        Real djLeak     = - M_v2;
        Real djXFer     = - 1.0 / M_tauXFer;
        Real constINaCa = M_kNaCa * ( 1.0 / ( std::pow(M_KmNa, 3) + std::pow(M_NaO, 3) ) ) * ( 1.0 / ( M_KmCa + M_CaO ) )
                            * ( 1.0 / ( 1.0 + M_kSat * std::exp( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) ) );

		Real cCai0  ( cCa );

		Real jLeak  ( 0 );
		Real jXFer  ( 0 );
		Real jUp    ( 0 );
		Real djUp   ( 0 );
		Real jTRPN  ( 0 );
		Real djTRPN ( 0 );
		Real bi     ( 0 );
		Real dbi    ( 0 );
        Real potCaN ( 0 );
        Real ICab   ( 0 );
        Real dICab  ( 0 );
        Real INaCa  ( 0 );
        Real dINaCa ( 0 );
        Real IpCa   ( 0 );
        Real dIpCa  ( 0 );

		Real newtF  ( 0 );
		Real newtdF ( 0 );

		for ( int i(0); i < nitermax; ++i )
		{
			jLeak  = M_v2 * ( cCaNSR - cCa );

            jXFer  = ( cCaSS - cCa ) / M_tauXFer;

            jUp    =  M_v3 * ( ( cCa * cCa ) / ( M_KmUp * M_KmUp + cCa * cCa ) );
            djUp   =  M_v3 * 2 * ( cCa * M_KmUp * M_KmUp / ( ( M_KmUp * M_KmUp + cCa * cCa ) * ( M_KmUp * M_KmUp + cCa * cCa ) ) );

            jTRPN  = M_kPHtrpn * cCa * ( M_HtrpnTot - cHTrpnCa ) - M_kNHtrpn * cHTrpnCa
						+ M_kPLtrpn * cCa * ( M_LtrpnTot - cLTrpnCa ) - M_kNLtrpn * cLTrpnCa;
            djTRPN = M_kPHtrpn * ( M_HtrpnTot - cHTrpnCa ) + M_kPLtrpn * ( M_LtrpnTot - cLTrpnCa );

			bi     = 1.0 / ( 1.0 + M_CmdnTot * M_KmCmdn / ( ( M_KmCmdn + cCa ) * ( M_KmCmdn + cCa ) ) );
			dbi    = 2.0 * std::pow ( 1.0 + M_CmdnTot * M_KmCmdn / ( ( M_KmCmdn + cCa ) * ( M_KmCmdn + cCa ) ), -2) * M_CmdnTot * M_KmCmdn / std::pow(M_KmCmdn + cCa, 3);

            potCaN = ( ( M_R * M_T ) / ( 2 * M_F ) ) * std::log( M_CaO / cCa );
            ICab   = M_GCab * ( V - potCaN );
            dICab  = M_GCab * ( ( M_R * M_T ) / ( 2 * M_F ) ) * 1.0 / cCa;

            INaCa  = constINaCa * ( std::exp( M_eta * ( V * M_F ) / ( M_R * M_T ) ) * std::pow(cNa, 3) * M_CaO -
						std::exp( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) * std::pow(M_NaO, 3) * cCa );
            dINaCa = constINaCa * ( - std::exp( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) * std::pow(M_NaO, 3) );

            IpCa   = M_IpCa * cCa / ( M_KmPCa + cCa );
            dIpCa  = M_IpCa * M_KmPCa / ( ( M_KmPCa + cCa ) * ( M_KmPCa + cCa ) );

			newtF = cCa - cCai0 - dt * bi *
								( jLeak + jXFer - jUp - jTRPN - ( ICab - 2 * INaCa + IpCa ) * M_ACap / ( 2 * M_VMyo * M_F ) );
			newtdF = 1.0 - dt * ( dbi * ( jLeak + jXFer - jUp - jTRPN - ( ICab - 2 * INaCa + IpCa ) * M_ACap / ( 2 * M_VMyo * M_F ) )
								+ bi * ( djLeak + djXFer - djUp - djTRPN - ( dICab - 2 * dINaCa + dIpCa ) * M_ACap / ( 2 * M_VMyo * M_F ) ) );

			cCa = cCa - newtF / newtdF;
		}

		return cCa;
	}

    Real IonicJafriRiceWinslow::computeNewtonCaNSR( const std::vector<Real>& v, const Real& dt, const int& nitermax )
    {
       	Real cCa      ( v[8] );
       	Real cCaNSR   ( v[9] );
       	Real cCaJSR   ( v[11] );

       	Real djUp   ( 0 );
        Real djLeak = M_v2;
        Real djTr   = 1.0 / M_tauTr;

        Real jUp    =  M_v3 * ( ( cCa * cCa ) / ( M_KmUp * M_KmUp + cCa * cCa ) );

      	Real cCaNSR0 ( cCaNSR );

    	Real jLeak ( 0 );
    	Real jTr   ( 0 );

    	Real newtF  ( 0 );
       	Real newtdF ( 0 );

       	for ( int i(0); i < nitermax; ++i )
       	{
       		jLeak  = M_v2 * ( cCaNSR - cCa );

       		jTr  = ( cCaNSR - cCaJSR ) / M_tauTr;

        	newtF = cCaNSR - cCaNSR0 - dt * ( ( jUp - jLeak ) * ( M_VMyo / M_VNsr ) - jTr * ( M_VJsr / M_VNsr ) );
        	newtdF = 1.0 - dt * ( ( djUp - djLeak ) * ( M_VMyo / M_VNsr ) - djTr * ( M_VJsr / M_VNsr ) );

        	cCaNSR = cCaNSR - newtF / newtdF;
        }

        return cCaNSR;
    }

    Real IonicJafriRiceWinslow::computeNewtonCaSS( const std::vector<Real>& v, const Real& dt, const int& nitermax )
	{
		Real V 	      ( v[0] );
		Real cCa      ( v[8] );
		Real cCaSS    ( v[10] );
		Real cCaJSR   ( v[11] );
		Real fracPO1  ( v[13] );
		Real fracPO2  ( v[14] );
		Real o        ( v[21] );
		Real oCa      ( v[27] );
		Real y        ( v[28] );


        Real djRel  = - M_v1 * ( fracPO1 + fracPO2 );
        Real djXFer = 1.0 / M_tauXFer;
		Real iCaMax = M_PCa * 4 * (V * M_F * M_F ) / ( M_R * M_T ) * ( 0.001 * std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 0.341 * M_CaO )
								/ ( std::exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 1.0 );
		Real iCa    = y * ( o + oCa ) * iCaMax;


		Real cCaSS0 ( cCaSS );

		Real jRel   ( 0 );
		Real jXFer  ( 0 );
		Real bSS    ( 0 );
		Real dbSS   ( 0 );

		Real newtF  ( 0 );
		Real newtdF ( 0 );

		for ( int i(0); i < nitermax; ++i )
		{
			jRel   = M_v1 * ( fracPO1 + fracPO2 ) * ( cCaJSR - cCaSS );

			jXFer  = ( cCaSS - cCa ) / M_tauXFer;

			bSS    = 1 / ( 1 + M_CmdnTot * M_KmCmdn / ( ( M_KmCmdn + cCaSS ) * ( M_KmCmdn + cCaSS ) ) );
			dbSS   = 2 * std::pow (1 + M_CmdnTot * M_KmCmdn / ( ( M_KmCmdn + cCaSS ) * ( M_KmCmdn + cCaSS ) ), -2) * M_CmdnTot * M_KmCmdn / std::pow(M_KmCmdn + cCaSS, 3);

			newtF = cCaSS - cCaSS0 - dt * bSS *
								( jRel * M_VJsr / M_VSs - jXFer * M_VMyo / M_VSs - iCa * M_ACap / ( 2 * M_VSs * M_F ) );
			newtdF = 1.0 - dt * ( dbSS * ( jRel * M_VJsr / M_VSs - jXFer * M_VMyo / M_VSs - iCa * M_ACap / ( 2 * M_VSs * M_F ) )
								+ bSS * ( djRel * M_VJsr / M_VSs - djXFer * M_VMyo / M_VSs ) );

			cCaSS = cCaSS - newtF / newtdF;
		}

		return cCaSS;
	}

    Real IonicJafriRiceWinslow::computeNewtonCaJSR( const std::vector<Real>& v, const Real& dt, const int& nitermax )
    {
    	Real cCaNSR   ( v[9] );
    	Real cCaSS    ( v[10] );
    	Real cCaJSR   ( v[11] );
    	Real fracPO1  ( v[13] );
    	Real fracPO2  ( v[14] );

        Real djRel  = M_v1 * ( fracPO1 + fracPO2 );
        Real djTr   = - 1.0 / M_tauTr;

    	Real cCaJSR0 ( cCaJSR );

   		Real jRel  ( 0 );
   		Real jTr   ( 0 );
   		Real bJSR  ( 0 );
   		Real dbJSR ( 0 );

   		Real newtF  ( 0 );
    	Real newtdF ( 0 );

    	for ( int i(0); i < nitermax; ++i )
    	{
    		jRel   = M_v1 * ( fracPO1 + fracPO2 ) * ( cCaJSR - cCaSS );

    		jTr  = ( cCaNSR - cCaJSR ) / M_tauTr;

    		bJSR    = 1 / ( 1 + M_CsqnTot * M_KmCsqn / ( ( M_KmCsqn + cCaJSR ) * ( M_KmCsqn + cCaJSR ) ) );
    		dbJSR   = 2 * std::pow (1 + M_CsqnTot * M_KmCsqn / ( ( M_KmCsqn + cCaJSR ) * ( M_KmCsqn + cCaJSR ) ), -2) * M_CsqnTot * M_KmCsqn / std::pow(M_KmCsqn + cCaJSR, 3);

    		newtF = cCaJSR - cCaJSR0 - dt * bJSR * ( jTr - jRel );
    		newtdF = 1.0 - dt * ( dbJSR * ( jTr - jRel ) + bJSR * ( djTr - djRel ) );

    		cCaJSR = cCaJSR - newtF / newtdF;
    	}

    	return cCaJSR;
    }


//    std::vector<Real> IonicJafriRiceWinslow::gateInf( const std::vector<Real>& v )
//    {
//
//		Real V ( v[0] );
//
//		std::vector<Real> gateInf (4);
//
//		Real alpha_m = 0.32 * ( V + 47.13 ) / ( 1.0 - std::exp( - 0.1 * ( V + 47.13 ) ) );
//		Real beta_m  = 0.08 * std::exp(- V / 11.0 );
//
//		Real alpha_h ( 0 );
//		Real alpha_j ( 0 );
//		Real beta_h  ( 0 );
//		Real beta_j  ( 0 );
//
//		if (V >= -40)
//		{
//			alpha_h = 0.0;
//			alpha_j = 0.0;
//			beta_h  = 1.0 / ( 0.13 * ( 1.0 + std::exp( ( V + 10.66 ) / -11.1 ) ) );
//			beta_j  = 0.3 * std::exp( -2.535e-7 * V) / ( 1.0 + std::exp( -0.1 * ( V + 32.0 ) ) );
//		}
//		else
//		{
//			alpha_h = 0.135 * std::exp( ( 80 + V ) / -6.8 );
//			alpha_j = ( -127140 * std::exp( 0.2444 * V) - 3.474e-5 * std::exp( -0.04391 * V ) ) * ( V + 37.78 ) / ( 1.0 + std::exp ( 0.311 * ( V + 79.23 ) ) );
//			beta_h  = 3.56 * std::exp( 0.079 * V ) + 3.1e5 * std::exp( 0.35 * V ) ;
//			beta_j  = 0.1212 * std::exp( -0.01052 * V ) / ( 1.0 + std::exp( -0.1378 * ( V + 40.14 ) ) );
//		}
//
//		Real alpha_x = 7.19e-5 * ( V + 30.0 ) / ( 1.0 - std::exp( -0.148 * ( V + 30.0 ) ) );
//		Real beta_x  = 1.31e-4 * ( V + 30.0 ) / ( -1.0 + std::exp( 0.0687 * ( V + 30.0 ) ) );
//
//	    gateInf[0] = - alpha_m / ( alpha_m + beta_m );
//	    gateInf[1] = - alpha_h / ( alpha_h + beta_h );
//	    gateInf[2] = - alpha_j / ( alpha_j + beta_j );
//	    gateInf[3] = - alpha_x / ( alpha_x + beta_x );
//
//	    return gateInf;
//    }
//
//
//    std::vector<Real> IonicJafriRiceWinslow::otherVarInf( const std::vector<Real>& v )
//    {
//
//		Real V 	      ( v[0] );
//		Real cCa      ( v[8] );
//		Real cCaSS    ( v[10] );
//		Real fracPC1  ( v[12] );
//		Real fracPO1  ( v[13] );
//		Real fracPO2  ( v[14] );
//		Real fracPC2  ( v[15] );
//		Real c0       ( v[16] );
//		Real c1 	  ( v[17] );
//		Real c2       ( v[18] );
//		Real c3       ( v[19] );
//		Real c4       ( v[20] );
//		Real o        ( v[21] );
//		Real cCa0     ( v[22] );
//		Real cCa1     ( v[23] );
//		Real cCa2     ( v[24] );
//		Real cCa3     ( v[25] );
//		Real cCa4     ( v[26] );
//		Real oCa      ( v[27] );
//
//		std::vector<Real> otherVarInf (19);
//
//		Real alpha      = 0.4 * std::exp( ( V + 12.0 ) / 10.0 );
//		Real beta       = 0.05 * std::exp( - ( V + 12 ) / 13.0 );
//		Real alphaprime = M_a * alpha;
//		Real betaprime  = beta / M_b;
//		Real gamma      = 0.1875 * cCaSS;
//
//
//		otherVarInf[0]  = ( M_kan * fracPO1 ) / ( M_kap * std::pow(cCaSS, M_n) );
//		otherVarInf[1]  = ( M_kap * std::pow(cCaSS, M_n) * fracPC1 + M_kbn * fracPO2 + M_kcn * fracPC2 ) / ( M_kbp * std::pow(cCaSS, M_m) + M_kcp + M_kan );
//		otherVarInf[2]  = ( M_kbp * std::pow(cCaSS, M_m) * fracPO1 ) / ( M_kbn );
//		otherVarInf[3]  = ( M_kbp * std::pow(cCaSS, M_m) * fracPO1 ) / ( M_kcn );
//		otherVarInf[4]  = ( beta * c1 + M_omega * cCa0 ) / ( 4 * alpha + gamma );
//		otherVarInf[5]  = ( 4 * alpha * c0 + 2 * beta * c2 + cCa1 * M_omega / M_b ) / ( beta + 3 * alpha + gamma * M_a );
//		otherVarInf[6]  = ( 3 * alpha * c1  + 3 * beta * c3 + cCa2 * M_omega * ( M_b * M_b ) ) / ( 2 * beta + 2 * alpha + gamma * ( M_a * M_a ) );
//		otherVarInf[7]  = ( 2 * alpha * c2 + 4 * beta * c4 + cCa3 * M_omega / std::pow(M_b, 3) ) / ( 3 * beta + alpha + gamma * std::pow(M_a, 3) );
// 		otherVarInf[8]  = ( alpha *c3 + M_g * o + cCa4 * M_omega / std::pow(M_b, 4) ) / ( 4 * beta + M_f + gamma * std::pow(M_a, 4) );
//		otherVarInf[9]  = ( M_f * c4 ) / M_g;
//		otherVarInf[10] = ( betaprime * cCa1 + gamma * c0 ) / ( 4 * alphaprime + M_omega );
//		otherVarInf[11] = ( 4 * alphaprime * cCa0 + 2 * betaprime * cCa2 + gamma * M_a * c1 ) / ( betaprime + 3 * alphaprime + M_omega / M_b );
//		otherVarInf[12] = ( 3 * alphaprime * cCa1 +  3 * betaprime * cCa3 + gamma * ( M_a * M_a ) * c2 ) / ( 2* betaprime + 2 * alphaprime + M_omega / ( M_b * M_b ) );
//		otherVarInf[13] = ( 2 * alphaprime * cCa2 + 4 * betaprime * cCa4 + gamma * std::pow(M_a, 3) * c3 ) / ( 3 * betaprime + alphaprime + M_omega / std::pow(M_b, 3) );
//		otherVarInf[14] = ( alphaprime * cCa3 + M_gprime * oCa + gamma * std::pow(M_a, 4) * c4 ) / ( 4 * betaprime + M_fprime + M_omega / std::pow(M_b, 4) );
//		otherVarInf[15] = ( M_fprime * cCa4 ) / M_gprime;
//		otherVarInf[16] = 1.0 / ( 1 + std::exp( ( V + 55 ) / 7.5 ) ) + 0.1 / ( 1 + std::exp( -V + 21 ) / 6 );
//		otherVarInf[17] = ( M_kPLtrpn * cCa * M_LtrpnTot ) / ( M_kNLtrpn + M_kPLtrpn * cCa );
//		otherVarInf[18] = ( M_kPHtrpn * cCa * M_HtrpnTot ) / ( M_kNHtrpn + M_kPHtrpn * cCa );
//
//		return otherVarInf;
//
//    }

//	std::vector<Real> IonicJafriRiceWinslow::computeYParameters( const std::vector<Real>& v)
//	{
//		std::vector<Real> intParam(2);
//
//		intParam[0] = 1 / ( 1 + std::exp( ( v[0] + 55 ) / 7.5 ) ) + 0.1 / ( 1 + std::exp( -v[0] + 21 ) / 6 );
//		intParam[1] = 20 + 600 / ( 1 + std::exp( v[0]+ 30 ) / 9.5 );
//
//		return intParam;
//	}

//	std::vector<Real> IonicJafriRiceWinslow::channelLCaMatrix( const std::vector<Real>& v )
//	{
//		V     = v[0];
//		cCaSS = v[10];
//
//		Real alpha      = 0.4 * std::exp( ( V + 12 ) / 10);
//		Real beta       = 0.05 * std::exp( -( V + 12 ) / 13);
//		Real alphaprime = M_a * alpha;
//		Real betaprime  = beta / M_b;
//		Real gamma      = 0.1875 * cCaSS;
//
//		Real lCaMatrix[12][12];
//
//		for(int i(0); i < 12; ++i)
//		{
//			for(int j(0); j < 12; ++j)
//			{
//				if ( j == i )
//				{
//					if ( i < 5 )
//					{
//						if ( i == 4 )
//							lCaMatrix[i][j] = - ( M_f + i * beta + gamma * std::pow(M_a, i) );
//						else
//							lCaMatrix[i][j] = - ( ( 4 - i ) * alpha + i * beta + gamma * std::pow(M_a, i) );
//					}
//					else if( i == 5 )
//					{
//						lCaMatrix[i][j] = - M_g;
//
//					}
//					else if ( i > 5 && i != 11 )
//					{
//						if ( i == 10 )
//							lCaMatrix[i][j] = - ( M_fprime + i * betaprime + M_omega * std::pow(1 / M_b, i) );
//						else
//							lCaMatrix[i][j] = - ( ( 4 - ( i - 6 ) ) * alphaprime + i * betaprime + M_omega * std::pow(1 / M_b, i) );
//					}
//					else
//					{
//						lCaMatrix[i][j] = - M_gprime;
//					}
//				}
//
//				else if ( j == i + 1 )
//				{
//					if ( i < 4 )
//					{
//						lCaMatrix[i][j] = i * beta;
//					}
//					else if( i == 4 )
//					{
//						lCaMatrix[i][j] = M_g;
//
//					}
//					else if ( i > 5 && i != 10 )
//					{
//						lCaMatrix[i][j] = i * betaprime;
//					}
//					else if ( i == 10 )
//					{
//						lCaMatrix[i][j] = M_gprime;
//					}
//					else
//						lCaMatrix[i][j] = 0;
//				}
//
//				else if ( j == i-1 )
//				{
//					if ( j < 4 )
//					{
//						lCaMatrix[i][j] = ( 4 - j ) * alpha;
//					}
//					else if( j == 4 )
//					{
//						lCaMatrix[i][j] = M_f;
//					}
//					else if ( i > 6 && i != 11 )
//					{
//						lCaMatrix[i][j] = ( 4 - ( j - 6 ) ) * alphaprime;
//					}
//					else if ( i == 11 )
//					{
//						lCaMatrix[i][j] = M_fprime;
//					}
//					else
//						lCaMatrix[i][j] = 0;
//				}
//
//				else if ( j == i + 6 )
//				{
//					if( i <= 4 )
//					{
//						lCaMatrix[i][j] = M_omega * std::pow(1 / M_b, i);
//					}
//					else
//						lCaMatrix[i][j] = 0;
//				}
//
//				else if ( j == i - 6 )
//				{
//					if( j <= 4 )
//					{
//						lCaMatrix[i][j] = gamma * std::pow(M_a, j);
//					}
//					else
//						lCaMatrix[i][j] = 0;
//				}
//
//				else
//					lCaMatrix[i][j] = 0;
//			}
//		}
//	}

//	std::vector<Real> IonicJafriRiceWinslow::channelRyRMatrix( const std::vector<Real>& v )
//	{
//		cCaSS = v[10];
//
//		Real lRyRMatrix[4][4];
//
//		lRyRMatrix[0][0] =  - M_kap * std::pow(cCaSS, M_n);
//		lRyRMatrix[0][1] = M_kan;
//		lRyRMatrix[0][2] = 0;
//		lRyRMatrix[0][3] = 0;
//		lRyRMatrix[1][0] = M_kap * std::pow(cCaSS, M_n);
//		lRyRMatrix[1][1] = - M_kan - M_kbp * std::pow(cCaSS, M_m) - M_kcp;
//		lRyRMatrix[1][2] = M_kbp;
//		lRyRMatrix[1][3] = M_kcn;
//		lRyRMatrix[2][0] = 0;
//		lRyRMatrix[2][1] = M_kbp * std::pow(cCaSS, M_m);
//		lRyRMatrix[2][2] = - M_kbp;
//		lRyRMatrix[2][3] = 0;
//		lRyRMatrix[3][0] = 0;
//		lRyRMatrix[3][1] = M_kcp;
//		lRyRMatrix[3][2] = 0;
//		lRyRMatrix[3][3] = - M_kcn;
//	}

	void IonicJafriRiceWinslow::showMe()
	{
		std::cout << "\n\n\t\tIonicJafriRiceWinslow Informations\n\n";
		std::cout << "number of unkowns: "  << this->Size() << std::endl;

		std::cout << "\n\t\tList of model parameters:\n\n";
		std::cout << "areaCap: " << this->areaCap() << std::endl;
		std::cout << "volMyo: " << this->volMyo() << std::endl;
		std::cout << "volJSR: " << this->volJSR() << std::endl;
		std::cout << "volNSR: " << this->volNSR() << std::endl;
		std::cout << "volSS: " << this->volSS() << std::endl;
		std::cout << "concNa0: " << this->concNa0() << std::endl;
		std::cout << "concCa0: " << this->concCa0() << std::endl;
		std::cout << "lTrpnTot: " << this->lTrpnTot() << std::endl;
		std::cout << "hTrpnTot: " << this->hTrpnTot() << std::endl;
		std::cout << "kpHtrpn: " << this->kpHtrpn() << std::endl;
		std::cout << "knHtrpn: " << this->knHtrpn() << std::endl;
		std::cout << "kpLtrpn: " << this->kpLtrpn() << std::endl;
		std::cout << "knLtrpn: " << this->knLtrpn() << std::endl;
		std::cout << "cmdnTot: " << this->cmdnTot() << std::endl;
		std::cout << "csqnTot: " << this->csqnTot() << std::endl;
		std::cout << "constmCmdn: " << this->constmCmdn() << std::endl;
		std::cout << "constmCsqn: " << this->constmCsqn() << std::endl;
		std::cout << "capMem: " << this->capMem() << std::endl;
		std::cout << "farad: " << this->farad() << std::endl;
		std::cout << "temp: " << this->temp() << std::endl;
		std::cout << "gasConst: " << this->gasConst() << std::endl;
		std::cout << "maxCondNa: " << this->maxCondNa() << std::endl;
		std::cout << "maxCondKp: " << this->maxCondKp() << std::endl;
		std::cout << "permNaK: " << this->permNaK() << std::endl;
		std::cout << "KNaCa: " << this->kNaCa() << std::endl;
		std::cout << "constmNa: " << this->constmNa() << std::endl;
		std::cout << "constmCa: " << this->constmCa() << std::endl;
		std::cout << "kSat: " << this->kSat() << std::endl;
		std::cout << "eta: " << this->eta() << std::endl;
		std::cout << "courNaK: " << this->courNaK() << std::endl;
		std::cout << "constmNai: " << this->constmNai() << std::endl;
		std::cout << "constmK0: " << this->constmK0() << std::endl;
		std::cout << "permNsCa: " << this->permNsCa() << std::endl;
		std::cout << "constmNsCa: " << this->constmNsCa() << std::endl;
		std::cout << "courpCa: " << this->courpCa() << std::endl;
		std::cout << "constmpCa: " << this->constmpCa() << std::endl;
		std::cout << "maxCondCab: " << this->maxCondCab() << std::endl;
		std::cout << "maxCondNab: " << this->maxCondNab() << std::endl;
		std::cout << "maxRyRPerm: " << this->maxRyRPerm() << std::endl;
		std::cout << "leakRateConst: " << this->leakRateConst() << std::endl;
		std::cout << "pumpRateATPase: " << this->pumpRateATPase() << std::endl;
		std::cout << "constmUp: " << this->constmUp() << std::endl;
		std::cout << "timeConstNsrJsr: " << this->timeConstNsrJsr() << std::endl;
		std::cout << "timeConstSubMyo: " << this->timeConstSubMyo() << std::endl;
		std::cout << "kAPlus: " << this->kAPlus() << std::endl;
		std::cout << "kANeg: " << this->kANeg() << std::endl;
		std::cout << "kBPlus: " << this->kBPlus() << std::endl;
		std::cout << "kBNeg: " << this->kBNeg() << std::endl;
		std::cout << "kCPlus: " << this->kCPlus() << std::endl;
		std::cout << "kCNeg: " << this->kCNeg() << std::endl;
		std::cout << "coopParamN: " << this->coopParamN() << std::endl;
		std::cout << "coopParamM: " << this->coopParamM() << std::endl;
		std::cout << "intoOpenSt: " << this->intoOpenSt() << std::endl;
		std::cout << "outOpenSt: " << this->outOpenSt() << std::endl;
		std::cout << "intoOpenStca: " << this->intoOpenStCa() << std::endl;
		std::cout << "outOpentSt2: " << this->outOpenSt2() << std::endl;
		std::cout << "modeTParamA: " << this->modeTParamA() << std::endl;
		std::cout << "modeTParamB: " << this->modeTParamB() << std::endl;
		std::cout << "modeTParamO: " << this->modeTParamO() << std::endl;
		std::cout << "permCa: " << this->permCa() << std::endl;
		std::cout << "permK: " << this->permK() << std::endl;
		std::cout << "courCaHalf: " << this->courCaHalf() << std::endl;


		std::cout << "\n\t\t End of IonicJafriRiceWinslow Informations\n\n\n";
	}


	}
