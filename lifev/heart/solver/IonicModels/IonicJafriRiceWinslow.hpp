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

#include <lifev/heart/solver/IonicModels/HeartIonicModel.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <cmath>


namespace LifeV
{
//! XbModel - This class implements a mean field model.

class IonicJafriRiceWinslow : public virtual HeartIonicModel
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

    inline const Real& concK0() const
    {
    	return M_K0;
    }
    inline void setK0(const Real& concK0)
    {
    	this->M_K0 = concK0;
    }

    inline const Real& concNa0() const
    {
    	return M_Na0;
    }
    inline void setNa0(const Real& concNa0)
    {
    	this->M_Na0 = concNa0;
    }

    inline const Real& concCa0() const
    {
    	return M_Ca0;
    }
    inline void setCa0(const Real& concCa0)
    {
    	this->M_Ca0 = concCa0;
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
    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

    void computeRhs ( const std::vector<Real>& v, const Real& Iapp, std::vector<Real>& rhs );


    //Compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp );

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

    Real noSpecInsCa( const std::vector<Real>& v );

    Real pumpIpCa( const std::vector<Real>& v );

    Real backICab( const std::vector<Real>& v );

    Real backINab( const std::vector<Real>& v);


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

    Real M_K0;
	Real M_Na0;
	Real M_Ca0;

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

//! Scientific notation calculus
inline double base10( int const& a )
{
	return pow(10, a);
}


// ===================================================
//! Constructors
// ===================================================
IonicJafriRiceWinslow::IonicJafriRiceWinslow()    :
    super       ( 9 ),
    M_ACap      ( 1.534 * base10( -4 ) ),
    M_VMyo      ( 25.84 * base10( -6 ) ),
    M_VJsr      ( 0.12 * base10( -6 ) ),
    M_VNsr      ( 2.098 * base10( -6 ) ),
    M_VSs       ( 1.485 * base10( -9 ) ),
    M_K0        ( 5.4 ),
    M_Na0       ( 140.0 ),
    M_Ca0       ( 1.8 ),
    M_LtrpnTot  ( 70.0 ),
    M_HtrpnTot  ( 140.0 ),
    M_kPHtrpn   ( 2.0 * base10( 7 ) ),
    M_kNHtrpn   ( 0.066 ),
    M_kPLtrpn   ( 4.0 * base10( 7 ) ),
    M_kNLtrpn   ( 40.0 ),
    M_CmdnTot   ( 50.0 ),
    M_CsqnTot   ( 15.0 ),
    M_KmCmdn    ( 2.38 ),
    M_KmCsqn    ( 0.8 ),
    M_F         ( 96.5 ),
    M_T         ( 310 ),
    M_R         ( 8.314 ),
    M_GNa       ( 12.8 ),
    M_GKp       ( 0.00828 ),
    M_PNaK      ( 0.01833 ),
    M_kNaCa     ( 5000.0 ),
    M_KmNa      ( 87.5 ),
    M_KmCa      ( 1.38 ),
    M_kSat      ( 0.1 ),
    M_eta       ( 0.35 ),
    M_INaK      ( 1.73 ),
    M_KmNai     ( 10.0 ),
    M_KmK0      ( 1.5 ),
    M_PnsCa     ( 0.0 ),
    M_KmNsCa    ( 1.2 ),
    M_IpCa      ( 1.15 ),
    M_KmPCa     ( 0.5 ),
    M_GCab      ( 0.00201 ),
    M_GNab      ( 0.00141 ),
	M_v1       	( 1.8 ),
	M_v2        ( 5.80 * base10( -5 ) ),
	M_v3        ( 1.8 ),
	M_KmUp      ( 0.5 ),
	M_tauTr     ( 34.48 ),
	M_tauXFer   ( 3.125 ),
	M_kap       ( 0.01215 ),
	M_kan       ( 0.1425 ),
	M_kbp       ( 0.00405 ),
	M_kbn       ( 1.930 ),
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
	M_PCa       ( 33.75 * base10( -4 ) ),
	M_PK        ( 62.5 * base10( -7 ) ),
	M_ICaHalf   ( -0.458 )
{}

IonicJafriRiceWinslow::IonicJafriRiceWinslow ( Teuchos::ParameterList& parameterList ) :
    super       ( 9 )
{
    M_ACap      = parameterList.get ( "areaCap", 1.534 * base10( -4 ) );
    M_VMyo      = parameterList.get ( "volMyo", 25.84 * base10( -6 ) );
    M_VJsr      = parameterList.get ( "volJSR", 0.12 * base10( -6 ) );
    M_VNsr      = parameterList.get ( "volNSR", 2.098 * base10( -6 ) );
    M_VSs       = parameterList.get ( "volSS", 1.485 * base10( -9 ) );
    M_K0        = parameterList.get ( "concK0", 5.4 );
    M_Na0       = parameterList.get ( "concNa0", 140.0 );
    M_Ca0       = parameterList.get ( "concCa0", 1.8 );
    M_LtrpnTot  = parameterList.get ( "lTrpnTot", 70.0 );
    M_HtrpnTot  = parameterList.get ( "hTrpnTot", 140.0 );
    M_kPHtrpn   = parameterList.get ( "kpHtrpn", 2.0 * base10( 7 ) );
    M_kNHtrpn   = parameterList.get ( "knHtrpn", 0.066 );
    M_kPLtrpn   = parameterList.get ( "kpLtrpn", 4.0 * base10( 7 ) );
    M_kNLtrpn   = parameterList.get ( "knLtrpn", 40.0 );
    M_CmdnTot   = parameterList.get ( "cmdnTot", 50.0 );
    M_CsqnTot   = parameterList.get ( "csqnTot", 15.0 );
    M_KmCmdn    = parameterList.get ( "constmCmdn", 2.38 );
    M_KmCsqn    = parameterList.get ( "constmCsqn", 0.8 );
    M_F         = parameterList.get ( "farad", 96.5 );
    M_T         = parameterList.get ( "temp", 310 );
    M_R         = parameterList.get ( "gasConst", 8.314 );
    M_GNa       = parameterList.get ( "maxCondNa", 12.8 );
    M_GKp       = parameterList.get ( "maxCondKp", 0.00828 );
    M_PNaK      = parameterList.get ( "permNaK", 0.01833 );
    M_kNaCa     = parameterList.get ( "kNaCA", 5000.0 );
    M_KmNa      = parameterList.get ( "constmNa", 87.5 );
    M_KmCa      = parameterList.get ( "constmCa", 1.38 );
    M_kSat      = parameterList.get ( "kSat", 0.1 );
    M_eta       = parameterList.get ( "eta", 0.35 );
    M_INaK      = parameterList.get ( "courNaK", 1.73 );
    M_KmNai     = parameterList.get ( "constmNai", 10.0 );
    M_KmK0      = parameterList.get ( "constmK0", 1.5 );
    M_PnsCa     = parameterList.get ( "permNsCa", 0.0 );
    M_KmNsCa    = parameterList.get ( "constmNsCa", 1.2 );
    M_IpCa      = parameterList.get ( "courpCa", 1.15 );
    M_KmPCa     = parameterList.get ( "constmpCa", 0.5 );
    M_GCab      = parameterList.get ( "maxCondCab", 0.00201 );
    M_GNab      = parameterList.get ( "maxCondNab", 0.00141 );
    M_v1       	= parameterList.get ( "maxRyRPerm", 1.8 );
    M_v2        = parameterList.get ( "leakRateConst", 5.80 * base10( -5 ) );
    M_v3        = parameterList.get ( "pumpRateATPase", 1.8 );
    M_KmUp      = parameterList.get ( "constmUp", 0.5 );
    M_tauTr     = parameterList.get ( "timeConstNsrJsr", 34.48 );
    M_tauXFer   = parameterList.get ( "timeConstSubMyo", 3.125 );
    M_kap       = parameterList.get ( "kAPlus", 0.01215 );
    M_kan       = parameterList.get ( "kANeg", 0.1425 );
    M_kbp       = parameterList.get ( "kBPlus", 0.00405 );
    M_kbn       = parameterList.get ( "kBNeg", 1.930 );
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
    M_PCa       = parameterList.get ( "permCa", 33.75 * base10( -4 ) );
    M_PK        = parameterList.get ( "permK", 62.5 * base10( -7 ) );
    M_ICaHalf   = parameterList.get ( "courCaHalf", -0.458 );
}

IonicJafriRiceWinslow::IonicJafriRiceWinslow ( const IonicJafriRiceWinslow& model )
{
    M_ACap		= model.M_ACap;
    M_VMyo		= model.M_VMyo;
    M_VJsr		= model.M_VJsr;
    M_VNsr		= model.M_VNsr;
    M_VSs		= model.M_VSs;
    M_K0		= model.M_K0;
    M_Na0		= model.M_Na0;
    M_Ca0 		= model.M_Ca0;
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
    M_K0        = model.M_K0;
    M_Na0       = model.M_Na0;
    M_Ca0       = model.M_Ca0;
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
	std::vector<Real> gatingRhs ( computeLocalGatingRhs(v) );
	std::vector<Real> concRhs ( computeLocalConcRhs(v) );
	std::vector<Real> subSysCaRhs ( computeLocalSubSysCaRhs(v) );
	std::vector<Real> channelRyrRhs ( computeLocalChannelRyrRhs(v) );

	std::copy(gatingRhs.begin(), gatingRhs.end(), rhs.begin());
	std::copy(concRhs.begin(), concRhs.end() - 2, rhs.begin() + gatingRhs.size());
	std::copy(subSysCaRhs.begin() + 2, subSysCaRhs.begin() + 5 , rhs.begin() + ( concRhs.size() - 2 ) + gatingRhs.size());
	std::copy(channelRyrRhs.begin(), channelRyrRhs.end(), rhs.begin() + ( concRhs.size() - 2 ) + gatingRhs.size() + 4);
	std::copy(subSysCaRhs.begin() + 6, subSysCaRhs.end(), rhs.begin() + ( concRhs.size() - 2 ) + gatingRhs.size() + 4 + channelRyrRhs.size());
	std::copy(concRhs.end() - 2, concRhs.end(), rhs.begin() + ( concRhs.size() - 2 ) + gatingRhs.size() + ( subSysCaRhs.size() - 4 ) );

}

//Potential and gating variables
void IonicJafriRiceWinslow::computeRhs (const   std::vector<Real>&  v,
                                         const   Real&           Iapp,
                                         std::vector<Real>& rhs )
{
	std::vector<Real> gatingRhs ( computeLocalGatingRhs(v) );
	std::vector<Real> concRhs ( computeLocalConcRhs(v) );
	std::vector<Real> subSysCaRhs ( computeLocalSubSysCaRhs(v) );
	std::vector<Real> channelRyrRhs ( computeLocalChannelRyrRhs(v) );

	rhs[0] = computeLocalPotentialRhs(v, Iapp);
	std::copy(gatingRhs.begin(), gatingRhs.end(), rhs.begin() + 1 );
	std::copy(concRhs.begin(), concRhs.end() - 2, rhs.begin() + 1 + gatingRhs.size());
	std::copy(subSysCaRhs.begin() + 2, subSysCaRhs.begin() + 5 , rhs.begin() + 1 + ( concRhs.size() - 2 ) + gatingRhs.size());
	std::copy(channelRyrRhs.begin(), channelRyrRhs.end(), rhs.begin() + 1 + ( concRhs.size() - 2 ) + gatingRhs.size() + 4);
	std::copy(subSysCaRhs.begin() + 6, subSysCaRhs.end(), rhs.begin() + 1 + ( concRhs.size() - 2 ) + gatingRhs.size() + 4 + channelRyrRhs.size());
	std::copy(concRhs.end() - 2, concRhs.end(), rhs.begin() + 1 + ( concRhs.size() - 2 ) + gatingRhs.size() + ( subSysCaRhs.size() - 4 ) );

}


Real IonicJafriRiceWinslow::computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp )
{
	std::vector<Real>param1 ( computeLocalSubSysCaRhs(v) );
	std::vector<Real>param2 ( fastINa(v) );
	std::vector<Real>param3 ( timeDIK(v) );

	return ( - ( param2[0] + param1[0] + param3[0] + timeIIK1(v) + plaIKp(v) +
			exINaCa(v) + pumpINaK(v) + noSpecInsCa(v) + pumpIpCa(v) + param1[1] + backICab(v) + backINab(v) )  + Iapp );
}

std::vector<Real> IonicJafriRiceWinslow::computeLocalGatingRhs ( const std::vector<Real>& v )
{
	Real m = v[1];
	Real h = v[2];
	Real j = v[3];
	Real x = v[4];

	std::vector<Real> gatingRhs(4);
	std::vector<Real>param1 ( fastINa(v) );
	std::vector<Real>param2 ( timeDIK(v) );

	gatingRhs[0] = param1[1] * ( 1 - m ) - param1[2] * m;
	gatingRhs[1] = param1[3] * ( 1 - h ) - param1[5] * h;
 	gatingRhs[2] = param1[4] * ( 1 - j ) - param1[6] * j;
	gatingRhs[3] = param2[1] * ( 1 - x ) - param2[2] * x;

	return gatingRhs;
}

std::vector<Real> IonicJafriRiceWinslow::computeLocalConcRhs ( const std::vector<Real>& v )
{

	std::vector<Real>param1 ( computeLocalSubSysCaRhs(v) );
	std::vector<Real>param2 ( fastINa(v) );
	std::vector<Real>param3 ( timeDIK(v) );

	std::vector<Real> concRhs(4);

	concRhs[0] = - ( param2[0] + backINab(v) + noSpecInsCa(v) + 3 * exINaCa(v) + 3 * pumpINaK(v) ) * M_ACap / ( M_VMyo * M_F );
	concRhs[1] = - ( param3[0] + timeIIK1(v) + plaIKp(v) + noSpecInsCa(v) - 2 * pumpINaK(v) + param1[1] ) * M_ACap / ( M_VMyo * M_F );
	concRhs[2] = M_kPHtrpn * v[7] * ( M_HtrpnTot - v[29] ) - M_kNHtrpn * v[29];
	concRhs[3] = M_kPLtrpn * v[7] * ( M_LtrpnTot - v[28] ) - M_kNLtrpn * v[28];

	return concRhs;
}



//! Ca2+ Subsystem

std::vector<Real> IonicJafriRiceWinslow::computeLocalSubSysCaRhs( const std::vector<Real>& v )
{
	Real V = v[0];
	Real cK = v[6];
	Real cCa = v[7];
	Real cCaNSR = v[8];
	Real cCaSS = v[9];
	Real cCaJSR = v[10];

	Real fracPO1 = v[12];
	Real fracPO2 = v[13];

	Real c0 = v[15];
	Real c1 = v[16];
	Real c2 = v[17];
	Real c3 = v[18];
	Real c4 = v[19];
	Real o = v[20];
	Real cCa0 = v[21];
	Real cCa1 = v[22];
	Real cCa2 = v[23];
	Real cCa3 = v[24];
	Real cCa4 = v[25];
	Real oCa = v[26];
	Real y = v[27];
	Real cLTRPNCa = v[28];
	Real cHTRPNCa = v[29];

	std::vector<Real> subSysCaRHS(19);


	// Internal Parameters


	Real jRel = M_v1 * ( fracPO1 + fracPO2 ) * ( cCaJSR - cCaSS );
	Real jLeak = M_v2 * ( cCaNSR - cCa );
	Real jUp = M_v3 * ( pow(cCa, 2) / ( pow(M_KmUp, 2) + pow(cCa, 2) ) );
	Real jTr = ( cCaNSR - cCaJSR ) / M_tauTr;
	Real jXFer = ( cCaSS - cCa ) / M_tauXFer;
	Real jTRPN = M_kPHtrpn * cCa * ( M_HtrpnTot - cHTRPNCa ) - M_kNHtrpn * cHTRPNCa
			+ M_kPLtrpn * ( M_LtrpnTot - cLTRPNCa ) - M_kNLtrpn * cLTRPNCa;

	Real bI = 1 / ( 1 + M_CmdnTot * M_KmCmdn / pow(M_KmCmdn + cCa, 2) );
	Real bSS = 1 / ( 1 + M_CmdnTot * M_KmCmdn / pow(M_KmCmdn + cCaSS, 2) );
	Real bJSR = 1 / ( 1 + M_CsqnTot * M_KmCsqn / pow(M_KmCsqn + cCaJSR, 2) );

	Real alpha = 0.4 * exp( ( V + 12 ) / 10);
	Real beta = 0.05 * exp( -( V + 12 ) / 13);
	Real alphaprime = M_a * alpha;
	Real betaprime = beta / M_b;
	Real gamma = 0.1875 * cCaSS;

	Real yinf = 1 / ( 1 + exp( ( V + 55 ) / 7.5 ) ) + 0.1 / ( 1 + exp( -V + 21 ) / 6 );
	Real tauY = 20 + 600 / ( 1 + exp( V + 30 ) / 9.5 );
	Real iCaMax = M_PCa * 4 * (V * pow(M_F, 2) ) / ( M_R * M_T ) * ( 0.001 * exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 0.341 * M_Ca0 )
			/ ( exp( 2 * ( V * M_F ) / ( M_R * M_T ) ) - 1);
	Real iCa = y * ( o + oCa ) * iCaMax;
	Real pK = M_PK / ( 1 + iCaMax / M_ICaHalf );
	Real iCaK = pK * y * (o + oCa ) * (V * pow(M_F, 2) ) / ( M_R * M_T ) * ( cK * exp( ( V * M_F ) / ( M_R * M_T ) ) - M_K0 )
			/ ( exp( ( V * M_F ) / ( M_R * M_T ) ) - 1 );


	// RHS of the Ca2+ subsystem

	subSysCaRHS[0] = iCa;
	subSysCaRHS[1] = iCaK;
	subSysCaRHS[2] = bI * ( jLeak + jXFer - jUp - jTRPN
			- ( backICab(v) - 2 * exINaCa(v) + pumpIpCa(v) ) * M_ACap / ( 2 * M_VMyo * M_F ) );
	subSysCaRHS[3] = bSS * ( jRel * M_VJsr / M_VSs - jXFer * M_VMyo / M_VSs - iCa * M_ACap / ( 2 * M_VSs * M_F ) );
	subSysCaRHS[4] = bJSR * ( jTr - jRel );
	subSysCaRHS[5] = ( jUp - jLeak ) * M_VMyo / M_VNsr - jTr * M_VJsr / M_VNsr;
	subSysCaRHS[6] = beta * c1 + M_omega * cCa0 - ( 4 * alpha + gamma ) * c0;
	subSysCaRHS[7] = 4 * alpha * c0 + 2 * beta * c2 + cCa1 * M_omega / M_b - ( beta + 3 * alpha + gamma * M_a ) * c1;
	subSysCaRHS[8] = 3 * alpha * c1  + 3 * beta * c3 + cCa2 * M_omega * pow(M_b, 2) - ( 2 * beta + 2* alpha + gamma * pow(M_a, 2) ) * c2;
	subSysCaRHS[9] = 2 * alpha * c2 + 4 * beta * c4 + cCa3 * M_omega / pow(M_b, 3) - ( 3 * beta + alpha + gamma * pow(M_a, 3) ) * c3;
	subSysCaRHS[10] = alpha *c3 + M_g * o + cCa4 * M_omega / pow(M_b, 4) - ( 4 * beta + M_f + gamma * pow(M_a, 4) ) * c4;
	subSysCaRHS[11] = M_f * c4 - M_g * o;
	subSysCaRHS[12] = betaprime * cCa1 + gamma * c0 - ( 4 * alphaprime + M_omega ) * cCa0;
	subSysCaRHS[13] = 4 * alphaprime * cCa0 + 2 * betaprime * cCa2 + gamma * M_a * c1 - ( betaprime + 3 * alphaprime + M_omega / M_b ) * cCa1;
	subSysCaRHS[14] = 3 * alphaprime * cCa1 +  3 * betaprime * cCa3 + gamma * pow(M_a, 2) * c2 - ( 2* betaprime + 2 * alphaprime + M_omega / pow(M_b, 2) ) * cCa2;
	subSysCaRHS[15] = 2 * alphaprime * cCa2 + 4 * betaprime * cCa4 + gamma * pow(M_a, 3) * c3  - ( 3 * betaprime + alphaprime + M_omega / pow(M_b, 3) ) * cCa3;
	subSysCaRHS[16] = alphaprime * cCa3 + M_gprime * oCa + gamma * pow(M_a, 4) * c4 -( 4 * betaprime + M_fprime + M_omega / pow(M_b, 4) ) * cCa4;
	subSysCaRHS[17] = M_fprime * cCa4 -M_gprime * oCa;
	subSysCaRHS[18] = ( yinf - y ) / tauY ;


	return subSysCaRHS;
}


//! Ionic Currents (Luo and Rudy)

// Fast Na+ Current INa
std::vector<Real> IonicJafriRiceWinslow::fastINa( const std::vector<Real>& v )
{
	std::vector<Real> fastNa(7);
	Real V = v[0];
	Real m = v[1];
	Real h = v[2];
	Real j = v[3];
	Real cNa = v[5];

	Real potNa = ( M_R * M_T / M_F ) * log( M_Na0 / cNa );

	fastNa[0] = M_GNa * pow( m, 3 )* h * j * ( V - potNa );
	fastNa[1] = 0.32 * ( V + 47.13 ) / ( 1- exp( 0.1 * ( V + 47.13 ) ) );
	fastNa[2] = 0.08 * exp(-V / 11.0 );

	if (V >= -40)
	{
		fastNa[3] = 0.0;
		fastNa[4] = 0.0;
		fastNa[5] = 1 / ( 0.13 * ( 1 + exp( ( V + 10.66 ) / -11.1 ) ) );
		fastNa[6] = 0.3 * exp( -2.535 * base10( -7 ) * V) / ( 1 + exp( -0.1 *( V + 32.0 ) ) );
	}
	else
	{
		fastNa[3] = 0.135 * exp( ( 80 + V ) / -6.8 );
		fastNa[4] = ( -127140 * exp( 0.2444 * V) - 3.474 * base10(-5) * exp( -0.04391 * V ) ) * ( V + 37.78 ) / ( 1 + exp ( 0.311 * ( V + 79.23 ) ) );
		fastNa[5] = 3.56 * exp( 0.079 * V ) + 3.1 * base10(5) * exp( 0.35 * V ) ;
		fastNa[6] = 0.1212 * exp( -0.01052 * V ) / ( 1 + exp( -0.1378 * (V + 40.14 ) ) );
	}

	return fastNa;
}

// Time dependent K+ currrent IK
std::vector<Real> IonicJafriRiceWinslow::timeDIK( const std::vector<Real>& v )
{
	std::vector<Real> timeDK(3);
	Real V = v[0];
	Real x = v[4];
	Real cNa = v[5];
	Real cK = v[6];


	Real potK = M_R * M_T / M_F * log( ( M_K0 + M_PNaK * M_Na0 ) / ( cK + M_PNaK * cNa) );
	Real maxCondK = 0.1128 * sqrt( M_K0 / 5.4);
	Real xi = 1 / ( 1 + exp( V - 56.26) / 32.1 );

	timeDK[0] = maxCondK * xi * pow(x, 2) * ( V - potK );
	timeDK[1] = 7.19 * base10(-5) * ( V + 30 ) / ( 1 - exp( -0.148 * ( V + 30 ) ) );
	timeDK[2] = 1.31 * base10(-4) * ( V + 30 ) / ( -1 + exp( 0.0687 * ( V + 30 ) ) );

	return timeDK;
}

// Time-independent K+ Current IK1
Real IonicJafriRiceWinslow::timeIIK1( const std::vector<Real>& v )
{
	Real V = v[0];
	Real cK = v[6];

	Real potK1 = M_R * M_T / M_F * log( M_K0 / cK );
	Real maxCondK1 = 0.75 * sqrt( M_K0 / 5.4);
	Real alphaK1 = 1.02 / ( 1 + exp( 0.2385 * ( V - potK1 -59.215 ) ) );
	Real betaK1 = ( 0.4912 * exp( 0.08032 * ( V - potK1 + 5.476 ) ) + exp( 0.06175 * (V - potK1 - 594.31 ) ) ) / ( 1 + exp( -0.5143 * (V - potK1 + 4.753 ) ) );
	Real constK1inf = alphaK1 / ( alphaK1 + betaK1 );

	return  maxCondK1 * constK1inf * ( V - potK1 );
}

// Plateau K+ current IKp
Real IonicJafriRiceWinslow::plaIKp( const std::vector<Real>& v )
{
	Real potKp = M_R * M_T / M_F * log( M_K0 / v[6] );
	Real constKp = 1 / ( 1 + exp( ( 7.488 - v[0] ) / 5.98 ) );

	return M_GKp * constKp *( v[0] - potKp );
}

// Na+/Ca2+ exchanger current INaCa
Real IonicJafriRiceWinslow::exINaCa( const std::vector<Real>& v )
{

	Real V = v[0];
	Real cNa = v[5];
	Real cCa = v[7];


	 return M_kNaCa * ( 1 / ( pow(M_KmNa, 3) + pow(M_Na0, 3) ) ) * ( 1 / ( M_KmCa + M_Ca0 ) )
			* (1 / ( 1 + M_kSat * exp( ( M_eta - 1) * ( V * M_F ) / ( M_R * M_T ) ) ) )
			* ( exp( M_eta * ( V * M_F ) / ( M_R * M_T ) ) * pow(cNa, 3) * M_Ca0 -
					exp( ( M_eta - 1 ) * ( V * M_F ) / ( M_R * M_T ) ) * pow(M_Na0, 3) * cCa );
}

// Na+/K+ pump INaK
Real IonicJafriRiceWinslow::pumpINaK( const std::vector<Real>& v )
{
	Real V = v[0];
	Real cNa = v[5];

	Real sigma = ( 1 / 7 ) * ( exp( M_Na0 / 67.3) - 1);
	Real fNak = 1 / ( 1 + 0.1245 * exp( -0.1 * ( V * M_F ) / ( M_R * M_T ) ) + 0.0365 * sigma * exp( -( V * M_F ) / ( M_R * M_T ) ) );

	return M_INaK * fNak * ( 1 / ( 1 + pow( M_KmNai / cNa, 1.5) ) ) * ( M_K0 / ( M_K0 + M_KmK0 ) );

}

// Nonspecific Ca2+ activated current InsCa
Real IonicJafriRiceWinslow::noSpecInsCa( const std::vector<Real>& v )
{

	Real V = v[0];
	Real cNa = v[5];
	Real cK = v[6];
	Real cCa = v[7];

	Real maxInsNa = M_PnsCa * ( ( V * pow(M_F, 2) ) / ( M_R * M_T ) ) *
			( ( 0.75 * cNa * exp( ( V * M_F ) / ( M_R * M_T ) ) - 0.75 * M_Na0 ) / ( exp( ( V * M_F ) / ( M_R * M_T ) ) - 1 ) );

	Real maxInsK = M_PnsCa * ( ( V * pow(M_F, 2) ) / ( M_R * M_T ) ) *
			( ( 0.75 * cK * exp( ( V * M_F ) / ( M_R * M_T ) ) - 0.75 * M_K0 ) / ( exp( ( V * M_F ) / ( M_R * M_T ) ) - 1 ) );

	Real courInsNa = maxInsNa * ( 1 / ( 1 + pow( M_KmNsCa /cCa, 3) ) );
	Real courInsK = maxInsK * ( 1 / ( 1 + pow( M_KmNsCa /cCa, 3) ) );

	return courInsK + courInsNa;

}

// Sarcolemmal Ca2+ pump current IpCa
Real IonicJafriRiceWinslow::pumpIpCa( const std::vector<Real>& v )
{
	return M_IpCa * v[7] / ( M_KmPCa * v[7]);
}

// Ca2+ background current ICab
Real IonicJafriRiceWinslow::backICab( const std::vector<Real>& v )
{
	Real V = v[0];
	Real cCa = v[7];

	Real potCaN = ( ( M_R * M_T ) / ( 2 * M_F ) ) * log( M_Ca0 / cCa);

	return M_GCab * ( V - potCaN );
}

// Na+ background current INab
Real IonicJafriRiceWinslow::backINab( const std::vector<Real>& v)
{
	Real V = v[0];
	Real cNa = v[5];

	Real potNaN = M_R * M_T / M_F * log( M_Na0 / cNa );

	return M_GNab * ( V - potNaN );
}


//! RyR Channel States (Keizer and Levine)

std::vector<Real> IonicJafriRiceWinslow::computeLocalChannelRyrRhs( const std::vector<Real>& v )
{
	Real cCaSS= v[9];
	Real fracPC1 = v[11];
	Real fracPO1 = v[12];
	Real fracPO2 = v[13];
	Real fracPC2 = v[14];

	std::vector<Real> channelRyrRHS(4);

	channelRyrRHS[0] = -M_kap * pow(cCaSS, M_n) * fracPC1 + M_kan * fracPO1;
	channelRyrRHS[1] = M_kap * pow(cCaSS, M_n) * fracPC1 - M_kan * fracPO1 -
			M_kbp * pow(cCaSS, M_m) * fracPO1 + M_kbn * fracPO2 - M_kcp * fracPO1 + M_kcn * fracPC2;
	channelRyrRHS[2] = M_kbp * pow(cCaSS, M_m) * fracPO1 - M_kbn * fracPO2;
	channelRyrRHS[3] = M_kcp * fracPO1 - M_kcn * fracPC2;

	return channelRyrRHS;
}

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
    std::cout << "concK0: " << this->concK0() << std::endl;
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
    std::cout << "constmpCa: " << this->constmCa() << std::endl;
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

// ===================================================
//! Additional Methods
// ===================================================








#endif // IonicJafriRiceWinslow_HPP_INCLUDED
