#include <cmath>

#include "physiological_bc.hpp"



Real pRescaleFactor(1.);
Real uRescaleFactor(1.);
Real flowrateCorrection(1.);
Real PinClosure(1.e4);




Real zeroFct( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/)
{
    return 0.0;
}


Real aortaVelIn(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
	Real vBar = -linearFluxIn(t) * uRescaleFactor / 300.655 / unitsFactor * flowrateCorrection; // velocità = flusso / area INLET (mm*T)

	std::vector<Real> inletNormal(3,0);
	inletNormal[0] = 0.00927153;
	inletNormal[1] = 0.095069;
	inletNormal[2] = -0.995428;

	switch (i) // proiezione delle tre componenti della velocità sulla sezione di ingresso
	{
		case 0:
			return vBar * inletNormal[0];
			break;
		case 1:
			return vBar * inletNormal[1];
			break;
		case 2:
			return vBar * inletNormal[2];
			break;
		default:
			ERROR_MSG("This entry is not allowed");
			break;
	}

}



Real linearFluxIn(Real  t)
{
    Real const T = ( t/0.8 - std::floor(t/0.8) ) * 0.8 + 0.075;

    Real const a0(1.059e+05), a1(-1.287e+04), b1(8.786e+04), a2(-1019), b2(2.226e+04), a3(4385), b3(5228), a4(2132), b4(1046);
    Real const a5(755), b5(133.3), a6(-160.4), b6(313.8), a7(-627.1), b7(618.9), a8(-454.5), b8(282.4), w(19.73);

    Real const flux =
            a0 + a1*std::cos(T*w) + b1*std::sin(T*w) +
            a2*std::cos(2*T*w) + b2*std::sin(2*T*w) + a3*std::cos(3*T*w) + b3*std::sin(3*T*w) + 
            a4*std::cos(4*T*w) + b4*std::sin(4*T*w) + a5*std::cos(5*T*w) + b5*std::sin(5*T*w) + 
            a6*std::cos(6*T*w) + b6*std::sin(6*T*w) + a7*std::cos(7*T*w) + b7*std::sin(7*T*w) + 
            a8*std::cos(8*T*w) + b8*std::sin(8*T*w);

    if ( T > 0.27 && T < 0.59)
        return std::max(0.0, flux);
    return 0.;
}



Real inletFixedPressureFct( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/ )
{
    return -1.*pRescaleFactor;
}



/* Imposed pressure at the inlet */
Real inletPressureFct( const Real& t, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/)
{
    Real const T = ( t/0.8 - std::floor(t/0.8) ) * 0.8;

    Real const c(1500);
    Real const xv2(0.15), yv2(c), xp2(0.2), yp2(8264.53);
    Real const a0(11860), a1(3188), b1(-256.9), a2(-689.7), b2(-834), a3(122.3), b3(345.9), a4(-66.48), b4(-187.7);
    Real const a5(83.66), b5(87.61), a6(-96.29), b6(-46.74), a7(68.79), b7(67.38), a8(-6.933), b8(-55.44), w(1.471);
    Real const T3 = (T - 0.3535) / 0.08906;
    Real const xv4(0.59), yv4(c), xp4(0.508), yp4(10214.1); // yp4(PinClosure/pRescaleFactor)

    Real const pressure1 = c;
    Real const pressure2 = yv2 + (yp2 - yv2) / (xp2 - xv2) / (xp2 - xv2) * (T - xv2) * (T - xv2);
    Real const pressure3 =
            a0 + a1*std::cos(T3*w) + b1*std::sin(T3*w) +
            a2*std::cos(2*T3*w) + b2*std::sin(2*T3*w) + a3*std::cos(3*T3*w) + b3*std::sin(3*T3*w) + 
            a4*std::cos(4*T3*w) + b4*std::sin(4*T3*w) + a5*std::cos(5*T3*w) + b5*std::sin(5*T3*w) + 
            a6*std::cos(6*T3*w) + b6*std::sin(6*T3*w) + a7*std::cos(7*T3*w) + b7*std::sin(7*T3*w) + 
            a8*std::cos(8*T3*w) + b8*std::sin(8*T3*w);
    Real const pressure4 = yv4 + (yp4 - yv4) / (xp4 - xv4) / (xp4 - xv4) * (T - xv4) * (T - xv4);

    if ( T < 0.15 )
        return - pressure1 * unitsFactor * pRescaleFactor;
    if ( T < 0.2 )
        return - pressure2 * unitsFactor * pRescaleFactor;
    if ( T < 0.508 )
        return - pressure3 * unitsFactor * pRescaleFactor;
    if ( T < 0.59 )
        return - pressure4 * unitsFactor * pRescaleFactor;
    return - pressure1 * unitsFactor * pRescaleFactor;
}


Real outletFixedPressureFct( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/ )
{
    return -2.*pRescaleFactor;
}


/* TRIAL: Pressure at the OUTLET */
Real outletPressureFct( const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real T = ( t/0.8 - std::floor(t/0.8) ) * 0.8;

    Real const y1(10480), y2(8264.53), x1(0.508), x2(1);
    Real const a0(1.071e+04), a1(1424), b1(1546), a2(336.4), b2(-509.8), a3(-126.9), b3(21.74), a4(32.84), b4(-20.45);
    Real const a5(-45.44), b5(28.01), a6(40.2), b6(3.565), a7(1.129), b7(-23.74), a8(-9.05), b8(-0.5837), w(15.57);

    Real const pressure1 = y1 + (y2 - y1) / (x2 - x1) * (T - x1 + 0.8);
    Real const pressure3 = y1 + (y2 - y1) / (x2 - x1) * (T - x1);
    Real const T2 = (T + 0.086723684210526) / 1.013157894736842;
    Real const pressure2 =
            a0 + a1*std::cos(T2*w) + b1*std::sin(T2*w) +
            a2*std::cos(2*T2*w) + b2*std::sin(2*T2*w) + a3*std::cos(3*T2*w) + b3*std::sin(3*T2*w) +
            a4*std::cos(4*T2*w) + b4*std::sin(4*T2*w) + a5*std::cos(5*T2*w) + b5*std::sin(5*T2*w) +
            a6*std::cos(6*T2*w) + b6*std::sin(6*T2*w) + a7*std::cos(7*T2*w) + b7*std::sin(7*T2*w) +
            a8*std::cos(8*T2*w) + b8*std::sin(8*T2*w);

    if ( T < 0.2 )
        return - pressure1 * unitsFactor * pRescaleFactor;
    if ( T < 0.508 )
        return - pressure2 * unitsFactor * pRescaleFactor;
    return - pressure3 * unitsFactor * pRescaleFactor;
}


Real initPressureFct( const Real& t, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    if ( phiGlobalFct(t,x,y,z,1) < 0 && psiGlobalFct(t,x,y,z,1) < 0 )
        return -inletPressureFct(t,x,y,z,1);
    return -outletPressureFct(t,x,y,z,1);
}
