#include "levelset_valve.hpp"


/*
    open leaflet defined with a 5th degree polynomial and divided per its gradient module to make it similar
    to a distance function near its zero level
*/
#define OPEN_FIFTH_DEGREE_POLY                                                                                      \
            (                                                                                                       \
                Zr - q00 - q10*Xr - q01*Yr - q20*Xr*Xr - q11*Xr*Yr - q02*Yr*Yr -                                    \
                q30*Xr*Xr*Xr - q21*Xr*Xr*Yr - q12*Xr*Yr*Yr - q03*Yr*Yr*Yr -                                         \
                q40*Xr*Xr*Xr*Xr - q31*Xr*Xr*Xr*Yr - q22*Xr*Xr*Yr*Yr - q13*Xr*Yr*Yr*Yr - q04*Yr*Yr*Yr*Yr -           \
                q50*Xr*Xr*Xr*Xr*Xr - q41*Xr*Xr*Xr*Xr*Yr - q32*Xr*Xr*Xr*Yr*Yr -                                      \
                q23*Xr*Xr*Yr*Yr*Yr - q14*Xr*Yr*Yr*Yr*Yr - q05*Yr*Yr*Yr*Yr*Yr                                        \
            ) /                                                                                                     \
            (                                                                                                       \
                unitsFactor * std::sqrt                                                                             \
                (                                                                                                   \
                    (                                                                                               \
                        q10 + 2*q20*Xr + q11*Yr + 3*q30*Xr*Xr + 2*q21*Yr*Xr + q12*Yr*Yr + 4*q40*Xr*Xr*Xr +          \
                        3*q31*Yr*Xr*Xr + 2*q22*Yr*Yr*Xr + q13*Yr*Yr*Yr + 5*q50*Xr*Xr*Xr*Xr + 4*q41*Yr*Xr*Xr*Xr +    \
                        3*q32*Yr*Yr*Xr*Xr + 2*q23*Yr*Yr*Yr*Xr + q14*Yr*Yr*Yr*Yr                                     \
                    ) *                                                                                             \
                    (                                                                                               \
                        q10 + 2*q20*Xr + q11*Yr + 3*q30*Xr*Xr + 2*q21*Yr*Xr + q12*Yr*Yr + 4*q40*Xr*Xr*Xr +          \
                        3*q31*Yr*Xr*Xr + 2*q22*Yr*Yr*Xr + q13*Yr*Yr*Yr + 5*q50*Xr*Xr*Xr*Xr + 4*q41*Yr*Xr*Xr*Xr +    \
                        3*q32*Yr*Yr*Xr*Xr + 2*q23*Yr*Yr*Yr*Xr + q14*Yr*Yr*Yr*Yr                                     \
                    ) +                                                                                             \
                    (                                                                                               \
                        q01 + 2*q02*Yr + q11*Xr + 3*q03*Yr*Yr + 2*q12*Xr*Yr + q21*Xr*Xr + 4*q04*Yr*Yr*Yr +          \
                        3*q13*Xr*Yr*Yr + 2*q22*Xr*Xr*Yr + q31*Xr*Xr*Xr + 5*q05*Yr*Yr*Yr*Yr + 4*q14*Xr*Yr*Yr*Yr +    \
                        3*q23*Xr*Xr*Yr*Yr + 2*q32*Xr*Xr*Xr*Yr + q41*Xr*Xr*Xr*Xr                                     \
                    ) *                                                                                             \
                    (                                                                                               \
                        q01 + 2*q02*Yr + q11*Xr + 3*q03*Yr*Yr + 2*q12*Xr*Yr + q21*Xr*Xr + 4*q04*Yr*Yr*Yr +          \
                        3*q13*Xr*Yr*Yr + 2*q22*Xr*Xr*Yr + q31*Xr*Xr*Xr + 5*q05*Yr*Yr*Yr*Yr + 4*q14*Xr*Yr*Yr*Yr +    \
                        3*q23*Xr*Xr*Yr*Yr + 2*q32*Xr*Xr*Xr*Yr + q41*Xr*Xr*Xr*Xr                                     \
                    )                                                                                               \
                    + 1.0                                                                                           \
                )                                                                                                   \
            )


#define CLOSED_SECOND_DEGREE_POLY                                                                                   \
            ( Z - p00 - p10*X - p01*Y - p20*X*X - p11*X*Y - p02*Y*Y) /                                              \
            ( unitsFactor * std::sqrt( (p10 + 2*p20*X + p11*Y)*(p10 + 2*p20*X + p11*Y) +                            \
                                       (p01 + 2*p02*Y + p11*X)*(p01 + 2*p02*Y + p11*X) +                            \
                                        1.0 ) )


#define OPEN_SECOND_DEGREE_POLY                                                                                     \
            ( Z - q00 - q10*X - q01*Y - q20*X*X - q11*X*Y - q02*Y*Y ) /                                             \
            ( unitsFactor * std::sqrt( (q10 + 2*q20*X + q11*Y)*(q10 + 2*q20*X + q11*Y) +                            \
                                       (q01 + 2*q02*Y + q11*X)*(q01 + 2*q02*Y + q11*X) +                            \
                                        1.0 ) )



Real angleCoeff(0.); //corresponds to angle=angleMin initially

// Factor to resize valve equation depending on the units of measure used during their computation
Real unitsFactor(10.); // equation in mm, mesh in cm

Real epsilon(0.1);
bool useSmoothHeaviside(false);


// Closed Valve Planes Coefficients
Real const nr00(-0.8508), nr10(-3.436), nr01(3.268);
Real const lr00(2.83), lr10(0.1644), lr01(4.224);
Real const ln00(-5.448), ln10(-3.292), ln01(-0.5809);



Real phiNonFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    Real const X = unitsFactor * x;
    Real const Y = unitsFactor * y;
    Real const Z = unitsFactor * z;

    Real const p00(-2.449), p10(2.235), p01(-0.5299), p20(0.1401), p11(-0.008372), p02(0.1007);

    Real closedPhiNon = CLOSED_SECOND_DEGREE_POLY;
    Real PlaneNonRight = (nr00 + nr10*X + nr01*Y - Z) / ( unitsFactor * std::sqrt(nr10*nr10 + nr01*nr01 + 1.0) );
    Real PlaneLeftNon = (ln00 + ln10*X + ln01*Y - Z) / ( unitsFactor * std::sqrt(ln10*ln10 + ln01*ln01 + 1.0) );

    closedPhiNon = std::min( closedPhiNon, std::min(PlaneLeftNon, PlaneNonRight) );

    // open valve template
    //Real const q00(), q10(), q01(), q20(), q11(), q02();
    //Real const q30(), q21(), q12(), q03();
    //Real const q40(), q31(), q22(), q13(), q04();
    //Real const q50(), q41(), q32(), q23(), q14(), q05();
    //Real const Xr = *X + *Y + *Z;
    //Real const Yr = *X + *Y + *Z;
    //Real const Zr = *X + *Y + *Z;
    Real const q00(11.27)/*11.72*/, q10(-0.0141), q01(0.9265), q20(-0.04777), q11(0.06885), q02(-0.1371);
    Real const q30(-0.0003492), q21(-0.01029), q12(0.007471), q03(-0.03105);
    Real const q40(0.0002975), q31(7.934e-6), q22(0.001547), q13(-0.001339), q04(0.001099);
    Real const q50(-2.48e-6), q41(3.747e-5), q32(-2.5e-5), q23(1.121e-5), q14(-4.677e-5), q05(0.0001805);
    Real const Xr = 0.00279372*X - 0.0527818*Y + 0.998602*Z;
    Real const Yr = -0.0527818*X + 0.997206*Y + 0.0528557*Z;
    Real const Zr = -0.998602*X - 0.0528557*Y + 0.0*Z;

    Real openPhiNon = OPEN_FIFTH_DEGREE_POLY;
    //if ( z < 0.7 )
    openPhiNon = std::min(closedPhiNon, openPhiNon);

    return angleCoeff*openPhiNon + (1.-angleCoeff)*closedPhiNon;
}


Real psiNonFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    Real const X = unitsFactor * x;
    Real const Y = unitsFactor * y;
    Real const Z = unitsFactor * z;

    Real const p00(1.788), p10(-0.356), p01(-0.1753), p20(-0.01966), p11(0.02025), p02(0.02268);

    Real closedPsiNon = CLOSED_SECOND_DEGREE_POLY;

    Real const q00(3.3404), q10(0), q01(0), q20(0), q11(0), q02(0);

    Real openPsiNon = OPEN_SECOND_DEGREE_POLY;

    return angleCoeff*openPsiNon + (1.-angleCoeff)*closedPsiNon;
}


Real phiRightFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    Real const X = unitsFactor * x;
    Real const Y = unitsFactor * y;
    Real const Z = unitsFactor * z;

    Real const p00(-0.4305), p10(-0.7846), p01(1.767), p20(0.0778), p11(-0.05065), p02(0.1136);

    Real closedPhiRight = CLOSED_SECOND_DEGREE_POLY;
    Real PlaneNonRight = (Z - nr00 - nr10*X - nr01*Y) / ( unitsFactor * std::sqrt( nr10*nr10 + nr01*nr01 + 1.0 ) );
    Real PlaneLeftRight = (Z - lr00 - lr10*X - lr01*Y) / ( unitsFactor * std::sqrt( lr10*lr10 + lr01*lr01 + 1.0 ) );

    closedPhiRight = std::min( closedPhiRight, std::min(PlaneLeftRight, PlaneNonRight) );

    // works only with mesh in cm, usefull (but not necessary) onlywith epsilon > 0.15
//    Real auxPlane = -( -0.42*(z-0.64) + 0.85*(y+0.69) - 0.31*(x-0.33) ) / ( std::sqrt( 0.31*0.31 + 0.85*0.85 + 0.35*0.35 ) );
//    closedPhiRight = std::max( closedPhiRight, auxPlane );

    Real const q00(11.2), q10(-0.2314), q01(-0.4394), q20(-0.1233), q11(-0.08756), q02(-0.1611);
    Real const q30(0.003853), q21(0.01493), q12(0.02175), q03(0.008018);
    Real const q40(0.000743), q31(0.0009229), q22(0.002858), q13(0.001287), q04(0.002998);
    Real const q50(-5.268e-06), q41(-0.0001251), q32(-0.0002307), q23(-0.0003675), q14(-0.0002329), q05(-0.0002151);
    Real const Xr = 0.868781*X + 0.33764*Y - 0.362242*Z;
    Real const Yr = 0.33764*X + 0.131219*Y + 0.932084*Z;
    Real const Zr = 0.362242*X - 0.932084*Y + 2.22045e-16*Z;

    Real openPhiRight = OPEN_FIFTH_DEGREE_POLY;
    //if ( z < 0.7 )
    openPhiRight = std::min(closedPhiRight, openPhiRight);

    return angleCoeff*openPhiRight + (1.-angleCoeff)*closedPhiRight;
}


Real psiRightFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    Real const X = unitsFactor * x;
    Real const Y = unitsFactor * y;
    Real const Z = unitsFactor * z;

    Real const p00(1.31), p10(0.1024), p01(-1.12), p20(0.001014), p11(0.009993), p02(-0.04911);//2nd commit

    Real closedPsiRight = CLOSED_SECOND_DEGREE_POLY;

    Real const q00(6.07), q10(0.2668), q01(-0.4473), q20(-0.03993), q11(0.03803), q02(-0.03829);

    Real openPsiRight = OPEN_SECOND_DEGREE_POLY;

    return angleCoeff*openPsiRight + (1.-angleCoeff)*closedPsiRight;
}


Real phiLeftFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    Real const X = unitsFactor * x;
    Real const Y = unitsFactor * y;
    Real const Z = unitsFactor * z;

    Real const p00(-3.333), p10(-1.02), p01(-0.6707), p20(0.09525), p11(-0.03106), p02(0.04808);

    Real closedPhiLeft = CLOSED_SECOND_DEGREE_POLY;
    Real PlaneLeftRight = (lr00 + lr10*X + lr01*Y - Z) / ( unitsFactor * std::sqrt( lr10*lr10 + lr01*lr01 + 1.0 ) );
    Real PlaneLeftNon = (Z - ln00 - ln10*X - ln01*Y) / ( unitsFactor * std::sqrt( ln10*ln10 + ln01*ln01 + 1.0 ) );

    closedPhiLeft = std::min(closedPhiLeft, std::min(PlaneLeftRight, PlaneLeftNon) );


    Real const q00(7.312), q10(0.07152), q01(0.1407), q20(-0.09866), q11(0.1557), q02(-0.07457);
    Real const q30(0.003669), q21(0.0009942), q12(-0.003091), q03(-0.00318);
    Real const q40(0.001082), q31(-0.001874), q22(0.002587), q13(-0.001961), q04(0.0009986);
    Real const q50(-7.379e-6), q41(5.793e-5), q32(-0.0001357), q23(0.0001255), q14(-2.977e-5), q05(8.082e-5);
    Real const Xr = 0.711314*X - 0.453152*Y - 0.537295*Z;
    Real const Yr = -0.453152*X + 0.288686*Y - 0.843394*Z;
    Real const Zr = 0.537295*X + 0.843394*Y + 5.55112e-17*Z;

    Real openPhiLeft = OPEN_FIFTH_DEGREE_POLY;
    //if ( z < 0.7 )
    openPhiLeft = std::min(closedPhiLeft, openPhiLeft);

    return angleCoeff*openPhiLeft + (1.-angleCoeff)*closedPhiLeft;
}


Real psiLeftFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    Real const X = unitsFactor * x;
    Real const Y = unitsFactor * y;
    Real const Z = unitsFactor * z;

    Real const p00(1.644), p10(0.1344), p01(0.1294), p20( 0.008303), p11(-0.02842), p02(-0.01664);

    Real closedPsiLeft = CLOSED_SECOND_DEGREE_POLY;

    Real const q00(1.671), q10(0.3334), q01(0.1243), q20(-0.01749), q11(-0.01958), q02(0.);

    Real openPsiLeft = OPEN_SECOND_DEGREE_POLY;

    return angleCoeff*openPsiLeft + (1.-angleCoeff)*closedPsiLeft;
}



Real phiGlobalFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
//    Real const X = unitsFactor * x;
//    Real const Y = unitsFactor * y;
//    Real const Z = unitsFactor * z;

//    Real closedPhiGlobal, planeLeftRight, planeLeftNon, planeNonRight, openPhiGlobal;

//    {
//        //  closed phi-non
//        Real const p00(-2.449), p10(2.235), p01(-0.5299), p20(0.1401), p11(-0.008372), p02(0.1007);
//        closedPhiGlobal = CLOSED_SECOND_DEGREE_POLY;
//    }
//    {
//        //  closed phi-left
//        Real const p00(-3.333), p10(-1.02), p01(-0.6707), p20(0.09525), p11(-0.03106), p02(0.04808);
//        closedPhiGlobal = std::max( CLOSED_SECOND_DEGREE_POLY, closedPhiGlobal );
//    }
//    {
//        //  closed phi-right
//        Real const p00(-0.4305), p10(-0.7846), p01(1.767), p20(0.0778), p11(-0.05065), p02(0.1136);
//        closedPhiGlobal = std::max( CLOSED_SECOND_DEGREE_POLY, closedPhiGlobal );
//    }

//    planeLeftRight = std::abs( Z - lr00 - lr10*X - lr01*Y ) / ( unitsFactor * std::sqrt( lr10*lr10 + lr01*lr01 + 1.0 ) );
//    planeLeftNon = std::abs( Z - ln00 - ln10*X - ln01*Y ) / ( unitsFactor * std::sqrt( ln10*ln10 + ln01*ln01 + 1.0 ) );
//    planeNonRight = std::abs( Z - nr00 - nr10*X - nr01*Y ) / ( unitsFactor * std::sqrt( nr10*nr10 + nr01*nr01 + 1.0 ) );

//    closedPhiGlobal = std::min( closedPhiGlobal, std::min( planeLeftRight, std::min( planeLeftNon, planeNonRight ) ) );

//    openPhiGlobal = 0;

//    return angleCoeff * openPhiGlobal + (1.-angleCoeff) * closedPhiGlobal;
    return std::max( phiNonFct(0.,x,y,z,0), std::max( phiLeftFct(0.,x,y,z,0), phiRightFct(0.,x,y,z,0) ) );
}



Real psiGlobalFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    Real const X = unitsFactor * x;
    Real const Y = unitsFactor * y;
    Real const Z = unitsFactor * z;

    Real const p00(2.085), p10(-0.07436), p01(-0.315), p20(0.01783), p11(-0.01123), p02(0.02592);

    Real closedPsiGlobal = CLOSED_SECOND_DEGREE_POLY;

//    Real const q00(3.721), q10(0.03491), q01(-0.4287), q20(-0.01058), q11(-0.03066), q02(0.04796);
//    Real const q30(2.657e-5), q21(0.007418), q12(0.001337), q03(-0.0004671);
//    Real const q40(1.799e-5), q31(7.965e-5), q22(1.914e-5), q13(0.0004656), q04(-0.0004828);
//    Real const q50(7.847e-7), q41(-2.637e-5), q32(1.056e-5), q23(-2.567e-5), q14(-3.188e-5), q05(2.8e-5);
    Real const q00(4.248), q10(0.01216), q01(-0.3981), q20(-0.01619), q11(-0.03216), q02(0.01407);
    Real const q30(0.0005644), q21(0.00431), q12(0.0006411), q03(0.0002822);
    Real const q40(2.563e-5), q31(0.0001611), q22(0.0002301), q13(0.0003876), q04(-7.486e-5);
    Real const q50(-1.496e-6), q41(-8.491e-6), q32(-2.574e-6), q23(1.216e-5), q14(3.871e-6), q05(8.227e-6);
    Real const Xr = X;
    Real const Yr = Y;
    Real const Zr = Z;
    
    Real openPsiGlobal = OPEN_FIFTH_DEGREE_POLY;

    openPsiGlobal = std::max( (Z - 8.0)/unitsFactor, openPsiGlobal );

    return angleCoeff*openPsiGlobal + (1.-angleCoeff)*closedPsiGlobal;
}



Real valveInterfaceFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    if(useSmoothHeaviside)
    {
    // to be implemented
//            Real frac(value/epsilon*unitsFactor);
//            if (value >= epsilon/unitsFactor)
//                return 0.0;
//            else if (value <= -epsilon/unitsFactor)
//                return -1.0;
//            else
//                return -0.5*( 1. + frac + std::sin(M_PI*frac)/M_PI );
    }

    if ( psiGlobalFct(0, x, y, z, 1) <= 0 )
    {
        Real const frac( phiGlobalFct(0, x, y, z, 1) / epsilon * unitsFactor );
        if ( frac < -1 || frac > 1)
            return 0.0;
        else
            return 0.5 * (1 + std::cos(M_PI*frac) ) / epsilon * unitsFactor;
    }
    return 0.0;
}

