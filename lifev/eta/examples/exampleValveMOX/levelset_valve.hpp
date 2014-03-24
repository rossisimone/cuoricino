#ifndef _LEVELSETVALVE_HPP
#define _LEVELSETVALVE_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorSmall.hpp>

using namespace LifeV;


extern Real angleCoeff;

extern Real unitsFactor;

extern Real epsilon;

extern bool useSmoothHeaviside;



Real phiNonFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real psiNonFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real phiRightFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real psiRightFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real phiLeftFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real psiLeftFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real phiGlobalFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real psiGlobalFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

Real valveInterfaceFct( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/);


class SmoothHeavisideFct
{
public:
    typedef Real return_Type;

    return_Type operator()(const Real& value)
    {
        if (useSmoothHeaviside)
        {
            Real frac(value/epsilon*unitsFactor);

            if (value >= epsilon/unitsFactor)
                return 1.0;
            else if (value <= -epsilon/unitsFactor)
                return 0.0;
            else
                return 0.5*( 1. + frac + std::sin(M_PI*frac)/M_PI );
        }
        if (value >= 0)
            return 1.0;
        return 0.0;
    }

    SmoothHeavisideFct(){}
    SmoothHeavisideFct(const SmoothHeavisideFct&){}
    ~SmoothHeavisideFct(){}
};


class SmoothDeltaFct
{
public:
    typedef Real return_Type;

    return_Type operator()(const Real& value)
    {
        Real frac(value/epsilon*unitsFactor);

        if (frac<-1)
            return 0.0;
        else if (frac > 1)
            return 0.0;
        else
            return 0.5 * (1 + std::cos(M_PI*frac) ) / epsilon * unitsFactor;
    }

    SmoothDeltaFct(){}
    SmoothDeltaFct(const SmoothDeltaFct&){}
    ~SmoothDeltaFct(){}
};


class ValveInterfaceFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return valveInterfaceFct( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2], 0  ) ;
    }

    ValveInterfaceFunctor() {}
    ValveInterfaceFunctor (const ValveInterfaceFunctor&) {}
    ~ValveInterfaceFunctor() {}
};

#endif
