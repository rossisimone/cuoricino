#ifndef _PHYSIOLOGICALBC_HPP
#define _PHYSIOLOGICALBC_HPP

#include <lifev/core/LifeV.hpp>

#include "levelset_valve.hpp"

using namespace LifeV;

extern Real pRescaleFactor;
extern Real uRescaleFactor;
extern Real flowrateCorrection;
extern Real PinClosure;


Real zeroFct( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/);

Real aortaVelIn(const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i);

Real linearFluxIn(Real  t);

Real inletFixedPressureFct( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/);

Real inletPressureFct( const Real& t, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/);

Real outletFixedPressureFct( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/);

Real outletPressureFct( const Real& t, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/);

Real initPressureFct( const Real& t, const Real& x , const Real& y, const Real& z , const ID& /*i*/);

#endif
