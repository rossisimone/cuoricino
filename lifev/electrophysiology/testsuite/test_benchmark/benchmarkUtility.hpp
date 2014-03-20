/*
 * benchmarkUtility.hpp
 *
 *  Created on: Mar 19, 2014
 *      Author: srossi
 */

#ifndef BENCHMARKUTILITY_HPP_
#define BENCHMARKUTILITY_HPP_

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"




#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicLuoRudyI.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher06.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicHodgkinHuxley.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicNoblePurkinje.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicFox.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>

#include <lifev/core/LifeV.hpp>

namespace LifeV{

namespace BenchmarkUtility{


typedef ElectroIonicModel                                        ionicModel_Type;
typedef boost::shared_ptr<ionicModel_Type>                       ionicModelPtr_Type;
typedef boost::function < Real (const Real& t,
                                 const Real& x,
                                 const Real& y,
                                 const Real& z,
                                 const ID&   i ) >   function_Type;

Real chooseIonicModel(ionicModelPtr_Type& model, std::string ionic_model, Epetra_Comm& Comm )
{
	Real activationThreshold(0.95);

    if ( Comm.MyPID() == 0 )
    {
        std::cout << "Constructing the ionic model ... !"; // << std::endl;
    }

    if ( ionic_model == "AlievPanfilov" )
    {
        model.reset ( new IonicAlievPanfilov() );
    }
    if ( ionic_model == "LuoRudyI" )
    {
        model.reset ( new IonicLuoRudyI() );
    }
    if ( ionic_model == "TenTusscher06")
    {
        model.reset (new IonicTenTusscher06() );
    }
    if ( ionic_model == "HodgkinHuxley")
    {
        model.reset (new IonicHodgkinHuxley() );
        activationThreshold = 10.0;
    }
    if ( ionic_model == "NoblePurkinje")
    {
        model.reset (new IonicNoblePurkinje() );
    }
    if ( ionic_model == "MinimalModel")
    {
        model.reset ( new IonicMinimalModel() );
    }
    if ( ionic_model == "Fox")
    {
//        model.reset ( new IonicFox() );
        if ( Comm.MyPID() == 0 )
        {
            assert(0 && "Fox model is not properly working in 3D."); //TO DO: Fix It!
        }
    }

    if ( Comm.MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    model -> showMe();

    return activationThreshold;
}




Real PacingProtocolMM ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 10;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

Real PacingProtocolHH ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 500.;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

Real PacingProtocol ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 40.0;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}


void setStimulus(function_Type& f, std::string ionic_model)
{
    if (ionic_model == "MinimalModel" || ionic_model == "AlievPanfilov")
    {
        f = &BenchmarkUtility::PacingProtocolMM;
    }
    else if (ionic_model == "HodgkinHuxley" )
    {
        f = &BenchmarkUtility::PacingProtocolHH;
    }
    else
    {
        f = &BenchmarkUtility::PacingProtocol;
    }
}

}//BenchmarkUtility

} //LifeV


#endif /* BENCHMARKUTILITY_HPP_ */
