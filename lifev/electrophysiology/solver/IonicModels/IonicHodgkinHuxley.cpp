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
  @brief Mean-field cross-bridge model by Negroni and Lascano 1996.
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#include <lifev/electrophysiology/solver/IonicModels/IonicHodgkinHuxley.hpp>


namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
//IonicMinimalModel::IonicMinimalModel()    :
//      super       ( 4     ),
//      M_uo        ( 0.    ),
//      M_uu        ( 1.61  ),
//      M_tetav     ( 0.3   ),
//      M_tetaw     ( 0.13  ),
//      M_tetavm    ( 0.1   ),
//      M_tetao     ( 0.005 ),
//      M_tauv1     ( 80.0  ),
//      M_tauv2     ( 1.4506),
//      M_tauvp     ( 1.4506),
//      M_tauw1     ( 70.0  ),
//      M_tauw2     ( 8.0   ),
//      M_kw        ( 200.0 ),
//      M_uw        ( 0.016 ),
//      M_tauwp     ( 280.0 ),
//      M_taufi     ( 0.078 ),
//      M_tauo1     ( 410.0 ),
//      M_tauo2     ( 7.0   ),
//      M_tauso1    ( 91.0  ),
//      M_tauso2    ( 0.8   ),
//      M_kso       ( 2.1   ),
//      M_uso       ( 0.6   ),
//      M_taus1     ( 2.7342),
//      M_taus2     ( 4.0   ),
//      M_ks        ( 2.0994),
//      M_us        ( 0.9087),
//      M_tausi     ( 3.3849),
//      M_tauwinf   ( 0.01  ),
//      M_winfstar  ( 0.5   )
//{
//}


IonicHodgkinHuxley::IonicHodgkinHuxley()  :
    super       ( 4, 3 ),
    M_gNa (120.0),
    M_gK (36.0),
    M_gL (0.3),
    M_vNa (115.0),
    M_vK (-12.0),
    M_vL (10.6)
{
    M_restingConditions.at (0) = 0.0;
    M_restingConditions.at (1) = 0.052932485257250;
    M_restingConditions.at (2) = 0.317676914060697;
    M_restingConditions.at (3) = 0.596120753508460;

}

IonicHodgkinHuxley::IonicHodgkinHuxley ( Teuchos::ParameterList& parameterList     )   :
    super       ( 4, 3 )
{
    M_gNa        =  parameterList.get ("gNa",      120.0     );
    M_gK        =  parameterList.get ("gK",     36.0     );
    M_gL        =  parameterList.get ("gL",      0.3     );
    M_vNa        =  parameterList.get ("vNa",      115.0     );
    M_vK        =  parameterList.get ("vK",      -12.0     );
    M_vL        =  parameterList.get ("vL",      10.6     );

}

IonicHodgkinHuxley::IonicHodgkinHuxley ( const IonicHodgkinHuxley& model )
{

    M_gNa        =  model.M_gNa;
    M_gK        =  model.M_gK;
    M_gL        =  model.M_gL;
    M_vNa        =  model.M_vNa;
    M_vK        =  model.M_vK;
    M_vL        =  model.M_vL;


    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;

}

// ===================================================
//! Operator
// ===================================================
IonicHodgkinHuxley& IonicHodgkinHuxley::operator= ( const IonicHodgkinHuxley& model )
{
    M_gNa        =  model.M_gNa;
    M_gK        =  model.M_gK;
    M_gL        =  model.M_gL;
    M_vNa        =  model.M_vNa;
    M_vK        =  model.M_vK;
    M_vL        =  model.M_vL;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;


    return      *this;
}


// ===================================================
//! Methods
// ===================================================
void IonicHodgkinHuxley::computeGatingRhs ( const   std::vector<Real>&  v,
                                            std::vector<Real>& rhs )
{

    Real V = v[0];
    Real M = v[1];
    Real N = v[2];
    Real H = v[3];

    Real alpham = 0.1 * (25. - V) / (std::exp ( (25. - V) / 10.) - 1.);
    Real betam = 4.*std::exp (-V / 18.0);
    Real alphah = 0.07 * std::exp (-V / 20.);
    Real betah = 1.0 / (std::exp ( (30. - V) / 10.) + 1.);
    Real alphan = 0.01 * (10. - V) / (std::exp ( (10. - V) / 10.) - 1.);
    Real betan = 0.125 * std::exp (-V / 80.0);

    rhs[0] = alpham * (1 - M) - betam * M;
    rhs[1] = alphan * (1 - N) - betan * N;
    rhs[2] = alphah * (1 - H) - betah * H;
}

void IonicHodgkinHuxley::computeRhs ( const   std::vector<Real>&  v,
                                      std::vector<Real>& rhs )
{

    Real V = v[0];
    Real M = v[1];
    Real N = v[2];
    Real H = v[3];

    Real alpham = 0.1 * (25. - V) / (std::exp ( (25. - V) / 10.) - 1.);
    Real betam = 4.*std::exp (-V / 18.0);
    Real alphah = 0.07 * std::exp (-V / 20.);
    Real betah = 1.0 / (std::exp ( (30. - V) / 10.) + 1.);
    Real alphan = 0.01 * (10. - V) / (std::exp ( (10. - V) / 10.) - 1.);
    Real betan = 0.125 * std::exp (-V / 80.0);

    rhs[0] = -M_gK * N * N * N * N * (V - M_vK) - M_gNa * M * M * M * H * (V - M_vNa) - M_gL * (V - M_vL); //+M_appliedCurrent;
    rhs[1] = alpham * (1 - M) - betam * M;
    rhs[2] = alphan * (1 - N) - betan * N;
    rhs[3] = alphah * (1 - H) - betah * H;
}


Real IonicHodgkinHuxley::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    Real dPotential (0.0);

    Real V = v[0];
    Real M = v[1];
    Real N = v[2];
    Real H = v[3];


    dPotential = -M_gK * N * N * N * N * (V - M_vK) - M_gNa * M * M * M * H * (V - M_vNa) - M_gL * (V - M_vL); //+M_appliedCurrent;

    return dPotential;
}

void IonicHodgkinHuxley::computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt )
{
    Real V = v[0];
    Real M = v[1];
    Real N = v[2];
    Real H = v[3];

    Real alpham = 0.1 * (25. - V) / (std::exp ( (25. - V) / 10.) - 1.);
    Real betam = 4.*std::exp (-V / 18.0);
    Real alphah = 0.07 * std::exp (-V / 20.);
    Real betah = 1.0 / (std::exp ( (30. - V) / 10.) + 1.);
    Real alphan = 0.01 * (10. - V) / (std::exp ( (10. - V) / 10.) - 1.);
    Real betan = 0.125 * std::exp (-V / 80.0);

    Real taum = alpham + betam;
    Real taun = alphan + betan;
    Real tauh = alphah + betah;

    Real mInf = alpham / (taum);
    Real nInf = alphan / (taun);
    Real hInf = alphah / (tauh);

    v[1] = mInf + (M - mInf) * std::exp (-dt * taum);
    v[2] = nInf + (N - nInf) * std::exp (-dt * taun);
    v[3] = hInf + (H - hInf) * std::exp (-dt * tauh);

}

void IonicHodgkinHuxley::showMe()
{

    std::cout << "\n\tHi, I'm the Hodgkin Huxley model for neurons. This is the first model implemented by Simone Palamara!!!\n\t See you soon\n\n";
}


}

