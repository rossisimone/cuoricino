/*
 * LuoRudy.hpp
 *
 *  Created on: 1-ago-2012
 *      Author: srossi
 */

#ifndef ELECTROLUORUDY_HPP_
#define ELECTROLUORUDY_HPP_


#include <lifev/electrophysiology/solver/ElectroIonicSolver.hpp>


namespace LifeV
{

template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class ElectroLuoRudy : public virtual ElectroIonicSolver<Mesh, SolverType>
{
public:
    typedef typename ElectroIonicSolver<Mesh, SolverType>::data_Type data_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::vector_Type vector_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::function_Type function_Type;

    ElectroLuoRudy ( const data_Type&          dataType,
                     const Mesh&          mesh,
                     FESpace<Mesh, MapEpetra>& uFEspace,
                     Epetra_Comm&              comm );

    virtual ~ElectroLuoRudy() {}

    void updateRepeated( );

    void updateElementSolution ( UInt eleID);

    void computeODECoefficients ( const Real& u_ig );

    void solveIonicModel ( const vector_Type& u, const Real timeStep );

    void computeIonicCurrent ( Real Capacitance,
                               VectorElemental& elvec,
                               VectorElemental& elvec_u,
                               FESpace<Mesh, MapEpetra>& uFESpace );

    void initialize ( );

    //! Returns the local solution vector for each field
    const vector_Type& solutionGatingH() const
    {
        return M_solutionGatingH;
    }

    const vector_Type& solutionGatingJ() const
    {
        return M_solutionGatingJ;
    }

    const vector_Type& solutionGatingM() const
    {
        return M_solutionGatingM;
    }

    const vector_Type& solutionGatingD() const
    {
        return M_solutionGatingD;
    }

    const vector_Type& solutionGatingF() const
    {
        return M_solutionGatingF;
    }

    const vector_Type& solutionGatingX() const
    {
        return M_solutionGatingX;
    }

    const vector_Type& solutionGatingCa() const
    {
        return M_solutionGatingCa;
    }

#ifdef ACTIVATED
    const vector_Type& solutionMechanicalActivation() const
    {
        return M_solutionMechanicalActivation;
    }
#endif

    Real M_K0, M_Ki, M_Na0, M_Nai, M_R, M_temperature, M_F, M_permeabilityRatio, M_c,
         M_Ena, M_Gk, M_Ek, M_Gk1, M_Ek1, M_Ekp, M_Esi, M_ah, M_bh, M_aj, M_bj, M_xii,
         M_am, M_bm, M_ad, M_bd, M_af, M_bf, M_aX, M_bX, M_ak1, M_bk1, M_Kp, M_K1inf,
         M_hinf, M_tauh, M_jinf, M_tauj, M_minf, M_taum, M_dinf, M_taud, M_finf, M_tauf,
         M_Xinf, M_tauX;

    //fast sodium current
    Real M_Ina;
    //slow inward current
    Real M_Islow;
    //time dependent potassium current
    Real M_Ik;
    //time independent potassium current
    Real M_Ik1;
    //plateau potassium current
    Real M_Ikp;
    //background current
    Real M_Iback;
    //Total time independent potassium current
    Real M_Ik1t;

    vector_Type M_vectorExponentialh;
    vector_Type M_vectorExponentialj;
    vector_Type M_vectorExponentialm;
    vector_Type M_vectorExponentiald;
    vector_Type M_vectorExponentialf;
    vector_Type M_vectorExponentialX;
    vector_Type M_vectorInfimumh;
    vector_Type M_vectorInfimumj;
    vector_Type M_vectorInfimumm;
    vector_Type M_vectorInfimumd;
    vector_Type M_vectorInfimumf;
    vector_Type M_vectorInfimumX;
    vector_Type M_vectorIonicChange;

#ifdef ACTIVATED
    vector_Type M_vectorActivationchange;
#endif

protected:
    //! Global solution h
    vector_Type                     M_solutionGatingH;
    //! Global solution j
    vector_Type                     M_solutionGatingJ;
    //! Global solution m
    vector_Type                     M_solutionGatingM;
    //! Global solution d
    vector_Type                     M_solutionGatingD;
    //! Global solution f
    vector_Type                     M_solutionGatingF;
    //! Global solution X
    vector_Type                     M_solutionGatingX;
    //! Global solution Ca_i
    vector_Type                     M_solutionGatingCa;

#ifdef ACTIVATED
    vector_Type                     M_solutionMechanicalActivation;
#endif

    vector_Type                     M_ionicCurrent;

    vector_Type                     M_ionicCurrentRepeated;

    VectorElemental                         M_elemVecIonicCurrent;

private:
};


//! Constructor
template<typename Mesh, typename SolverType>
ElectroLuoRudy<Mesh, SolverType>::ElectroLuoRudy ( const data_Type& dataType,
                                                   const Mesh& mesh,
                                                   FESpace<Mesh, MapEpetra>& uFEspace,
                                                   Epetra_Comm& comm ) :
    ElectroIonicSolver<Mesh, SolverType> ( dataType, mesh, uFEspace, comm),
    M_K0 (5.4),
    M_Ki (145.),
    M_Na0 (140.),
    M_Nai (18.),
    M_R (8.314472), //% Joules/(Kelvin*mole)
    M_temperature (307.7532), //%kelvins
    M_F (96485.33838), //% coulumbs/mole
    M_permeabilityRatio (0.01833), //% Na/K permeability ratio
    M_c (1.), //% membrane capacitance set as 1 mui-F/cm^2
    M_Ena (1000. * (M_R* M_temperature / M_F ) * log (M_Na0 / M_Nai) ),
    M_Gk (0.282 * sqrt (M_K0 / 5.4) ),
    M_Ek (1000. * (M_R* M_temperature / M_F) * log ( (M_K0 + M_permeabilityRatio* M_Na0)
                                                     / (M_Ki + M_permeabilityRatio* M_Nai) ) ),
    M_Gk1 (0.6047 * sqrt (M_K0 / 5.4) ),
    M_Ek1 (1000.* (M_R* M_temperature / M_F) * log (M_K0 / M_Ki) ),
    M_Ekp (M_Ek1),
    M_vectorExponentialh (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialj (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialm (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentiald (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialf (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialX (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumh (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumj (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumm (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumd (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumf (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumX (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorIonicChange (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_solutionGatingH                  ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingJ                  ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingM                  ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingD                  ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingF                  ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingX                  ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingCa                  ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
#ifdef ACTIVATED
    M_solutionMechanicalActivation ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorActivationchange ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
#endif
    M_ionicCurrent ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_ionicCurrentRepeated ( M_ionicCurrent, Repeated ),
    M_elemVecIonicCurrent ( ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 )
{
}

template<typename Mesh, typename SolverType>
void ElectroLuoRudy<Mesh, SolverType>::updateRepeated( )
{
    M_ionicCurrentRepeated = M_ionicCurrent;
}

template<typename Mesh, typename SolverType>
void ElectroLuoRudy<Mesh, SolverType>::updateElementSolution ( UInt eleID)
{
    M_elemVecIonicCurrent.zero();
    UInt ig;
    //! Filling local elvec with recovery variable values in the nodes
    for ( UInt iNode = 0 ; iNode < ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobalMap ( eleID, iNode );
        M_elemVecIonicCurrent.vec() [ iNode ] = M_ionicCurrentRepeated[ig];
    }

}


template<typename Mesh, typename SolverType>
void ElectroLuoRudy<Mesh, SolverType>::solveIonicModel ( const vector_Type& u, const Real timeStep )
{
    //! Solving dw/dt=eta2 (u/vp -  eta3 w)
    LifeChrono chronoionmodelsolve;
    chronoionmodelsolve.start();
    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

    for ( Int i = 0 ; i < u.epetraVector().MyLength() ; i++ )
    {
        Int ig = u.blockMap().MyGlobalElements() [i];
        Real u_ig = u[ig];
        computeODECoefficients (u_ig);
        M_Esi = 7.7 - 13.0287 * log (M_solutionGatingCa[ig]);
        //fast sodium current
        M_Ina = 23.* M_solutionGatingM[ig] * M_solutionGatingM[ig] * M_solutionGatingM[ig] * M_solutionGatingH[ig] * M_solutionGatingJ[ig] * (u_ig - M_Ena);
        //slow inward current
        M_Islow = 0.09 * M_solutionGatingD[ig] * M_solutionGatingF[ig] * (u_ig - M_Esi);

        //Mechanical activation change
#ifdef ACTIVATED
        M_vectorActivationchange.epetraVector().ReplaceGlobalValue (ig,
                                                                    0,
                                                                    -0.66 * M_solutionGatingCa[ig] - 0.011 * M_solutionMechanicalActivation[ig]);
#endif
        //change in ioniq concentration
        M_vectorIonicChange.epetraVector().ReplaceGlobalValue (ig,
                                                               0,
                                                               -1e-4 * M_Islow + 0.07 * (1e-4 - M_solutionGatingCa[ig]) );
        //time dependent potassium current
        M_Ik = M_Gk * M_solutionGatingX[ig] * M_xii * (u_ig - M_Ek);
        //time independent potassium current
        M_Ik1 = M_Gk1 * M_K1inf * (u_ig - M_Ek1);
        //plateau potassium current
        M_Ikp = 0.0183 * M_Kp * (u_ig - M_Ekp);
        //background current
        M_Iback = 0.03921 * (u_ig + 59.87);
        //Total time independent potassium current
        M_Ik1t = M_Ik1 + M_Ikp + M_Iback;
        // adding up the six ionic currents
        M_ionicCurrent.epetraVector().ReplaceGlobalValue (ig,
                                                          0,
                                                          M_Ina + M_Islow + M_Ik + M_Ik1t);
        M_vectorExponentialh.epetraVector().ReplaceGlobalValue (ig,
                                                                0,
                                                                exp (-timeStep / M_tauh) );
        M_vectorExponentialj.epetraVector().ReplaceGlobalValue (ig,
                                                                0,
                                                                exp (-timeStep / M_tauj) );
        M_vectorExponentialm.epetraVector().ReplaceGlobalValue (ig,
                                                                0,
                                                                exp (-timeStep / M_taum) );
        M_vectorExponentiald.epetraVector().ReplaceGlobalValue (ig,
                                                                0,
                                                                exp (-timeStep / M_taud) );
        M_vectorExponentialf.epetraVector().ReplaceGlobalValue (ig,
                                                                0,
                                                                exp (-timeStep / M_tauf) );
        M_vectorExponentialX.epetraVector().ReplaceGlobalValue (ig,
                                                                0,
                                                                exp (-timeStep / M_tauX) );
        M_vectorInfimumh.epetraVector().ReplaceGlobalValue (ig,
                                                            0,
                                                            M_hinf);
        M_vectorInfimumj.epetraVector().ReplaceGlobalValue (ig,
                                                            0,
                                                            M_jinf);
        M_vectorInfimumm.epetraVector().ReplaceGlobalValue (ig,
                                                            0,
                                                            M_minf);
        M_vectorInfimumd.epetraVector().ReplaceGlobalValue (ig,
                                                            0,
                                                            M_dinf);
        M_vectorInfimumf.epetraVector().ReplaceGlobalValue (ig,
                                                            0,
                                                            M_finf);
        M_vectorInfimumX.epetraVector().ReplaceGlobalValue (ig,
                                                            0,
                                                            M_Xinf);
    }
    M_vectorExponentialh.globalAssemble();
    M_vectorExponentialj.globalAssemble();
    M_vectorExponentialm.globalAssemble();
    M_vectorExponentiald.globalAssemble();
    M_vectorExponentialf.globalAssemble();
    M_vectorExponentialX.globalAssemble();
    M_vectorInfimumh.globalAssemble();
    M_vectorInfimumj.globalAssemble();
    M_vectorInfimumm.globalAssemble();
    M_vectorInfimumd.globalAssemble();
    M_vectorInfimumf.globalAssemble();
    M_vectorInfimumX.globalAssemble();
    M_ionicCurrent.globalAssemble();
    M_vectorIonicChange.globalAssemble();

#ifdef ACTIVATED
    M_vectorActivationchange.globalAssemble();
#endif

    M_solutionGatingH -= M_vectorInfimumh;
    M_solutionGatingH.epetraVector().Multiply (1.,
                                               M_solutionGatingH.epetraVector(),
                                               M_vectorExponentialh.epetraVector(),
                                               0.);
    M_solutionGatingH += M_vectorInfimumh;
    M_solutionGatingJ -= M_vectorInfimumj;

    M_solutionGatingJ.epetraVector().Multiply (1.,
                                               M_solutionGatingJ.epetraVector(),
                                               M_vectorExponentialj.epetraVector(),
                                               0.);
    M_solutionGatingJ += M_vectorInfimumj;
    M_solutionGatingM -= M_vectorInfimumm;
    M_solutionGatingM.epetraVector().Multiply (1.,
                                               M_solutionGatingM.epetraVector(),
                                               M_vectorExponentialm.epetraVector(),
                                               0.);
    M_solutionGatingM += M_vectorInfimumm;
    M_solutionGatingD -= M_vectorInfimumd;

    M_solutionGatingD.epetraVector().Multiply (1.,
                                               M_solutionGatingD.epetraVector(),
                                               M_vectorExponentiald.epetraVector(),
                                               0.);
    M_solutionGatingD += M_vectorInfimumd;
    M_solutionGatingF -= M_vectorInfimumf;

    M_solutionGatingF.epetraVector().Multiply (1.,
                                               M_solutionGatingF.epetraVector(),
                                               M_vectorExponentialf.epetraVector(),
                                               0.);
    M_solutionGatingF += M_vectorInfimumf;
    M_solutionGatingX -= M_vectorInfimumX;

    M_solutionGatingX.epetraVector().Multiply (1.,
                                               M_solutionGatingX.epetraVector(),
                                               M_vectorExponentialX.epetraVector(),
                                               0.);
    M_solutionGatingX += M_vectorInfimumX;
    M_solutionGatingCa += timeStep * M_vectorIonicChange;

#ifdef ACTIVATED
    M_solutionMechanicalActivation += timeStep * M_vectorActivationchange;
#endif
    M_solutionGatingH.globalAssemble();
    M_solutionGatingJ.globalAssemble();
    M_solutionGatingM.globalAssemble();
    M_solutionGatingD.globalAssemble();
    M_solutionGatingF.globalAssemble();
    M_solutionGatingX.globalAssemble();
    M_solutionGatingCa.globalAssemble();

#ifdef ACTIVATED
    M_solutionMechanicalActivation.globalAssemble();
#endif
    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

    chronoionmodelsolve.stop();
    if (ElectroIonicSolver<Mesh, SolverType>::M_comm->MyPID() == 0)
    {
        std::cout << "Total ionmodelsolve time " << chronoionmodelsolve.diff() << " s." << std::endl;
    }
}



template<typename Mesh, typename SolverType>
void ElectroLuoRudy<Mesh, SolverType>::computeODECoefficients ( const Real& u_ig )
{
    if (u_ig >= -40.)
    {
        M_ah = 0.;
        M_bh = 1. / (0.13 * (1. + exp ( (u_ig + 10.66) / (-11.1) ) ) );
        M_aj = 0.;
        M_bj = 0.3 * exp (-2.535e-7 * u_ig) / (1. + exp (-0.1 * (u_ig + 32.) ) );
    }
    else
    {
        M_ah = 0.135 * exp ( (80. + u_ig) / -6.8);
        M_bh = 3.56 * exp (0.079 * u_ig) + 3.1e5 * exp (0.35 * u_ig);
        M_aj = (-1.2714e5 * exp (0.2444 * u_ig) - 3.474e-5 * exp (-0.04391 * u_ig) ) *
               (u_ig + 37.78) / (1 + exp (0.311 * (u_ig + 79.23) ) );
        M_bj = 0.1212 * exp (-0.01052 * u_ig) / (1. + exp (-0.1378 * (u_ig + 40.14) ) );
    }
    M_am = 0.32 * (u_ig + 47.13) / (1. - exp (-0.1 * (u_ig + 47.13) ) );
    M_bm = 0.08 * exp (-u_ig / 11.);

    //slow inward current
    M_ad = 0.095 * exp (-0.01 * (u_ig - 5.) ) / (1. + exp (-0.072 * (u_ig - 5.) ) );
    M_bd = 0.07  * exp (-0.017 * (u_ig + 44.) ) / (1. + exp ( 0.05 * (u_ig + 44.) ) );
    M_af = 0.012 * exp (-0.008 * (u_ig + 28.) ) / (1. + exp ( 0.15 * (u_ig + 28.) ) );
    M_bf = 0.0065 * exp (-0.02 * (u_ig + 30.) ) / (1. + exp ( -0.2 * (u_ig + 30.) ) );

    //Time dependent potassium outward current
    M_aX = 0.0005 * exp (0.083 * (u_ig + 50.) ) / (1. + exp (0.057 * (u_ig + 50.) ) );
    M_bX = 0.0013 * exp (-0.06 * (u_ig + 20.) ) / (1. + exp (-0.04 * (u_ig + 20.) ) );

    if (u_ig <= -100)
    {
        M_xii = 1.;
    }
    else
    {
        M_xii = 2.837 * (exp (0.04 * (u_ig + 77.) ) - 1.) / ( (u_ig + 77.) * exp (0.04 * (u_ig + 35.) ) );
    }
    M_ak1 = 1.02 / (1. + exp (0.2385 * (u_ig - M_Ek1 - 59.215) ) );
    M_bk1 = (0.49124 * exp (0.08032 * (u_ig - M_Ek1 + 5.476) ) +
             exp (0.06175 * (u_ig - M_Ek1 - 594.31) ) ) / (1. + exp (-0.5143 * (u_ig - M_Ek1 + 4.753) ) );
    //Plateau potassium outward current
    M_Kp = 1. / (1. + exp ( (7.488 - u_ig) / 5.98) );

    M_K1inf = M_ak1 / (M_ak1 + M_bk1);
    M_hinf = M_ah   / (M_ah  + M_bh);
    M_tauh = 1.   / (M_ah + M_bh);
    M_jinf = M_aj / (M_aj + M_bj);
    M_tauj = 1.   / (M_aj + M_bj);
    M_minf = M_am / (M_am + M_bm);
    M_taum = 1.   / (M_am + M_bm);
    M_dinf = M_ad / (M_ad + M_bd);
    M_taud = 1.   / (M_ad + M_bd);
    M_finf = M_af / (M_af + M_bf);
    M_tauf = 1.   / (M_af + M_bf);
    M_Xinf = M_aX / (M_aX + M_bX);
    M_tauX = 1.   / (M_aX + M_bX);
}

template<typename Mesh, typename SolverType>
void ElectroLuoRudy<Mesh, SolverType>::computeIonicCurrent (  Real Capacitance,
                                                              VectorElemental& elvec,
                                                              VectorElemental& /*elvec_u*/,
                                                              FESpace<Mesh, MapEpetra>& uFESpace )
{
    Real Iion_ig;
    for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
    {
        Iion_ig = 0.;
        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            Iion_ig += M_elemVecIonicCurrent ( i ) * uFESpace.fe().phi ( i, ig );
        }
        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            // divide by 1000 to convert microA in mA
            elvec ( i ) -= Iion_ig * Capacitance *
                           uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
        }
    }
}

template<typename Mesh, typename SolverType>
void ElectroLuoRudy<Mesh, SolverType>::
initialize( )
{
    M_solutionGatingH.epetraVector().PutScalar (1.);
    M_solutionGatingJ.epetraVector().PutScalar (1.);
    M_solutionGatingM.epetraVector().PutScalar (0.);
    M_solutionGatingD.epetraVector().PutScalar (0.);
    M_solutionGatingF.epetraVector().PutScalar (1.);
    M_solutionGatingX.epetraVector().PutScalar (0.);
    M_solutionGatingCa.epetraVector().PutScalar (0.0002);

#ifdef ACTIVATED
    M_solutionMechanicalActivation.epetraVector().PutScalar (0.0);
#endif

    M_solutionGatingH.globalAssemble();
    M_solutionGatingJ.globalAssemble();
    M_solutionGatingM.globalAssemble();
    M_solutionGatingD.globalAssemble();
    M_solutionGatingF.globalAssemble();
    M_solutionGatingX.globalAssemble();
    M_solutionGatingCa.globalAssemble();

#ifdef ACTIVATED
    M_solutionMechanicalActivation.globalAssemble();
#endif
}


}

#endif /* LUORUDYI_HPP_ */
