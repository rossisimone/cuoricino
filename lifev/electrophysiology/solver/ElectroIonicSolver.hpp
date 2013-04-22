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
  @brief Class for choosing and solving ionic models in electrophysiology.
  @Mitchel-Schaeffer, Rogers-McCulloch, Luo-Rudy I, Bueno-Orovio et al., and Courtemanche et al.
  @All ventricular models include the coupling with the mechanical activation in the fibers direction

  @date 11-2007
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>

  @contributors J.Castelneau (INRIA), Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>, Simon Labarthe <simon.labarthe@u-bordeaux2.fr>
  @mantainer Ricardo Ruiz Baier <ricardo.ruiz@epfl.ch>
  @last update 02-2012
 */

#ifndef _IONICSOLVER_H_
#define _IONICSOLVER_H_

#define ACTIVATED

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/electrophysiology/solver/ElectroIonicData.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <boost/shared_ptr.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

namespace LifeV
{
//! IonicSolver - This class implements a ionic model solver.


template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class ElectroIonicSolver
{

public:
    //! @name Type definitions
    //@{

    typedef ElectroIonicData data_Type;

    typedef Real ( *function_Type ) ( const Real&,
                                      const Real&,
                                      const Real&,
                                      const Real&,
                                      const ID& );

    typedef boost::function < Real (const markerID_Type& ref,
                                    const Real& x,
                                    const Real& y,
                                    const Real& z,
                                    const ID& i) > functorTauClose_Type;

    typedef Mesh mesh_Type;

    typedef typename SolverType::matrix_type      matrix_Type;
    typedef typename SolverType::vector_type      vector_Type;

    typedef typename SolverType::prec_raw_type    precRaw_Type;
    typedef typename SolverType::prec_type        prec_Type;

    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param dataType
     * @param mesh
     * @param recovery FE space
     * @param Epetra communicator
     */
    ElectroIonicSolver ( const data_Type& dataType,
                       const Mesh& mesh,
                       FESpace<Mesh, MapEpetra>& uFEspace,
                       Epetra_Comm& comm );

    //! Destructor
    virtual ~ElectroIonicSolver() {}

    //@}

    //! @name Methods
    //@{


    //! Returns the recovery variable FE space
    FESpace<Mesh, MapEpetra>& recoveryFESpace()
    {
        return M_uFESpace;
    }

    //! Return maps
    Epetra_Map const& getRepeatedMapEpetra() const
    {
        return *M_localMap.map (Repeated);
    }

    MapEpetra const& getMap() const
    {
        return M_localMap;
    }

    virtual void updateRepeated( ) = 0;

    //! Update the ionic model elvecs
    virtual void updateElementSolution ( UInt eleIDw ) = 0;

    //! Solves the ionic model
    virtual void solveIonicModel ( const vector_Type& u,
                                   const Real timeStep ) = 0;

    //! Computes the term -1/ \Cm u^n (G (1-u^n/vp) (1-u^n/v_th) + eta_1 v^{n+1})
    //! for the PDE righthand side
    virtual void computeIonicCurrent ( Real Capacitance,
                                       VectorElemental& elvec,
                                       VectorElemental& elvec_u,
                                       FESpace<Mesh, MapEpetra>& uFESpace ) = 0;

    //! Initialize
    virtual void initialize( ) = 0;

    const data_Type&               M_data;

    const Mesh&                    M_mesh;

    // FE space
    FESpace<Mesh, MapEpetra>&      M_uFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;

    Int                            M_me;

    //! Map
    MapEpetra                      M_localMap;

    //! Boolean that indicates if output is sent to cout
    bool                           M_verbose;

    //@}

protected:

    UInt solutionUDimension() const
    {
        return M_uFESpace.dim();
    }

}; // class IonicSolver


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
ElectroIonicSolver<Mesh, SolverType>::
ElectroIonicSolver ( const data_Type& dataType,
                   const Mesh& mesh,
                   FESpace<Mesh, MapEpetra>& uFEspace,
                   Epetra_Comm& comm ) :
    M_data                   ( dataType ),
    M_mesh                   ( mesh ),
    M_uFESpace               ( uFEspace ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_localMap               ( M_uFESpace.map() ),
    M_verbose                ( M_me == 0)
{
}

template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class MitchellSchaeffer : public virtual ElectroIonicSolver<Mesh, SolverType>
{
public:
    typedef typename ElectroIonicSolver<Mesh, SolverType>::data_Type  data_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::vector_Type    vector_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::function_Type  function_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::functorTauClose_Type  functorTauClose_Type;

    MitchellSchaeffer ( const data_Type& dataType,
                        const Mesh& mesh,
                        FESpace<Mesh, MapEpetra>& uFEspace,
                        Epetra_Comm& comm );

    virtual ~MitchellSchaeffer();

    void updateRepeated( );

    void updateElementSolution ( UInt eleID );

    void solveIonicModel ( const vector_Type& u, const Real timeStep );

    void computeIonicCurrent ( Real Capacitance,
                               VectorElemental& elvec,
                               VectorElemental& elvec_u,
                               FESpace<Mesh, MapEpetra>& uFESpace );

    const vector_Type& solutionGatingW() const
    {
        return M_solutionGatingW;
    }

    void initialize( );

    void setHeteroTauClose (functorTauClose_Type);

    Real functorTauClose (const markerID_Type& ref,
                          const Real& x,
                          const Real& y,
                          const Real& z,
                          const ID& i) const;

protected:

    //! Global solution _w
    vector_Type             M_solutionGatingW;
    vector_Type             M_solutionGatingWRepeated;
    VectorElemental                 M_elvec;
    UInt                    M_BDForder;
    TimeAdvanceBDF<vector_Type>         M_BDFW;
    functorTauClose_Type    M_tauClose;

private:
};


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
MitchellSchaeffer<Mesh, SolverType>::
MitchellSchaeffer ( const data_Type& dataType,
                    const Mesh& mesh,
                    FESpace<Mesh, MapEpetra>& uFEspace,
                    Epetra_Comm& comm ) :
    ElectroIonicSolver<Mesh, SolverType> ( dataType, mesh, uFEspace, comm),
    M_solutionGatingW ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingWRepeated ( M_solutionGatingW, Repeated ),
    M_elvec ( ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 ),
    M_BDForder ( ElectroIonicSolver<Mesh, SolverType>::M_data.MSBDForder() ),
    M_BDFW( )
{
    M_BDFW.setup ( M_BDForder );
}

template<typename Mesh, typename SolverType>
MitchellSchaeffer<Mesh, SolverType>::
~MitchellSchaeffer()
{
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::updateRepeated( )
{
    M_solutionGatingWRepeated = M_solutionGatingW;
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::updateElementSolution ( UInt eleID )
{
    M_elvec.zero();
    UInt ig;
    //! Filling local elvec_w with recovery variable values in the nodes
    for ( UInt iNode = 0 ; iNode < ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobalMap ( eleID, iNode );
        M_elvec.vec() [ iNode ] = M_solutionGatingWRepeated[ig];
    }
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::setHeteroTauClose (functorTauClose_Type fct)
{
    M_tauClose = fct;
}

template<typename Mesh, typename SolverType>
Real MitchellSchaeffer<Mesh, SolverType>::functorTauClose (const markerID_Type& ref,
                                                           const Real& x,
                                                           const Real& y,
                                                           const Real& z, const ID& i) const
{
    return M_tauClose (ref, x, y, z, i);
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::solveIonicModel ( const vector_Type& u, const Real timeStep )
{
    //! Solving :
    //!           ((v_max-v_min)^{-2}  - w )/tau_open  if u < vcrit
    //! dw/dt ={
    //!            -w/tau_close   if u > vcrit

    Real aux1 = 1.0 / (M_BDFW.coefficientFirstDerivative (0) / timeStep +
                       1.0 / this->M_data.MSTauOpen() );
    Real aux = 1.0 / ( (this->M_data.MSPotentialMaximum() -
                        this->M_data.MSPotentialMinimum() ) *
                       (this->M_data.MSPotentialMaximum() -
                        this->M_data.MSPotentialMinimum() ) *
                       this->M_data.MSTauOpen() );
    Real aux2 = 1.0 / (M_BDFW.coefficientFirstDerivative (0) / timeStep +
                       1.0 / this->M_data.MSTauClose() );

    M_BDFW.updateRHSContribution (timeStep);
    vector_Type M_time_der = M_BDFW.rhsContributionFirstDerivative();

    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();
    UInt ID;
    for ( Int i = 0 ; i < u.epetraVector().MyLength() ; ++i )
    {
        Int ig = u.blockMap().MyGlobalElements() [i];
        ID  = ig;
        /*ref = this->M_mesh->point(ig).markerID();
        x   = this->M_mesh->point(ig).x();
        y   = this->M_mesh->point(ig).y();
        z   = this->M_mesh->point(ig).z();*/
        if (u[ig] < this->M_data.MSCriticalPotential() )
        {
            M_solutionGatingW[ig] = aux1 * (aux + M_time_der[ig]);
        }
        else if (this->M_data.MSHasHeterogeneousTauClose() )
            M_solutionGatingW[ig] = (1.0 / (M_BDFW.coefficientFirstDerivative (0) / timeStep  +
                                            1.0/*fct_Tau_Close(ref,x,y,z,ID)*/) ) *
                                    M_time_der[ig];//aux2 * M_time_der[ig];
        else
        {
            M_solutionGatingW[ig] = aux2 *  M_time_der[ig];
        }
    }
    M_BDFW.shiftRight (M_solutionGatingW);

    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::computeIonicCurrent (  Real,
                                                                 VectorElemental& elvec,
                                                                 VectorElemental& elvec_u,
                                                                 FESpace<Mesh, MapEpetra>& uFESpace )
{
    for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
    {
        elvec ( i ) =  this->M_data.MSReactionAmplitude() *
                       ( ( (M_elvec ( i ) / this->M_data.MSTauIn() ) *
                           (elvec_u ( i ) - this->M_data.MSPotentialMinimum() ) *
                           (elvec_u ( i ) - this->M_data.MSPotentialMinimum() ) *
                           (this->M_data.MSPotentialMaximum() - elvec_u ( i ) ) /
                           (this->M_data.MSPotentialMaximum() - this->M_data.MSPotentialMinimum() ) ) -
                         ( ( elvec_u ( i ) - this->M_data.MSPotentialMinimum() ) /
                           (  this->M_data.MSTauOut() * (this->M_data.MSPotentialMaximum() -
                                                         this->M_data.MSPotentialMinimum() ) ) ) ) ;
    }
}

template<typename Mesh, typename SolverType>
void MitchellSchaeffer<Mesh, SolverType>::
initialize( )
{
    M_solutionGatingW.epetraVector().PutScalar (1.0 /
                                                ( (this->M_data.MSPotentialMaximum() -
                                                   this->M_data.MSPotentialMinimum() ) *
                                                  (this->M_data.MSPotentialMaximum() -
                                                   this->M_data.MSPotentialMinimum() ) ) );
    M_BDFW.setInitialCondition (M_solutionGatingW);
    M_BDFW.showMe();
}

/////////////////////////////////////////////////////////////////////////
//Phenomenological model for the electrophysiology of the heart
//Parameters and methods as in Rogers-McCulloch 98.
////////////////////////////////////////////////////////////////////////


template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class RogersMcCulloch : public virtual ElectroIonicSolver<Mesh, SolverType>
{
public:
    typedef typename ElectroIonicSolver<Mesh, SolverType>::data_Type data_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::vector_Type vector_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::function_Type function_Type;

    RogersMcCulloch ( const data_Type& dataType,
                      const Mesh& mesh,
                      FESpace<Mesh, MapEpetra>& uFEspace,
                      Epetra_Comm& comm );
    virtual ~RogersMcCulloch();

    void updateRepeated( );

    void updateElementSolution ( UInt eleID );

    void solveIonicModel ( const vector_Type& u, const Real timeStep );

    void computeIonicCurrent ( Real Capacitance,
                               VectorElemental& elvec,
                               VectorElemental& elvec_u,
                               FESpace<Mesh, MapEpetra>& uFESpace );

    const vector_Type& solutionGatingW() const
    {
        return M_solutionGatingW;
    }

#ifdef ACTIVATED
    const vector_Type& solutionMechanicalActivation() const
    {
        return M_solutionMechanicalActivation;
    }
    vector_Type M_vectorActivationchange;
#endif

    void initialize( );



protected:

    //! Global solution _w
    vector_Type                       M_solutionGatingW;

#ifdef ACTIVATED
    vector_Type                       M_solutionMechanicalActivation;
#endif

    vector_Type                       M_solutionGatingWRepeated;

    VectorElemental                       M_elvec;
private:
};


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
RogersMcCulloch<Mesh, SolverType>::
RogersMcCulloch ( const data_Type& dataType,
                  const Mesh& mesh,
                  FESpace<Mesh, MapEpetra>& uFEspace,
                  Epetra_Comm& comm ) :
    ElectroIonicSolver<Mesh, SolverType> ( dataType, mesh, uFEspace, comm),
    M_solutionGatingW ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
#ifdef ACTIVATED
    M_solutionMechanicalActivation ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorActivationchange ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
#endif
    M_solutionGatingWRepeated ( M_solutionGatingW, Repeated ),
    M_elvec ( ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 )
{
}

template<typename Mesh, typename SolverType>
RogersMcCulloch<Mesh, SolverType>::
~RogersMcCulloch()
{
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::updateRepeated( )
{
    M_solutionGatingWRepeated = M_solutionGatingW;
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::updateElementSolution ( UInt eleID )
{
    M_elvec.zero();
    UInt ig;
    for ( UInt iNode = 0 ;
            iNode < ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobalMap ( eleID, iNode );
        M_elvec.vec() [ iNode ] = M_solutionGatingWRepeated[ig];
    }
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::solveIonicModel ( const vector_Type& u, const Real timeStep )
{
    //! Solving dw/dt= b/(A*T) (u - u0 - A d w)
    Real G = this->M_data.RMCParameterB() / this->M_data.RMCPotentialAmplitude() / this->M_data.RMCTimeUnit();

    Real alpha = 1 / timeStep + G * this -> M_data.RMCPotentialAmplitude() * this -> M_data.RMCParameterD() ;

    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

    M_solutionGatingW *= 1 / timeStep;
    vector_Type temp (u);
    temp.epetraVector().PutScalar ( G * this -> M_data.RMCRestPotential() );
    M_solutionGatingW += G * u;
    M_solutionGatingW -= temp;
    M_solutionGatingW *= 1 / alpha;
    M_solutionGatingW.globalAssemble();

    //Evolution of mechanical activation
#ifdef ACTIVATED
    for ( Int i = 0 ; i < u.epetraVector().MyLength() ; i++ )
    {
        Int ig = u.blockMap().MyGlobalElements() [i];
        M_vectorActivationchange.epetraVector().ReplaceGlobalValue (ig,
                                                                    0,
                                                                    -0.0075 * M_solutionGatingW[ig] - 0.0105 * M_solutionMechanicalActivation[ig]);
    }
    M_vectorActivationchange.globalAssemble();
    M_solutionMechanicalActivation += timeStep * M_vectorActivationchange;
    M_solutionMechanicalActivation.globalAssemble();
#endif
    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::computeIonicCurrent (  Real Capacitance,
                                                               VectorElemental& elvec,
                                                               VectorElemental& elvec_u,
                                                               FESpace<Mesh, MapEpetra>& uFESpace )
{
    Real u_ig, w_ig;

    Real G1 = this->M_data.RMCParameterC1() / this->M_data.RMCTimeUnit() / std::pow (this->M_data.RMCPotentialAmplitude(), 2.0);
    Real G2 = this->M_data.RMCParameterC2() / this->M_data.RMCTimeUnit();

    for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
    {
        u_ig = w_ig = 0.;
        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            u_ig += elvec_u ( i ) * uFESpace.fe().phi ( i, ig );
        }
        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            w_ig += M_elvec ( i ) * uFESpace.fe().phi ( i, ig );
        }

        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            elvec ( i ) -= Capacitance * ( G1 * ( u_ig - this -> M_data.RMCRestPotential() ) *
                                           (u_ig - this->M_data.RMCRestPotential() -
                                            this->M_data.RMCParameterA() * this->M_data.RMCPotentialAmplitude() ) *
                                           (u_ig - this->M_data.RMCRestPotential() -
                                            this->M_data.RMCPotentialAmplitude() ) + G2 *
                                           (u_ig - this->M_data.RMCRestPotential() ) * w_ig) *
                           uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
        }
    }
}

template<typename Mesh, typename SolverType>
void RogersMcCulloch<Mesh, SolverType>::
initialize( )
{
    M_solutionGatingW.epetraVector().PutScalar (0.);

#ifdef ACTIVATED
    M_solutionMechanicalActivation.epetraVector().PutScalar (0.0);
#endif
}



///////////////////////////////////////////////////////////////////////////////
///////////// MINIMAL MODEL: Bueno-Orovio et al. (2008) ///////////////////////
// Parameters and data as in Ruiz-Baier et al. (2012) IUTAM Proc.         /////
//            only epicardial parameters are considered                      //
//Includes the evolution of mechanical activation on the fibers direction    //
///////////////////////////////////////////////////////////////////////////////

template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class MinimalModel : public virtual ElectroIonicSolver<Mesh, SolverType>
{
public:
    typedef typename ElectroIonicSolver<Mesh, SolverType>::data_Type data_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::vector_Type vector_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::function_Type function_Type;

    MinimalModel ( const data_Type& dataType,
                   const Mesh& mesh,
                   FESpace<Mesh, MapEpetra>& uFEspace,
                   Epetra_Comm& comm );
    virtual ~MinimalModel();

    void updateRepeated( );

    void updateElementSolution ( UInt eleID );

    void computeODECoefficients ( const Real& u_ig );

    void solveIonicModel ( const vector_Type& u, const Real timeStep );

    void computeIonicCurrent (Real Capacitance,
                              VectorElemental& elvec,
                              VectorElemental& elvec_u,
                              FESpace<Mesh, MapEpetra>& uFESpace );


    void initialize( );

    const vector_Type& solutionGatingW1() const
    {
        return M_solutionGatingW1;
    }
    const vector_Type& solutionGatingW2() const
    {
        return M_solutionGatingW2;
    }
    const vector_Type& solutionGatingW3() const
    {
        return M_solutionGatingW3;
    }

#ifdef ACTIVATED
    const vector_Type& solutionMechanicalActivation() const
    {
        return M_solutionMechanicalActivation;
    }
#endif



    Real M_MinimalEpitau1minus, M_MinimalEpitau2minus, M_MinimalEpitau30,
         M_MinimalEpitau3, M_MinimalEpitau0, M_MinimalEpiw1inf, M_MinimalEpiw2inf;

    Real M_W1change;
    Real M_W2change;
    Real M_W3change;
    Real M_current1;
    Real M_current2;


#ifdef ACTIVATED
    vector_Type M_vectorActivationchange;
#endif


protected:

    //! Global solution
    vector_Type                     M_solutionGatingW1;
    vector_Type                     M_solutionGatingW2;
    vector_Type                     M_solutionGatingW3;


    vector_Type M_vectorW1change;
    vector_Type M_vectorW2change;
    vector_Type M_vectorW3change;


#ifdef ACTIVATED
    vector_Type                     M_solutionMechanicalActivation;
#endif

    vector_Type                     M_ionicCurrent;

    vector_Type                     M_ionicCurrentRepeated;

    VectorElemental                 M_elemvecIonicCurrent;


private:
};


//
// IMPLEMENTATION
//

//! Constructor
template<typename Mesh, typename SolverType>
MinimalModel<Mesh, SolverType>::
MinimalModel ( const data_Type& dataType,
               const Mesh& mesh,
               FESpace<Mesh, MapEpetra>& uFEspace,
               Epetra_Comm& comm ) :
    ElectroIonicSolver<Mesh, SolverType> ( dataType, mesh, uFEspace, comm),
    M_solutionGatingW1 ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingW2 ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingW3 ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorW1change ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorW2change ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorW3change ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),

#ifdef ACTIVATED
    M_solutionMechanicalActivation ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorActivationchange ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
#endif

    M_ionicCurrent ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_ionicCurrentRepeated ( M_ionicCurrent, Repeated ),
    M_elemvecIonicCurrent ( ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 )
{
}

template<typename Mesh, typename SolverType>
MinimalModel<Mesh, SolverType>::
~MinimalModel()
{
}

template<typename Mesh, typename SolverType>
void MinimalModel<Mesh, SolverType>::updateRepeated( )
{
    M_ionicCurrentRepeated = M_ionicCurrent;
}

template<typename Mesh, typename SolverType>
void MinimalModel<Mesh, SolverType>::updateElementSolution ( UInt eleID )
{
    M_elemvecIonicCurrent.zero();
    UInt ig;
    for ( UInt iNode = 0 ;
            iNode < ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobalMap ( eleID, iNode );
        M_elemvecIonicCurrent.vec() [ iNode ] = M_ionicCurrentRepeated[ig];
    }
}

template<typename Mesh, typename SolverType>
void MinimalModel<Mesh, SolverType>::solveIonicModel ( const vector_Type& u,
                                                       const Real timeStep )
{

    LifeChrono chronoionmodelsolve;
    chronoionmodelsolve.start();
    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

    for ( Int i = 0 ; i < u.epetraVector().MyLength() ; i++ )
    {
        Int ig = u.blockMap().MyGlobalElements() [i];
        Real u_ig = u[ig];
        computeODECoefficients (u_ig);

        if (u_ig >= this->M_data.MinimalEpitheta1() )
        {
            M_W1change = - M_solutionGatingW1[ig] / this->M_data.MinimalEpitau1plus();
            M_current1 = - M_solutionGatingW1[ig] * (u_ig - this->M_data.MinimalEpitheta1() ) *
                         (this->M_data.MinimalEpivv() - u_ig) / this->M_data.MinimalEpitauFi();
        }
        else
        {
            M_W1change = (M_MinimalEpiw1inf - M_solutionGatingW1[ig]) / M_MinimalEpitau1minus;
            M_current1 = 0.;
        }
        if (u_ig >= this->M_data.MinimalEpitheta2() )
        {
            M_W2change = - M_solutionGatingW2[ig] / this->M_data.MinimalEpitau2plus();
            M_current2 = 1. / M_MinimalEpitau30 - M_solutionGatingW2[ig] * M_solutionGatingW3[ig]
                         / this->M_data.MinimalEpitauSi();
        }
        else
        {
            M_W2change = (M_MinimalEpiw2inf - M_solutionGatingW2[ig]) / M_MinimalEpitau2minus;
            M_current2 = (u_ig - 0.) / M_MinimalEpitau0;
        }

        M_vectorW1change.epetraVector().ReplaceGlobalValue (ig, 0., M_W1change);
        M_vectorW2change.epetraVector().ReplaceGlobalValue (ig, 0., M_W2change);
        M_vectorW3change.epetraVector().ReplaceGlobalValue (ig, 0,
                                                            ( (1.0 + tanh (this->M_data.MinimalEpik3() * (u_ig
                                                                           - this->M_data.MinimalEpiv3() ) ) ) * 0.5
                                                              - M_solutionGatingW3[ig]) / M_MinimalEpitau3);

#ifdef ACTIVATED
        M_vectorActivationchange.epetraVector().ReplaceGlobalValue (ig,
                                                                    0,
                                                                    -0.02 * M_solutionGatingW3[ig] - 0.04 * M_solutionMechanicalActivation[ig]);
#endif
        //Compute Iion= total current
        M_ionicCurrent.epetraVector().ReplaceGlobalValue (ig, 0, M_current1 + M_current2);
    }

    M_ionicCurrent.globalAssemble();
    M_vectorW1change.globalAssemble();
    M_vectorW2change.globalAssemble();
    M_vectorW3change.globalAssemble();

    M_solutionGatingW1 += timeStep * M_vectorW1change;
    M_solutionGatingW2 += timeStep * M_vectorW2change;
    M_solutionGatingW3 += timeStep * M_vectorW3change;

    M_solutionGatingW1.globalAssemble();
    M_solutionGatingW2.globalAssemble();
    M_solutionGatingW3.globalAssemble();


#ifdef ACTIVATED
    M_vectorActivationchange.globalAssemble();
    M_solutionMechanicalActivation += timeStep * M_vectorActivationchange;
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
void MinimalModel<Mesh, SolverType>::computeODECoefficients ( const Real& u_ig )
{
    if (u_ig >= this->M_data.MinimalEpitheta1minus() )
    {
        M_MinimalEpitau1minus = this->M_data.MinimalEpitau12minus();
        M_MinimalEpiw1inf = 0.;
    }
    else
    {
        M_MinimalEpitau1minus = this->M_data.MinimalEpitau11minus();
        M_MinimalEpiw1inf = 1.;
    }
    if (u_ig >= this->M_data.MinimalEpitheta2() )
    {
        M_MinimalEpitau3 = this->M_data.MinimalEpitau32();
    }
    else
    {
        M_MinimalEpitau3 = this->M_data.MinimalEpitau31();
    }
    if (u_ig >= this->M_data.MinimalEpitheta0() )
    {
        M_MinimalEpitau0 = this->M_data.MinimalEpitau02();
        M_MinimalEpiw2inf = this->M_data.MinimalEpiwStar2inf();
    }
    else
    {
        M_MinimalEpitau0 = this->M_data.MinimalEpitau01();
        M_MinimalEpiw2inf = 1.0 - u_ig / this->M_data.MinimalEpitau2inf();
    }

    M_MinimalEpitau2minus = this->M_data.MinimalEpitau21minus() +
                            (this->M_data.MinimalEpitau22minus() - this->M_data.MinimalEpitau21minus() ) *
                            0.5 * (1 + tanh (this->M_data.MinimalEpik2minus() * (u_ig -
                                                                                 this->M_data.MinimalEpiv2minus() ) ) );

    M_MinimalEpitau30 = this->M_data.MinimalEpitau301() +
                        (this->M_data.MinimalEpitau302() - this->M_data.MinimalEpitau301() ) *
                        0.5 * (1 + tanh (this->M_data.MinimalEpik30() * (u_ig -
                                                                         this->M_data.MinimalEpiv30() ) ) );
}

template<typename Mesh, typename SolverType>
void MinimalModel<Mesh, SolverType>::computeIonicCurrent (  Real Capacitance,
                                                            VectorElemental& elvec,
                                                            VectorElemental&,
                                                            FESpace<Mesh, MapEpetra>& uFESpace )
{
    Real Iion_ig ;

    for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
    {
        Iion_ig = 0.;
        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            Iion_ig += M_elemvecIonicCurrent ( i ) * uFESpace.fe().phi ( i, ig );
        }
        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            elvec ( i ) -= Iion_ig * Capacitance *
                           uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
        }
    }
}


template<typename Mesh, typename SolverType>
void MinimalModel<Mesh, SolverType>::
initialize( )
{
    M_solutionGatingW1.epetraVector().PutScalar (0.0);
    M_solutionGatingW2.epetraVector().PutScalar (1.0);
    M_solutionGatingW3.epetraVector().PutScalar (0.0);

#ifdef ACTIVATED
    M_solutionMechanicalActivation.epetraVector().PutScalar (0.0);
    M_solutionMechanicalActivation.globalAssemble();
#endif
    M_solutionGatingW1.globalAssemble();
    M_solutionGatingW2.globalAssemble();
    M_solutionGatingW3.globalAssemble();

}




/////////////////////////////////////////////////////////////////////////
//Model for the electrophysiology in the left atrium.
//Parameters and methods as in Courtemanche-Ramirez-Nattel 98.
////////////////////////////////////////////////////////////////////////


template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class CourtemancheRamirezNattel : public virtual ElectroIonicSolver<Mesh, SolverType>
{
public:
    typedef typename ElectroIonicSolver<Mesh, SolverType>::data_Type data_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::vector_Type vector_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::function_Type function_Type;

    CourtemancheRamirezNattel ( const data_Type&          dataType,
                                const Mesh&          mesh,
                                FESpace<Mesh, MapEpetra>& uFEspace,
                                Epetra_Comm&              comm );

    virtual ~CourtemancheRamirezNattel() {}

    void updateRepeated( );

    void updateElementSolution ( UInt eleID);

    void computeODECoefficients ( const Real& u_ig, const Int& ig );

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
    const vector_Type& solutionGatingAA() const
    {
        return M_solutionGatingAA;
    }
    const vector_Type& solutionGatingAI() const
    {
        return M_solutionGatingAI;
    }
    const vector_Type& solutionGatingUA() const
    {
        return M_solutionGatingUA;
    }
    const vector_Type& solutionGatingUI() const
    {
        return M_solutionGatingUI;
    }
    const vector_Type& solutionGatingXR() const
    {
        return M_solutionGatingXR;
    }
    const vector_Type& solutionGatingXS() const
    {
        return M_solutionGatingXS;
    }
    const vector_Type& solutionGatingD() const
    {
        return M_solutionGatingD;
    }
    const vector_Type& solutionGatingF() const
    {
        return M_solutionGatingF;
    }
    const vector_Type& solutionGatingFCa() const
    {
        return M_solutionGatingFCa;
    }
    const vector_Type& solutionGatingPU() const
    {
        return M_solutionGatingPU;
    }
    const vector_Type& solutionGatingPV() const
    {
        return M_solutionGatingPV;
    }
    const vector_Type& solutionGatingPW() const
    {
        return M_solutionGatingPW;
    }
    const vector_Type& vectorConcentrationNa() const
    {
        return M_vectorConcentrationNa;
    }
    const vector_Type& vectorConcentrationK() const
    {
        return M_vectorConcentrationK;
    }
    const vector_Type& vectorConcentrationCa() const
    {
        return M_vectorConcentrationCa;
    }
    const vector_Type& vectorConcentrationCaRel() const
    {
        return M_vectorConcentrationCaRel;
    }
    const vector_Type& vectorConcentrationCaUp() const
    {
        return M_vectorConcentrationCaUp;
    }

    Real M_R, M_temperature, M_F, M_permeabilityRatio, M_C, M_Ca0, M_Na0, M_K0,
         M_Vi, M_Vup, M_Vrel, M_KmNai, M_KmK0, M_KmNa, M_KmCa, M_Kup, M_consCaUpMax,
         M_consTrpnMax, M_consCmdnMax, M_consCsqnMax, M_KmTrpn, M_KmCmdn,
         M_KmCsqn,  M_gNa, M_gKl, M_gTo, M_gKr, M_gKs, M_gCal, M_gbCa, M_gbNa,
         M_INaKmax, M_INaCamax, M_IpCamax, M_Iupmax, M_krel, M_Irel,
         M_ah, M_aj, M_bh, M_bj, M_am, M_bm, M_tauh, M_tauj, M_taum,
         M_hinf, M_jinf, M_minf, M_akl, M_bkl, M_aaa, M_baa, M_aai, M_bai,
         M_tauai, M_aiinf, M_tauaa, M_aainf, M_aua, M_bua, M_aui, M_bui, M_tauui,
         M_uiinf, M_tauua, M_uainf, M_axr, M_bxr, M_tauxr, M_xrinf, M_axs, M_bxs,
         M_tauxs, M_xsinf, M_dinf, M_finf, M_taud, M_tauf, M_fCainf, M_taufca, M_taupu,
         M_puinf, M_taupv, M_pvinf, M_taupw, M_pwinf, M_EK, M_EKl, M_ENa, M_Klinf, M_ECan, M_gKur;

    // Na current
    Real M_INa;
    // Potassium current
    Real M_IKl;
    // K transient current
    Real M_Ito;
    // K ultra rapid current
    Real M_IKur;
    // K rapid current
    Real M_IKr;
    // K slow current
    Real M_IKs;
    // Ca current
    Real M_ICal;
    // Na-K exchange
    Real M_INaK;
    // Na-Ca exchange
    Real M_INaCa;
    // Ca backward
    Real M_ICab;
    // Na backward
    Real M_INab;
    // Ca pump
    Real M_IpCa;

    // Ca transfert into the reticulum
    Real M_Itr;
    // Ca transfert into the reticulum from the cell
    Real M_Iup;
    // leak of Ca
    Real M_Iupleak;

    // Total Na currents
    Real M_INatot;
    // Total K currents
    Real M_IKtot;
    // Total Ca currents
    Real M_ICatot;


    vector_Type M_vectorExponentialh;
    vector_Type M_vectorExponentialj;
    vector_Type M_vectorExponentialm;
    vector_Type M_vectorExponentialaa;
    vector_Type M_vectorExponentialai;
    vector_Type M_vectorExponentialua;
    vector_Type M_vectorExponentialui;
    vector_Type M_vectorExponentialxr;
    vector_Type M_vectorExponentialxs;
    vector_Type M_vectorExponentiald;
    vector_Type M_vectorExponentialf;
    vector_Type M_vectorExponentialfca;
    vector_Type M_vectorExponentialpu;
    vector_Type M_vectorExponentialpv;
    vector_Type M_vectorExponentialpw;
    vector_Type M_vectorInfimumh;
    vector_Type M_vectorInfimumj;
    vector_Type M_vectorInfimumm;
    vector_Type M_vectorInfimumaa;
    vector_Type M_vectorInfimumai;
    vector_Type M_vectorInfimumua;
    vector_Type M_vectorInfimumui;
    vector_Type M_vectorInfimumxr;
    vector_Type M_vectorInfimumxs;
    vector_Type M_vectorInfimumd;
    vector_Type M_vectorInfimumf;
    vector_Type M_vectorInfimumfca;
    vector_Type M_vectorInfimumpu;
    vector_Type M_vectorInfimumpv;
    vector_Type M_vectorInfimumpw;
    vector_Type M_vectorFn;
    vector_Type M_vectorIonicChange;

protected:
    //! Global solution h
    vector_Type                     M_solutionGatingH;
    //! Global solution j
    vector_Type                     M_solutionGatingJ;
    //! Global solution m
    vector_Type                     M_solutionGatingM;
    //!Global solution AA
    vector_Type                     M_solutionGatingAA;
    //!Global solution AI
    vector_Type                     M_solutionGatingAI;
    //!Global solution UA
    vector_Type                     M_solutionGatingUA;
    //!Global solution UI
    vector_Type                     M_solutionGatingUI;
    //!Global solution XR
    vector_Type                     M_solutionGatingXR;
    //!Global solution XS
    vector_Type                     M_solutionGatingXS;
    //! Global solution d
    vector_Type                     M_solutionGatingD;
    //! Global solution f
    vector_Type                     M_solutionGatingF;
    //!Global solution FCa
    vector_Type                     M_solutionGatingFCa;
    //!Global solution PU
    vector_Type                     M_solutionGatingPU;
    //!Global solution PV
    vector_Type                     M_solutionGatingPV;
    //!Global solution PW
    vector_Type                     M_solutionGatingPW;
    //!Global solution Na_i
    vector_Type                     M_vectorConcentrationNa;
    //!Global solution K_i
    vector_Type                     M_vectorConcentrationK;
    //! Global solution Ca_i
    vector_Type                     M_vectorConcentrationCa;
    //!Global solution Ca_rel
    vector_Type                     M_vectorConcentrationCaRel;
    //!Global solution Ca_up
    vector_Type                     M_vectorConcentrationCaUp;

    vector_Type         M_vectorIonicChangeCa;
    vector_Type         M_vectorIonicChangeNa;
    vector_Type         M_vectorIonicChangeK;
    vector_Type         M_vectorIonicChangeCaUp;
    vector_Type         M_vectorIonicChangeCaRel;
    // Release of Ca from the reticulum
    vector_Type         M_vectorCurrentIrel;
    vector_Type         M_ionicCurrent;
    vector_Type         M_ionicCurrentRepeated;
    VectorElemental     M_elemVecIonicCurrent;

private:
};


//
// IMPLEMENTATION
//


//! Constructor
template<typename Mesh, typename SolverType>
CourtemancheRamirezNattel<Mesh, SolverType>::
CourtemancheRamirezNattel ( const data_Type& dataType,
                            const Mesh& mesh,
                            FESpace<Mesh, MapEpetra>& uFEspace,
                            Epetra_Comm& comm ) :
    ElectroIonicSolver<Mesh, SolverType> ( dataType, mesh, uFEspace, comm),
    M_R (8.314472), //% Joules/(Kelvin*mole)
    M_temperature (307.7532), //%kelvins
    M_F (96.48533838), //% coulumbs/mmole
    M_permeabilityRatio (0.01833), //% Na/K permeability ratio
    M_C (100.), //% membrane capacitance set as 1 p-F/cm^2
    // Concentration constants
    M_Ca0 (1.8),
    M_Na0 (140.),
    M_K0 (5.4),
    // volume constants
    M_Vi (13668),
    M_Vup (1109.52),
    M_Vrel (96.48),

    M_KmNai (10.),
    M_KmK0 (1.5),
    M_KmNa (87.5),
    M_KmCa (1.38),
    M_Kup (0.00092),

    M_consCaUpMax (15),
    M_consTrpnMax (0.07),
    M_consCmdnMax (0.05),
    M_consCsqnMax (10),

    M_KmTrpn (0.0005),
    M_KmCmdn (0.00238),
    M_KmCsqn (0.8),

    // model constants
    M_gNa (7.8),
    M_gKl (0.09),
    M_gTo (0.1652),
    M_gKr (0.0294),
    M_gKs (0.129),
    M_gCal (0.1238),
    M_gbCa (0.00113),
    M_gbNa (0.000674),
    M_INaKmax (0.6),
    M_INaCamax (1600),
    M_IpCamax (0.275),
    M_Iupmax (0.005),
    M_krel (30),

    // gating variables
    M_vectorExponentialh (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialj (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialm (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialaa (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialai (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialua (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialui (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialxr (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialxs (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentiald (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialf (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialfca (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialpu (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialpv (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorExponentialpw (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumh (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumj (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumm (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumaa (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumai (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumua (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumui (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumxr (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumxs (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumd (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumf (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumfca (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumpu (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumpv (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorInfimumpw (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorFn (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_vectorIonicChange (ElectroIonicSolver<Mesh, SolverType>::M_localMap),
    M_solutionGatingH ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingJ ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingM ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingAA ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingAI ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingUA ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingUI ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingXR ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingXS ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingD ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingF ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingFCa ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingPU ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingPV ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_solutionGatingPW ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),

    // Concentration
    M_vectorConcentrationNa ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorConcentrationK ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorConcentrationCa ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorConcentrationCaRel ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorConcentrationCaUp ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorIonicChangeCa ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorIonicChangeNa ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorIonicChangeK ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorIonicChangeCaUp ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_vectorIonicChangeCaRel ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),

    //Currents
    M_vectorCurrentIrel ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_ionicCurrent ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_ionicCurrentRepeated ( M_ionicCurrent, Repeated ),
    M_elemVecIonicCurrent ( ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 )
{
}

template<typename Mesh, typename SolverType>
void CourtemancheRamirezNattel<Mesh, SolverType>::updateRepeated( )
{

    M_ionicCurrentRepeated = M_ionicCurrent;
}


template<typename Mesh, typename SolverType>
void CourtemancheRamirezNattel<Mesh, SolverType>::updateElementSolution ( UInt eleID)
{
    M_elemVecIonicCurrent.zero();
    UInt ig;
    for ( UInt iNode = 0 ; iNode < ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobalMap ( eleID, iNode );
        M_elemVecIonicCurrent.vec() [ iNode ] = M_ionicCurrentRepeated[ig];
    }
}




template<typename Mesh, typename SolverType>
void CourtemancheRamirezNattel<Mesh, SolverType>::solveIonicModel ( const vector_Type& u, const Real timeStep )
{
    LifeChrono chronoionmodelsolve;
    chronoionmodelsolve.start();
    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

    for ( Int i = 0 ; i < u.epetraVector().MyLength() ; i++ )
    {
        Int ig = u.blockMap().MyGlobalElements() [i];
        Real u_ig = u[ig];
        computeODECoefficients (u_ig, ig);

        M_ENa = (M_R * M_temperature / M_F ) * std::log (M_Na0 / M_vectorConcentrationNa[ig]);
        // fast Na current
        M_INa = M_gNa * M_solutionGatingM[ig] * M_solutionGatingM[ig] * M_solutionGatingM[ig] * M_solutionGatingH[ig] * M_solutionGatingJ[ig] * (u_ig - M_ENa);

        // time independant K current
        M_EKl = (M_R * M_temperature / M_F ) * std::log (M_K0 / M_vectorConcentrationK[ig]);
        M_Klinf = 1.0 / (1.0 + std::exp (0.07 * (u_ig + 80.0) ) );
        M_IKl = M_gKl * M_Klinf * (u_ig - M_EKl);

        // transient outward K current
        M_EK = (M_R * M_temperature / M_F ) * std::log (M_K0 / M_vectorConcentrationK[ig]);
        M_Ito = M_gTo * M_solutionGatingAA[ig] * M_solutionGatingAA[ig] * M_solutionGatingAA[ig] * M_solutionGatingAI[ig] * (u_ig - M_EK);

        // ultra rapid  delayed rectifier K current
        M_gKur = 0.005 + 0.05 / (1.0 + std::exp (- (u_ig - 15.0) / 13.0) );
        M_IKur = M_gKur * M_solutionGatingUA[ig] * M_solutionGatingUA[ig] * M_solutionGatingUA[ig] * M_solutionGatingUI[ig] * (u_ig - M_EK);

        // rapid  delayed rectifier outward K current
        M_IKr = M_gKr * M_solutionGatingXR[ig] * (u_ig - M_EK) / (1.0 + std::exp ( (u_ig + 15) / 22.4) );

        // slow  delayed rectifier outward K current
        M_IKs = M_gKs * M_solutionGatingXS[ig] * M_solutionGatingXS[ig] * (u_ig - M_EK);

        // L-type Ca current
        M_ICal = M_gCal * M_solutionGatingD[ig] * M_solutionGatingF[ig] * M_solutionGatingFCa[ig] * (u_ig - 65.0);

        Real satiga;
        satiga =
            // Na-K pump
            M_INaK = M_INaKmax * (1.0 / (1.0 + 0.1245
                                         * std::exp (-0.1 * u_ig * M_F / (M_R * M_temperature) )
                                         + 0.0365 / 7.0 * (std::exp (M_Na0 / 67.3) - 1.0)
                                         * std::exp (-u_ig * M_F / (M_R * M_temperature) ) ) )
                     * (1.0 / (1.0 + std::pow (M_KmNai / M_vectorConcentrationNa[ig], 1.5) ) ) * (M_K0 / (M_K0 + M_KmK0) );

        // Na-Ca exchanger
        M_INaCa = M_INaCamax / ( (std::pow (M_KmNa, 3.0) + std::pow (M_Na0, 3.0) ) * (M_KmCa + M_Ca0) * (1.0 + 0.1 * std::exp ( (0.35 - 1.0) * u_ig * M_F / (M_R * M_temperature) ) )
                                 * (std::exp (0.35 * u_ig * M_F / (M_R * M_temperature) ) * std::pow (M_vectorConcentrationNa[ig], 3.0) * M_Ca0 - std::exp ( (0.35 - 1.0) * u_ig * M_F / (M_R * M_temperature) ) * std::pow (M_Na0, 3.0) * M_vectorConcentrationCa[ig]) );

        // background current
        M_ECan = 0.5 * (M_R * M_temperature) / M_F * std::log (M_Ca0 / M_vectorConcentrationCa[ig]);
        M_ICab = M_gbCa * (u_ig - M_ECan);
        M_ENa = (M_R * M_temperature) * std::log (M_Na0 / M_vectorConcentrationNa[ig]);
        M_INab = M_gbNa * (u_ig - M_ENa);

        // Ca pump
        M_IpCa = M_IpCamax * (M_vectorConcentrationCa[ig]) / (0.0005 + M_vectorConcentrationCa[ig]);

        // Sarcoplasmic currents

        // Ca release current from JSR
        M_vectorFn.epetraVector().ReplaceGlobalValue (ig, 0,
                                                      M_Vrel * M_vectorCurrentIrel[ig] - ( (5.0 * std::pow (10.0, 2.0) / M_F) * (0.5 * M_ICal) - 0.2 * M_INaCa) );
        M_vectorCurrentIrel.epetraVector().ReplaceGlobalValue (ig, 0,
                                                               M_krel * M_solutionGatingPU[ig]*M_solutionGatingPU[ig]*M_solutionGatingPV[ig]*M_solutionGatingPW[ig] * (M_vectorConcentrationCaRel[ig] - M_vectorConcentrationCa[ig]) );

        // transfert current from NSR to JSR
        M_Itr = (M_vectorConcentrationCaUp[ig] - M_vectorConcentrationCaRel[ig]) / 180.0;

        // Ca uptake current by the NSR
        M_Iup = M_Iupmax * M_vectorConcentrationCa[ig] / (M_vectorConcentrationCaUp[ig] + M_Kup);

        // Ca leak current by the NSR
        M_Iupleak = M_Iupmax * M_vectorConcentrationCaUp[ig] / M_consCaUpMax;


        ////////////////////////////////////////////////////////////////////////////////
        // The original model from the CRN paper contains buffer modeling. They were not
        // implemented here because they have no incidence on transmembrane current.
        ////////////////////////////////////////////////////////////////////////////////

        M_INatot = 3.0 * M_INaK + 3.0 * M_INaCa + M_INab + M_INa;
        M_IKtot = -2.0 * M_INaK + M_IKl + M_Ito + M_IKur + M_IKr + M_IKs;
        M_ICatot = -2.0 * M_INaCa + M_IpCa + M_ICal + M_ICab;

        // for homogeneity reasons the transmembrane currents IKtot, INatot and ICatot
        // are multiplied by M_C (this multiplication doesn't appear in the original
        // paper), but the sarcoplasmic currents are not multiplied

        M_vectorIonicChangeK.epetraVector().ReplaceGlobalValue (ig, 0,
                                                                -M_C * M_IKtot / (M_F * M_Vi) );
        M_vectorIonicChangeNa.epetraVector().ReplaceGlobalValue (ig, 0,
                                                                 -M_C * M_INatot / (M_F * M_Vi) );
        M_vectorIonicChangeCa.epetraVector().ReplaceGlobalValue (ig, 0,
                                                                 (-M_C * M_ICatot / (2.0 * M_F * M_Vi) + (M_Vup * (M_Iupleak - M_Iup)
                                                                         + M_Irel * M_Vrel) / M_Vi) / (1.0 + (M_consTrpnMax * M_KmTrpn / (std::pow ( M_vectorConcentrationCa[ig] + M_KmTrpn, 2.0) ) ) + (M_consCmdnMax * M_KmCmdn / (std::pow (M_vectorConcentrationCa[ig] + M_KmCmdn, 2.0) ) ) ) );
        M_vectorIonicChangeCaUp.epetraVector().ReplaceGlobalValue (ig, 0,
                                                                   M_Iup - M_Iupleak - M_Itr * M_Vrel / M_Vup);
        M_vectorIonicChangeCaRel.epetraVector().ReplaceGlobalValue (ig, 0,
                                                                    (M_Itr - M_Irel) / (1.0 + M_consCsqnMax * M_KmCsqn / std::pow (M_vectorConcentrationCaRel[ig] + M_KmCsqn, 2.0) ) );
        M_ionicCurrent.epetraVector().ReplaceGlobalValue (ig, 0,
                                                          M_INatot  + M_IKtot + M_ICatot);

        //evolution of gating variables.
        M_vectorExponentialh.epetraVector().ReplaceGlobalValue (ig, 0, std::exp (-timeStep / M_tauh) );
        M_vectorExponentialj.epetraVector().ReplaceGlobalValue (ig, 0, std::exp (-timeStep / M_tauj) );
        M_vectorExponentialm.epetraVector().ReplaceGlobalValue (ig, 0, std::exp (-timeStep / M_taum) );
        M_vectorExponentialaa.epetraVector().ReplaceGlobalValue (ig, 0, std::exp (-timeStep / M_tauaa) );
        M_vectorExponentialai.epetraVector().ReplaceGlobalValue (ig, 0, std::exp (-timeStep / M_tauai) );
        M_vectorExponentialua.epetraVector().ReplaceGlobalValue (ig, 0, std::exp (-timeStep / M_tauua) );
        M_vectorExponentialui.epetraVector().ReplaceGlobalValue (ig, 0, std::exp (-timeStep / M_tauui) );
        M_vectorExponentialxr.epetraVector().ReplaceGlobalValue (ig,  0,  std::exp (-timeStep / M_tauxr) );
        M_vectorExponentialxs.epetraVector().ReplaceGlobalValue (ig,  0, std::exp (-timeStep / M_tauxs) );
        M_vectorExponentiald.epetraVector().ReplaceGlobalValue (ig,  0, std::exp (-timeStep / M_taud) );
        M_vectorExponentialf.epetraVector().ReplaceGlobalValue (ig, 0,  std::exp (-timeStep / M_tauf) );
        M_vectorExponentialfca.epetraVector().ReplaceGlobalValue (ig,  0,  std::exp (-timeStep / M_taufca) );
        M_vectorExponentialpu.epetraVector().ReplaceGlobalValue (ig,  0,  std::exp (-timeStep / M_taupu) );
        M_vectorExponentialpv.epetraVector().ReplaceGlobalValue (ig,  0,  std::exp (-timeStep / M_taupv) );
        M_vectorExponentialpw.epetraVector().ReplaceGlobalValue (ig, 0,  std::exp (-timeStep / M_taupw) );
        M_vectorInfimumh.epetraVector().ReplaceGlobalValue (ig, 0, M_hinf);
        M_vectorInfimumj.epetraVector().ReplaceGlobalValue (ig,  0,  M_jinf);
        M_vectorInfimumm.epetraVector().ReplaceGlobalValue (ig, 0, M_minf);
        M_vectorInfimumaa.epetraVector().ReplaceGlobalValue (ig, 0, M_aainf);
        M_vectorInfimumai.epetraVector().ReplaceGlobalValue (ig,  0, M_aiinf);
        M_vectorInfimumua.epetraVector().ReplaceGlobalValue (ig,  0,  M_uainf);
        M_vectorInfimumui.epetraVector().ReplaceGlobalValue (ig,  0,  M_uiinf);
        M_vectorInfimumxr.epetraVector().ReplaceGlobalValue (ig, 0,  M_xrinf);
        M_vectorInfimumxs.epetraVector().ReplaceGlobalValue (ig,  0, M_xsinf);
        M_vectorInfimumd.epetraVector().ReplaceGlobalValue (ig, 0,  M_dinf);
        M_vectorInfimumf.epetraVector().ReplaceGlobalValue (ig, 0, M_finf);
        M_vectorInfimumfca.epetraVector().ReplaceGlobalValue (ig, 0, M_fCainf);
        M_vectorInfimumpu.epetraVector().ReplaceGlobalValue (ig, 0, M_puinf);
        M_vectorInfimumpv.epetraVector().ReplaceGlobalValue (ig, 0, M_pvinf);
        M_vectorInfimumpw.epetraVector().ReplaceGlobalValue (ig, 0, M_pwinf);
    }
    M_vectorExponentialh.globalAssemble();
    M_vectorExponentialj.globalAssemble();
    M_vectorExponentialm.globalAssemble();
    M_vectorExponentialaa.globalAssemble();
    M_vectorExponentialai.globalAssemble();
    M_vectorExponentialua.globalAssemble();
    M_vectorExponentialui.globalAssemble();
    M_vectorExponentialxr.globalAssemble();
    M_vectorExponentialxs.globalAssemble();
    M_vectorExponentiald.globalAssemble();
    M_vectorExponentialf.globalAssemble();
    M_vectorExponentialfca.globalAssemble();
    M_vectorExponentialpu.globalAssemble();
    M_vectorExponentialpv.globalAssemble();
    M_vectorExponentialpw.globalAssemble();

    M_vectorInfimumh.globalAssemble();
    M_vectorInfimumj.globalAssemble();
    M_vectorInfimumm.globalAssemble();
    M_vectorInfimumaa.globalAssemble();
    M_vectorInfimumai.globalAssemble();
    M_vectorInfimumua.globalAssemble();
    M_vectorInfimumui.globalAssemble();
    M_vectorInfimumxr.globalAssemble();
    M_vectorInfimumxs.globalAssemble();
    M_vectorInfimumd.globalAssemble();
    M_vectorInfimumf.globalAssemble();
    M_vectorInfimumfca.globalAssemble();
    M_vectorInfimumpu.globalAssemble();
    M_vectorInfimumpv.globalAssemble();
    M_vectorInfimumpw.globalAssemble();

    M_vectorFn.globalAssemble();
    M_vectorCurrentIrel.globalAssemble();
    M_vectorIonicChangeK.globalAssemble();
    M_vectorIonicChangeNa.globalAssemble();
    M_vectorIonicChangeCa.globalAssemble();
    M_vectorIonicChangeCaUp.globalAssemble();
    M_vectorIonicChangeCaRel.globalAssemble();
    M_ionicCurrent.globalAssemble();


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
    M_solutionGatingAA -= M_vectorInfimumaa;
    M_solutionGatingAA.epetraVector().Multiply (1.,
                                                M_solutionGatingAA.epetraVector(),
                                                M_vectorExponentialaa.epetraVector(),
                                                0.);
    M_solutionGatingAA += M_vectorInfimumaa;
    M_solutionGatingAI -= M_vectorInfimumai;
    M_solutionGatingAI.epetraVector().Multiply (1.,
                                                M_solutionGatingAI.epetraVector(),
                                                M_vectorExponentialai.epetraVector(),
                                                0.);
    M_solutionGatingAI += M_vectorInfimumai;
    M_solutionGatingUA -= M_vectorInfimumua;
    M_solutionGatingUA.epetraVector().Multiply (1.,
                                                M_solutionGatingUA.epetraVector(),
                                                M_vectorExponentialua.epetraVector(),
                                                0.);
    M_solutionGatingUA += M_vectorInfimumua;
    M_solutionGatingUI -= M_vectorInfimumui;
    M_solutionGatingUI.epetraVector().Multiply (1.,
                                                M_solutionGatingUI.epetraVector(),
                                                M_vectorExponentialui.epetraVector(),
                                                0.);
    M_solutionGatingUI += M_vectorInfimumui;
    M_solutionGatingXR -= M_vectorInfimumxr;
    M_solutionGatingXR.epetraVector().Multiply (1.,
                                                M_solutionGatingXR.epetraVector(),
                                                M_vectorExponentialxr.epetraVector(),
                                                0.);
    M_solutionGatingXR += M_vectorInfimumxr;
    M_solutionGatingXS -= M_vectorInfimumxs;
    M_solutionGatingXS.epetraVector().Multiply (1.,
                                                M_solutionGatingXS.epetraVector(),
                                                M_vectorExponentialxs.epetraVector(),
                                                0.);
    M_solutionGatingXS += M_vectorInfimumxs;
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
    M_solutionGatingFCa -= M_vectorInfimumfca;
    M_solutionGatingFCa.epetraVector().Multiply (1.,
                                                 M_solutionGatingFCa.epetraVector(),
                                                 M_vectorExponentialfca.epetraVector(),
                                                 0.);
    M_solutionGatingFCa += M_vectorInfimumfca;
    M_solutionGatingPU -= M_vectorInfimumpu;
    M_solutionGatingPU.epetraVector().Multiply (1.,
                                                M_solutionGatingPU.epetraVector(),
                                                M_vectorExponentialpu.epetraVector(),
                                                0.);
    M_solutionGatingPU += M_vectorInfimumpu;
    M_solutionGatingPV -= M_vectorInfimumpv;
    M_solutionGatingPV.epetraVector().Multiply (1.,
                                                M_solutionGatingPV.epetraVector(),
                                                M_vectorExponentialpv.epetraVector(),
                                                0.);
    M_solutionGatingPV += M_vectorInfimumpv;
    M_solutionGatingPW -= M_vectorInfimumpw;
    M_solutionGatingPW.epetraVector().Multiply (1.,
                                                M_solutionGatingPW.epetraVector(),
                                                M_vectorExponentialpw.epetraVector(),
                                                0.);
    M_solutionGatingPW += M_vectorInfimumpw;
    M_vectorConcentrationNa += timeStep * M_vectorIonicChangeNa;
    M_vectorConcentrationK += timeStep * M_vectorIonicChangeK;
    M_vectorConcentrationCa += timeStep * M_vectorIonicChangeCa;
    M_vectorConcentrationCaUp += timeStep * M_vectorIonicChangeCaUp;
    M_vectorConcentrationCaRel += timeStep * M_vectorIonicChangeCaRel;


    M_solutionGatingH.globalAssemble();
    M_solutionGatingJ.globalAssemble();
    M_solutionGatingM.globalAssemble();
    M_solutionGatingAA.globalAssemble();
    M_solutionGatingAI.globalAssemble();
    M_solutionGatingUA.globalAssemble();
    M_solutionGatingUI.globalAssemble();
    M_solutionGatingXR.globalAssemble();
    M_solutionGatingXS.globalAssemble();
    M_solutionGatingD.globalAssemble();
    M_solutionGatingF.globalAssemble();
    M_solutionGatingFCa.globalAssemble();
    M_solutionGatingPU.globalAssemble();
    M_solutionGatingPV.globalAssemble();
    M_solutionGatingPW.globalAssemble();
    M_vectorConcentrationNa.globalAssemble();
    M_vectorConcentrationK.globalAssemble();
    M_vectorConcentrationCa.globalAssemble();
    M_vectorConcentrationCaUp.globalAssemble();
    M_vectorConcentrationCaRel.globalAssemble();

    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

    chronoionmodelsolve.stop();
    if (ElectroIonicSolver<Mesh, SolverType>::M_comm->MyPID() == 0)
    {
        std::cout << "Total ionmodelsolve time " << chronoionmodelsolve.diff() << " s." << std::endl;
    }
}



template<typename Mesh, typename SolverType>
void CourtemancheRamirezNattel<Mesh, SolverType>::computeODECoefficients ( const Real& u_ig,
                                                                           const Int& ig //index of the current point in the mapping
                                                                         )
{

    // Na current
    if (u_ig >= -40.)
    {
        M_ah = 0.;
        M_bh = 1. / (0.13 * (1. + std::exp ( (u_ig + 10.66) / (-11.1) ) ) );
        M_aj = 0.;
        M_bj = 0.3 * std::exp (-2.535e-7 * u_ig) / (1. + std::exp (-0.1 * (u_ig + 32.) ) );
    }
    else
    {
        M_ah = 0.135 * std::exp ( (80. + u_ig) / -6.8);
        M_bh = 3.56 * std::exp (0.079 * u_ig) + 3.1e5 * std::exp (0.35 * u_ig);
        M_aj = (-1.2714e5 * std::exp (0.2444 * u_ig) - 3.474e-5 * std::exp (-0.04391 * u_ig) ) *
               (u_ig + 37.78) / (1 + std::exp (0.311 * (u_ig + 79.23) ) );
        M_bj = 0.1212 * std::exp (-0.01052 * u_ig) / (1. + std::exp (-0.1378 * (u_ig + 40.14) ) );
    }
    M_am = 0.32 * (u_ig + 47.13) / (1. - std::exp (-0.1 * (u_ig + 47.13) ) );
    M_bm = 0.08 * std::exp (-u_ig / 11.);

    M_hinf = M_ah   / (M_ah  + M_bh);
    M_tauh = 1.   / (M_ah + M_bh);
    M_jinf = M_aj / (M_aj + M_bj);
    M_tauj = 1.   / (M_aj + M_bj);
    M_minf = M_am / (M_am + M_bm);
    M_taum = 1.   / (M_am + M_bm);

    // transient outward K current

    M_aaa = 0.65 * 1.0 / (std::exp (- (u_ig + 10.0) / 8.5) + std::exp (- (u_ig - 30.0) / 59.0) );
    M_baa = 0.65 * 1.0 / (2.5 + std::exp ( (u_ig + 82.0) / 17.0) );

    M_aai = 1.0 / (18.53 + std::exp ( (u_ig + 113.7) / 10.95) );
    M_bai = 1.0 / (35.56 + std::exp (- (u_ig + 1.26) / 7.44) );

    M_tauaa = 1.0 / ( (M_aaa + M_baa) * 3.0);
    M_aainf = 1.0 / (1.0 + std::exp (- (u_ig + 20.47) / 17.54) );

    M_tauai = 1.0 / ( (M_aai + M_bai) * 3.0);
    M_aiinf = 1.0 / (1.0 + std::exp ( (u_ig + 43.1) / 5.3) );


    // ultra rapid delayed rectifier K current

    M_aua = 0.65 * 1.0 / (std::exp (- (u_ig + 10.0) / 8.5) + std::exp (- (u_ig - 30.0) / 59.0) );
    M_bua = 0.65 * 1.0 / (2.5 + std::exp ( (u_ig + 82.0) / 17.0) );

    M_aui = 1.0 / (21.0 + std::exp (- (u_ig - 185.0) / 28.0) );
    M_bui = std::exp ( (u_ig - 158.0) / 16.0);

    M_tauua = (1.0 / (M_aua + M_bua) ) / 3.0;
    M_uainf = 1.0 / (1.0 + std::exp (- (u_ig + 30.3) / 9.6) );

    M_tauui = (1.0 / (M_aui + M_bui) ) / 3.0;
    M_uiinf = 1.0 / (1.0 + std::exp ( (u_ig - 99.45) / 27.48) );


    // rapid delayed rectifier outward K current
    M_axr = 0.0003 * (u_ig + 14.1) / (1.0 - std::exp (- (u_ig + 14.1) / 5.0) );
    M_bxr = 7.3898e-5 * (u_ig - 3.3328) / (std::exp ( (u_ig - 3.3328) / 5.1237) - 1.0);

    M_tauxr = 1.0 / (M_axr + M_bxr);
    M_xrinf = 1.0 / (1.0 + std::exp (- (u_ig + 14.1) / 6.5) );


    // slow  delayed rectifier outward K current
    M_axs = 4e-5 * (u_ig - 19.9) / (1.0 - std::exp (- (u_ig - 19.9) / 17.0) );
    M_bxs = 3.5e-5 * (u_ig - 19.9) / (std::exp ( (u_ig - 19.9) / 9.0) - 1.0);

    M_tauxs = 0.5 / (M_axs + M_bxs);
    M_xsinf = std::pow ( (1.0 + std::exp (- (u_ig - 19.9) / 12.7) ), (-1.0 / 2.0) );


    // L-type Ca current

    M_taud = (1.0 - std::exp (- (u_ig + 10.0) / 6.24) ) / (0.035 * (u_ig + 10.0) * (1.0 + std::exp (- (u_ig + 10.0) / 6.24) ) );
    M_dinf = 1.0 / (1.0 + std::exp (- (u_ig + 10.0) / 8.0) );

    M_tauf = 9.0 * 1.0 / (0.0197 * std::exp (-std::pow (0.0337, 2.0) * std::pow ( (u_ig + 10.0), 2.0) ) + 0.02);
    M_finf = 1.0 / (1.0 + std::exp ( (u_ig + 28.0) / 6.9) );

    M_taufca = 2.0;
    M_fCainf = 1.0 / (1.0 + (M_vectorConcentrationCa[ig] / 0.00035) );

    // Ca release current from JSR
    //Fn expression has been simplified, being equivalent as the initial expression (simplification of the power of 10, in the
    //following expressions using Fn). Furthermore, a modification has been done : instead of 10^-13, we put 10^-10

    M_taupu = 8.0;
    M_puinf = 1.0 / (1.0 + std::exp (- (M_vectorFn[ig] * std::pow (10, 4) - 3.4175 * std::pow (10, 3) ) / (13.67) ) );

    M_taupv = 1.91 + 2.09 / (1.0 + std::exp (- (M_vectorFn[ig] * std::pow (10, 4.0) - 3.4175 * std::pow (10, 3) ) / (13.67) ) );
    M_pvinf = 1.0 - 1.0 / (1.0 + std::exp (- (M_vectorFn[ig] * std::pow (10, 4.0) - 6.835 * std::pow (10, 2) ) / (13.67) ) );


    M_taupw = 6.0 * (1.0 - std::exp (- (u_ig - 7.9) / 5.0) ) / ( (1.0 + 0.3 * std::exp (- (u_ig - 7.9) / 5.0) ) * (u_ig - 7.9) );
    M_pwinf = 1.0 - 1.0 / (1.0 + std::exp (- (u_ig - 40.0) / 17.0) );
}

template<typename Mesh, typename SolverType>
void CourtemancheRamirezNattel<Mesh, SolverType>::computeIonicCurrent (  Real Capacitance,
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
            elvec ( i ) -= Iion_ig * Capacitance *
                           uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
        }
    }
}



//////////////////////////////////////////////////////////////////////////////
///////// N.B. The values are taken from the CRN paper and correspond to the
// asymptotic values of the model at resting potential. If this method is used,
// the corresponding value of u should be resting potential, i.e. -81.2 mV.
//////////////////////////////////////////////////////////////////////////////

template<typename Mesh, typename SolverType>
void CourtemancheRamirezNattel<Mesh, SolverType>::
initialize( )
{

    //initial values of gating variables (asymptotic values of the model at resting potential)
    M_solutionGatingH.epetraVector().PutScalar (0.965);
    M_solutionGatingD.epetraVector().PutScalar (1.37e-4);
    M_solutionGatingXR.epetraVector().PutScalar (3.29e-5);
    M_solutionGatingAI.epetraVector().PutScalar (0.999);
    M_solutionGatingUI.epetraVector().PutScalar (0.999);
    M_solutionGatingM.epetraVector().PutScalar (2.91e-3);
    M_solutionGatingJ.epetraVector().PutScalar (0.978);
    M_solutionGatingF.epetraVector().PutScalar (0.999);
    M_solutionGatingXS.epetraVector().PutScalar (1.87e-2);
    M_solutionGatingAA.epetraVector().PutScalar (3.04e-2);
    M_solutionGatingUA.epetraVector().PutScalar (4.96e-3);
    M_solutionGatingFCa.epetraVector().PutScalar (0.775);
    M_solutionGatingPV.epetraVector().PutScalar (1.0);
    M_solutionGatingPU.epetraVector().PutScalar (0.0);
    M_solutionGatingPW.epetraVector().PutScalar (0.999);



    //initial concentrations.
    M_vectorConcentrationNa.epetraVector().PutScalar (11.2);
    M_vectorConcentrationK.epetraVector().PutScalar (139.0);
    M_vectorConcentrationCa.epetraVector().PutScalar (1.02e-4);
    M_vectorConcentrationCaUp.epetraVector().PutScalar (1.49);
    M_vectorConcentrationCaRel.epetraVector().PutScalar (1.49);

    M_solutionGatingH.globalAssemble();
    M_solutionGatingD.globalAssemble();
    M_solutionGatingXR.globalAssemble();
    M_solutionGatingAI.globalAssemble();
    M_solutionGatingUI.globalAssemble();
    M_solutionGatingM.globalAssemble();
    M_solutionGatingJ.globalAssemble();
    M_solutionGatingF.globalAssemble();
    M_solutionGatingXS.globalAssemble();
    M_solutionGatingAA.globalAssemble();
    M_solutionGatingUA.globalAssemble();
    M_solutionGatingFCa.globalAssemble();
    M_solutionGatingPV.globalAssemble();
    M_solutionGatingPU.globalAssemble();
    M_solutionGatingPW.globalAssemble();
    M_vectorConcentrationNa.globalAssemble();
    M_vectorConcentrationK.globalAssemble();
    M_vectorConcentrationCa.globalAssemble();
    M_vectorConcentrationCaUp.globalAssemble();
    M_vectorConcentrationCaRel.globalAssemble();
}






} // namespace LifeV


#endif //_IONICSOLVER_H_
