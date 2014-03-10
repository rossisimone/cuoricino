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
  @brief Implementation of the Aliev Panfilov ionic model.
  @See Nash MP, Panfilov AV bla bla bla
  @

  @date 07-2012
  @author Aymen Laadhari <aymen.laadhari@epfl.ch> and Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Aymen Laadhari <aymen.laadhari@epfl.ch> and Simone Rossi <simone.rossi@epfl.ch>
  @last update 07-2012
 */


#ifndef ELECTROALIEVPANFILOV_HPP_
#define ELECTROALIEVPANFILOV_HPP_


#include <lifev/electrophysiology/solver/ElectroIonicSolver.hpp>



using namespace LifeV;

//////////////////////////////////////////////////////////////////////////////
// Aliev Panfilov Model
//////////////////////////////////////////////////////////////////////////////


template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class ElectroAlievPanfilov : public virtual ElectroIonicSolver<Mesh, SolverType>
{
public:

    //! @name Type definitions
    //@{

    typedef typename ElectroIonicSolver<Mesh, SolverType>::data_Type  data_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::vector_Type     vector_Type;
    typedef typename ElectroIonicSolver<Mesh, SolverType>::function_Type  function_Type;

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
    ElectroAlievPanfilov ( const data_Type& dataType,
                           const Mesh& mesh,
                           FESpace<Mesh, MapEpetra>& uFEspace,
                           Epetra_Comm& comm );

    //! Destructor
    virtual ~ElectroAlievPanfilov();

    //@}



    //! @name Methods
    //@{

    inline const vector_Type& getPotential() const
    {
        return M_potential;
    }
    inline const vector_Type& getRecoveryVariable() const
    {
        return M_recoveryVariable;
    }


    void updateRepeated( );

    //! Update the ionic model elvecs
    void updateElementSolution ( UInt eleID );

    //! Solves the ionic model
    void solveIonicModel ( const vector_Type& u, const Real timeStep );


    //! Solves the ionic model
    void solveIonicModel ( ElectroFunctors& data, const Real timeStep, Real t );

    inline Real applyStimulus()
    {
        return M_stimulus;
    };
    inline void setStimulus (Real stim)
    {
        M_stimulus = stim;
    }

    //! Computes the ionic currents
    //! for the PDE righthand side
    void computeIonicCurrent ( Real Capacitance,
                               VectorElemental& elvec,
                               VectorElemental& elvec_u,
                               FESpace<Mesh, MapEpetra>& uFESpace );

    //const vector_Type& solutionGatingW() const {return M_solutionGatingW;}

    //! Initialize
    void initialize( );

    void updatePotential ( const VectorEpetra& V );

    //@}


    //recovery variable r of the Aliev-Panfilov model
    vector_Type                 M_recoveryVariable;
    //intermediate potential of the Aliev-Panfilov model
    vector_Type                 M_potential;
    //copy for parallel computations
    vector_Type                 M_recoveryVariableRepeated;
    //recovery variable r of the Aliev-Panfilov model
    vector_Type                 M_potentialRepeated;
    //Vector on the element
    VectorElemental             M_elvec;
    //order for time integration
    UInt                        M_BDForder;
    //dunno!
    TimeAdvanceBDF<vector_Type> M_BDFr;
    TimeAdvanceBDF<vector_Type> M_BDFV;
    //set stimulus
    Real                        M_stimulus;

};


// ===================================================
//! Constructors
// ===================================================
template<typename Mesh, typename SolverType>
ElectroAlievPanfilov<Mesh, SolverType>::
ElectroAlievPanfilov ( const data_Type& dataType,
                       const Mesh& mesh,
                       FESpace<Mesh, MapEpetra>& uFEspace,
                       Epetra_Comm& comm ) :
    ElectroIonicSolver<Mesh, SolverType> ( dataType, mesh, uFEspace, comm),
    M_recoveryVariable ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_potential ( ElectroIonicSolver<Mesh, SolverType>::M_localMap ),
    M_recoveryVariableRepeated ( M_recoveryVariable, Repeated ),
    M_potentialRepeated ( M_potential, Repeated ),
    M_elvec ( ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof(), 1 ),
    M_BDForder ( ElectroIonicSolver<Mesh, SolverType>::M_data.MSBDForder() ),
    M_stimulus ( 10.0 ) //,
{
    M_BDFr.setup ( M_BDForder );
    M_BDFV.setup ( M_BDForder );
}

template<typename Mesh, typename SolverType>
ElectroAlievPanfilov<Mesh, SolverType>::
~ElectroAlievPanfilov()
{
}

// ===================================================
//! Methods
// ===================================================

template<typename Mesh, typename SolverType>
void ElectroAlievPanfilov<Mesh, SolverType>::updateRepeated( )
{
    M_recoveryVariableRepeated = M_recoveryVariable;
    M_potentialRepeated = M_potential;
}


template<typename Mesh, typename SolverType>
void ElectroAlievPanfilov<Mesh, SolverType>::updateElementSolution ( UInt eleID )
{
    M_elvec.zero();
    UInt ig;
    //! Filling local elvec_r with recovery variable values in the nodes
    for ( UInt iNode = 0 ; iNode < ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.fe().nbFEDof() ; iNode++ )
    {
        ig = ElectroIonicSolver<Mesh, SolverType>::M_uFESpace.dof().localToGlobalMap ( eleID, iNode );
        M_elvec.vec() [ iNode ] = M_recoveryVariableRepeated[ig];
        M_elvec.vec() [ iNode ] = M_potentialRepeated[ig];
    }
}

template<typename Mesh, typename SolverType>
void ElectroAlievPanfilov<Mesh, SolverType>::solveIonicModel ( const vector_Type& u, const Real timeStep )
{

    Real epsilon = 0.01;
    Real mu1     = 0.12;
    Real mu2     = 0.30;
    Real k       = 8.00;
    Real b       = 0.10;

    Real aux     = 0.0;
    Real dr      = 0.0;

    M_BDFr.updateRHSContribution (timeStep);
    //vector_Type M_time_der=M_BDFr.rhsContributionFirstDerivative();

    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();
    for ( Int i = 0 ; i < u.epetraVector().MyLength() ; ++i )
    {
        Int ig = u.blockMap().MyGlobalElements() [i];


        aux = epsilon + mu1 * M_BDFr.solution() [ig] / ( mu2 + u[ig] );
        dr  = - aux * ( M_BDFr.solution() [ig] + k * u[ig] * ( u[ig] - b - 1.0 ) );
        M_recoveryVariable[ig] = M_BDFr.solution() [ig] + timeStep * dr;
    }
    M_BDFr.shiftRight (M_recoveryVariable);

    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

}


template<typename Mesh, typename SolverType>
void ElectroAlievPanfilov<Mesh, SolverType>::solveIonicModel ( ElectroFunctors& data, const Real timeStep, Real t )
{

    //Real dt = timeStep / iter;

    Real epsilon = 0.01;
    Real mu1     = 0.12;
    Real mu2     = 0.30;
    Real k       = 8.00;
    Real b       = 0.10;

    Real auxr     = 0.0;
    Real dr      = 0.0;

    Real dV      = 0.0;

    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();


    for ( Int i = 0 ; i < M_potential.epetraVector().MyLength() ; ++i )
    {
        Int ig = M_potential.blockMap().MyGlobalElements() [i];

        auxr = epsilon + mu1 * M_recoveryVariable[ig] / ( mu2 + M_potential[ig] );
        dr   = - auxr * ( M_recoveryVariable[ig] + k * M_potential[ig] * ( M_potential[ig] - b - 1.0 ) );
        M_recoveryVariable[ig] = M_recoveryVariable[ig] + timeStep * dr;

        dV = - k * M_potential[ig] * ( M_potential[ig] - b ) * ( M_potential[ig] - 1.0 ) - M_recoveryVariable[ig] * M_potential[ig];

        if ( t >= data.M_stimulusStart1 && t <= data.M_stimulusStop1 && ig == 1 )
        {
            dV += data.M_stimulusValue1;
        }
        if ( t >= data.M_stimulusStart1 && t <= data.M_stimulusStop1 && ig == 4 )
        {
            dV += data.M_stimulusValue1;
        }
        if ( t >= data.M_stimulusStart1 && t <= data.M_stimulusStop1 && ig == 63 )
        {
            dV += data.M_stimulusValue1;
        }
        if ( t >= data.M_stimulusStart1 && t <= data.M_stimulusStop1 && ig == 274)
        {
            dV += data.M_stimulusValue1;
        }

        M_potential[ig] = M_potential[ig] + timeStep * dV;
    }


    ElectroIonicSolver<Mesh, SolverType>::M_comm->Barrier();

}


template<typename Mesh, typename SolverType>
void ElectroAlievPanfilov<Mesh, SolverType>::computeIonicCurrent (  Real,
                                                                    VectorElemental& elvec,
                                                                    VectorElemental& elvec_u,
                                                                    FESpace<Mesh, MapEpetra>& uFESpace )
{

    Real k       = 8.00;
    Real b       = 0.10;
    for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
    {
        for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
        {
            elvec ( i ) += ( - k * elvec_u ( i ) * ( elvec_u ( i ) - b ) * ( elvec_u ( i ) - 1.0 )
                             - M_elvec ( i ) * elvec_u ( i ) ) *
                           uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );
        }
    }
}


template<typename Mesh, typename SolverType>
void ElectroAlievPanfilov<Mesh, SolverType>::
initialize( )
{
    M_recoveryVariable.epetraVector().PutScalar (0.0);
    M_potential.epetraVector().PutScalar (0.0);
    M_BDFr.setInitialCondition (M_recoveryVariable);
    M_BDFV.setInitialCondition (M_potential);
    M_BDFr.showMe();
    M_BDFV.showMe();
}


template<typename Mesh, typename SolverType>
void ElectroAlievPanfilov<Mesh, SolverType>::
updatePotential ( const VectorEpetra& V )
{
    this->M_potential = V;
}



#endif /* ALIEVPANFILOV_HPP_ */

