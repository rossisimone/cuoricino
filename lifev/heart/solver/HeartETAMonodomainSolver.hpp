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
  @brief Class for solving the Monodomain equations in electrophysiology.

  @date 02-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @last update 02-2013

  This class provides interfaces to solve the monodomain equation
  ( reaction diffusion equation ) using the ETA framework.
  The solution can be performed using three different methods:
  -operator splitting method (at this point available only with forward Euler
      for the reaction step and backward Euler for the diffusion step. );
  -Ionic Currents Interpolation (at this point only forward Euler);
  -State Variable interpolation (at this point only forward Euler).
 */

#ifndef _HEARTETAMONODOMAINSOLVER_H_
#define _HEARTETAMONODOMAINSOLVER_H_



#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/heart/solver/IonicModels/HeartIonicModel.hpp>

#include <lifev/core/util/LifeChrono.hpp>
#include <boost/shared_ptr.hpp>
#include <lifev/core/fem/FESpace.hpp>


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>


namespace LifeV
{

//! monodomainSolver - Class featuring the usual solver for monodomain equations

template< typename Mesh, typename IonicModel >
class HeartETAMonodomainSolver
{


    //!Monodomain Solver
    /*!
       The monodomain equation reads
       \f \Chi

    */


public:


    //! @name Type definitions
    //@{

    typedef Mesh                                   mesh_Type;
    typedef boost::shared_ptr< mesh_Type >         meshPtr_Type;

    typedef VectorEpetra                           vector_Type;
    typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;
    typedef std::vector<vectorPtr_Type>            vectorOfPtr_Type;

    typedef MatrixEpetra<Real>                     matrix_Type;
    typedef boost::shared_ptr<matrix_Type>     matrixPtr_Type;

    typedef boost::shared_ptr<Epetra_Comm>     commPtr_Type;

    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >                        ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > >    ETFESpacePtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >                                    feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>                                    feSpacePtr_Type;

    typedef boost::shared_ptr<LinearSolver>                        linearSolverPtr_Type;

    typedef ExporterHDF5< mesh_Type >          exporter_Type;
    typedef boost::shared_ptr<exporter_Type>                       exporterPtr_Type;

    typedef LifeV::Preconditioner               basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack     prec_Type;
    typedef boost::shared_ptr<prec_Type>        precPtr_Type;


    typedef HeartIonicModel                     superIonicModel;
    typedef boost::shared_ptr<IonicModel>   ionicModelPtr_Type;

    typedef Teuchos::ParameterList              list_Type;

    typedef boost::function < Real (const Real& t,
                                    const Real& x,
                                    const Real& y,
                                    const Real& z,
                                    const ID&   i ) >   function_Type;


    //@}



    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    HeartETAMonodomainSolver();

    //! Constructor
    /*!
     * @param Teuchos::ParameterList parameter list
     * @param GetPot  datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     */

    HeartETAMonodomainSolver ( list_Type             list,
                               GetPot&            dataFile,
                               ionicModelPtr_Type model   );

    //! Constructor
    /*!
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel>  chosen ionic model pointer
     * @param boost::shared_ptr<Mesh> Pointer to the partitioned mesh
     */
    HeartETAMonodomainSolver ( GetPot&           dataFile,
                               ionicModelPtr_Type model,
                               meshPtr_Type       meshPtr );
    //! Constructor
    /*!
     * @param Teuchos::ParameterList parameter list
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Epetra_Comm> Epetra communicator
     */

    HeartETAMonodomainSolver ( list_Type             list,
                               GetPot&            dataFile,
                               ionicModelPtr_Type model,
                               commPtr_Type       comm    );

    //! Constructor
    /*!
     * @param Teuchos::ParameterList parameter list
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Mesh> Pointer to the partitioned mesh
     */

    HeartETAMonodomainSolver ( list_Type             list,
                               GetPot&            dataFile,
                               ionicModelPtr_Type model,
                               meshPtr_Type       meshPtr );


    //! Constructor
    /*!
     * @param string file name of the mesh
     * @param string path to the mesh
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     */

    HeartETAMonodomainSolver ( std::string           meshName,
                               std::string        meshPath,
                               GetPot&            dataFile,
                               ionicModelPtr_Type model   );

    //! Constructor
    /*!
     * @param string file name of the mesh
     * @param string path to the mesh
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Epetra_Comm> Epetra communicator
     */
    HeartETAMonodomainSolver ( std::string       meshName,
                               std::string        meshPath,
                               GetPot&            dataFile,
                               ionicModelPtr_Type model,
                               commPtr_Type       comm    );



    //! Copy Constructor
    /*!
     * @param HeartETAmonodomainSolver object
     */
    HeartETAMonodomainSolver ( const HeartETAMonodomainSolver& solver    );

    //!Operator=()
    /*!
     * @param HeartETAmonodomainSolver object
     */
    HeartETAMonodomainSolver<Mesh,  IonicModel>& operator= (const HeartETAMonodomainSolver& solver   );

    //! Destructor
    virtual ~HeartETAMonodomainSolver() {}

    //@}

    //! @name Get Methods
    //@{

    //! get the surface to volume ratio
    /*!
     * Not used in the code ( implicit definition inside the diffusion tensor)
     */
    inline const Real& surfaceVolumeRatio() const
    {
        return M_surfaceVolumeRatio;
    }

    //! get the membrane capacitance
    /*!
     * Not used in the code (  Usually equal 1 )
     */
    inline const Real& membraneCapacitance()    const
    {
        return M_membraneCapacitance;
    }

    //! get the initial time (by default 0)
    inline const Real& initialTime()            const
    {
        return M_initialTime;
    }
    //! get the final time
    inline const Real& timeStep()               const
    {
        return M_timeStep;
    }
    //! get the time step
    inline const Real& endTime()                const
    {
        return M_endTime;
    }
    //! get the diagonal diffusion tensor
    inline const VectorSmall<3>& diffusionTensor()      const
    {
        return M_diffusionTensor;
    }
    //! get the order of the elements
    inline const std::string elementsOrder()    const
    {
        return M_elementsOrder;
    }
    //! get the pointer to the ionic model
    inline const ionicModelPtr_Type     ionicModelPtr       ()  const
    {
        return M_ionicModelPtr;
    }
    //! get the pointer to the Epetra communicator
    inline const commPtr_Type           commPtr             ()  const
    {
        return M_commPtr;
    }
    //! get the pointer to the partitioned mesh
    inline const meshPtr_Type           localMeshPtr             ()  const
    {
        return M_localMeshPtr;
    }
    //! get the pointer to the partitioned mesh
    inline const meshPtr_Type           fullMeshPtr             ()  const
    {
        return M_fullMeshPtr;
    }
    //! get the pointer to the ETA finite element space
    inline const ETFESpacePtr_Type  ETFESpacePtr        ()  const
    {
        return M_ETFESpacePtr;
    }
    //! get the pointer to the usual finite element space
    inline const feSpacePtr_Type        feSpacePtr          ()  const
    {
        return M_feSpacePtr;
    }
    //! get the pointer to the mass matrix
    inline const matrixPtr_Type     massMatrixPtr       ()  const
    {
        return M_massMatrixPtr;
    }
    //! get the pointer to the stiffness matrix
    inline const matrixPtr_Type     stiffnessMatrixPtr()    const
    {
        return M_stiffnessMatrixPtr;
    }
    //! get the pointer to the global matrix
    /*!
     *  \f[
     *  A = \frac{M}{\Delta t} + K(\mathbf{f})
     *  \f]
     */
    inline const matrixPtr_Type     globalMatrixPtr ()  const
    {
        return M_globalMatrixPtr;
    }
    //! get the pointer to the right hand side
    inline const vectorPtr_Type     rhsPtr              ()  const
    {
        return M_rhsPtr;
    }
    //! get the pointer to the unique version of the right hand side
    inline const vectorPtr_Type     rhsPtrUnique        ()  const
    {
        return M_rhsPtrUnique;
    }
    //! get the pointer to the transmembrane potential
    inline const vectorPtr_Type     potentialPtr        ()  const
    {
        return M_potentialPtr;
    }
    //! get the pointer to the fiber vector
    inline const vectorPtr_Type     fiberPtr            ()  const
    {
        return M_fiberPtr;
    }
    //! get the pointer to the applied current vector
    inline const vectorPtr_Type     appliedCurrentPtr   ()  const
    {
        return M_appliedCurrentPtr;
    }
    //! get the pointer to the linear solver
    inline const linearSolverPtr_Type   linearSolverPtr ()  const
    {
        return M_linearSolverPtr;
    }
    //! get the pointer to the vector of pointers containing the transmembrane potential (at 0) and the gating variables
    inline const vectorOfPtr_Type&      globalSolution      ()  const
    {
        return M_globalSolution;
    }
    //! get the pointer to the vector of pointers containing the rhs for transmembrane potential (at 0) and the gating variables
    inline const vectorOfPtr_Type&      globalRhs           ()  const
    {
        return M_globalRhs;
    }

    //@}

    //! @name Set Methods
    //@{

    //! set the surface to volume ratio (NOT USED IN THE CODE)
    /*!
        @param Real surface to volume ratio
     */
    inline void setSurfaceVolumeRatio   ( const Real& p )
    {
        this -> M_surfaceVolumeRatio = p;
    }
    //! set the membrane capacitance (NOT USED IN THE CODE, SET TO 1)
    /*!
        @param Real membrane capacitance
     */
    inline void setMembraneCapacitance  ( const Real& p )
    {
        this -> M_membraneCapacitance = p;
    }
    //! set the starting time
    /*!
        @param Real initial time
     */
    inline void setInitialTime          ( const Real& p )
    {
        this -> M_initialTime = p;
    }
    //! set the ending time
    /*!
        @param Real ending time
     */
    inline void setTimeStep                 ( const Real& p )
    {
        this -> M_timeStep = p;
    }
    //! set the time step
    /*!
        @param Real time step
     */
    inline void setEndTime              ( const Real& p )
    {
        this -> M_endTime = p;
    }
    //! set the diagonal diffusion tensor
    /*!
        @param  VectorSmall<3> diagonal diffusion tensor
     */
    inline void setDiffusionTensor      ( const VectorSmall<3>& p )
    {
        this -> M_diffusionTensor = p;
    }
    //! set the pointer to the ionic model
    /*!
        @param boost::shared_ptr<IonicModel> pointer to the ionic model
     */
    inline void setIonicModelPtr        ( const ionicModelPtr_Type  p       )
    {
        this -> M_ionicModelPtr = p ;
    }
    //! set the pointer to the Epetra communicator
    /*!
        @param boost::shared_ptr<Epetra_Comm> pointer to the Epetra communicator
     */

    inline void setCommPtr              ( const commPtr_Type        p       )
    {
        this -> M_commPtr = p ;
    }
    //! set the pointer to the partitioned mesh
    /*!
        @param boost::shared_ptr<Mesh> pointer to the partitioned mesh
     */
    inline void setLocalMeshPtr              ( const meshPtr_Type        p       )
    {
        this -> M_localMeshPtr = p ;
    }
    //! set the pointer to the partitioned mesh
    /*!
        @param boost::shared_ptr<Mesh> pointer to the partitioned mesh
     */
    inline void setFullMeshPtr              ( const meshPtr_Type        p       )
    {
        this -> M_fullMeshPtr = p ;
    }
    //! set the pointer to the ETA fe space
    /*!
        @param  boost::shared_ptr<ETFESpace<Mesh,MapEpetra,3,1>> pointer to the ETA fe space
     */
    inline void setETFESpacePtr             ( const ETFESpacePtr_Type   p       )
    {
        this -> M_ETFESpacePtr = p ;
    }
    //! set the pointer to the usual fe space
    /*!
        @param boost::shared_ptr<IFESpace<Mesh,MapEpetra>> pointer to the usual fe space
     */
    inline void setFeSpacePtr           ( const feSpacePtr_Type     p       )
    {
        this -> M_feSpacePtr = p ;
    }

    //! set the pointer to the  mass matrix
    /*!
        @param boost::shared_ptr<MatrixEpetra<Real>> pointer to the mass matrix
     */
    inline void setMassMatrixPtr            ( const matrixPtr_Type      p       )
    {
        this -> M_massMatrixPtr = p ;
    }
    //! set the pointer to the stiffness matrix
    /*!
        @param boost::shared_ptr<MatrixEpetra<Real>> pointer to the stiffness matrix
     */
    inline void setStiffnessMatrixPtr   ( const matrixPtr_Type      p       )
    {
        this -> M_stiffnessMatrixPtr = p ;
    }
    //! set the pointer to the global matrix
    /*!
        @param boost::shared_ptr<MatrixEpetra<Real>> pointer to the global matrix
     */
    inline void setGlobalMatrixPtr      ( const matrixPtr_Type      p       )
    {
        this -> M_globalMatrixPtr = p ;
    }
    //! set the pointer to the right hand side
    /*!
        @param boost::shared_ptr<VectorEpetra> pointer to the right hand side
     */
    inline void setRhsPtr                   ( const vectorPtr_Type      p       )
    {
        this -> M_rhsPtr = p ;
    }
    //! set the pointer to the unique version of the right hand side
    /*!
        @param boost::shared_ptr<VectorEpetra>  pointer to the  unique version of the right hand side
     */
    inline void setRhsPtrUnique         ( const vectorPtr_Type      p       )
    {
        this -> M_rhsPtrUnique = p ;
        //! set the pointer to the transmembrane potential
        /*!
            @param boost::shared_ptr<VectorEpetra>  pointer to the transmembrane potential
         */
        this -> M_globalRhs .at (0) = M_rhsPtrUnique;
    }
    inline void setPotentialPtr         ( const vectorPtr_Type      p       )
    {
        this -> M_potentialPtr = p ;
        this -> M_globalSolution.at (0) = M_potentialPtr;
    }
    //! set the pointer to the applied current vector
    /*!
        @param boost::shared_ptr<VectorEpetra>  pointer to the applied current vector
     */
    inline void setAppliedCurrentPtr        ( const vectorPtr_Type      p       )
    {
        this -> M_appliedCurrentPtr = p ;
    }
    //! set the pointer to the linear solver
    /*!
        @param boost::shared_ptr<LinearSolver> pointer to the linear solver
     */
    inline void setLinearSolverPtr      ( const linearSolverPtr_Type    p   )
    {
        this -> M_linearSolverPtr = p ;
    }
    //! set the vector of pointers containing the transmembrane potential (at 0) and the gating variables
    /*!
        @param std::vector<boost::shared_ptr<VectorEpetra>> vector of pointers containing the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalSolution       ( const vectorOfPtr_Type&   p       )
    {
        this -> M_globalSolution = p;
    }
    //! set the vector of pointers containing the rhs for the transmembrane potential (at 0) and the gating variables
    /*!
        @param std::vector<boost::shared_ptr<VectorEpetra>> vector of pointers containing the rhs for the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalRhs                ( const vectorOfPtr_Type&   p       )
    {
        this -> M_globalRhs      = p;
    }
    //! set the pointer to the fiber direction vector
    /*!
        @param boost::shared_ptr<VectorEpetra> pointer to the fiber direction vector
     */
    inline void setFiberPtr         ( const vectorPtr_Type      p       )
    {
        this -> M_fiberPtr = p ;
    }

    //@}

    //! @name Copy Methods
    //@{
    inline void copyIonicModel      ( const ionicModelPtr_Type  p   )
    {
        (* (M_ionicModelPtr) ) = *p ;
    }
    inline void copyComm                ( const commPtr_Type        p   )
    {
        (* (M_commPtr) ) = *p ;
    }
    inline void copyLocalMesh                ( const meshPtr_Type        p   )
    {
        (* (M_localMeshPtr) ) = *p ;
    }
    inline void copyFullMesh                ( const meshPtr_Type        p   )
    {
        (* (M_fullMeshPtr) ) = *p ;
    }
    inline void copyETFESpace       ( const ETFESpacePtr_Type   p   )
    {
        (* (M_ETFESpacePtr) ) = *p ;
    }
    inline void copyFeSpace         ( const feSpacePtr_Type     p   )
    {
        (* (M_feSpacePtr) ) = *p ;
    }
    inline void copyMassMatrix      ( const matrixPtr_Type      p   )
    {
        (* (M_massMatrixPtr) ) = *p ;
    }
    inline void copyStiffnessMatrix ( const matrixPtr_Type      p   )
    {
        (* (M_stiffnessMatrixPtr) ) = *p ;
    }
    inline void copyGlobalMatrix        ( const matrixPtr_Type      p   )
    {
        (* (M_globalMatrixPtr) ) = *p ;
    }
    inline void copyRhs             ( const vectorPtr_Type      p   )
    {
        (* (M_rhsPtr) ) = *p ;
    }
    inline void copyRhsUnique       ( const vectorPtr_Type      p   )
    {
        (* (M_rhsPtrUnique) ) = *p ;
    }
    inline void copyPotential       ( const vectorPtr_Type      p   )
    {
        (* (M_potentialPtr) ) = *p ;
    }
    inline void copyFiber               ( const vectorPtr_Type      p   )
    {
        (* (M_fiberPtr) ) = *p ;
    }
    inline void copyAppliedCurrent  ( const vectorPtr_Type      p   )
    {
        (* (M_appliedCurrentPtr) ) = *p ;
    }
    inline void copyLinearSolver        ( const linearSolverPtr_Type p  )
    {
        (* (M_linearSolverPtr) ) = *p ;
    }
    inline void copyGlobalSolution  ( const vectorOfPtr_Type&   p   )
    {
        for (int j = 0; j < M_ionicModelPtr -> Size(); j++ )
        {
            ( * ( M_globalSolution.at (j) ) )   = (* (p.at (j) ) );
        }
    }
    inline void copyGlobalRhs       ( const vectorOfPtr_Type&   p   )
    {
        for (int j = 0; j < M_ionicModelPtr -> Size(); j++ )
        {
            ( * ( M_globalRhs.at (j) ) )    = (* (p.at (j) ) );
        }
    }

    //@}

    //! @name Methods
    //@{

    void setup ( GetPot& dataFile, short int ionicSize);

    void setup (std::string meshName, std::string meshPath, GetPot& dataFile, short int ionicSize);

    //! create mass matrix
    void setupMassMatrix();
    //! create mass matrix
    void setupLumpedMassMatrix();
    //! create stiffness matrix
    void setupStiffnessMatrix();
    //! create stiffness matrix given a diagonal diffusion tensor
    void setupStiffnessMatrix (VectorSmall<3> diffusion);
    //! create stiffness matrix given the fiber direction and a diagonal diffusion tensor
    void setupStiffnessMatrix (VectorEpetra& fiber, VectorSmall<3> diffusion);
    //! setup the total matrix
    /*!
     *  \f[
     *  A = \frac{M}{\Delta t} + K(\mathbf{f})
     *  \f]
     */
    void setupGlobalMatrix();
    //! setup the linear solver
    /*!
     * A file named MonodomainSolverParamList.xml must be in the execution folder
     * with the parameters to set the linear solver
     */
    void setupLinearSolver ( GetPot dataFile );
    //! setup the linear solver
    void setupLinearSolver ( GetPot dataFile, list_Type list );
    //! Initialize the potential to zero
    void inline initializePotential()
    {
        (*M_potentialPtr)     *= 0;
    }
    //! Initialize the potential to the value k
    void inline initializePotential (Real k)
    {
        (*M_potentialPtr)      = k;
    }
    //! Initialize the applied current to zero
    void inline initializeAppliedCurrent()
    {
        (*M_appliedCurrentPtr) *= 0;
    }
    //! Initialize the applied current to the value k
    void inline initializeAppliedCurrent (Real k)
    {
        (*M_appliedCurrentPtr) = k;
    }
    //! creates a vector of pointers to store the solution
    /*!
     * The first pointer points to the vector of the transmembrane potential,
     * while the others point to the gating variables
     */
    void setupGlobalSolution (short int ionicSize);
    //! creates a vector of pointers to store the rhs
    /*!
     * The first pointer points to the rhs of the transmembrane potential,
     * while the others point to the rhs of the gating variables
     */
    void setupGlobalRhs (short int ionicSize);
    //! Set parameters from an xml file
    void setParameters (list_Type    list);
    //! partition the mesh
    void inline partitionMesh ( std::string  meshName, std::string   meshPath)
    {
        MeshUtility::fillWithFullMesh ( M_fullMeshPtr, M_localMeshPtr, meshName, meshPath );
    }
    //! given an boost function initialize the potential
    void inline setPotentialFromFunction ( function_Type f )
    {
        M_feSpacePtr -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *M_potentialPtr , 0);
    }
    //! given an boost function initialize the applied current
    void inline setAppliedCurrentFromFunction ( function_Type f )
    {
        M_feSpacePtr -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *M_appliedCurrentPtr , 0);
    }
    //! Solves one reaction step using the forward Euler scheme
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^n + \Delta t I_{ion}(\mathbf{V}^n).
     * \f]
     */
    void solveOneReactionStepFE();

    //! Update the rhs
    /*!
     * \f[
     * rhs \leftarrow \frac{M}{\Delta t} \mathbf{V}^n
     * \f]
     */
    void inline updateRhs()
    {
        (*M_rhsPtrUnique) += (*M_massMatrixPtr) * (*M_potentialPtr) * ( 1.0 / M_timeStep );
    }
    //! Solves one diffusion step using the backward Euler scheme
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^*.
     * \f]
     */
    void solveOneDiffusionStepBE();

    //!Solve one full step with operator splitting
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^n + \Delta t I_{ion}(\mathbf{V}^n).
     * \f]
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^*.
     * \f]
     */
    void solveOneSplittingStep();
    //!Solve the system with operator splitting from M_initialTime to the M_endTime with time step M_timeStep
    void solveSplitting();
    //!Solve one full step with operator splitting and export the solution
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^n + \Delta t I_{ion}(\mathbf{V}^n).
     * \f]
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^*.
     * \f]
     */
    void solveOneSplittingStep ( exporter_Type& exporter, Real t );
    //!Solve the system with operator splitting from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    void solveSplitting ( exporter_Type& exporter );
    //! add to a given exporter the pointer to the potential
    void setupPotentialExporter (exporter_Type& exporter);
    //! add to a given exporter the pointer to the potential and to the gating variables
    void setupExporter (exporter_Type& exporter);

    //! add to a given exporter the pointer to the potential saved with name fileName
    void setupPotentialExporter (exporter_Type& exporter, std::string fileName);
    //! add to a given exporter the pointer to the potential and to the gating variables saved with name fileName
    void setupExporter (exporter_Type& exporter, std::string fileName);

    //! Generates a default fiber direction (0,1,0)
    void setupFibers();
    //! Generates the fiber direction given the three component of the vector (F_x,F_y,F_z)
    void setupFibers (VectorSmall<3> fibers);
    //! Generates the fiber direction from a  given file
    void setupFibers (std::string fibersFile);
    //! Solves the gating variables with forward Euler
    void solveOneStepGatingVariablesFE();
    //! Compute the rhs using state variable interpolation
    void computeRhsSVI();
    //! Compute the rhs using ionic current interpolation
    void computeRhsICI();
    //!Solve one full step with ionic current interpolation
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+M\mathbf{I},
     * \f]
     * where $\mathbf{I}$ is the vector of the ionic currents $I_j = I_{ion}(V_j^n)$
     */
    void solveOneICIStep();
    //!Solve one full step with ionic current interpolation
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+\mathbf{I}_{ion}(\mathbf{V}^n).
     * \f]
     */
    void solveOneSVIStep();
    //! solve system using ICI from M_initialTime to the M_endTime with time step M_timeStep
    void solveICI();
    //! solve system using SVI from M_initialTime to the M_endTime with time step M_timeStep
    void solveSVI();
    //!Solve one full step with ionic current interpolation  and export the solution
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+M\mathbf{I},
     * \f]
     * where $\mathbf{I}$ is the vector of the ionic currents $I_j = I_{ion}(V_j^n)$
     */
    void solveOneICIStep (exporter_Type& exporter, Real t);
    //!Solve one full step with ionic current interpolation  and export the solution
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+\mathbf{I}_{ion}(\mathbf{V}^n).
     * \f]
     */
    void solveOneSVIStep (exporter_Type& exporter, Real t);

    //! solve system using ICI from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    void solveICI (exporter_Type& exporter);
    //! solve system using SVI from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    void solveSVI (exporter_Type& exporter);
    //! Generates a file where the fiber direction is saved
    void exportFiberDirection();
    //! save the fiber direction into the given exporter
    void exportFiberDirection (exporter_Type& exporter);

    //! Save the solution in the exporter
    void inline exportSolution (exporter_Type& exporter, Real t)
    {
        exporter.postProcess (t);
    }

    //@}

private:

    void setParameters();
    void init();
    void init ( commPtr_Type comm );
    void init ( meshPtr_Type meshPtr );
    void init ( ionicModelPtr_Type model);
    void init ( commPtr_Type comm, ionicModelPtr_Type model);
    void init ( meshPtr_Type meshPtr, ionicModelPtr_Type model );

    Real                M_surfaceVolumeRatio;
    Real                M_membraneCapacitance;

    ionicModelPtr_Type  M_ionicModelPtr;

    commPtr_Type        M_commPtr;
    meshPtr_Type        M_localMeshPtr;
    meshPtr_Type        M_fullMeshPtr;
    ETFESpacePtr_Type   M_ETFESpacePtr;
    feSpacePtr_Type     M_feSpacePtr;
    matrixPtr_Type      M_massMatrixPtr;
    matrixPtr_Type      M_stiffnessMatrixPtr;
    matrixPtr_Type      M_globalMatrixPtr;

    Real                M_initialTime;
    Real                M_endTime;
    Real                M_timeStep;

    VectorSmall<3>      M_diffusionTensor;

    vectorPtr_Type          M_rhsPtr;
    vectorPtr_Type          M_rhsPtrUnique;
    vectorPtr_Type          M_potentialPtr;
    vectorPtr_Type          M_appliedCurrentPtr;

    linearSolverPtr_Type    M_linearSolverPtr;


    vectorOfPtr_Type        M_globalSolution;
    vectorOfPtr_Type        M_globalRhs;

    std::string             M_elementsOrder;


    vectorPtr_Type          M_fiberPtr;

}; // class MonodomainSolver



//
// IMPLEMENTATION
//
// ===================================================
//! Constructors
// ===================================================
template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver()
{
    setParameters();
    init();
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver (   list_Type list,
                             GetPot& dataFile,
                             ionicModelPtr_Type model)
{
    init (model);
    setParameters (list);
    setup (list.get ("meshName", "lid16.mesh"), list.get ("meshPath", "./"), dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver (   list_Type           list,
                             GetPot&             dataFile,
                             ionicModelPtr_Type  model,
                             commPtr_Type        comm    )
{
    setParameters (list);
    init (comm);
    setup (list.get ("meshName", "lid16.mesh"), list.get ("meshPath", "./"), dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh, IonicModel>::
HeartETAMonodomainSolver (   list_Type           list,
                             GetPot&             dataFile,
                             ionicModelPtr_Type  model,
                             meshPtr_Type        meshPtr )
{
    setParameters (list);
    init (meshPtr);
    setup (dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::
HeartETAMonodomainSolver (   std::string         meshName,
                             std::string         meshPath,
                             GetPot&             dataFile,
                             ionicModelPtr_Type  model   )
{
    setParameters();
    init (model);
    setup (meshName, meshPath, dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::
HeartETAMonodomainSolver (   std::string         meshName,
                             std::string         meshPath,
                             GetPot&             dataFile,
                             ionicModelPtr_Type  model,
                             commPtr_Type        comm    ) :
    M_ionicModelPtr (model)
{
    setParameters();
    init (comm);
    setup (meshName, meshPath, dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::
HeartETAMonodomainSolver (   GetPot&             dataFile,
                             ionicModelPtr_Type  model,
                             meshPtr_Type        meshPtr ) :
    M_ionicModelPtr (model)
{
    setParameters();
    init (meshPtr);
    setup (dataFile, M_ionicModelPtr -> Size() );
}

template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>::
HeartETAMonodomainSolver ( const HeartETAMonodomainSolver& solver    ) :
    M_surfaceVolumeRatio    ( solver.M_surfaceVolumeRatio ),
    M_membraneCapacitance   ( solver.M_membraneCapacitance ),
    M_ionicModelPtr         ( solver.M_ionicModelPtr ),
    M_commPtr               ( solver.M_commPtr ),
    M_localMeshPtr          ( solver.M_localMeshPtr ),
    M_fullMeshPtr           ( solver.M_fullMeshPtr ),
    M_ETFESpacePtr          ( solver.M_ETFESpacePtr ),
    M_feSpacePtr            (  solver.M_feSpacePtr ),
    M_massMatrixPtr         ( new matrix_Type ( * (solver.M_massMatrixPtr) ) ),
    M_stiffnessMatrixPtr    ( new matrix_Type ( * (solver.M_stiffnessMatrixPtr) ) ),
    M_globalMatrixPtr       ( new matrix_Type ( * (solver.M_globalMatrixPtr) ) ),
    M_initialTime           ( solver.M_initialTime ),
    M_endTime               ( solver.M_endTime ),
    M_timeStep              ( solver.M_timeStep ),
    M_diffusionTensor       ( solver.M_diffusionTensor ),
    M_rhsPtr                ( new vector_Type ( * (solver.M_rhsPtr ) ) ),
    M_rhsPtrUnique          ( new vector_Type ( * (M_rhsPtr), Unique) ),
    M_potentialPtr          ( new vector_Type ( solver.M_ETFESpacePtr->map() ) ),
    //  M_potentialPtr          ( new vector_Type( *(solver.M_potentialPtr     ) ) ),
    M_appliedCurrentPtr     ( new vector_Type ( * (solver.M_appliedCurrentPtr) ) ),
    M_linearSolverPtr       ( new LinearSolver ( * (solver.M_linearSolverPtr  ) ) ),
    M_elementsOrder         ( solver.M_elementsOrder ),
    M_fiberPtr              (  new vector_Type ( * (solver.M_fiberPtr) ) )
{
    std::cout << "new potential: " << M_potentialPtr << endl;
    std::cout << "given potential: " << solver.M_potentialPtr << endl;

    setupGlobalSolution (M_ionicModelPtr -> Size() );
    std::cout << "new potential: " << M_potentialPtr << endl;
    std::cout << "given potential: " << solver.M_potentialPtr << endl;
    copyGlobalSolution (solver.M_globalSolution);
    setupGlobalRhs (M_ionicModelPtr -> Size() );
    copyGlobalRhs (solver.M_globalRhs);
}


template<typename Mesh, typename IonicModel>
HeartETAMonodomainSolver<Mesh,  IonicModel>& HeartETAMonodomainSolver<Mesh,  IonicModel>::
operator= (const HeartETAMonodomainSolver& solver    )
{
    M_surfaceVolumeRatio    = solver.M_surfaceVolumeRatio;
    M_membraneCapacitance   = solver.M_membraneCapacitance;
    copyIonicModel ( solver.M_ionicModelPtr );
    M_commPtr               = solver.M_commPtr ;
    M_localMeshPtr               = solver.M_localMeshPtr ;
    M_fullMeshPtr               = solver.M_fullMeshPtr ;
    copyETFESpace ( solver.M_ETFESpacePtr );
    copyFeSpace ( solver.M_feSpacePtr );
    copyMassMatrix ( solver.M_massMatrixPtr );
    copyStiffnessMatrix (solver.M_stiffnessMatrixPtr );
    copyGlobalMatrix ( solver.M_globalMatrixPtr );
    M_initialTime           = solver.M_initialTime ;
    M_endTime               = solver.M_endTime ;
    M_timeStep              = solver.M_timeStep ;
    M_diffusionTensor       = solver.M_diffusionTensor ;
    copyRhs (solver.M_rhsPtr);
    copyRhsUnique (solver.M_rhsPtrUnique);
    std::cout << "new potential: " << M_potentialPtr << endl;
    std::cout << "given potential: " << solver.M_potentialPtr << endl;
    copyPotential (solver.M_potentialPtr);
    std::cout << "new potential: " << M_potentialPtr << endl;
    std::cout << "given potential: " << solver.M_potentialPtr << endl;

    copyAppliedCurrent (solver.M_appliedCurrentPtr);
    copyLinearSolver (solver.M_linearSolverPtr);
    copyGlobalSolution (solver.M_globalSolution);
    copyGlobalRhs (solver.M_globalRhs);
    M_elementsOrder         = solver.M_elementsOrder ;
    copyFiber ( solver.M_fiberPtr );
    return      *this;
}


/********* SETUP METHODS *///////

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupFibers()
{

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( M_localMeshPtr, M_elementsOrder, 3, M_commPtr) );

    M_fiberPtr.reset ( new vector_Type ( Space3D -> map() ) );

    int d1 = (*M_fiberPtr).epetraVector().MyLength() / 3;
    (*M_fiberPtr) *= 0;
    int j (0);
    for ( int k (0); k < d1; k++)
    {
        j = (*M_fiberPtr).blockMap().GID (k + d1);
        (*M_fiberPtr) [j] = 1.0;
    }
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupFibers (VectorSmall<3> fibers)
{
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( M_localMeshPtr, M_elementsOrder, 3, M_commPtr) );

    M_fiberPtr.reset ( new vector_Type ( Space3D -> map() ) );

    int n1 = (*M_fiberPtr).epetraVector().MyLength();
    int d1 = n1 / 3;
    (*M_fiberPtr) *= 0;
    int i (0);
    int j (0);
    int k (0);
    for ( int l (0); l < d1; l++)
    {

        i = (*M_fiberPtr).blockMap().GID (l);
        j = (*M_fiberPtr).blockMap().GID (l + d1);
        k = (*M_fiberPtr).blockMap().GID (l + 2 * d1);
        (*M_fiberPtr) [i] = fibers[0];
        (*M_fiberPtr) [j] = fibers[1];
        (*M_fiberPtr) [k] = fibers[2];

    }
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::setupFibers (std::string fibersFile)
{
//    std::stringstream MyPID;
//    ifstream fibers (fibersFile.c_str() );
//
//    std::cout << "fiber_file: " <<  fibersFile.c_str() << std::endl;
//    UInt NumGlobalElements =  M_feSpacePtr.map (Repeated)->NumGlobalElements();
//    std::vector<Real> fiber_global_vector (NumGlobalElements);
//
//    for ( UInt i = 0; i < NumGlobalElements; ++i)
//    {
//        fibers >> fiber_global_vector[i];
//    }
//    UInt NumMyElements = M_localMapVector.map (Repeated)->NumMyElements();
//    for (UInt j = 0; j < NumMyElements; ++j)
//    {
//        UInt ig = M_localMapVector.map (Repeated)->MyGlobalElements() [j];
//        (*M_fiberPtr)[ig] = fiber_global_vector[ig];
//    }
//    std::cout << std::endl;
//    fiber_global_vector.clear();

}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
exportFiberDirection()
{
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( M_localMeshPtr, M_elementsOrder, 3, M_commPtr) );

    ExporterHDF5<mesh_Type> exp;
    exp.setMeshProcId ( M_localMeshPtr, M_commPtr -> MyPID() );
    exp.setPrefix ("FiberDirection");
    exp.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", Space3D, M_fiberPtr, UInt (0) );
    exp.postProcess (0);
    exp.closeFile();
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
exportFiberDirection (exporter_Type& exporter)
{
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( M_localMeshPtr, M_elementsOrder, 3, M_commPtr) );

    exporter.addVariable ( ExporterData<mesh_Type>::VectorField,  "fibers", Space3D, M_fiberPtr, UInt (0) );
    exporter.postProcess (0);
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setup (GetPot& dataFile, short int ionicSize)
{

    M_ETFESpacePtr.reset ( new ETFESpace_Type ( M_localMeshPtr, &feTetraP1, M_commPtr) );
    M_feSpacePtr.reset (new FESpace< mesh_Type, MapEpetra > (M_localMeshPtr, M_elementsOrder, 1, M_commPtr) );
    M_feSpacePtr.reset ( new feSpace_Type (M_localMeshPtr, M_elementsOrder, 1, M_commPtr)  );
    M_massMatrixPtr.reset ( new matrix_Type ( M_ETFESpacePtr->map() ) );
    M_stiffnessMatrixPtr.reset ( new matrix_Type ( M_ETFESpacePtr->map() ) );
    M_globalMatrixPtr.reset ( new matrix_Type ( M_ETFESpacePtr->map() ) );
    M_rhsPtr.reset ( new vector_Type ( M_ETFESpacePtr->map(), Repeated ) );
    M_rhsPtrUnique.reset ( new vector_Type ( * (M_rhsPtr), Unique ) );
    M_potentialPtr.reset ( new vector_Type ( M_ETFESpacePtr->map() ) );
    M_appliedCurrentPtr.reset ( new vector_Type ( M_ETFESpacePtr->map() ) );

    //***********************//
    //  Setup Linear Solver  //
    //***********************//
    setupLinearSolver ( dataFile );


    //**************************//
    //  Setup Initial condition //
    //**************************//
    initializePotential (0.0);
    initializeAppliedCurrent (0.0);
    setupGlobalSolution (ionicSize);
    setupGlobalRhs (ionicSize);
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setup ( std::string meshName, std::string meshPath, GetPot& dataFile, short int ionicSize)
{
    //partitioning the mesh
    partitionMesh ( meshName, meshPath );
    setup (dataFile, ionicSize);
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupMassMatrix()
{
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_localMeshPtr  ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   phi_i * phi_j
			   )
			   >> M_massMatrixPtr;

	}
	M_massMatrixPtr -> globalAssemble();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupLumpedMassMatrix()
{
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_localMeshPtr  ),
				   quadRuleTetra4ptNodal,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   phi_i * phi_j
			   )
			   >> M_massMatrixPtr;

	}
	M_massMatrixPtr -> globalAssemble();
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupStiffnessMatrix()
{
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_localMeshPtr ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   dot(  rotate( M_ETFESpacePtr, *M_fiberPtr, M_diffusionTensor ) * grad(phi_i) , grad(phi_j) )
		   )
		   >> M_stiffnessMatrixPtr;

	}
	M_stiffnessMatrixPtr -> globalAssemble();
}



template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupStiffnessMatrix (VectorSmall<3> diffusion)
{
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_localMeshPtr  ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   dot(  rotate( M_ETFESpacePtr, *M_fiberPtr, diffusion ) * grad(phi_i) , grad(phi_j) )
		   )
		   >> M_stiffnessMatrixPtr;

	}
	M_stiffnessMatrixPtr -> globalAssemble();
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupStiffnessMatrix (VectorEpetra& fiber, VectorSmall<3> diffusion)
{
	{
	   using namespace ExpressionAssembly;

	   integrate(  elements( M_localMeshPtr  ),
				   quadRuleTetra4pt,
				   M_ETFESpacePtr,
				   M_ETFESpacePtr,
				   dot( rotate( M_ETFESpacePtr, fiber, diffusion ) * grad(phi_i) , grad(phi_j) )
		   )
		   >> M_stiffnessMatrixPtr;

	}
	M_stiffnessMatrixPtr -> globalAssemble();
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupGlobalMatrix()
{
    (*M_globalMatrixPtr) *= 0;
    (*M_globalMatrixPtr) = (*M_stiffnessMatrixPtr);
    (*M_globalMatrixPtr) += ( (*M_massMatrixPtr) * ( 1.0 / M_timeStep ) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupLinearSolver ( GetPot dataFile )
{
    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );
    Teuchos::RCP< Teuchos::ParameterList > solverParamList = Teuchos::rcp ( new Teuchos::ParameterList );
    solverParamList = Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" );

    M_linearSolverPtr -> setCommunicator ( M_commPtr );
    M_linearSolverPtr -> setParameters ( *solverParamList );
    M_linearSolverPtr -> setPreconditioner ( precPtr );
    M_linearSolverPtr -> setOperator ( M_globalMatrixPtr );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupLinearSolver ( GetPot dataFile, list_Type list )
{
    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    M_linearSolverPtr -> setCommunicator ( M_commPtr );
    M_linearSolverPtr -> setParameters ( list );
    M_linearSolverPtr -> setPreconditioner ( precPtr );
    M_linearSolverPtr -> setOperator ( M_globalMatrixPtr );
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupGlobalSolution (short int ionicSize)
{
    M_globalSolution.push_back ( M_potentialPtr );
    for (int k = 1; k <  ionicSize; ++k )
    {
        M_globalSolution.push_back ( * (new vectorPtr_Type ( new VectorEpetra ( M_ETFESpacePtr -> map() ) ) ) );
    }
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh,  IonicModel>::
setupGlobalRhs (short int ionicSize)
{
    M_globalRhs.push_back ( M_rhsPtrUnique );
    for (int k = 1; k < ionicSize; ++k )
    {
        M_globalRhs.push_back ( * (new vectorPtr_Type ( new VectorEpetra ( M_ETFESpacePtr -> map() ) ) ) );
    }
}


/************** EXPORTER *///////////////

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupPotentialExporter (exporter_Type& exporter)
{
    exporter.setMeshProcId ( M_localMeshPtr, M_commPtr -> MyPID() );
    exporter.setPrefix ("Potential");
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "Potential", M_feSpacePtr, M_potentialPtr, UInt (0) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupExporter (exporter_Type& exporter)
{
    exporter.setMeshProcId ( M_localMeshPtr, M_commPtr -> MyPID() );
    exporter.setPrefix ("Solution");
    std::string variableName;
    for ( int i = 0; i < M_ionicModelPtr -> Size() ; i++ )
    {
        variableName = "Variable" + boost::lexical_cast<std::string> ( i );
        exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  variableName, M_feSpacePtr, M_globalSolution.at (i), UInt (0) );
    }
}



template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupPotentialExporter (exporter_Type& exporter, std::string fileName)
{
    exporter.setMeshProcId ( M_localMeshPtr, M_commPtr -> MyPID() );
    exporter.setPrefix (fileName);
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "Potential", M_feSpacePtr, M_potentialPtr, UInt (0) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setupExporter (exporter_Type& exporter, std::string fileName)
{
    exporter.setMeshProcId ( M_localMeshPtr, M_commPtr -> MyPID() );
    exporter.setPrefix (fileName);
    std::string variableName;
    for ( int i = 0; i < M_ionicModelPtr -> Size() ; i++ )
    {
        variableName = "Variable" + boost::lexical_cast<std::string> ( i );
        exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  variableName, M_feSpacePtr, M_globalSolution.at (i), UInt (0) );
    }
}


/********* SOLVING METHODS */////////////////////////

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneReactionStepFE()
{
    M_ionicModelPtr -> superIonicModel::computeRhs ( M_globalSolution, *M_appliedCurrentPtr, M_globalRhs);

    for ( UInt i = 0; i < M_ionicModelPtr -> Size() ; i++ )
    {
        * ( M_globalSolution.at (i) ) = * ( M_globalSolution.at (i) ) + M_timeStep * ( * ( M_globalRhs.at (i) ) );
    }

}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneDiffusionStepBE()
{
    M_linearSolverPtr -> setRightHandSide ( M_rhsPtrUnique );
    M_linearSolverPtr -> solve ( M_potentialPtr );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSplittingStep()
{
    solveOneReactionStepFE();
    (*M_rhsPtrUnique) *= 0;
    updateRhs();
    solveOneDiffusionStepBE();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSplittingStep (exporter_Type& exporter, Real t)
{
    solveOneSplittingStep();
    exportSolution (exporter, t);
}



template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSplitting()
{
    for ( Real t = M_initialTime; t < M_endTime; )
    {
        solveOneSplittingStep();
        t = t + M_timeStep;
    }
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSplitting (exporter_Type& exporter)
{
    for ( Real t = M_initialTime; t < M_endTime; )
    {
        t = t + M_timeStep;
        solveOneSplittingStep (exporter, t);

    }
}



template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneStepGatingVariablesFE()
{
    M_ionicModelPtr -> superIonicModel::computeRhs ( M_globalSolution, M_globalRhs);

    for ( UInt i = 1; i < M_ionicModelPtr -> Size() ; i++ )
    {
        * ( M_globalSolution.at (i) ) = * ( M_globalSolution.at (i) ) + M_timeStep * ( * ( M_globalRhs.at (i) ) );
    }

}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
computeRhsICI()
{
    M_ionicModelPtr -> superIonicModel::computePotentialRhsICI ( M_globalSolution, (*M_appliedCurrentPtr), M_globalRhs, (*M_massMatrixPtr) );
    updateRhs();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
computeRhsSVI()
{
    M_ionicModelPtr -> superIonicModel::computePotentialRhsSVI ( M_globalSolution, (*M_appliedCurrentPtr), M_globalRhs, (*M_feSpacePtr) );
    updateRhs();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneICIStep()
{
    computeRhsICI();
    M_linearSolverPtr -> setRightHandSide ( M_rhsPtrUnique );
    M_linearSolverPtr -> solve ( M_potentialPtr );
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSVIStep()
{
    computeRhsSVI();
    M_linearSolverPtr -> setRightHandSide ( M_rhsPtrUnique );
    M_linearSolverPtr -> solve ( M_potentialPtr );
}




template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveICI()
{
    for ( Real t = M_initialTime; t < M_endTime; )
    {
        t = t + M_timeStep;
        solveOneStepGatingVariablesFE();
        solveOneICIStep();
    }
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSVI()
{
    for ( Real t = M_initialTime; t < M_endTime; )
    {
        t = t + M_timeStep;
        solveOneStepGatingVariablesFE();
        solveOneSVIStep();
    }
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneICIStep (exporter_Type& exporter, Real t)
{
    solveOneICIStep();
    exportSolution (exporter, t);
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveOneSVIStep (exporter_Type& exporter, Real t)
{
    solveOneSVIStep();
    exportSolution (exporter, t);
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveICI (exporter_Type& exporter)
{
    for ( Real t = M_initialTime; t < M_endTime; )
    {
        t = t + M_timeStep;
        solveOneStepGatingVariablesFE();
        solveOneICIStep ( exporter, t);
    }
}


template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
solveSVI (exporter_Type& exporter)
{
    for ( Real t = M_initialTime; t < M_endTime; )
    {
        t = t + M_timeStep;
        solveOneStepGatingVariablesFE();
        solveOneSVIStep ( exporter, t);
    }
}



/********   INITIALIZITION FOR CONSTRUCTOR ****///////

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init()
{
    M_linearSolverPtr.reset ( new LinearSolver() );
    M_globalSolution = * (new vectorOfPtr_Type() ) ;
    M_globalRhs =  * (new vectorOfPtr_Type() ) ;
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init (ionicModelPtr_Type model)
{
    init();
    M_commPtr.reset ( new Epetra_MpiComm (MPI_COMM_WORLD)  );
    M_localMeshPtr.reset ( new mesh_Type ( M_commPtr ) );
    M_fullMeshPtr.reset ( new mesh_Type ( M_commPtr ) );
    M_ionicModelPtr = model;

}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init ( commPtr_Type comm  )
{
    init();
    M_commPtr = comm;
    M_localMeshPtr.reset ( new mesh_Type ( M_commPtr ) );
    M_fullMeshPtr.reset ( new mesh_Type ( M_commPtr ) );
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init (meshPtr_Type meshPtr)
{
    init();
    //TODO change the meshPtr to pass the fullMeshPtr
    M_localMeshPtr = meshPtr;
    M_fullMeshPtr.reset ( new mesh_Type ( M_commPtr ) );
    M_commPtr = meshPtr -> comm();
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init ( commPtr_Type comm, ionicModelPtr_Type model)
{
    init ( comm );
    M_ionicModelPtr = model;
}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
init ( meshPtr_Type meshPtr, ionicModelPtr_Type model )
{
    init ( meshPtr );
    M_ionicModelPtr = model;
}

/********* parameter initialization */////////
template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setParameters()
{
    M_surfaceVolumeRatio = 2400.0;
    M_membraneCapacitance = 1.0;
    M_diffusionTensor[0] = 0.001;
    M_diffusionTensor[1] = 0.001;
    M_diffusionTensor[2] = 0.001;
    M_initialTime        = 0.0;
    M_endTime            = 100.0;
    M_timeStep           = 0.01;
    M_elementsOrder      = "P1";

}

template<typename Mesh, typename IonicModel>
void HeartETAMonodomainSolver<Mesh, IonicModel>::
setParameters (list_Type list)
{
    M_surfaceVolumeRatio = list.get ("surfaceVolumeRatio", 2400.0  );
    M_membraneCapacitance = list.get ("membraneCapacitance", 1.0 );
    M_diffusionTensor[0] = list.get ("longitudinalDiffusion", 0.001 );
    M_diffusionTensor[1] = list.get ("transversalDiffusion", 0.001 );
    M_diffusionTensor[2] = M_diffusionTensor[1];
    M_initialTime        = list.get ("initialTime", 0.0 );
    M_endTime            = list.get ("endTime", 100.0 );
    M_timeStep           = list.get ("timeStep", 0.01 );
    M_elementsOrder      = list.get ("elementsOrder", "P1" );

}


} // namespace LifeV


#endif //_MONODOMAINSOLVER_H_
