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

/**
   \file FSIMonolithicGI.hpp
   @brief Monolithic Geometry--Implicit FSI Solver
   \author crosetto <Paolo Crosetto>
   \date 18 Sep 2008

   This file implements the Monolithic Geometry--Implicit solver, see \cite CrosettoEtAl2009 for details

*/
#ifndef _MONOLITHICGI_HPP
#define _MONOLITHICGI_HPP

#include <life/lifesolver/FSIMonolithic.hpp>

namespace LifeV
{
#ifdef OBSOLETE
class Epetra_FullMonolithic;
#endif

   typedef FactorySingleton<Factory<FSIOperator, std::string> >                    FSIFactory_Type;

/**
   FSIMonolithic Geomitry-Implicit solver
 * Class handling the nonlinear monolithic solver for FSI problems. The (exact or inexact)
 * Newton algorithm is used to solve the nonlinearity.
 * The block structure of the jacobian matrix is
 *\f$\left(\begin{array}{ccc}
 C&B&S\\
 D&N&0\\
 0&E&H
 \end{array}\right)\f$
 * where \f$N\f$ represents the solid block, \f$C\f$ the fluid block, \f$H\f$ is the harmonic extension block,
 * while the extra
 * diagonal blocks represent the coupling. The implementation of the stress continuity coupling condition
 * is obtained by means of an augmented formulation.
 * Different possible preconditioners are implemented.


 Important parameters to set properly in the data file:
 - useShapeDerivatives: if true the shape derivatives block is added to the Jacobian matrix;
 - domainVelImplicit: if true the domain velocity w in the convective term is considered an unknown (at the time n+1);
 - convectiveTermDer: false if the convective term is linearized (\f$u^{n+1}\nabla(u^n-w^n)\f$),
 otherwise it can be either true (if we use the Newton method to solve the convective term nonlinearity) or false
 (fixed-point method);
 - semiImplicit:  if true only one iteration of the nonlinear solver is performed. Otherwise
 the nonlinear iterations continue up to the specified tolerance. Set it to true for the GCE;
 - method: can be either monolithicGE, monolithicGI if the geometry is treated respectively explicitly or implicitly,
 or exactJacobians, fixedPoint for partitioned strategies;
 - blockOper: specifies the matrix type to be used for the linear system: if AdditiveSchwarz, the matrix is the standard
 ine for GE; if AdditiveSchwarzRN the coupling blocks are of Robin type instead of Dirichlet and Neumann. The parameters
 for the Robin coupling are alphaf and alphas in the data file. NOTE: this method has currently been tested only for
 alphas=0.
 - DDBlockPrec: specifies the possible preconditioners to use. Can be: AdditiveSchwarz, MonolithicBlockComposedDN, MonolithicBlockComposedDN2,
 MonolithicBlockComposedNN, MonolithicBlockComposedDNND.
 */

class FSIMonolithicGI : public FSIMonolithic
{
public:

    typedef FSIMonolithic                                  super_Type;
    typedef Preconditioner                                 prec_Type;
    typedef boost::shared_ptr<prec_Type>                   prec_type;

    //!@name Constructor and Destructor
    //@{

    //! Empty Constructor
    FSIMonolithicGI();

    //! Destructor
    virtual ~FSIMonolithicGI() {}

    //@}


    //!@name Public Methods
    //@{

    //! Initializes the system with functions
    void initialize( fluidPtr_Type::value_type::function_Type const& u0,
                     fluidPtr_Type::value_type::function_Type const& p0,
                     solidPtr_Type::value_type::Function const& d0,
                     solidPtr_Type::value_type::Function const& /*w0*/,
                     fluidPtr_Type::value_type::function_Type const& /*df0*/ );

    //! Initializes the system with vectors
    void initialize( const vector_Type& un )
    {
        M_un.reset( new vector_Type( un ) );
        M_uk.reset( new vector_Type( un ) );
    }

    /**
       Sets the parameters read from data file
    */
    void setUp( const GetPot& dataFile );

    //! initializes the fluid and mesh problems, creates the map of the global matrix
    void setupFluidSolid( UInt const fluxes );

    //! builds the constant part of the fluid-structure-mesh motion matrix
    void buildSystem ();

    /**
       updates the solution, advances of a time step
    */
    void updateSystem();

    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: fluid domain displacement solution
       \param iter: current NonLinearRichardson (block Gauss Seidel for the tangent system) iteration
    */
    void evalResidual( vector_Type&  res, const vector_Type& sol, const UInt iter );

    //!Apply the boundary conditions to each block composing the monolithic problem
    /**
       Sets the vectors of: boundary conditions, FESpaces, couplings, offsets, and sets the blocks in the composed operator
       which constitutes the monolithic problem. Then calls the applyBoundaryConditions of the MonolithicBlockMatrix operator, passing
       also the right hand side.
     */
    void applyBoundaryConditions();
    //@}

    //!@name Set Methods
    //@{

    //! set the solution
    void setSolution( const vector_Type& solution ) { M_uk.reset( new vector_Type( solution ) ); }

    void setSolutionPtr( const vectorPtr_Type& sol) { M_uk = sol; }

    //!Builds an extrapolation of the solution to initialize the Newton scheme
    void couplingVariableExtrap( );

    //@}


    //!@name Get Methods
    //@{

    //! getter for the map of fluid-structure-interface (without the mesh motion)
    const MapEpetra& mapWithoutMesh() const { return *M_mapWithoutMesh; }

    //! getter for the global matrix of the system
    const matrixPtr_Type matrixPtr() const { return M_monolithicMatrix->matrix(); }

    //! getter for the pointer to the current iteration solution
    //const vectorPtr_Type  uk()  const      {return M_uk;}

    //! get the current solution vector.
    const vector_Type& solution() const { return *M_uk; }

    //! get the solution.
    vectorPtr_Type& solutionPtr() { return M_uk; }

    //@}

protected:

    //!@name Protected Methods
    //@{

    //! set the block preconditioner
    void setupBlockPrec();

    //@}

private:

    //! @name Private Methods
    //@{

    //! Factory method for the system matrix, of type MonolithicBlockBase
    void createOperator( std::string& operType )
    {
        M_monolithicMatrix.reset(MonolithicBlockMatrix::Factory_Type::instance().createObject( operType ));
        M_monolithicMatrix.reset(MonolithicBlockMatrix::Factory_Type::instance().createObject( operType ));
    }

    /**
       calculates the terms due to the shape derivatives given the mesh increment deltaDisp. The shape derivative block is assembled in a matrix
       (not in a right hand side representing the matrix-vector multiplication)
       \param sdMatrix: output. Shape derivatives block to be summed to the Jacobian matrix.
    */
    void shapeDerivatives( matrixPtr_Type sdMatrix );

    //! assembles the mesh motion matrix.
    /*!In Particular it diagonalize the part of the matrix corresponding to the
      Dirichlet condition expressing the coupling
      \param iter: current iteration: used as flag to distinguish the first nonlinear iteration from the others
     */
    void assembleMeshBlock(UInt iter);

    //@}


    //!@name Private Members
    //@{

    boost::shared_ptr<MapEpetra>         M_mapWithoutMesh;
    vectorPtr_Type                       M_uk;
    bool                                 M_domainVelImplicit;
    bool                                 M_convectiveTermDer;
    UInt                                 M_interface;
    matrixPtr_Type                       M_meshBlock;
    matrixPtr_Type                       M_shapeDerivativesBlock;
    matrixPtr_Type                       M_solidDerBlock;
    //std::vector<fluidBchandlerPtr_Type>    M_BChsLin;
    static bool                          S_register;
    //@}

    //! Factory method
    static FSIOperator* instantiate() { return new FSIMonolithicGI(); }

};

}
#endif

