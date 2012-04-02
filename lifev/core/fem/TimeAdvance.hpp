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
    @brief File containing a class to  deal the time advancing scheme

    @date
    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
*/

#ifndef TIMEADVANCE_H
#define TIMEADVANCE_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <boost/numeric/ublas/vector.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
typedef boost::numeric::ublas::vector<Real> ScalarVector;

//! timeAdvance_template - File containing a class to deal the time advancing scheme
/*!
  @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

  This class define an abstract method to build temporal discretization schemes.
  In particular we consider problems of the first order and the second  order in time;
  that after space and time discretization, and suitable linearitazion of non linear operator,
  we obtain a linear system to solve at each time step (or each iteration ):

  \f[ K U^{n+1} =F^{n+1}\f]

  where \f$K\f$ is an opportune matrix, \f$ U^{n+1}\f$ is the unknown vector and \f$F^{n+1}\f$ is the
  right hand side vector at the time \f$t^{n+1}\f$ .
  To determine \f$F^{n+1}\f$ we define the state vector \f$X^{n+1}\f$ that contained the informations
  about previous solutions.
  <ol>
  <li>   First order problems:

  \f[ M \frac{\partial u}{\partial t} + L( u) = f ,\f]

  where L is a generic nonlinear operator.
  We define U an approximation of u and the velocity vector \f$V\f$  an approximation of \f$\dot{u}\f$
  at the time step \f$n+1\f$ as

 \f[ V^{n+1}=\frac{\alpha_0}{\Delta t} U^{n+1}-  f_V^{n+1},  \f]

 where \f$\alpha_0\f$ is a suitable coefficient and  \f$f_V\f$  is a linear combination of the previous  solutions
 \f$X^{n+1}\f$ with coefficients \f$\alpha_i\f$, that  will be specified in the following.


 Then the time discrete version is
 \f[ \frac{\alpha_0}{\Delta t}M U^{n+1}  +A (U^{n+1})= f^{n+1}+ M  f_V^{n+1} \f],
 that can be solved by any non-linear iterative solver (fixed point, Newton, ....).
 This  class  provides also a suitable extrapolation \f$ U^*\f$ of  \f$U^{n+1}\f$, given by a linear
 combination of previous solutions and of order consistent with the time discretization formula.
 </li>
 <li> In this part we want to extend the previous approach to   second order
   problems in time.

  Second order problems:

 \f[ \ddot{ u}  + L( \dot{u}, u) = f \f]
 where \f$L\f$  is non linear operator.

  We consider the following semidiscrete problem

   \f[ M \frac{d^2 U}{d t^2} + D ( U,  \frac{d U}{d t} ) + A ( U ) = f  \f]

   where \f$M\f$ is the mass matrix, \f$f_V\f$ the forcing term
   and \f$U\f$ is the  vector of unknowns.

   We define the following quantities

   \f[ V :=\frac{d U}{d  t }\qquad  W := \frac{d^2  U} {d  t^2},  \f]
   where \f$V\f$ and \f$W\f$ are the velocity and the  acceleration vectors, respectively.

    At the time step \f$n+1\f$, we consider the following approssimations of \f$V\f$ and  \f$W\f$:

     \f[ V^{n+1}=\frac{\alpha_0}{\Delta t} U^{n+1}- f_V^{n+1},  \f]

     and

     \f[ W^{n+1}=\frac{\xi_0}{\Delta t^2} U^{n+1} - f_W^{n+1}, \f]

     where \f$f_V^{n+1}\f$ and \f$f_W^{n+1}\f$ are linear combinations of the previous  solutions with
     suitable coefficients \f$\alpha_i\f$ and  \f$\xi_i\f$, respectively.
     If  \f$A\f$ and \f$D\f$ depend on  \f$U\f$ and \f$V\f$ we can linearize the
     problem using suitable extrapolations \f$U^*\f$ \f$V^*\f$ obtained   by linear combinations
     of previous solutions with coefficients  \f$\beta_i\f$ and \f$\beta_i^V\f$ respectively.
     </li>
     </ol>
*/

template<typename feVectorType = VectorEpetra >

class TimeAdvance
{
public:

    //! @name Public Types
    //@{
    typedef ScalarVector                             container_Type;
    typedef std::vector<feVectorType>                feVectorContainer_Type;
    typedef std::vector<feVectorType*>               feVectorContainerPtr_Type;
    typedef typename feVectorContainerPtr_Type::iterator  feVectorContainerPtrIterate_Type;
    typedef std::vector<boost::shared_ptr<feVectorType> > feVectorSharedPtrContainer_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    TimeAdvance();

    //! Destructor
    ~TimeAdvance();

     //@}

    //! @name Methods
    //@{

     //!Update the vectors of the previous time steps by shifting on the right
     /*!
     Update the vectors of the previous time steps by shifting on the right
     the old values.
     @param solution  is  a (new) value of the state vector
     */
     virtual void shiftRight(const feVectorType& solution ) = 0;

    //!Update the right hand side
    /*
    update rhs contributions: \f$f_V\f$ and \$f_W\f$
    */
    void updateRHSContribution( const Real& timeStep);

     //! Update the right hand side \f$ f_V \f$ of the time derivative formula
     /*!
     Sets and Returns the right hand side \f$ f_V \f$ of the time derivative formula
     @param timeStep defined the  time step need to compute the
     @return  rhsV the first order Rhs
     */
     virtual  void updateRHSFirstDerivative(const Real& timeStep = 1 )  = 0;

     //! Update the right hand side \f$ f_W \f$ of the time derivative formula
     /*
     Sets and Returns the right hand side \f$ f_W \f$ of the time derivative formula
     @param timeStep defined the  time step need to compute the \f$ f_W \f$
     @return  rhsW the fsecond order Rhs
     */
     virtual void updateRHSSecondDerivative(const Real& timeStep = 1 ) = 0;

     //!Show the properties  of temporal scheme
     /*!
     Show the properties  of temporal scheme
     */
     virtual void showMe()  const = 0;

     //! Spy state vector
     /*!
     Spy of stateVector;
     this method saves  n-vectors  ( unknownsIJ.m)
     containing the each  element of state vector;
     the index J defines the J-element of StateVector;
     the index I defines the I-th time that this method is called;
     */
     void spyStateVector();

     //! Spy rhs vector
     /*!
     Spy  of rhsVector;
     this method saves  n-vectors  (rhsIJ.m) containing the each  element of rhs vector;
     the index J defines the J-element of rhsVector;
     the index I defines the I-th time that this method is called;
     */
     void spyRHS();

     //@}

     //! @name Set Methods
     //@{

     //!Initialize parameters of time advance scheme;
     /*!
     Initialize parameters of time advance scheme;
     @param  orderDerivative  define the maximum  order of derivate
     */
     inline void setup ( const  UInt& orderDerivative ) { M_orderDerivative = orderDerivative;}

     //! Initialize the parameters of time advance scheme
     /*
     @param  order define the order of BDF;
     @param  orderDerivatve  define the order of derivate;
     */
     virtual void setup ( const UInt& order,  const  UInt& orderDerivative ) = 0;

     //! Initialize the parameters of time advance scheme
     /*
     @param  coefficients define the TimeAdvanceNewmark's coefficients (\theta, \gamma);
     @param  orderDerivative  define the order of derivate;
     */
     virtual void setup ( const std::vector<Real>&  coefficients, const  UInt& orderDerivative ) = 0;

     //! Initialize the State Vector
     /*!
     Initialize all the entries of the unknown vector to be derived with the vector x0 (duplicated).
     this class is virtual because used in bdf;
     @param x0 is the initial unk;
     */
     virtual void setInitialCondition( const feVectorType& x0) = 0;


     //! initialize the State Vector
     /*!
     Initialize all the entries of the unknown vector to be derived with the vector x0, v0 (duplicated).
     this class is virtual because used in \f$\theta\f$-method scheme;
     @param x0 is the initial unk;
     @param v0 is the initial velocity
     */
     virtual void setInitialCondition( const feVectorType& x0, const feVectorType& v0) = 0;

     //! initialize the StateVector
     /*!
     Initialize all the entries of the unknown vector to be derived with the vector x0, v0,w0 (duplicated).
     this class is virtual because used in Newamrk scheme;
     @param x0 is the initial unk;
     @param v0 is the initial velocity
     @param w0 is the initial accelerate
     */
     virtual void setInitialCondition( const feVectorType& x0, const feVectorType& v0, const feVectorType& w0) = 0;

    //! initialize the StateVector
    /*!
    Initialize all the entries of the unknown vector to be derived with the vector x0.
    this class is virtual because used in TimeAdvanceNewmark scheme;
    @param x0 is a vector of feVectorType containing the state vector;
    */
    virtual void setInitialCondition( const feVectorSharedPtrContainer_Type& x0){}

  //!Initialize the RhsVector:
    /*!
    Initialize all the entries of the unknown vector to be derived with the vector x0.
    this class is virtual because used in Newamrk scheme;
    @param rhs0 is a vector of feVectorType containing the state vector;
    */
    void setInitialRHS(const feVectorType & rhs0 ) ;

    //! Set time step
    /*!
    @param timeStep is time step
    */
    inline void setTimeStep(const Real& timeStep) {M_timeStep = timeStep; }

    //@}

   //!@name Get Methods
   //@{

   //! Return the i-th coefficient of the first time derivative
   /*!
   @param i index of coefficient alpha
   @returns the i-th coefficient of the time derivative alpha_i
   */
   Real coefficientFirstDerivative(const UInt& i) const;

   //!Return the \f$i-\f$th coefficient of the second time derivative
   /*
   @param \f$i\f$ index of coefficient \f$xi\f$
   @returns the i-th coefficient of the second time derivative \f$xi_i\f$
   */
    Real coefficientSecondDerivative(const UInt& i) const;

    //!Return the\f$ i-\f$th coefficient of the solution's extrapolation
    /*!
     @param \f$i\f$ index of  extrapolation coefficient
     @returns the \f$i-\f$th coefficient of the extrapolation of the first order derivative
    */
    virtual Real coefficientExtrapolation(const UInt& i )  const = 0;

    //!Returns the \f$i-\f$th coefficient of the velocity's extrap
    /*!
    @param i index of velocity's extrapolation  coefficient
    @returns the \f$i-\f$th coefficient of the velocity's extrapolation
   */

    virtual Real coefficientExtrapolationFirstDerivative(const UInt& i ) const =0;

    //! Compute the polynomial extrapolation of solution
    /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
    @returns  extrap of state vector u^*
    */
    virtual void extrapolation(feVectorType& extrapolation) const =0;

    //! Compute the polynomial extrapolation of solution
    /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
    @returns  extrap of state vector u^*
    */
    virtual void extrapolationFirstDerivative(feVectorType& extrapolation) const =0;


    //! Return the state vector
     /*!
    @returns the state vector
    */
    // inline const feVectorContainerPtr_Type unknowns() const {return this->M_unknowns;}

    //! Return the \f$i-\f$th element of state vector
    /*!
    @param \f$i\f$ index of element
    @returns the i-th element of state vector
    */
    const  feVectorType& singleElement(const UInt& i)  const;

    //! Return the last solution (the first element of state vector)
    const  feVectorType& solution()  const;

    void setSolution( feVectorType& solution )
    {
        delete M_unknowns[0];
        M_unknowns[0]= new feVectorType(solution);
    }


    //! Return the current velocity
    virtual  feVectorType velocity() const = 0;

    //!Return the velocity
    /*!
    @param u unk  to compute the current velocity;
    @returns the velocity associated to \f$u\f$
    this method is used for example in FSI to return the value of solid
    in the internal loop
    */
    feVectorType  velocity(const  feVectorType& u);

    //!Return the current accelerate
    virtual feVectorType accelerate() const =0;

    //! Return the accelerate
    /*!
    @param \f$u\f$ is necessary to compute wnk;
    @return the accelerate associated to \f$u\f$;
    this method is used for example in FSI to return the value of solid in the internal loop;
     */
    feVectorType accelerate(const  feVectorType& u);

    //!Return velocity's right hand side
    /*!
    @returns velocity's right hand side
    */
    inline const feVectorType& rhsContributionFirstDerivative() {return *M_rhsContribution[0];}

    //! Return accelerate's right hand side
    /*!
    @return accelerate's right hand side
    */
    inline const feVectorType& rhsContributionSecondDerivative(){return *M_rhsContribution[1];}

    //! Return order of accuracy of the scheme
    /*!
    @returns the order of accuracy of the scheme
    */
    inline UInt order() const  {return M_order;}

    //! Returns size of the stencil used by time integration scheme
    /*!
    @returns the size of the stencil
    */
    inline UInt size() const  {return M_size;}

    /*!Returns a pointer to the vector of solutions stored in M_unknowns
     */
    feVectorContainerPtr_Type& stencil() { return M_unknowns; }
    //@}

protected:

    //! Order of the BDF derivative/extrapolation: the time-derivative
    //! coefficients vector has size \f$n+1\f$, the extrapolation vector has size \f$n\f$
    UInt M_order;

    //! Order of temporal derivate: the time-derivative
    //! coefficients vector has size \f$n+1\f$, the extrapolation vector has size \f$n\f$
    UInt M_orderDerivative;

    //! time step
    Real M_timeStep;

    //! Size of the unknown vector
    UInt M_size;

    //! Size for firstOrderDerivative loop (for bdf  equal M_order, for TimeAdvanceNewmark equal  M_size/2)
    UInt M_firstOrderDerivativeSize;

    //! Size for setSecondOrderDerivatve loop  (for bdf  equal M_order, for TimeAdvanceNewmark equal M_size/2)
    UInt M_secondOrderDerivativeSize;

    //!Size of coefficients (for bdf equal M_order + M_orderDerivative, for theta-method is 3, and TimeAdvanceNewmark is 4)
    UInt M_coefficientsSize;

    //! Coefficients \f$ \alpha_i \f$ of the time advance discretization
    container_Type M_xi;

    //! Coefficients \f$ \alpha_i \f$ of the time advance discretization
    container_Type M_alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    container_Type M_beta;

    //! Coefficients \f$ \beta^V_i \f$ of the velocity's extrapolation
    container_Type M_betaFirstDerivative;

    //! Last n state vectors
    feVectorContainerPtr_Type M_unknowns;

    //! Vector of rhs (rhsV and rhsW)
    feVectorContainerPtr_Type M_rhsContribution;
};

// ===================================================
// Constructors & Destructor
// ===================================================

//! Empty Constructor
template<typename feVectorType>
TimeAdvance<feVectorType>::TimeAdvance()
        :
        M_unknowns(),
        M_rhsContribution()
{
    M_unknowns.reserve( 1 );
    M_rhsContribution.reserve(2);
}

//! Destructor
template<typename feVectorType>
TimeAdvance<feVectorType>::~TimeAdvance()
{
    feVectorContainerPtrIterate_Type iter     = M_unknowns.begin();
    feVectorContainerPtrIterate_Type iter_end = M_unknowns.end();

    for ( ; iter != iter_end; iter++ )
        delete *iter;
}

// ===================================================
// Methods
// ===================================================

template<typename feVectorType>
void
TimeAdvance<feVectorType>::
updateRHSContribution(const Real& timeStep )
{
  //! update rhsContribution  of the first Derivative
  this->updateRHSFirstDerivative( timeStep);

  //! update rhsContribution  of the second Derivative
  if( M_orderDerivative == 2 )
    this->updateRHSSecondDerivative( timeStep );
}


template<typename feVectorType>
void
TimeAdvance<feVectorType>::
spyStateVector()
{
    static UInt saveUnknowns=0;
    std::string unknowns="unknowns";

    for ( UInt i=0 ; i< M_size ; i++ )
    {
        std::ostringstream j;
        j<<saveUnknowns;
        j<<i;

        unknowns+j.str();

        M_unknowns[i]->spy(unknowns+j.str());
    }
    saveUnknowns++;
}

template<typename feVectorType>
void
TimeAdvance<feVectorType>::
spyRHS()
{
    static UInt saveRhs=0;
    std::string rhs="rhs";
    for ( UInt i=0 ; i< 2 ; ++i )
    {
        std::ostringstream j;
        j<<saveRhs;
        j<<i;

        rhs+j.str();
        M_rhsContribution[i]->spy(rhs+j.str());
    }
    saveRhs++;
}

// ===================================================
// Set Methods
// ===================================================

template<typename feVectorType>
void
TimeAdvance<feVectorType>::setInitialRHS(const feVectorType& rhs )
{
    for (UInt i=0; i<2; ++i )
    {
      M_rhsContribution.push_back(new feVectorType(rhs));
    }
}

// ===================================================
// Get Methods
// ===================================================
template<typename feVectorType>
Real
TimeAdvance<feVectorType>::coefficientFirstDerivative(const UInt& i) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i < M_coefficientsSize,
            "Error in specification of the time derivative coefficient for the time scheme" );
    return M_alpha[ i ];
}

template<typename feVectorType>
Real
TimeAdvance<feVectorType>::coefficientSecondDerivative( const UInt& i )  const
{ // you should replace any call to coef_derOrder2() with a call to coefficientSecondDerivative()
    // Pay attention: i is c-based indexed
    ASSERT( i < M_coefficientsSize,
            "Error in specification of the time derivative coefficient for the time scheme" );
    return M_xi[ i ];
}

template<typename feVectorType>
const feVectorType&
TimeAdvance<feVectorType>::solution() const
{
   return *M_unknowns[0];
}


template<typename feVectorType>
const feVectorType&
TimeAdvance<feVectorType>::singleElement( const UInt& i) const
{
// Pay attention: i is c-based indexed
    ASSERT( i < M_size,
            "Error there isn't unk(i), i must be shorter than M_size" );

    return *M_unknowns[i];
}

template<typename feVectorType>
feVectorType
TimeAdvance<feVectorType>::velocity( const feVectorType& u )
{  // you should replace any call to vnk() with a call to vnk()
    feVectorType vel( u );
    vel  *= M_alpha[ 0 ] /M_timeStep;
    vel  -= (*this->M_rhsContribution[ 0 ]);
    return vel;
}

template<typename feVectorType>
feVectorType
TimeAdvance<feVectorType>::accelerate(const feVectorType& u)
{
    feVectorType accelerate(u);
    accelerate  *= M_xi[ 0 ] / (M_timeStep *M_timeStep );
    accelerate  -= (*this->M_rhsContribution[1]);
    return accelerate;
}
// ===================================================
// Macros
// ===================================================

//! create factory
typedef FactorySingleton< Factory < TimeAdvance<>,  std::string> > TimeAdvanceFactory;

}
#endif  /* TIMEADVANCE_H */
