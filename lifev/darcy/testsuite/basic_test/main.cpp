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

/* ========================================================

Simple Darcy test with Dirichlet, Neumann and Robin boundary conditions

Solve the problem

               div u - f = 0            in \Omega

               K^{-1} u + \nabla p = 0  in \Omega

*/


/**
   @file main.hpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2010-07-29
*/


// ===================================================
//! Includes
// ===================================================

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

#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include "darcy.hpp"


// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}

// ===================================================
//! Main
// ===================================================
int main(int argc, char** argv)
{

#ifdef HAVE_MPI

    MPI_Init( &argc, &argv );

    std::cout << "MPI Initialization" << std::endl;

#endif

    // Error known
    const LifeV::Real errorKnown( 0.2003822844278755 );

    // Tolerance between the error and the errorKnown
    const LifeV::Real tolerance( 1e-8 );

    darcy Darcy( argc, argv );

    // Error of the problem
    const LifeV::Real error = Darcy.run();
    bool unsuccess=std::fabs( error - errorKnown ) > tolerance;
    // For tribits handling of success/failure
    //! @todo Add verbose to avoid all processes printing this stuff
    if (unsuccess)
      std::cout << "End Result: TEST NOT PASSED" << std::endl;
    else
      std::cout << "End Result: TEST PASSED" << std::endl;
#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif
    
    if (unsuccess)
      {
        return ( EXIT_FAILURE );
      }
    else
      {
        return ( EXIT_SUCCESS );
      }
}

