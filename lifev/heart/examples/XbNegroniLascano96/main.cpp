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
    @brief main for the test_heart

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

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

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/heart/solver/XbModels/HeartXbModel.hpp>
#include <lifev/heart/solver/XbModels/XbNegroniLascano96.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace LifeV;

Int main( Int argc, char** argv )
{
	   //! Initializing Epetra communicator
	    MPI_Init(&argc, &argv);
	    Epetra_MpiComm Comm(MPI_COMM_WORLD);
	    if ( Comm.MyPID() == 0 )
	        cout << "% using MPI" << endl;


	std::cout << "Importing parameters list.. ";
    Teuchos::ParameterList NLParameterList = *( Teuchos::getParametersFromXmlFile( "NegroniLascano96Parameters.xml" ) );
	std::cout << " Done!" << endl;

	std::cout << NLParameterList.get("alpha1", 0.0) << endl;



	std::cout << "Empty constructor HeartXbSolver ... ";
	HeartXbModel Xb;
	std::cout << " Done!" << endl;

	std::cout << "Constructor HeartXbSolver ... ";
	HeartXbModel Xx( Comm );
	std::cout << " Done!" << endl;

	std::cout << "Copy constructor ... ";
	HeartXbModel Xd = Xx;
	std::cout << " Done!" << endl;

	std::cout << "0perator = ... ";
	Xb = Xx;
	std::cout << " Done!" << endl;

	std::cout << "New operator  ... ";
	HeartXbModel * Xc = new HeartXbModel( Comm );
	std::cout << " Done!" << endl;


	std::cout << "Empty constructor NL96 ... ";
	XbNegroniLascano96 * Xq = new XbNegroniLascano96();
	std::cout << " Done!" << endl;

	std::cout << "Comm Constructor NL96 ... ";
	XbNegroniLascano96 * Xw = new XbNegroniLascano96( Comm );
	std::cout << " Done!" << endl;

	std::cout << "Parameter ConstructorNL96 ... ";
	XbNegroniLascano96 * Xe = new XbNegroniLascano96( Comm, NLParameterList );
	std::cout << " Done!" << endl;

	std::cout << "Operator =  ... ";
	Xq = Xe;
	std::cout << " Done!" << endl;

	std::cout << "Copy Constructor NL96 ... ";
	XbNegroniLascano96 * Xr = Xq;
	std::cout << " Done!" << endl;

	std::cout << "alpha1: " << Xr->Alpha1() << endl;
	std::cout << "alpha2: " << Xr->Alpha2() << endl;
	std::cout << "alpha3: " << Xr->Alpha3() << endl;
	std::cout << "alpha4: " << Xr->Alpha4() << endl;
	std::cout << "alpha5: " << Xr->Alpha5() << endl;
	std::cout << "beta1: " << Xr->Beta1() << endl;
	std::cout << "beta2: " << Xr->Beta2() << endl;
	std::cout << "beta3: " << Xr->Beta3() << endl;



    //! Finalizing Epetra communicator
    MPI_Finalize();
    return( EXIT_SUCCESS );
}
