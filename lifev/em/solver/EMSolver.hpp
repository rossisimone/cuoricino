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

#ifndef _EMSOLVER_H_
#define _EMSOLVER_H_

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

#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

namespace LifeV {

//! EMSolver - Class featuring the solution of the electromechanical problem with monodomain equation

template<typename Mesh, typename IonicModel>
class EMSolver {

	//!Monodomain Solver
	/*!
	 The monodomain equation reads
	 \f \Chi

	 */

public:

	//! @name Type definitions
	//@{

	typedef Mesh mesh_Type;
	typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

	typedef VectorEpetra vector_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

	typedef std::vector<vectorPtr_Type> vectorOfPtr_Type;

	typedef MatrixEpetra<Real> matrix_Type;
	typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<comm_Type> commPtr_Type;

	typedef ETFESpace<mesh_Type, MapEpetra, 3, 1> ETFESpace_Type;
	typedef boost::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 1> > ETFESpacePtr_Type;

	typedef ETFESpace<mesh_Type, MapEpetra, 3, 3> ETFESpaceVectorial_Type;
	typedef boost::shared_ptr< ETFESpaceVectorial_Type > ETFESpaceVectorialPtr_Type;

	typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
	typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

	typedef LinearSolver linearSolver_Type;
	typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;

//    typedef ExporterHDF5< mesh_Type >          exporter_Type;
//    typedef boost::shared_ptr<exporter_Type>                       exporterPtr_Type;
	typedef Exporter<mesh_Type> exporter_Type;    //                IOFile_Type;
	typedef boost::shared_ptr<exporter_Type> exporterPtr_Type; //                IOFilePtr_Type;

	typedef LifeV::Preconditioner basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack prec_Type;
	typedef boost::shared_ptr<prec_Type> precPtr_Type;

	typedef IonicModel ionicModel_Type;
	typedef ElectroIonicModel superIonicModel;
	typedef boost::shared_ptr<ionicModel_Type> ionicModelPtr_Type;

	typedef Teuchos::ParameterList list_Type;

	typedef boost::function<
			Real(const Real& t, const Real& x, const Real& y, const Real& z,
					const ID& i)> function_Type;

	typedef MatrixSmall<3, 3>                          matrixSmall_Type;


	//@}

	//! @name Constructors & Destructor
	//@{

	//!Empty Constructor
	/*!
	 */


};

} // namespace LifeV

#endif //_MONODOMAINSOLVER_H_
