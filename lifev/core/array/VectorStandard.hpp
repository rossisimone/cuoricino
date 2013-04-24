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
    @brief Implemetation of a small vector used with MatrixStandard

    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
    @maintainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @date 01-04-2013
 */

#ifndef VectorStandard_HPP_
#define VectorStandard_HPP_

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include <lifev/core/LifeV.hpp>


namespace LifeV
{

class VectorStandard : public std::vector<Real>
{

public:
	VectorStandard(){};
	VectorStandard(UInt size, Real r);
	virtual ~VectorStandard(){};

	Real norm2() const;
	void disp() const;

	VectorStandard operator+ (const VectorStandard& w) const;
	VectorStandard& operator+=(const VectorStandard& w);
	VectorStandard operator- (const VectorStandard& w) const;
	VectorStandard& operator-=(const VectorStandard& w);
	VectorStandard operator* (const Real r) const;
	VectorStandard& operator*=(const Real r);
	VectorStandard operator/ (const Real r) const;
	VectorStandard& operator/=(const Real r);

};

}// namespace LifeV



#endif /* VectorStandard_HPP_ */
