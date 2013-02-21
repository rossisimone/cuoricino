/* -*- mode: c++ -*- */
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

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/BCManage.hpp>

#include <lifev/fsi/solver/MonolithicBlockComposedNN.hpp>

namespace LifeV
{

// ===================================================
//! Public Methods
// ===================================================

int MonolithicBlockComposedNN::solveSystem ( const vector_Type& rhs, vector_Type& step, solverPtr_Type& linearSolver )
{
    M_firstCompPrec .reset (new composed_prec (M_comm) );
    M_secondCompPrec.reset (new composed_prec (M_comm) );

    LifeChrono chrono;

    if ( !M_blockPrecs.get() )
    {
        M_blockPrecs.reset (new ComposedOperator< composed_prec > (M_comm) );
    }
    if ( !set() )
    {
        M_overlapLevel     = M_list.get ("overlap level", 2);
        M_precType = M_list.get ("prectype", "Amesos");

        //////////////// \todo Copy: should be avoided ////////////////////////
        M_matrixVector.push_back (*M_blocks[ (*M_blockReordering) [0]]);
        M_matrixVector.push_back (*M_blocks[ (*M_blockReordering) [1]]);
        M_matrixVector.push_back (*M_blocks[ (*M_blockReordering) [2]]);
        M_matrixVector.push_back (*M_blocks[ (*M_blockReordering) [3]]);
        /////////////////////////////////////////////////////////////

        for (ID k (0); k < M_blocks.size(); ++k)
        {
            M_blockPrecs->displayer().leaderPrint ("  M-  Computing double prec. factorization ...        ");
            chrono.start();
            //////////////// \todo Copy: should be avoided ////////////////////////
            if (!M_recompute[ (*M_blockReordering) [k]])
            {
                M_prec[k].reset (M_factory.Create (M_precType, M_matrixVector[k].matrixPtr().get(), M_overlapLevel) );
            }
            else
                //////////////////////////////////////////////////////////////

            {
                M_prec[k].reset (M_factory.Create (M_precType, M_blocks[ (*M_blockReordering) [k]]->matrixPtr().get(), M_overlapLevel) );
            }
            if ( !M_prec[k].get() )
            {
                ERROR_MSG ( "Preconditioner not set, something went wrong in its computation\n" );
            }
            M_prec[k]->SetParameters (M_list);
            M_prec[k]->Initialize();
            M_prec[k]->Compute();
            chrono.stop();
            M_blockPrecs->displayer().leaderPrintMax ("done in ", chrono.diff() );
        }
    }
    else
    {
        for (ID k (0); k < M_blocks.size(); ++k)
        {
            if (M_recompute[ (*M_blockReordering) [k]])
            {
                M_blockPrecs->displayer().leaderPrint ("  M-  Computing double prec. factorization ...        ");
                chrono.start();
                M_prec[k].reset (M_factory.Create (M_precType, M_blocks[ (*M_blockReordering) [k]]->matrixPtr().get(), M_overlapLevel) );
                if ( !M_prec[k].get() )
                {
                    ERROR_MSG ( "Preconditioner not set, something went wrong in its computation\n" );
                }
                M_prec[k]->SetParameters (M_list);
                M_prec[k]->Initialize();
                M_prec[k]->Compute();
                chrono.stop();
                M_blockPrecs->displayer().leaderPrintMax ("done in ", chrono.diff() );
            }
            else
            {
                M_blockPrecs->displayer().leaderPrint ("  M-  Reusing double prec. factorization ...        \n");
            }
        }
    }
    M_firstCompPrec->push_back (M_prec[0], false, false);
    M_firstCompPrec->push_back (M_prec[1], false, false);

    M_secondCompPrec->push_back (M_prec[2], false, false);
    M_secondCompPrec->push_back (M_prec[3], false, false);

    //M_blockPrecs->resetOperator();
    if (! (M_blockPrecs->number() ) )
    {
        M_blockPrecs->push_back (M_firstCompPrec, false, false, false);
        M_blockPrecs->push_back (M_secondCompPrec, false, false, true /*sum*/);
    }
    else
    {
        M_blockPrecs->replace (M_firstCompPrec, (UInt) 0, false, false);
        M_blockPrecs->replace (M_secondCompPrec, (UInt) 1, false, false);
    }

    return linearSolver->solveSystem (rhs, step, boost::static_pointer_cast<Epetra_Operator> (M_blockPrecs) );
}



void MonolithicBlockComposedNN::setDataFromGetPot (const GetPot& data, const std::string& section)
{
    PreconditionerIfpack::createIfpackList ( M_list, data, section, "ifpack");
}

void MonolithicBlockComposedNN::coupler (mapPtr_Type& map,
                                         const std::map<ID, ID>& locDofMap,
                                         const vectorPtr_Type& numerationInterface,
                                         const Real& timeStep,
                                         const Real& coefficient,
                                         const Real& rescaleFactor)
{
    UInt totalDofs = map->map (Unique)->NumGlobalElements();
    UInt fluidSolid = M_offset[0] + M_FESpace[0]->map().map (Unique)->NumGlobalElements();

    for (ID k = 0; k < 2; ++k)
    {
        M_blocks[k]->globalAssemble();
        matrixPtr_Type block (new matrix_Type (*M_blocks[k]) );
        M_blocks.push_back (block);
        M_bch.push_back (M_bch[k]);
        M_FESpace.push_back (M_FESpace[k]);
        M_offset.push_back (M_offset[k]);
        M_recompute[2 + k] = (M_recompute[k]);
    }

    matrixPtr_Type coupling (new matrix_Type (*map) );

    coupling.reset (new matrix_Type (*map, 0) );
    coupling->insertValueDiagonal (1., M_offset[fluid], M_offset[solid] );
    coupling->insertValueDiagonal (1.,  fluidSolid, totalDofs);
    couplingMatrix (coupling, (*M_couplingFlags) [0]/*8*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor);
    M_coupling.push_back (coupling);

    coupling.reset (new matrix_Type (*map, 0) );
    coupling->insertValueDiagonal ( 1., M_offset[0], fluidSolid );
    coupling->insertValueDiagonal ( 1., fluidSolid , totalDofs);
    couplingMatrix (coupling, (*M_couplingFlags) [1]/*4*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor);
    M_coupling.push_back (coupling);

    coupling.reset (new matrix_Type (*map, 0) );
    coupling->insertValueDiagonal (1., M_offset[fluid], M_offset[solid] );
    coupling->insertValueDiagonal (-1,  fluidSolid, totalDofs);
    couplingMatrix (coupling, (*M_couplingFlags) [2]/*1*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor);
    M_coupling.push_back (coupling);

    coupling.reset (new matrix_Type (*map, 0) );
    coupling->insertValueDiagonal ( 1., M_offset[0], fluidSolid );
    coupling->insertValueDiagonal ( 1., fluidSolid, totalDofs );
    couplingMatrix (coupling, (*M_couplingFlags) [3]/*2*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor);
    M_coupling.push_back (coupling);

    M_prec.resize (M_blocks.size() );
}

void MonolithicBlockComposedNN::applyBoundaryConditions (const Real& time, const UInt i)
{
    M_blocks[i]->openCrsMatrix();
    if ( !M_bch[i]->bcUpdateDone() )
    {
        M_bch[i]->bcUpdate ( *M_FESpace[i]->mesh(), M_FESpace[i]->feBd(), M_FESpace[i]->dof() );
        M_bch[i]->setOffset (M_offset[i]);
    }
    bcManageMatrix ( *M_blocks[i] , *M_FESpace[i]->mesh(), M_FESpace[i]->dof(), *M_bch[i], M_FESpace[i]->feBd(), 2., time);
}


void MonolithicBlockComposedNN::push_back_matrix (const matrixPtr_Type& Mat, const  bool recompute)
{
    Mat->globalAssemble();
    *Mat *= 2.;
    super_Type::push_back_matrix (Mat, recompute);
}



void MonolithicBlockComposedNN::replace_matrix ( const matrixPtr_Type& oper, UInt position )
{
    oper->globalAssemble();
    *oper *= 2.;
    M_blocks[position] = oper;
    M_blocks[ (position + 2) % 4] = oper;
}

} // Namespace LifeV
