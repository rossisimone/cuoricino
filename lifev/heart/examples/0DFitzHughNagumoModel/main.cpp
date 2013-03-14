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
    @brief 0D test with the Fitz-Hugh Nagumo model.

    @date 01−2013
    @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

    @contributor
    @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
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


#include <fstream>
#include <string>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/heart/solver/IonicModels/IonicFitzHughNagumo.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using std::cout;
using std::endl;
using namespace LifeV;

void EulerExplicit (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output);
void ChebychevStabilized (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output);

Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }


    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    std::cout << "Importing parameters list...";
    Teuchos::ParameterList FHNParameterList = * ( Teuchos::getParametersFromXmlFile ( "FitzHughNagumoParameters.xml" ) );
    std::cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Fitz-Hugh Nagumo. The           //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for Fitz-Hugh Nagumo Model with parameters ... ";
    IonicFitzHughNagumo  model ( FHNParameterList );
    std::cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    model.showMe();

    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//

    Real TF (FHNParameterList.get ("TF", 300.0) );
    Real dt (FHNParameterList.get ("dt", 0.01) );

    cout << "Time parameters : " << endl;
    cout << "TF = " << TF << endl;
    cout << "dt = " << dt << endl;


    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    string filename = "output.txt";
    std::ofstream output ("output.txt");


    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";

    if ( FHNParameterList.get ("Cheby", 1.0) == 0.0 )
    {
        EulerExplicit (dt, TF, model, FHNParameterList.get ("Iapp", 2000.0), output);
    }
    else
    {
        ChebychevStabilized (dt, TF, model, FHNParameterList.get ("Iapp", 2000.0), output);
    }

    std::cout << "\n...Time loop ends.\n";
    std::cout << "Solution written on file: " << filename << "\n";
    //********************************************//
    // Close exported file.                       //
    //********************************************//
    output.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}

void EulerExplicit (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output)
{
    std::vector<Real> unknowns ( model.Size(), 0.0);
    std::vector<Real> rhs ( model.Size(), 0.0);
    Real Iapp;

    cout << "Computing using Explicit Euler" << endl;

    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if ( t > 10.0 && t < 10.1 )
        {
            Iapp = I;
        }
        else
        {
            Iapp = 0.0;
        }

        std::cout << "\r " << t << " ms.       " << std::flush;

        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        model.computeRhs ( unknowns, Iapp, rhs);

        //********************************************//
        // Use forward Euler method to advance the    //
        // solution in time.                          //
        //********************************************//
        unknowns.at (0) = unknowns.at (0)  + dt * rhs.at (0);
        unknowns.at (1) = unknowns.at (1)  + dt * rhs.at (1);


        //********************************************//
        // Writes solution on file.                   //
        //********************************************//
        output << t << ", " << unknowns.at (0) << ", " << unknowns.at (1) << "\n";

        //********************************************//
        // Update the time.                           //
        //********************************************//
        t = t + dt;
    }
}

void ChebychevStabilized (Real& dt, const Real& TF, IonicFitzHughNagumo model, const Real& I, std::ofstream& output)
{
    Real x1, x2;
    Real y1, y2;
    Real r, u;
    Real a, b, c, d;
    Real s = 3.0;
    Real s2 = 9.0;
    std::vector<Real> g1 (2, 0.0);
    std::vector<Real> g2 (2, 0.0);
    std::vector<Real> g3 (2, 0.0);
    std::vector<Real> rhs (model.Size(), 0.0);
    std::vector<Real> unknowns ( model.Size(), 0.0);
    Real Iapp;

    cout << "Computing using Chebychev Stabilized" << endl;

    for ( Real t = 0.0; t < TF; )
    {
        //********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if ( t >= 10.0 && t < 10.1 )
        {
            Iapp = I;
        }
        else
        {
            Iapp = 0.0;
        }

        std::cout << "\r " << t << " ms.       " << std::flush;


        //Computing the Jacobian, J = [a b ; c d]
        model.computeJ (a, b, c, d, unknowns);

        //Iterations of PowerMethod to approximate the dominant eigenvalue of the Jacobian
        //Works fine, exactly the same values as MATLAB
        while (r > 0.0001)
        {
            // Computing y = J * x
            y1 = a * x1 + b * x2;
            y2 = c * x1 + d * x2;

            //Normalizing y and setting x = y
            r = sqrt (y1 * y1 + y2 * y2);
            x1 = y1 / r;
            x2 = y2 / r;

            //Approximating the dominant eigenvalue u with the Rayleight quotient
            u = x1 * (a * x1 + b * x2) + x2 * (c * x1 + d * x2);

            //Computing the residual J*x - u*x
            y1 = a * x1 + b * x2 - u * x1;
            y2 = c * x1 + d * x2 - u * x2;

            //Squared norm of the residual
            r = y1 * y1 + y2 * y2;
        }

        //Approximation of spectral radius
        u = abs (u);

        if ( u < s2 / 0.02 )
        {
            u = s2 / 0.02;
        }

        //New stepsize
        dt = s2 / u;

        if ( t < 10.0 && t + dt >= 10.1)
        {
            dt = 10.0 - t;
        }

        cout << "Iteration : " << t << endl;
        cout << "Spectral Radius : " << u << endl;
        cout << "Step Size : " << dt << endl << endl;

        g1 = unknowns;

        model.computeRhs (g1, Iapp, rhs);
        g2.at (0) = g1.at (0) + (dt / s2) * rhs.at (0);
        g2.at (1) = g1.at (1) + (dt / s2) * rhs.at (1);

        for (int j = 2; j <= s; j++)
        {
            model.computeRhs (g2, Iapp, rhs);
            g3.at (0) = (2.0 * dt / s2) * rhs.at (0) + 2.0 * g2.at (0) - g1.at (0) ;
            g3.at (1) = (2.0 * dt / s2) * rhs.at (1) + 2.0 * g2.at (1) - g1.at (1) ;

            g2 = g3;
            g1 = g2;
        }

        unknowns.at (0) = g3.at (0);
        unknowns.at (1) = g3.at (1);


        //********************************************//
        // Writes solution on file.                   //
        //********************************************//
        output << t << ", " << unknowns.at (0) << ", " << unknowns.at (1) << ", " << u << ", " << dt << "\n";

        //********************************************//
        // Update the time.                           //
        //********************************************//

        t = t + dt;
    }
}
