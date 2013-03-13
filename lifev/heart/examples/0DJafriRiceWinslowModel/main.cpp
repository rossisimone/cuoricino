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
    @brief 0D test with the Jafri Rice Winslow model

    @date 04−2013
    @author Luis Miguel de Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

    @contributor
    @mantainer Luis Miguel de Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
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

#include <lifev/heart/solver/IonicModels/IonicJafriRiceWinslow.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

using std::cout;
using std::endl;
using namespace LifeV;

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

    cout << "Importing parameters list...";
    Teuchos::ParameterList ParameterList = * ( Teuchos::getParametersFromXmlFile ( "JafriRiceWinslowParameters.xml" ) );
    cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    cout << "Building Constructor for JafriRiceWinslow Model with parameters ... ";
    IonicJafriRiceWinslow  model ( ParameterList );
    cout << " Done!" << endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    model.showMe();


    //********************************************//
    // Initialize the solution to 0. The model    //
    // consist of three state variables. Xe.Size()//
    // returns the number of state variables of   //
    // the model. rStates is the reference to the //
    // the vector states                          //
    //********************************************//
    cout << "Initializing solution vector...";
    std::vector<Real> unknowns (model.Size(), 0 );
    unknowns[0] = - 84.1638;
    unknowns[1] = 0.0328302;
    unknowns[2] = 0.988354;
    unknowns[3] = 0.992540;
    unknowns[4] = 0.000928836;
    unknowns[5] = 10.2042;
    unknowns[6] = 143.727;
    unknowns[7] = 140.0;
    unknowns[8] = 0.0000994893;
    unknowns[9] = 1.24891;
    unknowns[10] = 0.000136058;
    unknowns[11] = 1.17504;
    unknowns[12] = 0.762527;
    unknowns[13] = 0.00119168;
    unknowns[14] = 0.00000000630613;
    unknowns[15] = 0.236283;
    unknowns[16] = 0.997208;
    unknowns[17] = 0.0000638897;
    unknowns[18] = 0.00000000153500;
    unknowns[19] = 0.0000000000000163909;
    unknowns[20] = 0.0000000000000000000656337;
    unknowns[21] = 0.00000000000000000000984546;
    unknowns[22] = 0.00272826;
    unknowns[23] = 0.000000699215;
    unknowns[24] = 0.0000000000671989;
    unknowns[25] = 0.00000000000000287031;
    unknowns[26] = 0.0000000000000000000459752;
    unknowns[27] = 0.0;
    unknowns[28] = 0.998983;
    unknowns[29] = 0.6349973;
    unknowns[30] = 135.9813;
    cout << " Done!" << endl;


    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//
    cout << "Initializing rhs..." ;
    std::vector<Real> rhs (model.Size(), 0);
    cout << " Done! "  << endl;


    //********************************************//
    // The model needs as external informations   //
    // the contraction velocity and the Calcium   //
    // concentration.                             //
    //********************************************//
    Real Iapp (0.0);


    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//
    Real TF (0.04);
    Real dt (0.02);


    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    string filename = "output.txt";
    std::ofstream output ("output.txt");


    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    cout << "Time loop starts...\n";
    for ( Real t = 0; t < TF; )
    {

    	//********************************************//
        // Compute Calcium concentration. Here it is  //
        // given as a function of time.               //
        //********************************************//
        if ( t > 1 && t < 2 )
        {
        	Iapp = 0.1;
        }
        else
        {
        	Iapp = 0;
        }

        cout << "\r " << t << " ms.       " << std::flush;

         //********************************************//
         // Compute the rhs using the model equations  //
         //********************************************//
         model.computeRhs ( unknowns, Iapp, rhs );

         int i(0);

         for(std::vector<Real>::iterator it=rhs.begin(); it!=rhs.end(); ++it)
             {
        	 	 cout << i << "/ " << *it << ", ";
        	 	 i = i + 1;
             }
         cout << "\n" << rhs.size() << endl;


         //********************************************//
         // Use forward Euler method to advance the    //
         // solution in time.                          //
         //********************************************//
         unknowns.at (0) = unknowns.at (0)   + dt * rhs.at (0);
         unknowns.at (1) = unknowns.at (1)   + dt * rhs.at (1);
         unknowns.at (2) = unknowns.at (2)   + dt * rhs.at (2);
         unknowns.at (3) = unknowns.at (3)   + dt * rhs.at (3);
         unknowns.at (4) = unknowns.at (4)   + dt * rhs.at (4);
         unknowns.at (5) = unknowns.at (5)   + dt * rhs.at (5);
         unknowns.at (6) = unknowns.at (6)   + dt * rhs.at (6);
         unknowns.at (7) = unknowns.at (7)   + dt * rhs.at (7);
         unknowns.at (8) = unknowns.at (8)   + dt * rhs.at (8);
         unknowns.at (9) = unknowns.at (9)   + dt * rhs.at (9);
         unknowns.at (10) = unknowns.at (10) + dt * rhs.at (10);
         unknowns.at (11) = unknowns.at (11) + dt * rhs.at (11);
         unknowns.at (12) = unknowns.at (12) + dt * rhs.at (12);
         unknowns.at (13) = unknowns.at (13) + dt * rhs.at (13);
         unknowns.at (14) = unknowns.at (14) + dt * rhs.at (14);
         unknowns.at (15) = unknowns.at (15) + dt * rhs.at (15);
         unknowns.at (16) = unknowns.at (16) + dt * rhs.at (16);
         unknowns.at (17) = unknowns.at (17) + dt * rhs.at (17);
         unknowns.at (18) = unknowns.at (18) + dt * rhs.at (18);
         unknowns.at (19) = unknowns.at (19) + dt * rhs.at (19);
         unknowns.at (20) = unknowns.at (20) + dt * rhs.at (20);
         unknowns.at (21) = unknowns.at (21) + dt * rhs.at (21);
         unknowns.at (22) = unknowns.at (22) + dt * rhs.at (22);
         unknowns.at (23) = unknowns.at (23) + dt * rhs.at (23);
         unknowns.at (24) = unknowns.at (24) + dt * rhs.at (24);
         unknowns.at (25) = unknowns.at (25) + dt * rhs.at (25);
         unknowns.at (26) = unknowns.at (26) + dt * rhs.at (26);
         unknowns.at (27) = unknowns.at (27) + dt * rhs.at (27);
         unknowns.at (28) = unknowns.at (28) + dt * rhs.at (28);
         unknowns.at (29) = unknowns.at (29) + dt * rhs.at (29);
         unknowns.at (30) = unknowns.at (30) + dt * rhs.at (30);


         //********************************************//
         // Writes solution on file.                   //
         //********************************************//
         output << t << ", " << unknowns.at (0) << ", " << unknowns.at (1) << ", "
        		 << unknowns.at (2) << ", " << unknowns.at (3) << ", "
        		 << unknowns.at (4) << ", " << unknowns.at (5) << ", "
        		 << unknowns.at (6) << ", " << unknowns.at (7) << ", "
        		 << unknowns.at (8) << ", " << unknowns.at (9) << ", "
        		 << unknowns.at (10) << ", " << unknowns.at (11) << ", "
        		 << unknowns.at (12) << ", " << unknowns.at (13) << ", "
        		 << unknowns.at (14) << ", " << unknowns.at (15) << ", "
        		 << unknowns.at (16) << ", " << unknowns.at (17) << ", "
        		 << unknowns.at (18) << ", " << unknowns.at (19) << ", "
        		 << unknowns.at (20) << ", " << unknowns.at (21) << ", "
        		 << unknowns.at (22) << ", " << unknowns.at (23) << ", "
        		 << unknowns.at (24) << ", " << unknowns.at (25) << ", "
        		 << unknowns.at (26) << ", " << unknowns.at (27) << ", "
        		 << unknowns.at (28) << ", " << unknowns.at (29) << "\n";

         //********************************************//
         // Update the time.                           //
         //********************************************//
         t = t + dt;
       }

    cout << "\n...Time loop ends.\n";
    cout << "Solution written on file: " << filename << "\n";

    //********************************************//
    // Close exported file.                       //
    //********************************************//

    output.close();


    //! Finalizing Epetra communicator
    MPI_Finalize();
    return ( EXIT_SUCCESS );
   }
