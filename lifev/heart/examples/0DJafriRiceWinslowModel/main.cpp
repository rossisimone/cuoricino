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
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "JafriRiceWinslowParameters.xml" ) );
    cout << " Done!" << endl;


    //********************************************//
    // Creates a new model object representing the//
    // model from Negroni and Lascano 1996. The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    cout << "Building Constructor for JafriRiceWinslow Model with parameters ... ";
    IonicJafriRiceWinslow  model ( parameterList );
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
    unknowns[7] = 5.4;
    unknowns[8] = 9.94893e-11;
    unknowns[9] = 1.243891;
    unknowns[10] = 1.36058e-4;
    unknowns[11] = 1.17504;
    unknowns[12] = 0.762527;
    unknowns[13] = 1.19168e-3;
    unknowns[14] = 6.30613e-9;
    unknowns[15] = 0.236283;
    unknowns[16] = 0.997208;
    unknowns[17] = 6.38897e-5;
    unknowns[18] = 1.535e-9;
    unknowns[19] = 1.63909e-14;
    unknowns[20] = 6.56337e-20;
    unknowns[21] = 9.84546e-21;
    unknowns[22] = 2.72826e-3;
    unknowns[23] = 6.99215e-7;
    unknowns[24] = 6.71989e-11;
    unknowns[25] = 2.87031e-15;
    unknowns[26] = 4.59752e-20;
    unknowns[27] = 0.0;
    unknowns[28] = 0.998983;
    unknowns[29] = 0.00635;
    unknowns[30] = 0.13598;
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
    Real TF = parameterList.get      ( "endTime", 5.0 );
    Real dt = parameterList.get      ( "timeStep", 5.77e-5 );
    Real firstst = parameterList.get ( "firstStimuliTime", 1.0 );
    Real st = parameterList.get      ( "stimuliTime", 400.0 );

    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    string filename = "output.txt";
    std::ofstream output  ("output.txt");


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
        if ( t > firstst && t < firstst + 1 )
        {
        	Iapp    = 0.516289;
        	firstst = firstst + st;
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
             << unknowns.at (28) << ", " << unknowns.at (29) << ", "
             << unknowns.at (30) << "\n";



         //********************************************//
         // Use forward Euler method to advance the    //
         // solution in time.                          //
         //********************************************//

    	 for(int j(0); j <= 30; ++j)
         {
//    		 if( ( t > 21 && t < 26.5 ) || ( t > 421 && t < 426.5 ) ||
//    	        		( t > 821 && t < 826.5 ) || ( t > 1221 && t < 1226.5 ) || ( t > 1621 && t < 1626.5 ) ||
//    	        		( t > 2001 && t < 2036.5 ) || ( t > 2401 && t < 2436.5 ) || ( t > 2801 && t < 2866.5 ) )
//    		 {
				 if(j!= 10)
					 unknowns.at (j) = unknowns.at (j)   + dt * rhs.at (j);
				 else
				 {
//					 for( int k(0) ; k < 150; k++ )
//					 {
//						 unknowns.at (10) = unknowns.at (10)   + dt / 150 * rhs.at (10);
//						 model.computeRhs ( unknowns, Iapp, rhs );
//					 }
					 unknowns.at (10) = model.computeNewtonCaSS(unknowns, dt, 10);
				 }
//    		 }
//    		 else unknowns.at (j) = unknowns.at (j)   + dt * rhs.at (j);
         }

//         unknowns.at(1) = ( unknowns.at(1) / dt + model.fastINa(unknowns).at(1) )
//        		 / ( 1 / dt + model.fastINa(unknowns).at(1) + model.fastINa(unknowns).at(2) );
//      	 unknowns.at(2) = ( unknowns.at(2) / dt + model.fastINa(unknowns).at(3) )
//      			 / ( 1 / dt + model.fastINa(unknowns).at(3) + model.fastINa(unknowns).at(5) );
//    	 unknowns.at(3) = ( unknowns.at(3) / dt + model.fastINa(unknowns).at(4) )
//    			 / ( 1 / dt + model.fastINa(unknowns).at(4) + model.fastINa(unknowns).at(6) );
//    	 unknowns.at(4) = ( unknowns.at(4) / dt + model.timeDIK(unknowns).at(1) )
//    			 / ( 1 / dt + model.timeDIK(unknowns).at(1) + model.timeDIK(unknowns).at(2) );


//    	 unknowns.at (28) = ( unknowns.at(28) / dt + model.computeYParameters(unknowns).at(0) /  model.computeYParameters(unknowns).at(1) )
//    			 / ( 1 / dt + 1 /  model.computeYParameters(unknowns).at(1) );


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
