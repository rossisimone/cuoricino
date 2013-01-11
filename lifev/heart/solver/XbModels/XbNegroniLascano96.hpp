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
  @brief Mean-field cross-bridge model by Negroni and Lascano 1996.
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _XBNEGRONILASCANO96_H_
#define _XBNEGRONILASCANO96_H_

#include <lifev/heart/solver/XbModels/HeartXbModel.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! XbModel - This class implements a mean field model.


class XbNegroniLascano96 : public virtual HeartXbModel
{

public:
    //! @name Type definitions
    //@{

    typedef VectorEpetra vector_Type;

 //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    XbNegroniLascano96();
    /*!
     * @param Epetra communicator
     */
    XbNegroniLascano96( Epetra_Comm& comm );
    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    XbNegroniLascano96( Epetra_Comm& comm,
    					  Teuchos::ParameterList& parameterList );

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     * @param State variables
     */
//    XbNegroniLascano96(  Epetra_Comm& comm,
//    						const Teuchos::ParameterList& parameterList,
//                			const vector_Type& TCa,
//                			const vector_Type& TCas,
//                			const vector_Type& Ts);
    /*!
     * @param XbNegroniLascano96 object
     */
    XbNegroniLascano96( const XbNegroniLascano96 & Xb );
    //! Destructor
    virtual ~XbNegroniLascano96() {}

 //@}

    //! @name Overloads
    //@{

    XbNegroniLascano96& operator=( const XbNegroniLascano96 &Xb );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& Alpha1()	const { return M_alpha1; 	}
    inline const Real& Alpha2()	const { return M_alpha2; 	}
    inline const Real& Alpha3() 	const { return M_alpha3; 	}
    inline const Real& Alpha4() 	const { return M_alpha4; 	}
    inline const Real& Alpha5() 	const { return M_alpha5; 	}
    inline const Real& Beta1()  	const { return M_beta1; 	}
    inline const Real& Beta2()  	const { return M_beta2; 	}
    inline const Real& Beta3()  	const { return M_beta3; 	}
//
//    inline void setAlpha1	( const Real &alpha1 ) { this->M_alpha1 = alpha1; 	}
//    inline void setAlpha2	( const Real &alpha2 ) { this->M_alpha2 = alpha2; 	}
//    inline void setAlpha3	( const Real &alpha3 ) { this->M_alpha3 = alpha3; 	}
//    inline void setAlpha4	( const Real &alpha4 ) { this->M_alpha4 = alpha4; 	}
//    inline void setAlpha5	( const Real &alpha5 ) { this->M_alpha5 = alpha5; 	}
//    inline void setBeta1 	( const Real &beta1  ) { this->M_beta1 	=  beta1; 	}
//    inline void setBeta2 	( const Real &beta2  ) { this->M_beta2 	=  beta2; 	}
//    inline void setBeta3 	( const Real &beta3  ) { this->M_beta3 	=  beta3; 	}
//
//
//    //members getters and setters
//
//    inline const vector_Type T()  	  const { return M_T; 		}
//    inline const vector_Type TCa()	  const { return M_TCa; 	}
//    inline const vector_Type TCas() const { return M_TCas; 	}
//    inline const vector_Type Ts() 	  const { return M_Ts; 	}
//
//    inline const void setT	( const vector_Type &T 		) { this->M_T 		= T; 	}
//    inline const void setTCa	( const vector_Type &TCa  	) { this->M_TCa		= TCa; 	}
//    inline const void setTCas( const vector_Type &TCas	) { this->M_TCas	= TCas;	}
//    inline const void setTs	( const vector_Type &Ts 	) { this->M_Ts 		= Ts; 	}
//
//    //repeated members getters and setters
//
//    inline const vector_Type TRepeated()  	const { return M_TRepeated; 	}
//    inline const vector_Type TCaRepeated()	const { return M_TCaRepeated; 	}
//    inline const vector_Type TCasRepeated() 	const { return M_TCasRepeated;	}
//    inline const vector_Type TsRepeated() 	const { return M_TsRepeated; 	}
//
//    inline const void setTRepeated		( vector_Type const &TRepeated 		) { this->M_TRepeated 		= TRepeated; 	}
//    inline const void setTCaRepeated		( vector_Type const &TCaRepeated  	) { this->M_TCaRepeated 	= TCaRepeated;	}
//    inline const void setTCasRepeated	( vector_Type const &TCaRepeateds	) { this->M_TCasRepeated	= TCasRepeated;	}
//    inline const void setTsRepeated		( vector_Type const &TsRepeated 	) { this->M_TsRepeated 		= TsRepeated; 	}
//    //@}
//
//
//    //! @name Methods
//    //@{
//
//    virtual void updateRepeated( )=0;
//
//    //! Update the xb model elvecs
// //   virtual void updateElementSolution( UInt eleIDw )=0;
//
//    //! Solves the ionic model
//    virtual void solveXbModel( const vector_Type& Calcium,
//                         	 	 const Real timeStep )=0;


private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_alpha1;
    Real M_alpha2;
    Real M_alpha3;
    Real M_alpha4;
    Real M_alpha5;
    Real M_beta1;
    Real M_beta2;
    Real M_beta3;



    //! Xb states
//    vector_Type M_T;
//    vector_Type M_TCa;
//    vector_Type M_TCas;
//    vector_Type M_Ts;
//
//    vector_Type M_TRepeated;
//    vector_Type M_TCaRepeated;
//    vector_Type M_TCasRepeated;
//    vector_Type M_TsRepeated;

    //@}

}; // class XbNegroniLascano96

// ===================================================
//! Constructors
// ===================================================
XbNegroniLascano96::XbNegroniLascano96()
{
}

XbNegroniLascano96::XbNegroniLascano96( Epetra_Comm& comm )	:
    HeartXbModel( comm ),
    M_alpha1 ( 0.0 ),
    M_alpha2 ( 0.0 ),
    M_alpha3 ( 0.0 ),
    M_alpha4 ( 0.0 ),
    M_alpha5 ( 0.0 ),
    M_beta1  ( 0.0 ),
    M_beta2  ( 0.0 ),
    M_beta3  ( 0.0 )
//    M_T 		  ( HeartXbModel::M_localMap ),
//    M_TRepeated   ( M_T, Repeated ),
//    M_TCa 		  ( HeartXbModel::M_localMap ),
//    M_TCaRepeated ( M_Tca, Repeated ),
//    M_TCas 		  ( HeartXbModel::M_localMap ),
//    M_TCasRepeated( M_Tcas, Repeated ),
//    M_Ts 		  ( HeartXbModel::M_localMap ),
//    M_TsRepeated  ( M_Ts, Repeated )
{
}



XbNegroniLascano96::XbNegroniLascano96( Epetra_Comm& comm,
											 Teuchos::ParameterList& parameterList 	)	:
	HeartXbModel( comm )
////	M_T 		  ( HeartXbModel::M_localMap ),
////	M_TRepeated   ( M_T, Repeated ),
////	M_TCa 		  ( HeartXbModel::M_localMap ),
////	M_TCaRepeated ( M_Tca, Repeated ),
////	M_TCas 		  ( HeartXbModel::M_localMap ),
////	M_TCasRepeated( M_Tcas, Repeated ),
////	M_Ts 		  ( HeartXbModel::M_localMap ),
////	M_TsRepeated  ( M_Ts, Repeated )
{
	M_alpha1 =  parameterList.get("alpha1", 0.0);
	M_alpha2 =  parameterList.get("alpha2", 0.0);
	M_alpha3 =  parameterList.get("alpha3", 0.0);
	M_alpha4 =  parameterList.get("alpha4", 0.0);
	M_alpha5 =  parameterList.get("alpha5", 0.0);
	M_beta1  =  parameterList.get("beta1",  0.0);
	M_beta2  =  parameterList.get("beta1",  0.0);
	M_beta3  =  parameterList.get("beta1",  0.0);
}

//XbNegroniLascano96::XbNegroniLascano96( Epetra_Comm& comm,
//										const Teuchos::ParameterList& parameterList,
//										const vector_Type& TCa,
//										const vector_Type& TCas,
//										const vector_Type& Ts							)	:
//    HeartXbModel( comm )
//{
//	M_alpha1 =  parameterList.get("alpha1", 0.0);
//	M_alpha2 =  parameterList.get("alpha2", 0.0);
//	M_alpha3 =  parameterList.get("alpha3", 0.0);
//	M_alpha4 =  parameterList.get("alpha4", 0.0);
//	M_alpha5 =  parameterList.get("alpha5", 0.0);
//	M_beta1  =  parameterList.get("beta1",  0.0);
//	M_beta2  =  parameterList.get("beta1",  0.0);
//	M_beta3  =  parameterList.get("beta1",  0.0);
//
////	M_TCa 	 =	TCa;
////	M_TCas 	 =	TCas;
////	M_Ts 	 =	Ts;
////	M_T		 =	1.0 - TCa - TCas - Ts;
//}


XbNegroniLascano96::XbNegroniLascano96( const XbNegroniLascano96& Xb )
{
	M_comm 		= Xb.M_comm;
	M_me		= Xb.M_me;
	M_verbose	= Xb.M_verbose;

	M_alpha1 =  Xb.M_alpha1;
	M_alpha2 =  Xb.M_alpha2;
	M_alpha3 =  Xb.M_alpha3;
	M_alpha4 =  Xb.M_alpha4;
	M_alpha5 =  Xb.M_alpha5;
	M_beta1  =  Xb.M_beta1;
	M_beta2  =  Xb.M_beta2;
	M_beta3  =  Xb.M_beta3;

//	M_TCa 	 =	Xb.M_TCa;
//	M_TCas 	 =	Xb.M_TCas;
//	M_Ts 	 =	Xb.M_Ts;
//	M_T		 =	Xb.M_T;
}





// ===================================================
//! Operator
// ===================================================
XbNegroniLascano96& XbNegroniLascano96::operator=( const XbNegroniLascano96& Xb )
{
	M_comm 		= Xb.M_comm;
	M_me		= Xb.M_me;
	M_verbose	= Xb.M_verbose;

	M_alpha1 =  Xb.M_alpha1;
	M_alpha2 =  Xb.M_alpha2;
	M_alpha3 =  Xb.M_alpha3;
	M_alpha4 =  Xb.M_alpha4;
	M_alpha5 =  Xb.M_alpha5;
	M_beta1  =  Xb.M_beta1;
	M_beta2  =  Xb.M_beta2;
	M_beta3  =  Xb.M_beta3;

	return 		*this;
}


}

#endif
