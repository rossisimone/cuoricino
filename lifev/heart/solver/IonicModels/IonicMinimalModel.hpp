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


#ifndef _IONICMINIMALMODEL_H_
#define _IONICMINIMALMODEL_H_

#include <lifev/heart/solver/IonicModels/HeartIonicModel.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! XbModel - This class implements a mean field model.

class IonicMinimalModel : public virtual HeartIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef HeartIonicModel                         super;
    typedef boost::shared_ptr<VectorEpetra>         vectorPtr_Type;
    typedef boost::shared_ptr<VectorElemental>  elvecPtr_Type;
    typedef RegionMesh<LinearTetra>                 mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicMinimalModel();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicMinimalModel ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicMinimalModel object
     */
    IonicMinimalModel ( const IonicMinimalModel& model );
    //! Destructor
    virtual ~IonicMinimalModel() {}

    //@}

    //! @name Overloads
    //@{

    IonicMinimalModel& operator= ( const IonicMinimalModel& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& Uo()             const
    {
        return M_uo;
    }
    inline const Real& Uu()             const
    {
        return M_uu;
    }
    inline const Real& Tetav()      const
    {
        return M_tetav;
    }
    inline const Real& Tetaw()      const
    {
        return M_tetaw;
    }
    inline const Real& Tetavm()         const
    {
        return M_tetavm;
    }
    inline const Real& Tetao()      const
    {
        return M_tetao;
    }
    inline const Real& Tauv1()      const
    {
        return M_tauv1;
    }
    inline const Real& Tauv2()      const
    {
        return M_tauv2;
    }
    inline const Real& Tauvp()      const
    {
        return M_tauvp;
    }
    inline const Real& Tauw1()      const
    {
        return M_tauw1;
    }
    inline const Real& Tauw2()      const
    {
        return M_tauw2;
    }
    inline const Real& Kw()             const
    {
        return M_kw;
    }
    inline const Real& Uw()             const
    {
        return M_uw;
    }
    inline const Real& Tauwp()          const
    {
        return M_tauwp;
    }
    inline const Real& Taufi()      const
    {
        return M_taufi;
    }
    inline const Real& Tauo1()      const
    {
        return M_tauo1;
    }
    inline const Real& Tauo2()      const
    {
        return M_tauo2;
    }
    inline const Real& Tauso1()         const
    {
        return M_tauso1;
    }
    inline const Real& Tauso2()         const
    {
        return M_tauso2;
    }
    inline const Real& Kso()            const
    {
        return M_kso;
    }
    inline const Real& Uso()            const
    {
        return M_uso;
    }
    inline const Real& Taus1()      const
    {
        return M_taus1;
    }
    inline const Real& Taus2()      const
    {
        return M_taus2;
    }
    inline const Real& Ks()             const
    {
        return M_ks;
    }
    inline const Real& Us()         const
    {
        return M_us;
    }
    inline const Real& Tausi()      const
    {
        return M_tausi;
    }
    inline const Real& Tauinf()         const
    {
        return M_tauwinf;
    }
    inline const Real& Winfstart()  const
    {
        return M_winfstar;
    }


    inline void setUo           ( const Real& p )
    {
        this->M_uo        = p;
    }
    inline void setUu           ( const Real& p )
    {
        this->M_uu        = p;
    }
    inline void setTetav        ( const Real& p )
    {
        this->M_tetav     = p;
    }
    inline void setTetaw        ( const Real& p )
    {
        this->M_tetaw     = p;
    }
    inline void setTetavm       ( const Real& p )
    {
        this->M_tetavm    = p;
    }
    inline void setTetao        ( const Real& p )
    {
        this->M_tetao     = p;
    }
    inline void setTauv1        ( const Real& p )
    {
        this->M_tauv1     = p;
    }
    inline void setTauv2        ( const Real& p )
    {
        this->M_tauv2     = p;
    }
    inline void setTauvp        ( const Real& p )
    {
        this->M_tauvp     = p;
    }
    inline void setTauw1        ( const Real& p )
    {
        this->M_tauw1     = p;
    }
    inline void setTauw2        ( const Real& p )
    {
        this->M_tauw2     = p;
    }
    inline void setKw           ( const Real& p )
    {
        this->M_kw        = p;
    }
    inline void setUw           ( const Real& p )
    {
        this->M_uw        = p;
    }
    inline void setTauwp        ( const Real& p )
    {
        this->M_tauwp     = p;
    }
    inline void setTaufi        ( const Real& p )
    {
        this->M_taufi     = p;
    }
    inline void setTauo1        ( const Real& p )
    {
        this->M_tauo1     = p;
    }
    inline void setTauo2        ( const Real& p )
    {
        this->M_tauo2     = p;
    }
    inline void setTauso1       ( const Real& p )
    {
        this->M_tauso1    = p;
    }
    inline void setTauso2       ( const Real& p )
    {
        this->M_tauso2    = p;
    }
    inline void setKso      ( const Real& p )
    {
        this->M_kso       = p;
    }
    inline void setUso      ( const Real& p )
    {
        this->M_uso       = p;
    }
    inline void setTaus1        ( const Real& p )
    {
        this->M_taus1     = p;
    }
    inline void setTaus2        ( const Real& p )
    {
        this->M_taus2     = p;
    }
    inline void setKs           ( const Real& p )
    {
        this->M_ks        = p;
    }
    inline void setUs           ( const Real& p )
    {
        this->M_us        = p;
    }
    inline void setTausi        ( const Real& p )
    {
        this->M_tausi     = p;
    }
    inline void setTauinf       ( const Real& p )
    {
        this->M_tauwinf   = p;
    }
    inline void setWinfstart    ( const Real& p )
    {
        this->M_winfstar  = p;
    }



    //inline const short int& Size() const { return M_numberOfEquations; }
    //@}

    //! @name Methods
    //@{

    inline static Real Heaviside ( const Real& x)
    {
        if ( x > 0 )
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    //Compute the rhs on a single node or for the 0D case
    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeRhs ( const std::vector<Real>& v, const Real& Iapp, std::vector<Real>& rhs);

    //Compute the rhs on a mesh/ 3D case
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );
    //
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Iapp, std::vector<vectorPtr_Type>& rhs );

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v, const Real& Iapp);

    //    void computePotentialRhs(     const std::vector<vectorPtr_Type>& v,
    //                                  const VectorEpetra& Iapp,
    //                                  std::vector<vectorPtr_Type>& rhs,
    //                                  FESpace<mesh_Type, MapEpetra>& uFESpace );

    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_uo;
    Real M_uu;
    Real M_tetav;
    Real M_tetaw;
    Real M_tetavm;
    Real M_tetao;
    Real M_tauv1;
    Real M_tauv2;
    Real M_tauvp;
    Real M_tauw1;
    Real M_tauw2;
    Real M_kw;
    Real M_uw;
    Real M_tauwp;
    Real M_taufi;
    Real M_tauo1;
    Real M_tauo2;
    Real M_tauso1;
    Real M_tauso2;
    Real M_kso;
    Real M_uso;
    Real M_taus1;
    Real M_taus2;
    Real M_ks;
    Real M_us;
    Real M_tausi;
    Real M_tauwinf;
    Real M_winfstar;


    //! Xb states == equivalent to the number of equations
    //short int M_numberOfEquations;

    //@}

}; // class IonicMinimalModel

// ===================================================
//! Constructors
// ===================================================
//IonicMinimalModel::IonicMinimalModel()	:
//		super 	    ( 4     ),
//		M_uo    	( 0. 	),
//		M_uu    	( 1.61 	),
//		M_tetav 	( 0.3 	),
//		M_tetaw 	( 0.13 	),
//		M_tetavm	( 0.1 	),
//		M_tetao 	( 0.005 ),
//		M_tauv1 	( 80.0	),
//		M_tauv2 	( 1.4506),
//		M_tauvp 	( 1.4506),
//		M_tauw1 	( 70.0 	),
//		M_tauw2 	( 8.0 	),
//		M_kw    	( 200.0	),
//		M_uw    	( 0.016 ),
//		M_tauwp  	( 280.0	),
//		M_taufi 	( 0.078 ),
//		M_tauo1 	( 410.0	),
//		M_tauo2 	( 7.0 	),
//		M_tauso1	( 91.0 	),
//		M_tauso2	( 0.8 	),
//		M_kso   	( 2.1 	),
//		M_uso   	( 0.6 	),
//		M_taus1 	( 2.7342),
//		M_taus2 	( 4.0	),
//		M_ks    	( 2.0994),
//		M_us    	( 0.9087),
//		M_tausi 	( 3.3849),
//		M_tauwinf	( 0.01 	),
//		M_winfstar 	( 0.5 	)
//{
//}

IonicMinimalModel::IonicMinimalModel()	:
		super 	    ( 4     ),
		M_uo    	( 0. 	),
		M_uu    	( 1.58 	),
		M_tetav 	( 0.3 	),
		M_tetaw 	( 0.015	),
		M_tetavm	( 0.015	),
		M_tetao 	( 0.006 ),
		M_tauv1 	( 60.0	),
		M_tauv2 	( 1150.0),
		M_tauvp 	( 1.4506),
		M_tauw1 	( 70.0 	),
		M_tauw2 	( 20.0 	),
		M_kw    	( 65.0	),
		M_uw    	( 0.03  ),
		M_tauwp  	( 280.0	),
		M_taufi 	( 0.11  ),
		M_tauo1 	( 6.0	),
		M_tauo2 	( 6.0 	),
		M_tauso1	( 43.0 	),
		M_tauso2	( 0.2 	),
		M_kso   	( 2.0 	),
		M_uso   	( 0.65 	),
		M_taus1 	( 2.7342),
		M_taus2 	( 3.0	),
		M_ks    	( 2.0994),
		M_us    	( 0.9087),
		M_tausi 	( 2.8723),
		M_tauwinf	( 0.07 	),
		M_winfstar 	( 0.94 	)
{
}

IonicMinimalModel::IonicMinimalModel ( Teuchos::ParameterList& parameterList     )   :
    super       ( 4 )
{
    M_uo        =  parameterList.get ("uo",      0.0     );
    M_uu        =  parameterList.get ("uu",      1.61    );
    M_tetav     =  parameterList.get ("tetav",   0.30    );
    M_tetaw     =  parameterList.get ("tetaw",   0.130   );
    M_tetavm    =  parameterList.get ("tetavm",  0.10    );
    M_tetao     =  parameterList.get ("tetao",   0.005   );
    M_tauv1     =  parameterList.get ("tauv1",   80.0    );
    M_tauv2     =  parameterList.get ("tauv2",   1.4506  );
    M_tauvp     =  parameterList.get ("tauvp",   1.4506  );
    M_tauw1     =  parameterList.get ("tauw1",   70.0    );
    M_tauw2     =  parameterList.get ("tauw2",   8.0     );
    M_kw        =  parameterList.get ("kw",      200.0   );
    M_uw        =  parameterList.get ("uw",      0.016   );
    M_tauwp     =  parameterList.get ("tauwp",   280.0   );
    M_taufi     =  parameterList.get ("taufi",   0.078   );
    M_tauo1     =  parameterList.get ("tauo1",   410.0   );
    M_tauo2     =  parameterList.get ("tauo2",   7.0     );
    M_tauso1    =  parameterList.get ("tauso1",  91.0    );
    M_tauso2    =  parameterList.get ("tauso2",  0.8     );
    M_kso       =  parameterList.get ("kso",     2.1     );
    M_uso       =  parameterList.get ("uso",     0.6     );
    M_taus1     =  parameterList.get ("taus1",   2.7342  );
    M_taus2     =  parameterList.get ("taus2",   4.0     );
    M_ks        =  parameterList.get ("ks",      2.0994  );
    M_us        =  parameterList.get ("us",      0.9087  );
    M_tausi     =  parameterList.get ("tausi",   3.3849  );
    M_tauwinf   =  parameterList.get ("tauwinf", 0.01    );
    M_winfstar  =  parameterList.get ("winfstar", 0.5     );
}

IonicMinimalModel::IonicMinimalModel ( const IonicMinimalModel& model )
{

    M_uo        =  model.M_uo;
    M_uu        =  model.M_uu;
    M_tetav     =  model.M_tetav;
    M_tetaw     =  model.M_tetaw;
    M_tetavm    =  model.M_tetavm;
    M_tetao     =  model.M_tetao;
    M_tauv1     =  model.M_tauv1;
    M_tauv2     =  model.M_tauv2;
    M_tauvp     =  model.M_tauvp;
    M_tauw1     =  model.M_tauw1;
    M_tauw2     =  model.M_tauw2;
    M_kw        =  model.M_kw;
    M_uw        =  model.M_uw;
    M_tauwp     =  model.M_tauwp;
    M_taufi     =  model.M_taufi;
    M_tauo1     =  model.M_tauo1;
    M_tauo2     =  model.M_tauo2;
    M_tauso1    =  model.M_tauso1;
    M_tauso2    =  model.M_tauso2;
    M_kso       =  model.M_kso;
    M_uso       =  model.M_uso;
    M_taus1     =  model.M_taus1;
    M_taus2     =  model.M_taus2;
    M_ks        =  model.M_ks;
    M_us        =  model.M_us;
    M_tausi     =  model.M_tausi;
    M_tauwinf   =  model.M_tauwinf;
    M_winfstar  =  model.M_winfstar;

    M_numberOfEquations = model.M_numberOfEquations;
}

// ===================================================
//! Operator
// ===================================================
IonicMinimalModel& IonicMinimalModel::operator= ( const IonicMinimalModel& model )
{
    M_uo        =  model.M_uo;
    M_uu        =  model.M_uu;
    M_tetav     =  model.M_tetav;
    M_tetaw     =  model.M_tetaw;
    M_tetavm    =  model.M_tetavm;
    M_tetao     =  model.M_tetao;
    M_tauv1     =  model.M_tauv1;
    M_tauv2     =  model.M_tauv2;
    M_tauvp     =  model.M_tauvp;
    M_tauw1     =  model.M_tauw1;
    M_tauw2     =  model.M_tauw2;
    M_kw        =  model.M_kw;
    M_uw        =  model.M_uw;
    M_tauwp     =  model.M_tauwp;
    M_taufi     =  model.M_taufi;
    M_tauo1     =  model.M_tauo1;
    M_tauo2     =  model.M_tauo2;
    M_tauso1    =  model.M_tauso1;
    M_tauso2    =  model.M_tauso2;
    M_kso       =  model.M_kso;
    M_uso       =  model.M_uso;
    M_taus1     =  model.M_taus1;
    M_taus2     =  model.M_taus2;
    M_ks        =  model.M_ks;
    M_us        =  model.M_us;
    M_tausi     =  model.M_tausi;
    M_tauwinf   =  model.M_tauwinf;
    M_winfstar  =  model.M_winfstar;

    M_numberOfEquations = model.M_numberOfEquations;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
void IonicMinimalModel::computeRhs ( const   std::vector<Real>&  v,
                                     std::vector<Real>& rhs )
{

    Real U = v[0];
    Real V = v[1];
    Real W = v[2];
    Real S = v[3];

    Real tauvm = ( 1.0 - Heaviside ( U - M_tetavm ) ) * M_tauv1 + Heaviside ( U - M_tetavm ) * M_tauv2;
    Real tauwm = M_tauw1 + ( M_tauw2  - M_tauw1  ) * ( 1.0 + std::tanh ( M_kw  * ( U - M_uw  ) ) ) / 2.0;
    Real taus  = ( 1.0 - Heaviside ( U - M_tetaw ) ) * M_taus1 + Heaviside ( U - M_tetaw ) * M_taus2;

    Real vinf  = Heaviside ( M_tetavm - U );
    Real winf  = ( 1.0 - Heaviside ( U - M_tetao ) ) * ( 1.0 - U / M_tauwinf ) + Heaviside ( U - M_tetao ) * M_winfstar;

    rhs[0] = ( 1.0 - Heaviside ( U - M_tetav ) ) * ( vinf - V ) / tauvm - Heaviside ( U - M_tetav ) * V / M_tauvp;
    rhs[1] = ( 1.0 - Heaviside ( U - M_tetaw ) ) * ( winf - W ) / tauwm - Heaviside ( U - M_tetaw ) * W / M_tauwp;
    rhs[2] = ( ( 1.0 + std::tanh ( M_ks * ( U - M_us ) ) ) / 2.0 - S ) / taus;

}

void IonicMinimalModel::computeRhs ( const   std::vector<Real>&  v,
                                     const   Real&           Iapp,
                                     std::vector<Real>& rhs )
{

    Real U = v[0];
    Real V = v[1];
    Real W = v[2];
    Real S = v[3];

    Real tauvm = ( 1.0 - Heaviside ( U - M_tetavm ) ) * M_tauv1 + Heaviside ( U - M_tetavm ) * M_tauv2;
    Real tauwm = M_tauw1 + ( M_tauw2  - M_tauw1  ) * ( 1.0 + std::tanh ( M_kw  * ( U - M_uw  ) ) ) / 2.0;
    Real tauso = M_tauso1 + ( M_tauso2 - M_tauso1 ) * ( 1.0 + std::tanh ( M_kso * ( U - M_uso ) ) ) / 2.0;
    Real taus  = ( 1.0 - Heaviside ( U - M_tetaw ) ) * M_taus1 + Heaviside ( U - M_tetaw ) * M_taus2;
    Real tauo  = ( 1.0 - Heaviside ( U - M_tetao ) ) * M_tauo1 + Heaviside ( U - M_tetao ) * M_tauo2;

    Real vinf  = Heaviside ( M_tetavm - U );
    Real winf  = ( 1.0 - Heaviside ( U - M_tetao ) ) * ( 1.0 - U / M_tauwinf ) + Heaviside ( U - M_tetao ) * M_winfstar;

    Real Jfi   = - V * Heaviside ( U - M_tetav ) * ( U - M_tetav ) * ( M_uu - U ) / M_taufi;
    Real Jso   = ( U - M_uo ) * ( 1.0 - Heaviside ( U - M_tetaw )  ) / tauo + Heaviside ( U - M_tetaw ) / tauso;
    Real Jsi   = - Heaviside ( U - M_tetaw ) * W * W / M_tausi;

	rhs[0] = - ( Jfi + Jso + Jsi ) + Iapp;
	rhs[1] = ( 1.0 - Heaviside( U - M_tetav ) ) * ( vinf - V ) / tauvm - Heaviside( U - M_tetav ) * V / M_tauvp;
	rhs[2] = ( 1.0 - Heaviside( U - M_tetaw ) ) * ( winf - W ) / tauwm - Heaviside( U - M_tetaw ) * W / M_tauwp;
	rhs[3] = ( ( 1.0 + std::tanh( M_ks * ( U - M_us ) ) ) / 2.0 - S ) / taus;
}


Real IonicMinimalModel::computeLocalPotentialRhs( const std::vector<Real>& v, const Real& Iapp)
{
    Real dPotential (0.0);

	Real U = v[0];
	Real V = v[1];
	Real W = v[2];
	//Real S = v[3];

    Real tauso = M_tauso1 + ( M_tauso2 - M_tauso1 ) * ( 1.0 + std::tanh ( M_kso * ( U - M_uso ) ) ) / 2.0;
    Real tauo  = ( 1.0 - Heaviside ( U - M_tetao ) ) * M_tauo1 + Heaviside ( U - M_tetao ) * M_tauo2;

    Real Jfi   = - V * Heaviside ( U - M_tetav ) * ( U - M_tetav ) * ( M_uu - U ) / M_taufi;
    Real Jso   = ( U - M_uo ) * ( 1.0 - Heaviside ( U - M_tetaw )  ) / tauo + Heaviside ( U - M_tetaw ) / tauso;
    Real Jsi   = - Heaviside ( U - M_tetaw ) * W * W / M_tausi;

    dPotential = - ( Jfi + Jso + Jsi ) + Iapp;

    return dPotential;
}




void IonicMinimalModel::showMe()
{

    std::cout << "\n\tHi, I'm the minimal model\n\t See you soon\n\n";
}


}

#endif
