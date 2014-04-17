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
 *  @file
 *  @brief This file contains solvers for different active materials.
 *  @warning: This is the most important issue related with this class.
 *  At the moment, BC they are applied on the residual directly. This does
 *  not work for nonhomogeneus Dirichlet conditions!!
 *  Check the structural operator class in the structure module
 *
 *  @version 1.0
 *  @date 04-2014
 *  @author Paolo Tricerri
 *
 *  @maintainer  Simone Rossi <simone.rossi@epfl.ch>
*/

#ifndef _EMSTRUCTURALOPERATOR_H_
#define _EMSTRUCTURALOPERATOR_H_ 1

#include <lifev/em/solver/EMActiveStructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>


namespace LifeV
{

using namespace ExpressionAssembly;

/*!
  \class EMStructuralSolver
  \brief
  This class solves the elastodynamics equations for different active materials


*/
template <typename Mesh>
class EMStructuralOperator : public StructuralOperator<Mesh>
{
public:

    //!@name Type definitions
    //@{
    typedef EMActiveStructuralConstitutiveLaw<Mesh>     material_Type;

    typedef boost::shared_ptr<material_Type>			materialPtr_Type;

    typedef StructuralOperator<Mesh>					super;

    typedef typename super::data_Type					data_Type;

    typedef typename super::FESpacePtr_Type				FESpacePtr_Type;

    typedef typename super::ETFESpacePtr_Type			ETFESpacePtr_Type;

    typedef typename super::bcHandler_Type				bcHandler_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    EMStructuralOperator();

    virtual ~EMStructuralOperator() {};

    //@}

    //!@name Methods
    //@{

    //! Setup the created object of the class Venantkirchhof
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param BCh boundary conditions for the displacement
      \param comm the Epetra Comunicator
    */
    void setup ( boost::shared_ptr<data_Type>  data,
                 const FESpacePtr_Type&        dFESpace,
                 const ETFESpacePtr_Type&      dETFESpace,
                 const bcHandler_Type&       BCh,
                 boost::shared_ptr<Epetra_Comm>&     comm
               );

    /*! Get the offset parameter. It is taken into account when the boundary conditions
      are applied and the matrices are assembled.
    */
    const materialPtr_Type& activeMaterial() const
    {
        return M_activeMaterial;
    }

    //@}

protected:

    //! Material class
    materialPtr_Type                     M_activeMaterial;

};

//====================================
// Constructor
//=====================================

template <typename Mesh>
EMStructuralOperator<Mesh>::EMStructuralOperator( ) :
    super(),
    M_activeMaterial()
{
}

template <typename Mesh>
void
EMStructuralOperator<Mesh>::setup (boost::shared_ptr<data_Type>          data,
                                 const FESpacePtr_Type& dFESpace,
                                 const ETFESpacePtr_Type& dETFESpace,
                                 const bcHandler_Type&                    BCh,
                                 boost::shared_ptr<Epetra_Comm>&   comm)
{
    super::setup (data, dFESpace, dETFESpace, BCh, comm);
//    material_Type * tmpMaterial = dynamic_cast<material_Type *>(this->M_material.get());
    //    M_activeMaterial = boost::dynamic_pointer_cast<materialPtr_Type>(const_pointer_cast<typename super::materialPtr_Type>(super::M_material) );
   M_activeMaterial.reset(dynamic_cast<material_Type *>(this->M_material.get()));// boost::dynamic_cast< material_Type *>( super::M_material.get() ) );// ( boost::const_pointer_cast<typename super::materialPtr_Type>(super::M_material) );
}




}
#endif
