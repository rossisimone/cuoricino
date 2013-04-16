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
 *  @brief File containing the Multiscale Model FSI3D Activated with electrophysiology
 *
 *  @date 15-04-2013
 *  @author Simone Rossi <simone.rossi@epfl.ch>
 *
 *  @maintainer Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef MultiscaleModelFSI3DActivated_H
#define MultiscaleModelFSI3DActivated_H 1

#include <lifev/multiscale/models/MultiscaleModelFSI3D.hpp>
#include <lifev/heart/solver/HeartETAMonodomainSolver.hpp>



namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelFSI3DActivated - Multiscale model for 3D FSI simulations with activation coming from electrophysiology
/*!
 *  @author Simone Rossi

 */
class MultiscaleModelFSI3DActivated: public virtual MultiscaleModelFSI3D
{
public:

    //! @name Public Types
    //@{
	typedef MultiscaleModel									base;
    typedef MultiscaleModelFSI3D                             super;

    typedef IonicMinimalModel							  minimalModel_Type;
    typedef boost::shared_ptr< minimalModel_Type >        minimalModelPtr_Type;
    typedef HeartETAMonodomainSolver< mesh_Type, minimalModel_Type >        monoSolver_Type;
//    typedef Heart
    typedef boost::shared_ptr< monoSolver_Type >         monoSolverPtr_Type;

    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const LifeV::ID&   /*i*/ ) >   function_Type;

    typedef VectorEpetra				vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef Exporter< mesh_Type >                              IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >                   IOFilePtr_Type;
    typedef ExporterData< mesh_Type >                          IOData_Type;

    typedef ExporterEnsight< mesh_Type >                       ensightIOFile_Type;
#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type >                          hdf5IOFile_Type;
#endif

   //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelFSI3DActivated();

    //! Destructor
    virtual ~MultiscaleModelFSI3DActivated() {}

    //@}


    //! @name MultiscaleModel Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData ( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Update the solution.
    void updateSolution();

    //! Save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //! Return a specific scalar quantity to be used for a comparison with a reference value.
    /*!
     * This method is meant to be used for night checks.
     * @return reference quantity.
     */
    Real checkSolution() const;

    //@}


    //! @name MultiscaleInterface Methods
    //@{


    //@}


    //! @name Get Methods
    //@{

    //! Get the FSI3D operator
    /*!
     * @return FSI3D operator
     */
//    const FSIOperatorPtr_Type& solver() const
//    {
//        return M_FSIoperator;
//    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelFSI3DActivated ( const MultiscaleModelFSI3DActivated& model );

    MultiscaleModelFSI3DActivated& operator= ( const MultiscaleModelFSI3DActivated& model );

    //@}


    //! @name Private Methods
    //@{

    //! Setup the global data of the model.
    /*!
     * In particular, it replaces the default local values with the ones in the global container.
     * If a value is already specified in the data file, do not perform the replacement.
     *
     * @param fileName File name of the specific model.
     */
//    void setupGlobalData ( const std::string& fileName );
//
//    //! Initialize the FSI solution.
//    void initializeSolution();
//
//    void setupCommunicator();
//
//    void setupBC ( const std::string& fileName );
//    void updateBC();
//
//    void setupExporter ( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label = "" );
//    void setupImporter ( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label = "" );
//
//    void setExporterFluid ( const IOFilePtr_Type& exporter );
//    void setExporterSolid ( const IOFilePtr_Type& exporter );
//
//    void exportFluidSolution();
//    void exportSolidSolution();
//
//    //! Setup the linear model
//    void setupLinearModel();
//
//    //! Update the linear system matrix and vectors
//    void updateLinearModel();
//
//    //! Solve the linear problem
//    void solveLinearModel ( bool& solveLinearSystem );
//
//    //! Impose the coupling perturbation on the correct BC inside the BCHandler
//    void imposePerturbation();
//
//    //! Reset all the coupling perturbations imposed on the BCHandler
//    void resetPerturbation();

    //@}


    vectorPtr_Type							M_fiber;
    vectorPtr_Type							M_gammaf;
    // Operator
    monoSolverPtr_Type                M_monodomain;

    // Exporters
    IOFilePtr_Type                          M_exporterElectro;

    // Importers
    IOFilePtr_Type                          M_importerElectro;


};



//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelFSI3DActivated()
{
    return new MultiscaleModelFSI3DActivated();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelFSI3D_H */
