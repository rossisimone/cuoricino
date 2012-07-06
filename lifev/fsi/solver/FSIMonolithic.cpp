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

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <EpetraExt_Reindex_MultiVector.h>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/fsi/solver/FSIMonolithic.hpp>

namespace LifeV
{

// ===================================================
// Constructors, Destructor
// ===================================================
FSIMonolithic::FSIMonolithic():
        super_Type(),
        M_monolithicMap(),
        M_interfaceMap(),
        M_beta(),
        M_monolithicMatrix(),
        M_precPtr(),
        M_rhsFull(),
        M_BCh_flux(),
        M_BChWS(),
        M_offset(0),
        M_solidAndFluidDim(0),
        M_fluidBlock(),
        M_solidBlockPrec(),
        M_linearSolver(),
        M_numerationInterface(),
        M_BChs(),
        M_FESpaces(),
        M_diagonalScale(false),
        M_reusePrec(true),
        M_resetPrec(true),
        M_maxIterSolver(-1),
        M_restarts(false),
        //end of protected attributes
        M_preconditionedSymmetrizedMatrix(),
        M_stress(),
        M_fluxes(0)
#ifdef OBSOLETE
        ,M_rhsShapeDerivatives()
#endif
{}

FSIMonolithic::~FSIMonolithic()
{
}

// ===================================================
// Setup Methods
// ===================================================
void
FSIMonolithic::setupFEspace()
{
    super_Type::setupFEspace();

    // Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset( new FESpace<mesh_Type, MapEpetra>( M_solidMesh,
                                                         M_data->dataSolid()->order(),
                                                         nDimensions,
                                                         M_epetraComm));

}

void
FSIMonolithic::setupDOF( void )
{

    M_dofStructureToHarmonicExtension    .reset( new DOFInterface3Dto3D );
    M_dofStructureToFluid    .reset( new DOFInterface3Dto3D );

    M_dofStructureToHarmonicExtension->setup(   M_mmFESpace->refFE(), M_mmFESpace->dof(),
                                                M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToHarmonicExtension->update( *M_mmFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                               *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                               M_data->interfaceTolerance(),
                                               M_data->fluidInterfaceVertexFlag() );

    M_dofStructureToFluid->setup(   M_uFESpace->refFE(), M_uFESpace->dof(),
                                    M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToFluid->update( *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                   *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                   M_data->interfaceTolerance(),
                                   M_data->fluidInterfaceVertexFlag() );


    createInterfaceMaps(M_dofStructureToFluid/*HarmonicExtension*/->localDofMap());

    M_fluidMeshPart->releaseUnpartitionedMesh();
    M_solidMeshPart->releaseUnpartitionedMesh();
    M_fluidMesh.reset();
    M_solidMesh.reset();
}

#ifdef HAVE_HDF5
void
FSIMonolithic::setupDOF( meshFilter_Type& filterMesh )
{
    createInterfaceMaps(*filterMesh.getStoredInterface(0));
}
#endif

void
FSIMonolithic::setupSystem( )
{
    M_fluid->setUp( M_dataFile );
    setUp( M_dataFile );
}

void
FSIMonolithic::setUp( const GetPot& dataFile )
{

    M_linearSolver.reset(new solver_Type(M_epetraComm));

    M_linearSolver->setDataFromGetPot( dataFile, "linear_system/solver" );
    std::string prectype = dataFile("problem/DDBlockPrec", "PRECTYPE_UNDEFINED");
    std::string opertype = dataFile("problem/blockOper", "PRECTYPE_UNDEFINED");

    M_precPtr.reset(BlockPrecFactory::instance().createObject( prectype ));

    M_precPtr->setupSolver(*M_linearSolver, dataFile);
    std::string section("linear_system/prec");
    M_precPtr->setComm(M_epetraComm);
    M_precPtr->setDataFromGetPot(dataFile, section);//to avoid if we build the prec from a matrix.
    M_monolithicMatrix->setComm(M_epetraComm);
    M_monolithicMatrix->setDataFromGetPot(dataFile, section);//to avoid if we build the prec from a matrix.
    M_reusePrec     = dataFile( "linear_system/prec/reuse", true);
    M_maxIterSolver = dataFile( "linear_system/solver/max_iter", -1);
    M_diagonalScale    = dataFile( "linear_system/prec/diagonalScaling",  false );
    M_restarts         = dataFile( "exporter/start"  ,  0   );
}

void
FSIMonolithic::setupFluidSolid( )
{
    M_BCh_flux = M_BCh_u; // For the moment M_BCh_u contains only the fluxes.
    M_fluxes=M_BCh_u->size( );

    setupFluidSolid( M_fluxes );

    M_BCh_flux->setOffset(M_offset-M_fluxes);
    std::vector<BCBase>::iterator fluxIt = M_BCh_flux->begin( );
    for ( UInt i = 0; i < M_fluxes; ++i, ++fluxIt )
        fluxIt->setOffset( i );

}

void
FSIMonolithic::setupFluidSolid( UInt const fluxes )
{

    // Note: up to now it works only with matching grids (and poly order) on the interface
    assert(M_fluidInterfaceMap->map(Unique)->NumGlobalElements() == M_solidInterfaceMap->map(Unique)->NumGlobalElements());

    M_interfaceMap = M_solidInterfaceMap;

    //std::map<ID, ID> const& localDofMap = M_dofStructureToHarmonicExtension->localDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    M_monolithicMap.reset(new MapEpetra(M_uFESpace->map()));

    std::string opertype = M_dataFile("problem/blockOper", "AdditiveSchwarz");

    createOperator( opertype );

    *M_monolithicMap+= M_pFESpace->map();
    *M_monolithicMap+= fluxes;
    *M_monolithicMap+= M_dFESpace->map();

    M_monolithicMatrix->createInterfaceMap( *M_interfaceMap, M_dofStructureToFluid->localDofMap(), M_dFESpace->map().map(Unique)->NumGlobalElements()/nDimensions, M_epetraWorldComm );
    *M_monolithicMap += *M_monolithicMatrix->interfaceMap();

    //the map for the interface coupling matrices should be done with respect to the coarser mesh.
    M_monolithicMatrix->numerationInterface(M_numerationInterface);
    M_beta.reset  (new vector_Type(/*M_monolithicMap*/M_uFESpace->map()));

    M_offset = M_uFESpace->dof().numTotalDof()*nDimensions + fluxes +  M_pFESpace->dof().numTotalDof();
    M_solidAndFluidDim= M_offset + M_dFESpace->dof().numTotalDof()*nDimensions;
    M_BCh_d->setOffset(M_offset);
}

// ===================================================
// Public Methods
// ===================================================
void
FSIMonolithic::monolithicToInterface(vector_Type& lambdaSolid, const vector_Type& disp)
{
    if (disp.mapType() == Repeated)
    {
        vector_Type const  dispUnique(disp, Unique);
        monolithicToInterface(lambdaSolid, dispUnique);
        return;
    }
    if (lambdaSolid.mapType() == Repeated)
    {
        vector_Type  lambdaSolidUn(lambdaSolid.map(), Unique);
        monolithicToInterface( lambdaSolidUn, disp);
        lambdaSolid = lambdaSolidUn;
        return;
    }
    /* UInt MyOffset(M_uFESpace->map().getMap(Unique)->NumMyElements()+M_pFESpace->map().getMap(Unique)->NumMyElements());
       vector_Type subDisp(this->M_dFESpace->map(), Unique);
       subDisp.mySubset(disp, MyOffset);
       lambdaSolid=subDisp;*/

    MapEpetra subMap(*disp.map().map(Unique), M_offset,disp.map().map(Unique)->NumGlobalElements() );
    vector_Type subDisp(subMap, Unique);
    subDisp.subset(disp, M_offset);
    lambdaSolid=subDisp;
}

void
FSIMonolithic::monolithicToX(const vector_Type& disp, vector_Type& dispFluid, MapEpetra& map, UInt offset)
{
    if(disp.mapType()== Repeated)
    {
        vector_Type dispUnique(disp, Unique);
        monolithicToX(dispUnique, dispFluid, map, offset);
        dispFluid = dispUnique;
        return;
    }
    dispFluid.subset(disp, map, offset, offset);
}


void
FSIMonolithic::buildSystem()
{
    M_solid->buildSystem( M_solidTimeAdvance->coefficientSecondDerivative( 0 )/(M_data->dataSolid()->dataTime()->timeStep()*M_data->dataSolid()->dataTime()->timeStep()));
}

#ifdef HAVE_TRILINOS_ANASAZI
Real&
FSIMonolithic::computeMaxSingularValue( )
{
    typedef Epetra_Operator                                                operatorPtr_Type;

    M_preconditionedSymmetrizedMatrix.reset(new ComposedOperator<Epetra_Operator>(M_epetraComm));

    boost::shared_ptr<operatorPtr_Type>  ComposedPrecPtr(M_linearSolver->preconditioner()->preconditioner());

    //M_monolithicMatrix->getMatrixPtr()->OptimizeStorage();
    boost::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr(new Epetra_FECrsMatrix(*M_monolithicMatrix->matrix()->matrixPtr()));
    M_preconditionedSymmetrizedMatrix->push_back(boost::static_pointer_cast<operatorPtr_Type>(/*ComposedPrecPtr*/matrCrsPtr));
    M_preconditionedSymmetrizedMatrix->push_back((ComposedPrecPtr/*matrCrsPtr*/), true);
    M_preconditionedSymmetrizedMatrix->push_back((ComposedPrecPtr/*matrCrsPtr*/), true, true);
    M_preconditionedSymmetrizedMatrix->push_back(boost::static_pointer_cast<operatorPtr_Type>(/*ComposedPrecPtr*/matrCrsPtr), false, true);

    std::vector<LifeV::Real> real;
    std::vector<LifeV::Real> imaginary;

    boost::shared_ptr<EigenSolver> eig;

    UInt nev = M_dataFile("eigensolver/nevec", 10);//number of eigenvectors
    if(nev)
    {
        eig.reset(new EigenSolver(M_preconditionedSymmetrizedMatrix, M_preconditionedSymmetrizedMatrix->OperatorDomainMap(), nev));
        eig->setDataFromGetPot(M_dataFile, "eigensolver/");
        eig->solve();
        eig->eigenvalues(real, imaginary);
    }
    else
    {
        throw UNDEF_EIGENSOLVER_EXCEPTION();
    }
    for (UInt i=0; i<real.size(); ++i)
    {
        displayer().leaderPrint("\n real part ", real[i]);
        displayer().leaderPrint("\n imaginary part ", imaginary[i]);
    }
    return real[0];
}
#endif

void
FSIMonolithic::computeFluidNormals( vector_Type& normals)
{
    BCManageNormal<matrix_Type> normalManager;
    if ( !M_BChWS->bcUpdateDone() )//possibly to avoid
        M_BChWS->bcUpdate(*M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    normalManager.init((*M_BChWS)[0], 0.);
    normalManager.computeIntegratedNormals(M_uFESpace->dof(), M_uFESpace->feBd(), normals, *M_uFESpace->mesh());
}

void
FSIMonolithic::solveJac( vector_Type& step, const vector_Type& res, const Real /*linearRelTol*/ )
{
    setupBlockPrec( );

    checkIfChangedFluxBC( M_precPtr );

    M_precPtr->blockAssembling();
    M_precPtr->applyBoundaryConditions( dataFluid()->dataTime()->time() );
    M_precPtr->GlobalAssemble();

#ifdef HAVE_LIFEV_DEBUG
    M_solid->displayer().leaderPrint("  M-  Residual NormInf:                        ", res.normInf(), "\n");
#endif
    iterateMonolithic(res, step);
#ifdef HAVE_LIFEV_DEBUG
    M_solid->displayer().leaderPrint("  M-  Solution NormInf:                        ", step.normInf(), "\n");
#endif
}

void
FSIMonolithic::updateSystem()
{
    //M_solidBlock->spy("solid");

    // this->fluid().updateUn(*this->M_un);
    *M_rhs*=0;
    *M_rhsFull*=0;
    this->M_fluid->resetStabilization();
}


// ===================================================
// Protected Methods
// ===================================================
void
FSIMonolithic::iterateMonolithic(const vector_Type& rhs, vector_Type& step)
{
    LifeChrono chrono;

    displayer().leaderPrint("  M-  Solving the system ... \n" );

    //M_monolithicMatrix->GlobalAssemble();
    //necessary if we did not imposed Dirichlet b.c.

    M_linearSolver->setOperator(*M_monolithicMatrix->matrix()->matrixPtr());

    M_linearSolver->setReusePreconditioner( (M_reusePrec) && (!M_resetPrec) );

    int numIter = M_precPtr->solveSystem( rhs, step, M_linearSolver );

    if (numIter < 0)
    {
        chrono.start();

        M_solid->displayer().leaderPrint("   Iterative solver failed, numiter = ", -numIter );

        if (numIter <= -M_maxIterSolver)
            M_solid->displayer().leaderPrint("   ERROR: Iterative solver failed.\n");
    }

    //M_solid->getDisplayer().leaderPrint("  M-  System solved.\n" );
}

void
FSIMonolithic::couplingRhs(vectorPtr_Type rhs) // not working with non-matching grids
{
    std::map<ID, ID> const& localDofMap = M_dofStructureToFluid->localDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //    UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;

    vector_Type extrapolation(M_solidTimeAdvance->rhsContributionFirstDerivative()*M_solid->rescaleFactor(), Unique);
    vector_Type lambda(*M_interfaceMap, Unique);

    this->monolithicToInterface(lambda, extrapolation);

    UInt interface(M_monolithicMatrix->interface());
    //Real rescale(M_solid->rescaleFactor());
    UInt totalDofs(M_dFESpace->dof().numTotalDof());


    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = localDofMap.begin(); ITrow != localDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap->map(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                if(rhs.get())
                    (*rhs)[  (int)(*M_numerationInterface)[ITrow->second ] + dim*interface +M_solidAndFluidDim ] = -lambda( ITrow->second + dim*totalDofs );
            }
        }
    }
}

void
FSIMonolithic::
evalResidual( const vector_Type& sol, vectorPtr_Type& rhs, vector_Type& res, bool diagonalScaling)
{
    if( diagonalScaling )
        diagonalScale(*rhs, M_monolithicMatrix->matrix());
if(!(M_data->dataSolid()->solidType().compare("exponential") && M_data->dataSolid()->solidType().compare("neoHookean")) )
{
    M_solid->Apply(sol*M_solid->rescaleFactor(), res);
    M_fluidBlock->globalAssemble();

    res += ((*M_fluidBlock)*sol);

    res += *M_monolithicMatrix->coupling()*sol;
}
else
{
    res = *(M_monolithicMatrix->matrix())*sol;
    res -= *rhs; // Ax-b
}
}

void
FSIMonolithic::
updateSolidSystem( vectorPtr_Type & rhsFluidCoupling )
{
    Real coefficient ( M_data->dataSolid()->dataTime()->timeStep() * M_data->dataSolid()->dataTime()->timeStep() * M_solid->rescaleFactor() /  M_solidTimeAdvance->coefficientSecondDerivative( 0 ) );
    M_solidTimeAdvance->updateRHSContribution( M_data->dataSolid()->dataTime()->timeStep() );
    *rhsFluidCoupling += *M_solid->Mass() * ( M_solidTimeAdvance->rhsContributionSecondDerivative() * coefficient );
}

void
FSIMonolithic::
diagonalScale(vector_Type& rhs, matrixPtr_Type matrFull)
{
    Epetra_Vector diagonal(*rhs.map().map(Unique));
    //M_matrFull->matrixPtr()->InvRowSums(diagonal);
    //M_matrFull->matrixPtr()->InvRowMaxs(diagonal);
    //M_matrFull->matrixPtr()->InvColSums(diagonal);
    matrFull->matrixPtr()->InvColMaxs(diagonal);
    matrFull->matrixPtr()->LeftScale(diagonal);
    rhs.epetraVector().Multiply(1, rhs.epetraVector(), diagonal,0);
}

void
FSIMonolithic::variablesInit(const std::string& dOrder)
{
    M_dFESpace.reset(new FESpace<mesh_Type, MapEpetra>(*M_solidMeshPart,
                                                       dOrder,
                                                       3,
                                                       M_epetraComm));
    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset(new vector_Type(*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset(new vector_Type(*M_fluidInterfaceMap, Repeated) );
}

void FSIMonolithic::setupBlockPrec( )
{
#ifdef HAVE_NS_PREC
  std::string PCD("PCD");
  UInt fluidPosition = M_precPtr->whereIsBlock(MonolithicBlockComposed::fluid);
  std::string precType(M_precPtr->blockPrecs()[fluidPosition]->preconditionerType());
#endif
    if(!(M_precPtr->set()))
     {
         M_precPtr->push_back_matrix(M_solidBlockPrec, M_structureNonLinear);
         M_precPtr->push_back_matrix(M_fluidBlock, true);
         M_precPtr->setConditions(M_BChs);
         M_precPtr->setSpaces(M_FESpaces);
         M_precPtr->setOffsets(2, M_offset, 0);
         M_precPtr->coupler(M_monolithicMap, M_dofStructureToFluid->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), M_solidTimeAdvance->coefficientFirstDerivative( 0 ), M_solid->rescaleFactor());
     }
    else
    {
        M_precPtr->replace_matrix(M_fluidBlock, 1);
        M_precPtr->replace_matrix(M_solidBlockPrec, 0);
    }

#ifdef HAVE_NS_PREC
    if(M_precPtr->blockPrecs().size()>1)
      {

	if(!precType.compare(PCD))
	  {
	    Preconditioner* prec=(M_precPtr->blockPrecs()[fluidPosition].get());
	    PreconditionerPCD* prec_PCD = dynamic_cast<PreconditionerPCD*>(prec);
	    ASSERT(prec, "The preconditioner corresponding to the fluid block is probably not PCD. Check in the data file");
	    prec_PCD->setFESpace(M_uFESpace, M_pFESpace);
	    prec_PCD->setBCHandler(M_BCh_u);
	    prec_PCD->setTimestep(M_data->dataFluid()->dataTime()->timeStep());
	    prec_PCD->setViscosity(M_data->dataFluid()->viscosity());
	    prec_PCD->setDensity(M_data->dataFluid()->density());
	    prec_PCD->setCouplingMatrixView(M_precPtr->couplingVector()[MonolithicBlockComposed::fluid]);
	    prec_PCD->setMapStructure(&M_dFESpace->map());
	    prec_PCD->updateBeta(M_fluidTimeAdvance->singleElement(0));
	  }
      }
#endif
}

void
FSIMonolithic::assembleSolidBlock( UInt iter, const vector_Type& solution )
{
    if (iter == 0)
        updateSolidSystem(this->M_rhs);


if(M_data->dataSolid()->solidType().compare("exponential") && M_data->dataSolid()->solidType().compare("neoHookean"))
{
    M_solid->material()->computeStiffness(solution*M_solid->rescaleFactor()/**M_data->dataFluid()->dataTime()->timeStep()*/, 1./*M_solid->rescaleFactor()*/, M_data->dataSolid(), M_solid->displayerPtr());
    M_solidBlockPrec.reset(new matrix_Type(*M_monolithicMap, 1));
    *M_solidBlockPrec += *M_solid->Mass();
    *M_solidBlockPrec += *M_solid->material()->stiffMatrix();
    M_solidBlockPrec->globalAssemble();
    *M_solidBlockPrec *= M_solid->rescaleFactor();
}
else
{
    M_solid->material()->updateJacobianMatrix( solution*M_solid->rescaleFactor()/**M_data->dataFluid()->dataTime()->timeStep()*/, dataSolid(), M_solid->displayerPtr() ); // computing the derivatives if nonlinear (comment this for inexact Newton);
    M_solidBlockPrec.reset(new matrix_Type(*M_monolithicMap, 1));
    *M_solidBlockPrec += *M_solid->Mass();
    *M_solidBlockPrec += *M_solid->material()->jacobian(); //stiffMatrix();
    M_solidBlockPrec->globalAssemble();
    *M_solidBlockPrec *= M_solid->rescaleFactor();
}

//     M_solidBlockPrec.reset( new matrix_Type( *M_monolithicMap, 1 ) );
//     *M_solidBlockPrec += *M_solidBlock;
}

void
FSIMonolithic::assembleFluidBlock(UInt iter, const vector_Type& solution)
{
    M_fluidBlock.reset(new  FSIOperator::fluidPtr_Type::value_type::matrix_Type(*M_monolithicMap));

    Real alpha = M_fluidTimeAdvance->coefficientFirstDerivative( 0 )/M_data->dataFluid()->dataTime()->timeStep();//mesh velocity w
    if(!M_data->dataFluid()->conservativeFormulation())
      {
	M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, M_fluidBlock, solution );
      }
    else
      if (! M_fluid->matrixMassPtr().get() )
	M_fluid->buildSystem( );

    if (iter==0)
      {
        M_resetPrec=true;
        M_fluidTimeAdvance->updateRHSContribution( M_data->dataFluid()->dataTime()->timeStep() );
        *M_rhs += M_fluid->matrixMass() * ( M_fluidTimeAdvance->rhsContributionFirstDerivative() );
        couplingRhs(M_rhs);
      }
    if(M_data->dataFluid()->conservativeFormulation())
        M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, M_fluidBlock, solution );
    this->M_fluid->updateStabilization(*M_fluidBlock);
}

void FSIMonolithic::updateRHS()
{
    // Update fluid (iter == 0)
    M_fluidTimeAdvance->updateRHSContribution( M_data->dataFluid()->dataTime()->timeStep() );
    *M_rhs += M_fluid->matrixMass() * ( M_fluidTimeAdvance->rhsContributionFirstDerivative() );
    couplingRhs(M_rhs);

    // Update solid (iter == 0)
    updateSolidSystem(M_rhs);

    // Update RHS
    *M_rhsFull = *M_rhs;
}

namespace
{
static Real fZero(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{return 0.;}
static Real fOne(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{return 1.;}
}

void FSIMonolithic::enableStressComputation(UInt flag)
{
    M_BChWS.reset(new solidBchandler_Type());
    BCFunctionBase bcfZero(fZero);
    BCFunctionBase bcfOne(fOne);
    M_bcfWs.setFunctions_Robin(bcfOne,bcfOne);

    M_BChWS->addBC("WS", (UInt) flag, Robin, Full, M_bcfWs, 3);
}

FSIMonolithic::vectorPtr_Type FSIMonolithic::computeStress()
{
    vector_Type lambda(M_monolithicMatrix->interfaceMap());
    lambda.subset(M_fluidTimeAdvance->singleElement(0), M_solidAndFluidDim);

    M_boundaryMass.reset(new matrix_Type(*M_interfaceMap));
    if ( !M_BChWS->bcUpdateDone() )
        M_BChWS->bcUpdate(*M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    bcManageMatrix(*M_boundaryMass, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BChWS, M_dFESpace->feBd(), 1., dataSolid()->dataTime()->time() );
    M_boundaryMass->globalAssemble();

    solver_Type solverMass(M_epetraComm);
    solverMass.setDataFromGetPot( M_dataFile, "solid/solver" );
    solverMass.setupPreconditioner(M_dataFile, "solid/prec");//to avoid if we have already a prec.

    boost::shared_ptr<PreconditionerIfpack> P(new PreconditionerIfpack());

    vectorPtr_Type sol(new vector_Type(M_monolithicMatrix->interfaceMap()));
    solverMass.setMatrix(*M_boundaryMass);
    solverMass.setReusePreconditioner(false);
    solverMass.solveSystem( lambda, *sol, M_boundaryMass);

    EpetraExt::MultiVector_Reindex reindexMV(*M_interfaceMap->map(Unique));
    boost::shared_ptr<MapEpetra> newMap(new MapEpetra( *M_interfaceMap ));
    M_stress.reset(new vector_Type(reindexMV(lambda.epetraVector()), newMap, Unique));
    return M_stress;
}

void
FSIMonolithic::checkIfChangedFluxBC( precPtr_Type oper )
{
   UInt nfluxes(M_BChs[1]->numberOfBCWithType(Flux));
    if(M_fluxes != nfluxes)
    {
        //std::vector<bcName_Type> names = M_BChs[1]->findAllBCWithType(Flux);
        for (UInt i=0; i<M_fluxes; ++i)
        {
            const BCBase* bc (M_BChs[1]->findBCWithName(M_BCFluxNames[i]));
            if(bc->type() != Flux)
            {
                oper->addToCoupling(1., M_fluxOffset[i], M_fluxOffset[i],1 );
            }
        }
    }
}


}
