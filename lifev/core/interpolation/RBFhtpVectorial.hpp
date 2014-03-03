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
    @brief A short description of the file content

    @author Davide Forti <forti@mathicsepc48.epfl.ch>
    @date 13 Mar 2013

    A more detailed description of the file (if necessary)
 */

#ifndef RBFHTPVECTORIAL_H
#define RBFHTPVECTORIAL_H 1

#include <lifev/core/interpolation/RBFInterpolation.hpp>

namespace LifeV
{

	typedef MatrixEpetra<double>                                                  matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                                        matrixPtr_Type;

template <typename mesh_Type>
class RBFhtpVectorial: public RBFInterpolation<mesh_Type>
{
public:

    typedef boost::shared_ptr<mesh_Type>                                          meshPtr_Type;

    typedef VectorEpetra                                                          vector_Type;
    typedef boost::shared_ptr<vector_Type >                                       vectorPtr_Type;

    //typedef MatrixEpetra<double>                                                  matrix_Type;
    //typedef boost::shared_ptr<matrix_Type>                                        matrixPtr_Type;

    typedef std::vector<int>                                                      flagContainer_Type;

    typedef boost::unordered_set<ID>                                                          idContainer_Type;

    typedef MapEpetra                                                             map_Type;
    typedef boost::shared_ptr<MapEpetra>                                          mapPtr_Type;

    typedef GhostHandler<mesh_Type>                                               neighbors_Type;
    typedef boost::shared_ptr<neighbors_Type>                                     neighborsPtr_Type;

    typedef LifeV::Preconditioner                                                 basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>                                      basePrecPtr_Type;

    typedef LifeV::PreconditionerIfpack                                           prec_Type;
    typedef boost::shared_ptr<prec_Type>                                          precPtr_Type;

    typedef Teuchos::RCP< Teuchos::ParameterList >                                parameterList_Type;

    RBFhtpVectorial();

    virtual ~RBFhtpVectorial();

    void setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags );

    void setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField);

    void setUpNeighbors ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags );

    void setUpNeighborsRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField);

    void buildOperators();

    void interpolationOperator();

    void projectionOperator();

    void buildRhs();

    void identifyNodes (meshPtr_Type LocalMesh, boost::unordered_set<ID>& GID_nodes, vectorPtr_Type CheckVector);

    bool isInside (ID pointMarker, flagContainer_Type Flags);

    double computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID);

    double rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius);

    void interpolate();

    void solution (vectorPtr_Type& Solution);

    void updateRhs(vectorPtr_Type newRhs);

    void approximateInverse();

    void buildInterpolationMatrix();

    matrixPtr_Type interpolationMatrix(){ return M_RescaledRBFMatrix;}

    void setUpNeighborsScalarMap(const mapPtr_Type& ScalarMap, const bool& mapOrganization, const bool& expand);

    mapPtr_Type scalarMap () { return M_scalarMap; }

    void createIdentityBlockVelocityCoupling(const mapPtr_Type& mappa);

    matrixPtr_Type identityStressBlock() {return M_identityStressBlock; }

    matrixPtr_Type identityVelocityBlock() {return M_identityVelocityBlock; }

    // some new getters

    mapPtr_Type rhsMap( ) { return M_rhsMap;} // ritorna la mappa per l'rhs in indici globali non riordinati

    mapPtr_Type solutionMap( ) { return M_solutionMap;} // ritorna la mappa per il vettore soluzione in indici globali non riordinati

    mapPtr_Type ProjectionMap(){ return M_projectionOperatorMap;} // ritorna la mappa della righe della matrice di proiezione con indici riordinati

    mapPtr_Type InterpolationMap(){ return M_interpolationOperatorMap;} // ritorna la mappa della righe della matrice di interpolazione con indici riordinati

    void createIdentityBlockStressCoupling(const mapPtr_Type& mappa);

    void buildMap(const meshPtr_Type & fullMeshKnown, const meshPtr_Type & localMeshKnown, const vectorPtr_Type & vector, flagContainer_Type flags_part,
                  mapPtr_Type & GIDindexMap, mapPtr_Type & reorderedMap, std::map<ID,ID> & GIDindexToReorderedMap, vectorPtr_Type & markerIDs);

    void expandInterpolationMatrixByRows(const mapPtr_Type& scalarMap);

    void expandInterpolationMatrixByColumns(const mapPtr_Type& scalarMap);

    matrixPtr_Type rowsExpandedInterpolationMatrix(){ return M_RowExpandedRescaledRBFMatrix;};

    matrixPtr_Type columnsExpandedInterpolationMatrix(){ return M_ColumnExpandedRescaledRBFMatrix;};

    matrixPtr_Type filteredMatrix(){return M_filteredRescaledRBFMatrix;};

    void createInterpolationSquareBlock(double value);

    matrixPtr_Type interpolationSquareBlock(){return M_interpolationSquareBlock;};

    // aggiunta domenica blocco identita

    void fillIdentityCouplingBlock(matrixPtr_Type & matrice, const mapPtr_Type & matrixMap, const Real & value);

    void interpolateVectors(const vectorPtr_Type& rhs, vectorPtr_Type & result);

    void filterMatrix(double flagValue, const mapPtr_Type &scalarMap);

    vectorPtr_Type orderedSolution(){return M_orderedSolution;};

    void restrictOnGammaOrdered(const vectorPtr_Type& vec, vectorPtr_Type vec_reordered);

    VectorEpetra vectorOnGamma(const vectorPtr_Type& vec);

    VectorEpetra globalToOrderedNumeration(const vectorPtr_Type& vec);

    matrixPtr_Type computeRobinBlock();

private:

    meshPtr_Type        M_fullMeshKnown;
    meshPtr_Type        M_localMeshKnown;

    meshPtr_Type        M_fullMeshUnknown;
    meshPtr_Type        M_localMeshUnknown;

    flagContainer_Type  M_flags;

    vectorPtr_Type      M_knownField;
    vectorPtr_Type      M_unknownField;

    idContainer_Type    M_GIdsKnownMesh;
    idContainer_Type    M_GIdsUnknownMesh;

    matrixPtr_Type      M_RBFMatrix;
    matrixPtr_Type      M_multiplicative;
    matrixPtr_Type      M_approximatedInverse;
    matrixPtr_Type      M_RescaledRBFMatrix;

    matrixPtr_Type      M_interpolationOperator;
    matrixPtr_Type      M_projectionOperator;

    vectorPtr_Type      M_RhsF1;
    vectorPtr_Type      M_RhsF2;
    vectorPtr_Type      M_RhsF3;

    vectorPtr_Type      M_RhsOne;

    mapPtr_Type         M_interpolationOperatorMap;
    mapPtr_Type         M_projectionOperatorMap;

    neighborsPtr_Type   M_neighbors;

    bool                M_rowMap;
    mapPtr_Type         M_scalarMap;

    //new stuff for new maps

    std::map<ID,ID>     M_globalToReorderedIDsKnown;
    std::map<ID,ID>     M_globalToReorderedIDsUnknown;

    mapPtr_Type         M_rhsMap;          // mappa per pezzo vettore di interfaccia
    mapPtr_Type         M_solutionMap;     // mappa per pezzo vettore da spedire nella soluzione
    vectorPtr_Type      M_RhsF1_reordered;
    vectorPtr_Type      M_RhsF2_reordered;
    vectorPtr_Type      M_RhsF3_reordered;
    vectorPtr_Type      M_solution;
    vectorPtr_Type      M_rhs;
    std::vector<int>    M_reordered_projection;
    std::vector<int>    M_flag_reordered_projection;
    bool                M_expand;
    matrixPtr_Type      M_identityStressBlock;
    matrixPtr_Type      M_identityVelocityBlock;
    std::map<ID,ID>     M_projectionMapGlobaltoLocal;

    matrixPtr_Type      M_RowExpandedRescaledRBFMatrix;
    matrixPtr_Type      M_ColumnExpandedRescaledRBFMatrix;

    vectorPtr_Type      M_markerIDsInterpolator;
    vectorPtr_Type      M_markerIDsProjector;

    matrixPtr_Type      M_interpolationSquareBlock;

    matrixPtr_Type      M_filteredRescaledRBFMatrix;

    mapPtr_Type         M_mapOrderedSolution;
    vectorPtr_Type      M_orderedSolution;

};

template <typename mesh_Type>
RBFhtpVectorial<mesh_Type>::RBFhtpVectorial()
{}

template <typename mesh_Type>
RBFhtpVectorial<mesh_Type>::~RBFhtpVectorial()
{}

template <typename mesh_Type>
VectorEpetra RBFhtpVectorial<mesh_Type>::globalToOrderedNumeration(const vectorPtr_Type& vec)
{
    vector_Type x(*M_rhsMap);
    vector_Type y(*M_rhsMap);
    vector_Type z(*M_rhsMap);

    x.subset (*vec, *M_rhsMap, 0, 0);
    y.subset (*vec, *M_rhsMap, vec->size()/3, 0);
    z.subset (*vec, *M_rhsMap, vec->size()/3*2, 0);

    vectorPtr_Type x_reordered;
    vectorPtr_Type y_reordered;
    vectorPtr_Type z_reordered;

    x_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    y_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    z_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );

    for(int i = 0; i < M_rhsMap->map(Unique)->NumMyElements(); ++i)
    {
        (*x_reordered)[x_reordered->blockMap().GID (i)] = x[x.blockMap().GID (i)];
        (*y_reordered)[y_reordered->blockMap().GID (i)] = y[y.blockMap().GID (i)];
        (*z_reordered)[z_reordered->blockMap().GID (i)] = z[z.blockMap().GID (i)];
    }

    map_Type mapReordered(*M_interpolationOperatorMap);
    mapReordered += *M_interpolationOperatorMap;
    mapReordered += *M_interpolationOperatorMap;

    vector_Type vec_reordered(mapReordered, Unique);

    vec_reordered.subset (*x_reordered, *M_interpolationOperatorMap, 0, 0);
    vec_reordered.subset (*y_reordered, *M_interpolationOperatorMap, 0, vec_reordered.size()/3);
    vec_reordered.subset (*z_reordered, *M_interpolationOperatorMap, 0, vec_reordered.size()/3*2);

    return vec_reordered;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField)
{
    M_knownField   = KnownField;
    M_unknownField = UnknownField;
}

template <typename mesh_Type>
VectorEpetra RBFhtpVectorial<mesh_Type>::vectorOnGamma(const vectorPtr_Type& vec)
{
    vector_Type x_comp(*M_rhsMap, Repeated);
    vector_Type y_comp(*M_rhsMap, Repeated);
    vector_Type z_comp(*M_rhsMap, Repeated);

    //ASSERT(0!=0,"PP");

    x_comp.subset (*vec, *M_rhsMap, 0, 0);
    y_comp.subset (*vec, *M_rhsMap, vec->mapPtr()->map(Unique)->NumGlobalElements()/3, 0);
    z_comp.subset (*vec, *M_rhsMap, vec->mapPtr()->map(Unique)->NumGlobalElements()/3*2, 0);

    vectorPtr_Type x_reordered;
    vectorPtr_Type y_reordered;
    vectorPtr_Type z_reordered;

    x_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    y_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    z_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );

    for(int i = 0; i < M_rhsMap->map(Unique)->NumMyElements(); ++i)
    {
        (*x_reordered)[x_reordered->blockMap().GID (i)] = x_comp[x_comp.blockMap().GID (i)];
        (*y_reordered)[y_reordered->blockMap().GID (i)] = y_comp[y_comp.blockMap().GID (i)];
        (*z_reordered)[z_reordered->blockMap().GID (i)] = z_comp[z_comp.blockMap().GID (i)];
    }

    mapPtr_Type mapReordered;
    mapReordered.reset(new map_Type(*M_interpolationOperatorMap));
    *mapReordered += *M_interpolationOperatorMap;
    *mapReordered += *M_interpolationOperatorMap;

    vector_Type vec_reordered(*mapReordered, Unique);

    vec_reordered.subset (*x_reordered, *M_interpolationOperatorMap, 0, 0);
    vec_reordered.subset (*y_reordered, *M_interpolationOperatorMap, 0, vec_reordered.size()/3);
    vec_reordered.subset (*z_reordered, *M_interpolationOperatorMap, 0, vec_reordered.size()/3*2);

    return vec_reordered;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::setUpNeighborsScalarMap(const mapPtr_Type& ScalarMap, const bool& mapOrganization, const bool& expand)
{
    M_scalarMap = ScalarMap;
    M_rowMap    = mapOrganization;
    M_expand    = expand;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags )
{
    M_fullMeshKnown = fullMeshKnown;
    M_localMeshKnown = localMeshKnown;
    M_fullMeshUnknown = fullMeshUnknown;
    M_localMeshUnknown = localMeshUnknown;
    M_flags = flags;
    M_expand = false;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::setUpNeighbors ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags )
{
    M_fullMeshKnown = fullMeshKnown;
    M_localMeshKnown = localMeshKnown;
    M_fullMeshUnknown = fullMeshUnknown;
    M_localMeshUnknown = localMeshUnknown;
    M_flags = flags;
    M_expand = false;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::setUpNeighborsRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField)
{
    M_knownField.reset(new vector_Type(KnownField->map(), Unique));
    M_unknownField.reset(new vector_Type(UnknownField->map(), Unique));
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::createInterpolationSquareBlock (double value)
{
    M_interpolationSquareBlock.reset(new matrix_Type(*M_interpolationOperatorMap,50));
    for( int i = 0 ; i < M_interpolationOperatorMap->map(Unique)->NumMyElements(); ++i )
        M_interpolationSquareBlock->addToCoefficient(M_interpolationOperatorMap->map(Unique)->GID(i), M_interpolationOperatorMap->map(Unique)->GID(i), value);
                                                                                                  //, M_projectionOperatorMap->map(Unique)->GID(i)
    M_interpolationSquareBlock->globalAssemble();
}

template <typename Mesh>
void RBFhtpVectorial<Mesh>::filterMatrix(double flagValue, const mapPtr_Type &scalarMap)
{
    matrixPtr_Type filter;
    filter.reset(new matrix_Type(*scalarMap, 50));
    for( int i = 0 ; i < M_solutionMap->map(Unique)->NumMyElements(); ++i )
        if( M_markerIDsProjector->blockMap().LID( M_markerIDsProjector->blockMap().GID(i) )!= -1 )
            if((*M_markerIDsProjector)[M_markerIDsProjector->blockMap().GID(i)]==flagValue )
                filter->addToCoefficient(M_solutionMap->map(Unique)->GID(i), M_solutionMap->map(Unique)->GID(i), 1.0);

    filter->globalAssemble();
    M_filteredRescaledRBFMatrix.reset(new matrix_Type (*scalarMap, 50) );
    filter->multiply ( false, *M_RowExpandedRescaledRBFMatrix, false, *M_filteredRescaledRBFMatrix, false);
    M_filteredRescaledRBFMatrix->globalAssemble(M_interpolationOperatorMap, scalarMap);
}

template <typename Mesh>
void RBFhtpVectorial<Mesh>::buildOperators()
{

    if(M_knownField->comm().MyPID()==0)
        std::cout << "\n[Assembling Interpolation and Projection operators ] -----> ";

    LifeChrono TimeBuilding;
    TimeBuilding.start();


    this->interpolationOperator();
    this->projectionOperator();

    TimeBuilding.stop();
    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "done in " << TimeBuilding.diff() << " s \n\n";


    this->buildRhs();

    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "[Computing the RBF interpolation matrix ] -----> ";

    LifeChrono TimeBuildingMatrix;
    TimeBuildingMatrix.start();

    this->approximateInverse();

    this->buildInterpolationMatrix();

    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "done in " << TimeBuildingMatrix.diff() << " s \n\n";
}

template <typename Mesh>
void RBFhtpVectorial<Mesh>::approximateInverse()
{
    M_multiplicative.reset(new matrix_Type(*M_interpolationOperator));
    *M_approximatedInverse -= *M_interpolationOperator;
    matrixPtr_Type result;

    for (int n = 2; n < 25; ++n)
    {
        if(n>2)
            M_multiplicative.reset(new matrix_Type(*result));
        result.reset( new matrix_Type (*M_interpolationOperatorMap, 5000));
        M_multiplicative->globalAssemble();
        M_multiplicative->multiply(false, *M_interpolationOperator, false, *result, true);
        ( n%2==0 ) ? *M_approximatedInverse += *result : *M_approximatedInverse -= *result;
    }
    M_approximatedInverse->globalAssemble();
    //M_approximatedInverse->spy("approssimata");
    M_multiplicative.reset();
    result.reset();
}

template <typename Mesh>
void RBFhtpVectorial<Mesh>::buildInterpolationMatrix()
{
    M_RBFMatrix.reset ( new matrix_Type (*M_projectionOperatorMap, 50) );
    M_projectionOperator->multiply(false, *M_approximatedInverse, false, *M_RBFMatrix, false);
    M_RBFMatrix->globalAssemble (M_interpolationOperatorMap, M_projectionOperatorMap);

    //M_RBFMatrix->spy("M_RBFMatrix");

    vectorPtr_Type solutionOne;
    solutionOne.reset (new vector_Type (*M_projectionOperatorMap) );
    M_RBFMatrix->multiply (false, *M_RhsOne, *solutionOne);

    matrixPtr_Type invDiag;
    invDiag.reset ( new matrix_Type(*M_projectionOperatorMap, 50) );

    double Values;
    int * Index = new int[1];

    for ( int i = 0 ; i < M_projectionOperatorMap->map(Unique)->NumMyElements(); ++i )
    {
        if(solutionOne->blockMap().LID(solutionOne->blockMap().GID(i))!=-1)
        {
            Index[0] = M_projectionOperator->matrixPtr()->RowMap().GID(i);
            Values = 1/((*solutionOne)[Index[0]]);
            invDiag->matrixPtr()->InsertGlobalValues (Index[0], 1, &Values, Index);
        }
    }
    invDiag->globalAssemble ();
    //invDiag->spy("invDiag");


    M_RescaledRBFMatrix.reset ( new matrix_Type(*M_projectionOperatorMap, 50) );
    invDiag->multiply(false, *M_RBFMatrix, false, *M_RescaledRBFMatrix, false);
    M_RescaledRBFMatrix->globalAssemble(M_interpolationOperatorMap, M_projectionOperatorMap);
    M_RescaledRBFMatrix->spy("RBFInterpolationMatrixHTP");
}

// questa funzione prende in ingresso la matrice di interpolazione e ne fa l'espansione sulle righe, ovvero dove proietti
template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::expandInterpolationMatrixByRows(const mapPtr_Type& scalarMap)
{
    matrixPtr_Type A;
    A.reset(new matrix_Type (*scalarMap, 50 ));
    double value = 1.0;

    for ( int i = 0 ; i < M_solutionMap->map(Unique)->NumMyElements(); ++i )
        if(this->isInside (M_fullMeshUnknown->point (M_solutionMap->map(Unique)->GID(i)).markerID(), M_flags))
            A->addToCoefficient(M_solutionMap->map(Unique)->GID(i),M_projectionOperatorMap->map(Unique)->GID(i),value);

    A->globalAssemble(M_projectionOperatorMap, scalarMap);
    M_RowExpandedRescaledRBFMatrix.reset(new matrix_Type (*scalarMap, 50) );
    A->multiply ( false, *M_RescaledRBFMatrix, false, *M_RowExpandedRescaledRBFMatrix, false);
    M_RowExpandedRescaledRBFMatrix->globalAssemble(M_interpolationOperatorMap, scalarMap);
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::fillIdentityCouplingBlock(matrixPtr_Type & matrice, const mapPtr_Type & matrixMap, const Real & value)
{
    for(int i = 0 ; i < M_rhsMap->map(Unique)->NumMyElements(); ++i)
        matrice->addToCoefficient(M_rhsMap->map(Unique)->GID(i), M_interpolationOperatorMap->map(Unique)->GID(i), value);

    matrice->globalAssemble(M_interpolationOperatorMap, matrixMap);
    //matrice->spy("CouplingSforzi");

}

// questa funzione prende in ingresso la matrice di interpolazione e ne fa l'espansione sulle righe, ovvero dove proietti
template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::expandInterpolationMatrixByColumns(const mapPtr_Type& scalarMap)
{
    matrixPtr_Type A;
    A.reset(new matrix_Type (*M_interpolationOperatorMap, 50 ));
    double value = 1.0;

    for(int i = 0 ; i < M_interpolationOperatorMap->map(Unique)->NumMyElements(); ++i)
        A->addToCoefficient(M_interpolationOperatorMap->map(Unique)->GID(i),M_rhsMap->map(Unique)->GID(i),value);

    A->globalAssemble( scalarMap, M_interpolationOperatorMap );
    M_ColumnExpandedRescaledRBFMatrix.reset(new matrix_Type (*M_projectionOperatorMap, 50 ));
    M_RescaledRBFMatrix->multiply ( false, *A, false, *M_ColumnExpandedRescaledRBFMatrix, false);
    M_ColumnExpandedRescaledRBFMatrix->globalAssemble(scalarMap, M_projectionOperatorMap );
}

// QUI IL VETTORE DI INPUT E' UN VETTORE AVENTE 3 COMPONENTI
template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::interpolateVectors(const vectorPtr_Type& rhs, vectorPtr_Type & result)
{
    vectorPtr_Type solution1;
    solution1.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution2;
    solution2.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution3;
    solution3.reset (new vector_Type (*M_projectionOperatorMap) );

    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "[Interpolate between vectors ] -----> ";

    *M_RhsF1_reordered *= 0;
    *M_RhsF2_reordered *= 0;
    *M_RhsF3_reordered *= 0;

    M_RhsF1_reordered->subset (*rhs, *M_interpolationOperatorMap, 0, 0);
    M_RhsF2_reordered->subset (*rhs, *M_interpolationOperatorMap, rhs->size()/3, 0);
    M_RhsF3_reordered->subset (*rhs, *M_interpolationOperatorMap, rhs->size()*2/3, 0);

    M_RescaledRBFMatrix->multiply (false, *M_RhsF1_reordered, *solution1);
    M_RescaledRBFMatrix->multiply (false, *M_RhsF2_reordered, *solution2);
    M_RescaledRBFMatrix->multiply (false, *M_RhsF3_reordered, *solution3);

    result->subset (*solution1, *M_projectionOperatorMap, 0, 0);
    result->subset (*solution2, *M_projectionOperatorMap, 0, result->size()/3);
    result->subset (*solution3, *M_projectionOperatorMap, 0, result->size()/3*2);
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::interpolate()
{

    vectorPtr_Type solution1;
    solution1.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution2;
    solution2.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution3;
    solution3.reset (new vector_Type (*M_projectionOperatorMap) );

    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "[Interpolate ] -----> ";
    LifeChrono TimeInterpolate;
    TimeInterpolate.start();

    M_RescaledRBFMatrix->multiply (false, *M_RhsF1_reordered, *solution1);
    M_RescaledRBFMatrix->multiply (false, *M_RhsF2_reordered, *solution2);
    M_RescaledRBFMatrix->multiply (false, *M_RhsF3_reordered, *solution3);

    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "done in " << TimeInterpolate.diff() << " s \n\n";

    vectorPtr_Type solution1_reordered;
    vectorPtr_Type solution2_reordered;
    vectorPtr_Type solution3_reordered;

    solution1_reordered.reset (new vector_Type (*M_solutionMap) );
    solution2_reordered.reset (new vector_Type (*M_solutionMap) );
    solution3_reordered.reset (new vector_Type (*M_solutionMap) );

    /*
    M_mapOrderedSolution.reset(new map_Type(*M_projectionOperatorMap));
    *M_mapOrderedSolution += *M_projectionOperatorMap;
    *M_mapOrderedSolution += *M_projectionOperatorMap;

    M_orderedSolution.reset(new vector_Type(*M_mapOrderedSolution,Unique));
    M_orderedSolution->subset(*solution1,*M_projectionOperatorMap, 0, 0);
    M_orderedSolution->subset(*solution2,*M_projectionOperatorMap, 0, M_orderedSolution->size()/3);
    M_orderedSolution->subset(*solution3,*M_projectionOperatorMap, 0, M_orderedSolution->size()/3*2);
    */

    for(int i = 0; i < M_projectionOperatorMap->map(Unique)->NumMyElements(); ++i)
    {
        (*solution1_reordered)[solution1_reordered->blockMap().GID (i)] = (*solution1)[solution1->blockMap().GID (i)];
        (*solution2_reordered)[solution2_reordered->blockMap().GID (i)] = (*solution2)[solution2->blockMap().GID (i)];
        (*solution3_reordered)[solution3_reordered->blockMap().GID (i)] = (*solution3)[solution3->blockMap().GID (i)];
    }

    M_unknownField->subset (*solution1_reordered, *M_solutionMap, 0, 0);
    M_unknownField->subset (*solution2_reordered, *M_solutionMap, 0, M_unknownField->size()/3);
    M_unknownField->subset (*solution3_reordered, *M_solutionMap, 0, M_unknownField->size()/3*2);
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::restrictOnGammaOrdered(const vectorPtr_Type& vec, vectorPtr_Type vec_reordered)
{
    vector_Type x(*M_rhsMap);
    vector_Type y(*M_rhsMap);
    vector_Type z(*M_rhsMap);

    x.subset (*vec, *M_rhsMap, 0, 0);
    y.subset (*vec, *M_rhsMap, vec->size()/3, 0);
    z.subset (*vec, *M_rhsMap, vec->size()/3*2, 0);

    vectorPtr_Type x_reordered;
    vectorPtr_Type y_reordered;
    vectorPtr_Type z_reordered;

    x_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    y_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    z_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );

    for(int i = 0; i < M_rhsMap->map(Unique)->NumMyElements(); ++i)
    {
        (*x_reordered)[x_reordered->blockMap().GID (i)] = x[x.blockMap().GID (i)];
        (*y_reordered)[y_reordered->blockMap().GID (i)] = y[y.blockMap().GID (i)];
        (*z_reordered)[z_reordered->blockMap().GID (i)] = z[z.blockMap().GID (i)];
    }

    map_Type mapReordered(*M_interpolationOperatorMap);
    mapReordered += *M_interpolationOperatorMap;
    mapReordered += *M_interpolationOperatorMap;

    //vectorPtr_Type vec_reordered;
    vec_reordered.reset(new vector_Type(mapReordered, Unique));

    vec_reordered->subset (*x_reordered, *M_interpolationOperatorMap, 0, 0);
    vec_reordered->subset (*y_reordered, *M_interpolationOperatorMap, 0, vec_reordered->size()/3);
    vec_reordered->subset (*z_reordered, *M_interpolationOperatorMap, 0, vec_reordered->size()/3*2);
    vec_reordered->globalAssemble();
    //return vec_reordered;
}


template <typename Mesh>
void RBFhtpVectorial<Mesh>::buildMap(const meshPtr_Type & fullMesh, const meshPtr_Type & localMesh, const vectorPtr_Type & vector, flagContainer_Type flags,
                                     mapPtr_Type & GIDglobalMap, mapPtr_Type & GIDreorderedMap, std::map<ID,ID> & MapGIDindexToReordered, vectorPtr_Type & markerIDs)
{
    boost::unordered_set<ID> fullGIDs;
    boost::unordered_set<ID> localGIDs;
    boost::unordered_set<ID> fullGIDsReordered;
    boost::unordered_set<ID> localGIDsReordered;

    for ( UInt i = 0; i < fullMesh->numVertices(); ++i )
        if ( this->isInside (fullMesh->point (i).markerID(), flags) )
            fullGIDs.insert (fullMesh->point (i).id() );

    this->identifyNodes (localMesh, localGIDs, vector);

    int k = -1;
    for(boost::unordered_set<ID>::iterator itG = fullGIDs.begin(); itG != fullGIDs.end(); ++itG)
    {
        ++k;
        fullGIDsReordered.insert(k);
        MapGIDindexToReordered.insert(std::pair<ID,ID>(*itG,k));
    }

    k = -1;
    int kk = 0;
    int LocalNodesNumber = localGIDs.size();
    int* globalID = new int[LocalNodesNumber];

    for(boost::unordered_set<ID>::iterator itL = localGIDs.begin(); itL != localGIDs.end(); ++itL)
    {
        globalID[kk] = *itL;
        for(boost::unordered_set<ID>::iterator itG = fullGIDs.begin(); itG != fullGIDs.end(); ++itG)
        {
            ++k;
            if(*itL==*itG)
            {
                localGIDsReordered.insert(k);
            }
        }
        ++kk;
        k = -1;
    }

    int* reorderedID = new int[LocalNodesNumber];
    k = 0;

    for(boost::unordered_set<ID>::iterator it = localGIDsReordered.begin(); it != localGIDsReordered.end(); ++it)
    {
        reorderedID[k] = *it;
        ++k;
    }

    GIDglobalMap.reset (new map_Type (-1, LocalNodesNumber, globalID, vector->mapPtr()->commPtr() ) );
    GIDreorderedMap.reset (new map_Type (-1, LocalNodesNumber, reorderedID, vector->mapPtr()->commPtr() ) );

    markerIDs.reset(new vector_Type (*GIDreorderedMap, Unique));
    for(int i = 0 ; i < GIDglobalMap->map(Unique)->NumMyElements(); ++i)
        (*markerIDs)[GIDreorderedMap->map(Unique)->GID(i)] = fullMesh->point(GIDglobalMap->map(Unique)->GID(i)).markerID();

    delete globalID;
    delete reorderedID;
}


template <typename Mesh>
void RBFhtpVectorial<Mesh>::interpolationOperator()
{

    this->buildMap(M_fullMeshKnown, M_localMeshKnown, M_knownField, M_flags, M_rhsMap, M_interpolationOperatorMap, M_globalToReorderedIDsKnown, M_markerIDsInterpolator);

    M_neighbors.reset ( new neighbors_Type ( M_fullMeshKnown, M_localMeshKnown, M_knownField->mapPtr(), M_knownField->mapPtr()->commPtr() ) );

    if (M_flags[0] == -1)
    {
        M_neighbors->setUpNeighbors();
    }
    else
    {
        M_neighbors->createPointPointNeighborsList (M_flags);
    }

    std::vector<boost::unordered_set<ID> > MatrixGraph (M_interpolationOperatorMap->map(Unique)->NumMyElements());
    std::vector<boost::unordered_set<ID> > MatrixGraphReordered (M_interpolationOperatorMap->map(Unique)->NumMyElements());

    int LocalNodesNumber = M_interpolationOperatorMap->map(Unique)->NumMyElements();
    std::vector<double> RBF_radius (LocalNodesNumber);

    for (int i = 0 ; i < LocalNodesNumber; ++i)
    {
        MatrixGraph[i] = M_neighbors->pointPointNeighborsList() [M_fullMeshKnown->point(M_rhsMap->map(Unique)->GID(i)).id()];
        for(boost::unordered_set<ID>::iterator it_matrix_graph = MatrixGraph[i].begin(); it_matrix_graph != MatrixGraph[i].end(); ++it_matrix_graph)
        {
            for(std::map<ID,ID>::iterator it_map = M_globalToReorderedIDsKnown.begin(); it_map != M_globalToReorderedIDsKnown.end(); ++it_map)
            {
                if(it_map->first == *it_matrix_graph)
                    MatrixGraphReordered[i].insert(it_map->second);
            }
        }
        RBF_radius[i] = computeRBFradius ( M_fullMeshKnown, M_fullMeshKnown, MatrixGraph[i], M_rhsMap->map(Unique)->GID(i));
    }

    M_interpolationOperator.reset (new matrix_Type (*M_interpolationOperatorMap, 50) );
    M_approximatedInverse.reset (new matrix_Type (*M_interpolationOperatorMap, 1) );

    int* Indices = new int[M_interpolationOperatorMap->map(Unique)->NumGlobalElements()];
    double* Values = new double[M_interpolationOperatorMap->map(Unique)->NumGlobalElements()];
    double d = 1.0;
    int k;

    for ( int i = 0 ; i < M_interpolationOperatorMap->map(Unique)->NumMyElements(); ++i )
    {
        k = 0;
        for ( boost::unordered_set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
        {
            Values[k]  = rbf ( M_fullMeshKnown->point (M_rhsMap->map(Unique)->GID(i)).x(),
                               M_fullMeshKnown->point (M_rhsMap->map(Unique)->GID(i)).y(),
                               M_fullMeshKnown->point (M_rhsMap->map(Unique)->GID(i)).z(),
                               M_fullMeshKnown->point (*it).x(),
                               M_fullMeshKnown->point (*it).y(),
                               M_fullMeshKnown->point (*it).z(),
                               RBF_radius[i]);
            ++k;
        }
        k = 0;
        for ( boost::unordered_set<ID>::iterator it = MatrixGraphReordered[i].begin(); it != MatrixGraphReordered[i].end(); ++it)
        {
            Indices[k] = *it;
            ++k;
        }
        M_interpolationOperator->matrixPtr()->InsertGlobalValues (M_interpolationOperatorMap->map(Unique)->GID(i), k, Values, Indices);
        M_approximatedInverse->addToCoefficient (M_interpolationOperatorMap->map(Unique)->GID(i), M_interpolationOperatorMap->map(Unique)->GID(i), d);
    }

    M_interpolationOperator->globalAssemble();
    //M_interpolationOperator->spy("matrice1");

    delete Indices;
    delete Values;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::projectionOperator()
{
    this->buildMap(M_fullMeshUnknown, M_localMeshUnknown, M_unknownField, M_flags, M_solutionMap, M_projectionOperatorMap, M_globalToReorderedIDsUnknown,M_markerIDsProjector);

    int LocalNodesNumber = M_projectionOperatorMap->map(Unique)->NumMyElements();

    std::vector<double>        RBF_radius (LocalNodesNumber);
    std::vector<boost::unordered_set<ID> > MatrixGraph (LocalNodesNumber);
    std::vector<boost::unordered_set<ID> > MatrixGraphReordered (LocalNodesNumber);

    double d;
    double d_min;
    int    nearestPoint;

    for (int i = 0 ; i < LocalNodesNumber; ++i)
    {
        d_min = 100;
        for (int j = 0; j <  M_fullMeshKnown->numVertices(); ++j)
        {
            if ( M_flags[0] == -1 || this->isInside (M_fullMeshKnown->point (j).markerID(), M_flags) )
            {
                d = std::sqrt(std::pow (M_fullMeshKnown->point (j).x() - M_fullMeshUnknown->point (M_solutionMap->map(Unique)->GID(i)).x(), 2)
                            + std::pow (M_fullMeshKnown->point (j).y() - M_fullMeshUnknown->point (M_solutionMap->map(Unique)->GID(i)).y(), 2)
                            + std::pow (M_fullMeshKnown->point (j).z() - M_fullMeshUnknown->point (M_solutionMap->map(Unique)->GID(i)).z(), 2)
                             );
                if (d < d_min)
                {
                    d_min = d;
                    nearestPoint = M_fullMeshKnown->point (j).id();
                }
            }
        }
        MatrixGraph[i] = M_neighbors->pointPointNeighborsList() [nearestPoint];
        MatrixGraph[i].insert (nearestPoint);
        for(boost::unordered_set<ID>::iterator it_matrix_graph = MatrixGraph[i].begin(); it_matrix_graph != MatrixGraph[i].end(); ++it_matrix_graph)
        {
            for(std::map<ID,ID>::iterator it_map = M_globalToReorderedIDsKnown.begin(); it_map != M_globalToReorderedIDsKnown.end(); ++it_map)
            {
                if(it_map->first == *it_matrix_graph)
                    MatrixGraphReordered[i].insert(it_map->second);
            }
        }
        RBF_radius[i] = computeRBFradius ( M_fullMeshKnown, M_fullMeshUnknown, MatrixGraph[i], M_fullMeshUnknown->point(M_solutionMap->map(Unique)->GID(i)).id());
    }

    M_projectionOperator.reset (new matrix_Type (*M_projectionOperatorMap, 50) );

    int* Indices = new int[M_projectionOperatorMap->map(Unique)->NumGlobalElements()];
    double* Values = new double[M_projectionOperatorMap->map(Unique)->NumGlobalElements()];
    int k = 0;
    for ( int i = 0 ; i < LocalNodesNumber; ++i )
    {
        k = 0;
        for ( boost::unordered_set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
        {
            Values[k]  = rbf ( M_fullMeshUnknown->point (M_solutionMap->map(Unique)->GID(i)).x(),
                               M_fullMeshUnknown->point (M_solutionMap->map(Unique)->GID(i)).y(),
                               M_fullMeshUnknown->point (M_solutionMap->map(Unique)->GID(i)).z(),
                               M_fullMeshKnown->point (*it).x(),
                               M_fullMeshKnown->point (*it).y(),
                               M_fullMeshKnown->point (*it).z(),
                               RBF_radius[i]);
            ++k;
        }
        k = 0;
        for ( boost::unordered_set<ID>::iterator it = MatrixGraphReordered[i].begin(); it != MatrixGraphReordered[i].end(); ++it)
        {
            Indices[k] = *it;
            ++k;
        }
        M_projectionOperator->matrixPtr()->InsertGlobalValues (M_projectionOperatorMap->map(Unique)->GID(i), k, Values, Indices);
    }
    M_projectionOperator->globalAssemble (M_interpolationOperatorMap, M_projectionOperatorMap);
    //M_projectionOperator->spy("matrice2");

    delete Indices;
    delete Values;
}

template <typename mesh_Type>
double RBFhtpVectorial<mesh_Type>::computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID)
{
    double r = 0;
    double r_max = 0;
    for (idContainer_Type::iterator it = Neighbors.begin(); it != Neighbors.end(); ++it)
    {
        r = std::sqrt ( std::pow ( MeshGID->point ( GlobalID ).x() - MeshNeighbors->point ( *it ).x(), 2 )
                        + std::pow ( MeshGID->point ( GlobalID ).y() - MeshNeighbors->point ( *it ).y(), 2 )
                        + std::pow ( MeshGID->point ( GlobalID ).z() - MeshNeighbors->point ( *it ).z(), 2 ) );
        r_max = ( r > r_max ) ? r : r_max;
    }
    return r_max;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::buildRhs()
{
    M_RhsF1.reset (new vector_Type (*M_rhsMap) );
    M_RhsF2.reset (new vector_Type (*M_rhsMap) );
    M_RhsF3.reset (new vector_Type (*M_rhsMap) );
    M_RhsOne.reset (new vector_Type (*M_interpolationOperatorMap) );

    M_RhsF1->subset (*M_knownField, *M_rhsMap, 0, 0);
    M_RhsF2->subset (*M_knownField, *M_rhsMap, M_knownField->size()/3, 0);
    M_RhsF3->subset (*M_knownField, *M_rhsMap, M_knownField->size()/3*2, 0);

    M_RhsF1_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF2_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF3_reordered.reset (new vector_Type (*M_interpolationOperatorMap) );

    for(int i = 0; i < M_rhsMap->map(Unique)->NumMyElements(); ++i)
    {
        (*M_RhsF1_reordered)[M_RhsF1_reordered->blockMap().GID (i)] = (*M_RhsF1)[M_RhsF1->blockMap().GID (i)];
        (*M_RhsF2_reordered)[M_RhsF2_reordered->blockMap().GID (i)] = (*M_RhsF2)[M_RhsF2->blockMap().GID (i)];
        (*M_RhsF3_reordered)[M_RhsF3_reordered->blockMap().GID (i)] = (*M_RhsF3)[M_RhsF3->blockMap().GID (i)];
    }

    *M_RhsOne += 1;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::identifyNodes (meshPtr_Type LocalMesh, boost::unordered_set<ID>& GID_nodes, vectorPtr_Type CheckVector)
{
    if (M_flags[0] == -1)
    {
        for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
            if (CheckVector->blockMap().LID (LocalMesh->point (i).id() ) != -1)
            {
                GID_nodes.insert (LocalMesh->point (i).id() );
            }
    }
    else
    {
        for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
            if ( this->isInside (LocalMesh->point (i).markerID(), M_flags) )
                if (CheckVector->blockMap().LID (LocalMesh->point (i).id() ) != -1)
                {
                    GID_nodes.insert (LocalMesh->point (i).id() );
                }
    }
}

template <typename mesh_Type>
bool RBFhtpVectorial<mesh_Type>::isInside (ID pointMarker, flagContainer_Type flags)
{
    int check = 0;
    for (UInt i = 0; i < flags.size(); ++i)
        if (pointMarker == flags[i])
        {
            ++check;
        }
    return (check > 0) ? true : false;
}

template <typename mesh_Type>
double RBFhtpVectorial<mesh_Type>::rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius)
{
    double distance = std::sqrt ( std::pow (x1 - x2, 2) + std::pow (y1 - y2, 2) + std::pow (z1 - z2, 2) );
    return std::pow (1 - distance / radius, 4) * (4 * distance / radius + 1);
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::updateRhs(vectorPtr_Type newRhs)
{
    *M_knownField = *newRhs;
    this->buildRhs();

    //M_RhsF1_reordered->spy("M_RhsF1_reordered");
    //M_RhsF2_reordered->spy("M_RhsF2_reordered");
    //M_RhsF3_reordered->spy("M_RhsF3_reordered");
    /*
    *M_RhsF1 *= 0;
    M_RhsF1->subset (*newRhs, *M_rhsMap, 0, 0);
    *M_RhsF2 *= 0;
    M_RhsF2->subset (*newRhs, *M_rhsMap, newRhs->size()/3, 0);
    *M_RhsF3 *= 0;
    M_RhsF3->subset (*newRhs, *M_rhsMap, newRhs->size()/3*2, 0);
    /*
    // ASSERT(0!=0, "Stop programmato RBFhtpVectorial.hpp linea 796");

    *M_RhsF1_reordered *= 0;
    *M_RhsF2_reordered *= 0;
    *M_RhsF3_reordered *= 0;

    for(int i = 0; i < M_rhsMap->map(Unique)->NumMyElements(); ++i)
    {
        (*M_RhsF1_reordered)[M_RhsF1_reordered->blockMap().GID (i)] = (*M_RhsF1)[M_RhsF1->blockMap().GID (i)];
        (*M_RhsF2_reordered)[M_RhsF2_reordered->blockMap().GID (i)] = (*M_RhsF2)[M_RhsF2->blockMap().GID (i)];
        (*M_RhsF3_reordered)[M_RhsF3_reordered->blockMap().GID (i)] = (*M_RhsF3)[M_RhsF3->blockMap().GID (i)];
    }
    */

}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::solution (vectorPtr_Type& Solution)
{
    Solution = M_unknownField;
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::createIdentityBlockVelocityCoupling(const mapPtr_Type & map)
{
    M_identityVelocityBlock.reset (new matrix_Type (*(M_interpolationOperatorMap), 50) );
    double value = 1.0;
    for( int i = 0 ; i < M_interpolationOperatorMap->map(Unique)->NumMyElements(); ++i )
    {
        M_identityVelocityBlock->addToCoefficient(M_interpolationOperatorMap->map(Unique)->GID(i), M_rhsMap->map(Unique)->GID(i), value);
    }
    M_identityVelocityBlock->globalAssemble(map, M_interpolationOperatorMap);
}

template <typename mesh_Type>
void RBFhtpVectorial<mesh_Type>::createIdentityBlockStressCoupling(const mapPtr_Type & map)
{
    M_identityStressBlock.reset (new matrix_Type (*(map), 50) );
    double value = -1.0;

    for( int i = 0 ; i < M_rhsMap->map(Unique)->NumMyElements(); ++i )
        M_identityStressBlock->addToCoefficient(M_rhsMap->map(Unique)->GID(i), M_interpolationOperatorMap->map(Unique)->GID(i), value);

    M_identityStressBlock->globalAssemble(M_interpolationOperatorMap, map);
}

template <typename mesh_Type>
matrixPtr_Type RBFhtpVectorial<mesh_Type>::computeRobinBlock()
{
	matrixPtr_Type robinBlock;
	robinBlock.reset(new matrix_Type(*M_interpolationOperatorMap, 50));
	for( int i = 0 ; i < M_interpolationOperatorMap->map(Unique)->NumMyElements(); ++i )
		if( M_markerIDsInterpolator->blockMap().LID( M_markerIDsInterpolator->blockMap().GID(i) )!= -1 )
			robinBlock->addToCoefficient(M_interpolationOperatorMap->map(Unique)->GID(i), M_solutionMap->map(Unique)->GID(i), 1.0);

	robinBlock->globalAssemble(M_solutionMap, M_interpolationOperatorMap);
	return robinBlock;
}

//! Factory create function
template <typename mesh_Type>
inline RBFInterpolation<mesh_Type> * createRBFhtpVectorial()
{
    return new RBFhtpVectorial< mesh_Type > ();
}
namespace
{
static bool S_registerTriHTPV = RBFInterpolation<LifeV::RegionMesh<LinearTriangle > >::InterpolationFactory::instance().registerProduct ( "RBFhtpVectorial", &createRBFhtpVectorial<LifeV::RegionMesh<LinearTriangle > > );
static bool S_registerTetHTPV = RBFInterpolation<LifeV::RegionMesh<LinearTetra > >::InterpolationFactory::instance().registerProduct ( "RBFhtpVectorial", &createRBFhtpVectorial<LifeV::RegionMesh<LinearTetra > > );
}

} // Namespace LifeV

#endif /* RBFHTPVECTORIAL_H */
