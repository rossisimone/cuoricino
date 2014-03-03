#ifndef STEKLOVPOINCAREOPERATOR_H
#define STEKLOVPOINCAREOPERATOR_H 1

namespace LifeV
{

typedef RegionMesh<LinearTetra>      mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
typedef MeshPartitioner<mesh_Type >  meshPartitioner_Type;


typedef FESpace<mesh_Type, MapEpetra >  FESpace_Type;
typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;

typedef MapEpetra 					map_Type;
typedef boost::shared_ptr<map_Type> mapPtr_Type;

typedef VectorEpetra 				   vector_Type;
typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

typedef MatrixEpetra<Real>             matrix_Type;
typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

typedef boost::shared_ptr<DOFInterface3Dto3D> dofInterfacePtr_Type;

typedef boost::shared_ptr<Epetra_Comm> communicatorPtr_Type;

class SteklovPoincareOperator
{
public:

	SteklovPoincareOperator(const FESpacePtr_Type& uFESpace,
							const FESpacePtr_Type& dFESpace,
							const boost::shared_ptr<DOFInterface3Dto3D>& dofFluidToStructure,
							const boost::shared_ptr<DOFInterface3Dto3D>& dofStructureToSolid);

    ~SteklovPoincareOperator();

    boost::shared_ptr<MapEpetra> createInterfaceMaps ( std::map<ID, ID> const& locDofMap,
    												   DOF& dof,
    												   const boost::shared_ptr<Epetra_Comm>& comm,
    												   const int type);

    void buildTranferOperators(	const mapPtr_Type& fluidInterfaceMap,
								const mapPtr_Type& structureInterfaceMap,
								const std::map<ID,ID> dofMap);

    void transferStructureOnFluid (const vectorPtr_Type& _vec1, vectorPtr_Type& _vec2);

private:

    FESpacePtr_Type M_uFESpace;
    FESpacePtr_Type M_dFESpace;
    dofInterfacePtr_Type M_dofFluidToStructure;
    dofInterfacePtr_Type M_dofStructureToSolid;
    mapPtr_Type M_interfaceMap;
    mapPtr_Type M_interfaceStructureMapScalar;
    mapPtr_Type M_interfaceFluidMapScalar;

    // Transfer operators
    matrixPtr_Type M_structureToGammaFluid;
};

SteklovPoincareOperator::SteklovPoincareOperator(const FESpacePtr_Type& uFESpace,
												 const FESpacePtr_Type& dFESpace,
												 const boost::shared_ptr<DOFInterface3Dto3D>& dofFluidToStructure,
												 const boost::shared_ptr<DOFInterface3Dto3D>& dofStructureToSolid):
M_uFESpace( uFESpace ),
M_dFESpace( dFESpace ),
M_dofFluidToStructure ( dofFluidToStructure ),
M_dofStructureToSolid ( dofStructureToSolid )
{}

SteklovPoincareOperator::~SteklovPoincareOperator()
{}

void SteklovPoincareOperator::buildTranferOperators(const mapPtr_Type& fluidInterfaceMap,
													const mapPtr_Type& structureInterfaceMap,
													const std::map<ID,ID> dofMap)
{
	M_structureToGammaFluid.reset( new matrix_Type( *M_interfaceFluidMapScalar, 50 ) );

	Real Value = 1.0;
	int * IndexRow = new int[1];
	int * IndexCol = new int[1];

	for( int i_map = 0; i_map <  M_interfaceFluidMapScalar->map(Unique)->NumMyElements(); ++i_map)
	{
		IndexRow[0] = M_interfaceFluidMapScalar->map(Unique)->GID(i_map);
		for ( std::map<ID, ID>::const_iterator i = dofMap.begin(); i != dofMap.end(); ++i )
			if ( i->first == M_interfaceFluidMapScalar->map(Unique)->GID(i_map) )
				IndexCol[0] = i->second;

		M_structureToGammaFluid->addToCoefficient(IndexRow[0], IndexCol[0], Value);
	}

	M_structureToGammaFluid->globalAssemble(M_interfaceStructureMapScalar, M_interfaceFluidMapScalar);
	M_structureToGammaFluid->spy("M_structureToFluid");

}

boost::shared_ptr<MapEpetra>
SteklovPoincareOperator::createInterfaceMaps ( std::map<ID, ID> const& locDofMap,
											  DOF& dof,
											  const boost::shared_ptr<Epetra_Comm>& comm,
											  const int type)
{
    std::vector<int> dofInterface;
    std::vector<int> dofInterfaceScalar;

    UInt nDimensions = 3;

    dofInterface.reserve ( locDofMap.size() );
    dofInterfaceScalar.reserve ( locDofMap.size() );

    for (UInt dim = 0; dim < nDimensions; ++dim)
    	for ( std::map<ID, ID>::const_iterator i = locDofMap.begin(); i != locDofMap.end(); ++i )
    	{
    		dofInterface.push_back (i->second+ dim * dof.numTotalDof() );
    		if( dim == 0 )
    			dofInterfaceScalar.push_back ( i->second );
    	}


    int* pointerToDofs(0);
    if (dofInterface.size() > 0)
        pointerToDofs = &dofInterface[0];

    boost::shared_ptr<MapEpetra> interfaceMap;
    interfaceMap.reset ( new MapEpetra ( -1, static_cast<int> (dofInterface.size() ), pointerToDofs, comm ) );

    int* pointerToDofsScalar(0);
    if (dofInterfaceScalar.size() > 0)
    	pointerToDofsScalar = &dofInterfaceScalar[0];

    switch(type)
    {
    case 0: //fluid
    	M_interfaceFluidMapScalar.reset ( new MapEpetra ( -1, static_cast<int> (dofInterfaceScalar.size() ), pointerToDofsScalar, comm ) );
    	break;
    case 1: //structure
    	M_interfaceStructureMapScalar.reset ( new MapEpetra ( -1, static_cast<int> (dofInterfaceScalar.size() ), pointerToDofsScalar, comm ) );
    	break;
    }

    boost::shared_ptr<MapEpetra> M_interfaceMapScalar;
    interfaceMap.reset ( new MapEpetra ( -1, static_cast<int> (dofInterface.size() ), pointerToDofs, comm ) );

    return interfaceMap;
}

void
SteklovPoincareOperator::transferStructureOnFluid (const vectorPtr_Type& vecStructure, vectorPtr_Type& vecFluidGamma)
{
	vector_Type vecStructure_x(*M_interfaceStructureMapScalar, Repeated);
	vector_Type vecStructure_y(*M_interfaceStructureMapScalar, Repeated);
	vector_Type vecStructure_z(*M_interfaceStructureMapScalar, Repeated);

	vecStructure_x.subset (*vecStructure, *M_interfaceStructureMapScalar, 0, 0);
	vecStructure_y.subset (*vecStructure, *M_interfaceStructureMapScalar, vecStructure->mapPtr()->map(Unique)->NumGlobalElements()/3, 0);
	vecStructure_z.subset (*vecStructure, *M_interfaceStructureMapScalar, vecStructure->mapPtr()->map(Unique)->NumGlobalElements()/3*2, 0);

	vector_Type vecFluidGamma_x(*M_interfaceFluidMapScalar, Repeated);
	vector_Type vecFluidGamma_y(*M_interfaceFluidMapScalar, Repeated);
	vector_Type vecFluidGamma_z(*M_interfaceFluidMapScalar, Repeated);

	M_structureToGammaFluid->multiply (false, vecStructure_x, vecFluidGamma_x);
	M_structureToGammaFluid->multiply (false, vecStructure_y, vecFluidGamma_y);
	M_structureToGammaFluid->multiply (false, vecStructure_z, vecFluidGamma_z);

	vecFluidGamma_x.spy("vecFluidGamma_x");
	vecFluidGamma_y.spy("vecFluidGamma_y");
	vecFluidGamma_z.spy("vecFluidGamma_z");

	/*
	vecFluidGamma->subset(vecFluidGamma_x, vecFluidGamma_x.map(), 0, 0);
	//vecFluidGamma->subset(vecFluidGamma_y, *M_interfaceFluidMapScalar, 0, vecFluidGamma->size()/3);
	//vecFluidGamma->subset(vecFluidGamma_z, *M_interfaceFluidMapScalar, 0, vecFluidGamma->size()/3*2);
	*/
}

} // end namespace LifeV

#endif
