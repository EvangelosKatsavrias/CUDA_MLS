#include"App_MeshAdapter.h"
#include<typeinfo>
#include<limits>


template<class T>
app_meshAdapter<T>::app_meshAdapter ( mesh_t* initialMesh, meshAdapter_t* boundaryAdapter ): m_initialMesh(initialMesh), m_meshBoundaryAdapter(boundaryAdapter), m_adaptedMesh(0), m_memAllocatedMembers( std::vector<bool>(3, 0) )
{ }


template<class T>
app_meshAdapter<T>::~app_meshAdapter() 
{ 
	if ( m_memAllocatedMembers[0] ) delete m_initialMesh; 
	if ( m_memAllocatedMembers[1] ) delete m_adaptedMesh; 
	if ( m_memAllocatedMembers[2] ) delete m_meshBoundaryAdapter; 
}


template<class T>
int app_meshAdapter<T>::set_initialMesh ( mesh_t* mesh ) {
	if ( m_memAllocatedMembers[0] ) delete m_initialMesh; m_initialMesh = mesh; m_memAllocatedMembers[0] = 0;
	if ( m_meshBoundaryAdapter ) m_meshBoundaryAdapter->set_initialTopology ( m_initialMesh );
}


template<class T>
int app_meshAdapter<T>::set_initialMesh (  ) {

	if ( !m_memAllocatedMembers[0] ) {
		switch ( m_domainCardinality ) {
			case 2:
				m_initialMesh = new Mesh2D<T, 3, int>;
				break;
			case 3:
				m_initialMesh = new Mesh3D<T, 3, int, 4>;
				break;
		}
		m_memAllocatedMembers[0] = 1;
	}

	if ( m_meshBoundaryAdapter ) m_meshBoundaryAdapter->set_initialTopology ( m_initialMesh );

	std::ifstream newFileStream; newFileStream.exceptions(std::ios_base::badbit); 
	

	newFileStream.open(m_nodesFile);
		int numOfTotalNodes; newFileStream >> numOfTotalNodes;
		m_initialMesh->set_numOfNodes( numOfTotalNodes );
		m_initialMesh->set_nodeFlags( newFileStream );
		m_initialMesh->remove_allComponents();
		for ( size_t component = 0; component < m_domainCardinality; component++ )
			m_initialMesh->insert_components ( newFileStream );
	newFileStream.close();


	std::set<int> interiorFlags = std::set<int>( {0} );
	m_initialMesh->set_interiorFlagValues( interiorFlags );


	newFileStream.open(m_connectivitiesFile);
		int numOfConnectivities; newFileStream >> numOfConnectivities;
		m_initialMesh->set_numOfConnectivities(numOfConnectivities);
		m_initialMesh->set_connectivities(newFileStream, -1);
	newFileStream.close();

	m_initialMesh->find_metrics();

	return 1;
}


template<class T>
int app_meshAdapter<T>::set_adaptedMesh ( mesh_t* mesh ) {
	if ( m_memAllocatedMembers[1] ) delete m_adaptedMesh; m_adaptedMesh = mesh; m_memAllocatedMembers[1] = 0;
	if ( m_meshBoundaryAdapter ) m_meshBoundaryAdapter->set_adaptedTopology ( m_adaptedMesh );
}


template<class T>
int app_meshAdapter<T>::set_adaptedMesh (  ) {

	if ( !m_memAllocatedMembers[1] ) m_memAllocatedMembers[1] = 1;
	else delete m_adaptedMesh;

	switch ( m_domainCardinality ) {
		case 2:
			m_adaptedMesh = new Mesh2D<T, 3, int>( *m_initialMesh );
			break;
		case 3:
			m_adaptedMesh = new Mesh3D<T, 3, int, 4>( *m_initialMesh );
			break;
	}

	if ( m_meshBoundaryAdapter ) m_meshBoundaryAdapter->set_adaptedTopology ( m_adaptedMesh );
}


template<class T>
int app_meshAdapter<T>::set_meshAdapter ( ) 
{
	if ( !m_meshBoundaryAdapter ) {
		m_meshBoundaryAdapter = new meshAdapter_t( m_initialMesh, m_adaptedMesh); m_memAllocatedMembers[2] = 1;
		m_meshBoundaryAdapter->set_solutionMethod	( m_option_solutionMethod );
		m_meshBoundaryAdapter->set_weightingFunctionType( m_option_weightingFunctionType );
		m_meshBoundaryAdapter->set_computingMode	( m_option_computingMode );
		m_meshBoundaryAdapter->set_linearSolverType	( m_option_linearSolverType );
		m_meshBoundaryAdapter->set_polynomialDegree	( m_polynomialDegree );
		m_meshBoundaryAdapter->set_weightingFunctionSpan( m_weightingFunctionSpan );
		m_meshBoundaryAdapter->set_cutOffMultiplier	( m_cutOffMultiplier );
		m_meshBoundaryAdapter->set_interpolationFactor	( m_interpolationFactor );
	}


	m_meshBoundaryAdapter->set_initialTopology ( m_initialMesh );
	m_meshBoundaryAdapter->set_adaptedTopology ( m_adaptedMesh );

	m_meshBoundaryAdapter->set_perturbations( );
	m_meshBoundaryAdapter->set_allNodes( );


	size_t numOfComponents = get_initialMesh()->get_numOfComponents();
	TransformationMatrix<T> trMatrix( numOfComponents ), trMatrixFinal( numOfComponents );
	std::vector<T> vectorValues;


	for ( auto& i : m_perturbedNodeSets ) {

		auto sources 		= m_dataSources[i].begin();
		auto transformations 	= m_transformations[i].begin();
		auto translations 	= m_translations[i].begin();
		auto rotationsEulerAngles= m_rotationsEulerAngles[i].begin();
		auto rotationsAxisAngle	= m_rotationsAxisAngle[i].begin();
		auto rotationsQuaternion= m_rotationsQuaternion[i].begin();
		auto dataFile		= m_dataFiles[i].begin();
		std::set<int> flags	= std::set<int>{int(i)};

		trMatrixFinal.set_identity();

		for ( std::vector<displacementsSource>::iterator j = sources; j != m_dataSources[i].end(); j++ ) {
			switch ( *j ) {
				case 0:
					//(dataFile.first++)->second;
					break;//continue;
				case 1:
					trMatrix = *(transformations++);
					break;//continue;
				case 2:
					vectorValues = *(translations++);
					trMatrix.set_translation( vectorValues.data() );
					break;//continue;
				case 3:
					vectorValues = *(rotationsEulerAngles++);
					trMatrix.set_rotation( vectorValues[0], vectorValues[1], vectorValues[2], vectorValues[5], vectorValues[4], vectorValues[3] );
					break;//continue;
				case 4:
					vectorValues = *(rotationsAxisAngle++);
					trMatrix.set_rotation( vectorValues[0], vectorValues[1], vectorValues[2], vectorValues[3] );
					break;//continue;
				case 5:
					vectorValues = *(rotationsQuaternion++);
					trMatrix.set_rotation_quaternion( vectorValues[0], vectorValues[1], vectorValues[2], vectorValues[3] );
					break;//continue;
			}

//			std::cout << std::endl << trMatrix << std::endl;
			trMatrixFinal = trMatrix*trMatrixFinal;
		}
//		std::cout << std::endl << "final matrix for node set " << i << std::endl << trMatrixFinal << std::endl;
		m_meshBoundaryAdapter->insert_perturbationData( flags, trMatrixFinal.data() );
	}

	for ( auto& data: m_weightingData ) {
		m_meshBoundaryAdapter->set_weightingFunctionSpans( data.first, data.second[0] );
		m_meshBoundaryAdapter->set_weightingFunctionCutOffLimits( data.first, data.second[1] );
		m_meshBoundaryAdapter->set_weightingFunctionInterpolationFactors( data.first, data.second[2] );
	}


//	std::vector<T>& weights = m_meshBoundaryAdapter->get_weightingFunctionSpans( );
//	std::vector<T>& limits = m_meshBoundaryAdapter->get_weightingFunctionCutOffLimits( );
//	std::cout << limits<< std::endl;
}


template<class T>
int app_meshAdapter<T>::set_meshAdapter ( meshAdapter_t* adapter )
{
	if ( m_memAllocatedMembers[2] ) delete m_meshBoundaryAdapter; m_meshBoundaryAdapter = adapter; m_memAllocatedMembers[2] = 0;

	set_option_solutionMethod	( m_meshBoundaryAdapter->get_solutionMethod() 		);
	set_option_weightingFunctionType( m_meshBoundaryAdapter->get_weightingFunctionType() 	);
	set_option_computingMode	( m_meshBoundaryAdapter->get_computingMode() 		);
	set_option_linearSolverType 	( m_meshBoundaryAdapter->get_linearSolverType() 	);
	set_polynomialDegree		( m_meshBoundaryAdapter->get_polynomialDegree() 	);
	set_weightingFunctionSpan	( m_meshBoundaryAdapter->get_weightingFunctionSpan() 	);
	set_cutOffMultiplier		( m_meshBoundaryAdapter->get_cutOffMultiplier() 	);
	set_interpolationFactor		( m_meshBoundaryAdapter->get_interpolationFactor() 	);

	set_meshAdapter();
}


template<class T>
int app_meshAdapter<T>::adaptMesh ()
{
	m_meshBoundaryAdapter->adapt_nodes();
	m_adaptedMesh->find_metrics();
	m_meshBoundaryAdapter->find_badElements();
}


template<class T>
int app_meshAdapter<T>::write_adaptedNodesFile ( )
{
	std::ofstream newFileStreamOut; newFileStreamOut.exceptions(std::ios_base::badbit);
	newFileStreamOut.open(m_fileOfAdaptedMeshNodes);
/*
	size_t numOfNodes = m_adaptedMesh->get_numOfNodes();
	newFileStreamOut << numOfNodes << std::endl;
	newFileStreamOut << *(m_adaptedMesh->get_nodeFlags()) << std::endl;
	m_adaptedMesh->get_component(newFileStreamOut, size_t(0));
	m_adaptedMesh->get_component(newFileStreamOut, size_t(1));
*/	newFileStreamOut.close();
}


template<class T>
int app_meshAdapter<T>::write_invalidElementsFile ( )
{
	std::ofstream newFileStreamOut; newFileStreamOut.exceptions(std::ios_base::badbit);

	newFileStreamOut.open(m_projectName.substr(0, m_projectName.find('.')) +"_invalid_tri.dat");
	newFileStreamOut << "List of Invalid triangles :\n";
	newFileStreamOut << "Columns of: triangle ID, Old Area, New Area \n";

/*	std::vector<size_t>& badElements = m_meshBoundaryAdapter->get_badElements();
	ElementsArea<T, size_t, 3>& elementsAreaInitial = *(m_initialMesh->get_elementsArea());
	ElementsArea<T, size_t, 3>& elementsAreaAdapted = *(m_adaptedMesh->get_elementsArea());

	for (int elemIndex = 0; elemIndex < badElements.size(); elemIndex++)
		newFileStreamOut << badElements[elemIndex] << "\t\t" << elementsAreaInitial[badElements[elemIndex]] << "\t\t" << elementsAreaAdapted[badElements[elemIndex]] << std::endl;


	newFileStreamOut << m_meshBoundaryAdapter->get_badElements().size() << " Invalid Triangles, out of " << m_initialMesh->get_numOfElements() << std::endl;
*/	newFileStreamOut.close();

}



template<class T>
int app_meshAdapter<T>::write_boundaryNodesDisplacements (std::ostream& out)
{ m_meshBoundaryAdapter->write_nodeDisplacements( out ); }




template<class T>
int app_meshAdapter<T>::assure_results ( std::string referenceResults )
{

	std::ifstream newFileStream; newFileStream.exceptions(std::ios_base::badbit); newFileStream.open(referenceResults);

	int numOfTotalNodes; newFileStream >> numOfTotalNodes;
	T temp;
	for ( int node=0; node<numOfTotalNodes;node++ ) newFileStream >> temp;

	std::vector<T> x(numOfTotalNodes), y(numOfTotalNodes);
	for ( int node=0; node<numOfTotalNodes;node++ ) newFileStream >> x[node];
	for ( int node=0; node<numOfTotalNodes;node++ ) newFileStream >> y[node];

	newFileStream.close();


	T *x_n = m_adaptedMesh->get_component(0), *y_n = m_adaptedMesh->get_component(1);

	T norm_x(0), norm_y(0);
	for ( int node=0; node < m_initialMesh->get_numOfNodes(); node++ )
	{
		norm_x += pow( (x_n[node] - x[node]), 2);
		norm_y += pow( (y_n[node] - y[node]), 2);
	}
	
	norm_x = sqrt(norm_x); norm_y = sqrt(norm_y);
	
	std::cout << "L2 norm of the error in x displacements: " << norm_x << std::endl << "L2 norm of the error in y displacements: " << norm_y << std::endl; 

}

