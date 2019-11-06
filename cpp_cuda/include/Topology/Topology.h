#ifndef TOPOLOGYHEADERFILE
#define TOPOLOGYHEADERFILE

#include<stdexcept>
#include<exception>
#include<iostream>
#include<fstream>
#include<sstream>
#include<set>
#include<map>
#include"Nodes.h"
#include"ContainerBase.h"
#include<algorithm>
#include<functional>
#include"myAlgorithms.h"

template<class T>
class NodeFlags: public ContainerBase< T, std::vector<T> >
{

typedef ContainerBase< T, std::vector<T> > ContainerBase_t;

public:
	using ContainerBase_t::ContainerBase;
	using ContainerBase_t::resize;
	int set_numOfNodes(size_t size) { resize(size); }
};


template<class T, size_t nodesPerConn>
class Connectivities: public ContainerBase< T, std::vector<T> >
{

typedef ContainerBase< T, std::vector<T> > ContainerBase_t;
	
protected:
	size_t								m_numOfConnectivities;
	std::set< std::pair<size_t, size_t >, unorderLess<size_t> > 	m_uniqueLinks;
	std::vector< std::set<size_t> > 				m_nodeLinks;
	size_t								m_numOfNodes = 0;

public:
	Connectivities ( size_t numOfConnectivities = 1 ) { set_numOfConnectivities (numOfConnectivities); }
	Connectivities ( size_t numOfConnectivities, T* data )
	{ set_numOfConnectivities (numOfConnectivities); set_connectivities(data); }
	Connectivities ( Connectivities& referenceConnectivities ): ContainerBase_t( referenceConnectivities )
	{ set_numOfConnectivities(referenceConnectivities.get_numOfConnectivities() ); }
	~Connectivities () {}


	size_t		get_numOfConnectivities 	() { return m_numOfConnectivities; }
	size_t		get_numOfNodesPerConnectivity 	() { return nodesPerConn; }
	int		set_numOfConnectivities 	(size_t newSize )  { return resize( newSize ); }
	T*		operator[] 			(size_t connIndex) { return this->m_data+connIndex*nodesPerConn; }
	T*		at				(size_t connIndex) { return this->m_data+connIndex*nodesPerConn; }

	int		resize(size_t numOfConnectivities )  { m_numOfConnectivities = numOfConnectivities; return ContainerBase_t::resize( numOfConnectivities*nodesPerConn ); }
	int		reserve(size_t numOfConnectivities ) { m_numOfConnectivities = numOfConnectivities; return ContainerBase_t::reserve( numOfConnectivities*nodesPerConn );  }


	int		set_connectivities ( T* data ) { return this->set_content ( data );  }
	int		set_connectivities ( std::istream& streamIn ) { return this->set_content ( streamIn );  }
	T*		get_connectivities (  ) { return this->m_data; }
	int		get_connectivities ( std::ostream& streamOut ) { this->get_content ( streamOut ); }
	std::set< std::pair<size_t, size_t >, unorderLess<size_t> >& 	get_uniqueLinks() { if ( m_uniqueLinks.empty() ) find_uniqueLinks(); return m_uniqueLinks; }
	std::vector< std::set<size_t> >& 				get_nodeLinks() { if ( m_nodeLinks.empty() ) find_nodeLinks(); return m_nodeLinks; }
	size_t								get_numOfNodes() { if (!m_numOfNodes) m_numOfNodes = (*std::max_element( this->begin(), this->end() ) ) +1; return m_numOfNodes; }


	void find_uniqueLinks ( ) {
		if ( !m_uniqueLinks.empty() ) m_uniqueLinks.clear();
		T* nodes;
		for ( size_t conn = 0; conn < m_numOfConnectivities; conn++) {
			nodes = at( conn );
			for ( size_t node = 0; node < nodesPerConn-1; node++ )
				m_uniqueLinks.insert( std::make_pair( nodes[ node ], nodes[ node+1 ] ) );
			m_uniqueLinks.insert( std::make_pair( nodes[ nodesPerConn-1 ], nodes[0] ) );
		}
	}


	void find_nodeLinks ( ) {
		if ( m_uniqueLinks.empty() ) find_uniqueLinks();

		if ( !m_nodeLinks.empty() ) m_nodeLinks.clear();
		m_nodeLinks.resize( get_numOfNodes() );
		for ( auto& el: m_uniqueLinks ) m_nodeLinks.at(el.first).insert( el.second );
		for ( auto& el: m_uniqueLinks ) m_nodeLinks.at(el.second).insert( el.first );
	}

};




/*
class ConnectivitiesN
{

protected:
	std::vector<std::vector<size_t> > container;

public:
	ConnectivitiesN( size_t numOfConnectivities ) { container.resize(numOfConnectivities); }
	size_t insert_node(size_t connectivityIndex, size_t node) {  container.at(connectivityIndex).push_back( node ); }
	size_t set_nodes(size_t connectivityIndex, size_t* nodes, size_t numOfNodes) { for (size_t i=0; i<numOfNodes; i++) container.at(connectivityIndex).push_back( nodes[i] ); }


	std::vector<size_t>& get_nodes( size_t connectivityIndex ) { return container.at(connectivityIndex); }

	size_t* operator [] ( size_t index ) { container.at(index).data(); }
};
*/

template<class nodalComponentType = float, size_t numOfConnectivityNodes = 3, class flagType = int>
class Topology: public Nodes<nodalComponentType>
{

public:

typedef Nodes<nodalComponentType> 					nodes_t;
typedef Connectivities<size_t, numOfConnectivityNodes> 			connectivities_t;
typedef Topology<nodalComponentType, numOfConnectivityNodes, flagType> 	topology_t;
typedef NodeFlags< flagType > 						nodeFlags_t;
typedef std::set< flagType >						flagValues_t;
typedef std::map< flagType, std::vector<size_t> >			flagsMap_t;
typedef std::map<size_t, size_t>					nodesMap_t;


private:
	std::vector< bool >	m_memAllocatedMembers;

protected:
	connectivities_t*	m_connectivities;
	nodeFlags_t*		m_nodeFlags;
	flagValues_t		m_flagValues;
	flagValues_t		m_interiorFlagValues;
	flagValues_t		m_boundaryFlagValues;
	flagsMap_t		m_flagToNodesMap;
	nodesMap_t		m_boundaryNodes;
	nodesMap_t		m_interiorNodes;


public:
	Topology ( 	size_t 			numOfComponents = 2, 
			size_t 			numOfNodes = 0, 
			nodalComponentType** 	nodesData=0, 
			connectivities_t* 	connectivities=0, 
			nodeFlags_t* 		nodeFlags=0, 
			flagValues_t* 		interiorFlagValues = 0, 
			flagsMap_t* 		flagToNodesMap = 0 );
	Topology ( topology_t* referenceTopology );
	Topology ( topology_t& referenceTopology );
	~Topology ();


virtual	int			find_interiorNodes	();
virtual	int			find_boundaryNodes	();

	int			find_complementaryNodes ( size_t* subset_begin, size_t* subset_end, size_t* complementarySet_begin );
	flagValues_t	static	find_flagValues		( nodeFlags_t* nodeFlags );
	flagsMap_t 	static	find_flagToNodesMap 	( nodeFlags_t* nodeFlags, flagValues_t* flagValues=0 );
	flagValues_t	static  find_complementarySet	( flagValues_t& fullSet, flagValues_t& subset );
	int		static  find_complementarySet	( flagType* set_begin, flagType* set_end, flagType* subset_begin, flagType* subset_end, flagType* complementary );


	int 			set_numOfNodes 		( size_t );
	int 			set_numOfConnectivities ( size_t numOfConnectivities );
	int 			set_nodeFlags 		( std::istream& inStream );
	int 			set_nodeFlags 		( nodeFlags_t* data );
	int			set_nodeFlagWithCriterionOnComponentValue( size_t componentIndex, unary_criterion<nodalComponentType>& fun, const flagType& newFlag );
//	int			set_nodeFlagWithCriterionOnConnectivities( std::unary_function<nodalComponentType, bool> fun, const flagType& newFlag );
	int 			set_connectivities 	( std::istream& inStream );
	int 			set_connectivities 	( std::istream& inStream, int nodeNumberingShift );
	int 			set_connectivities 	( connectivities_t* data );
	int			set_interiorFlagValues 	( flagValues_t& newFlag );
	int			set_interiorFlagValues 	( flagType* newFlag, int numOfFlags );
	int			set_boundaryFlagValues 	( flagValues_t& newFlag );
	int			set_boundaryFlagValues 	( flagType* newFlag, int numOfFlags );
	int			set_flagToNodesMap	( flagsMap_t& data );


	std::string 		get_typename		() { return "Topology";}
	nodeFlags_t*		get_nodeFlags 		();
	connectivities_t*	get_connectivities	();
	size_t			get_numOfConnectivities ();
	flagValues_t&		get_flagValues		();
	flagValues_t&		get_interiorFlagValues 	();
	flagValues_t&		get_boundaryFlagValues	();
	flagsMap_t& 		get_flagToNodesMap	();
	int			get_nodesWithFlag 	( flagValues_t& flags, std::vector<size_t>* containerToPopulate );
	int			get_nodesWithoutFlag 	( flagValues_t& flags, std::vector<size_t>* containerToPopulate );
	size_t			get_numOfNodesWithFlag 	( flagValues_t& flags );
	size_t			get_numOfNodesWithoutFlag ( flagValues_t& flags );
	nodesMap_t&		get_boundaryNodes	( );
	nodesMap_t&		get_interiorNodes	( );
	size_t			get_numOfBoundaryNodes	( );
	size_t			get_numOfInteriorNodes	( );

};



template<class T, size_t T2, class T3>
Topology<T, T2, T3>::Topology (	size_t 			numOfComponents, 
				size_t 			numOfNodes, 
				T** 			nodesData, 
				connectivities_t* 	connectivities, 
				nodeFlags_t* 		nodeFlags, 
				flagValues_t* 		interiorFlagValues, 
				flagsMap_t* 		flagToNodesMap ):
Nodes<T>(numOfComponents, numOfNodes, nodesData), m_connectivities(connectivities), m_nodeFlags(nodeFlags), m_memAllocatedMembers( std::vector< bool >(3, 0) )
{
	if (!nodeFlags) 		{ m_nodeFlags 		= new nodeFlags_t; m_memAllocatedMembers[0]=1;}
	else set_nodeFlags( nodeFlags );
	if (!connectivities)		{ m_connectivities 	= new connectivities_t; m_memAllocatedMembers[1]=1; }
	else set_connectivities(connectivities);
	if ( interiorFlagValues )	m_interiorFlagValues 	= *interiorFlagValues;
	if ( flagToNodesMap ) 		m_flagToNodesMap 	= *flagToNodesMap;
}


template<class T, size_t T2, class T3>
Topology<T, T2, T3>::Topology ( Topology<T, T2, T3> *referenceTopology ): Nodes<T>(referenceTopology), m_memAllocatedMembers( std::vector< bool >(3, 0) )
{
	set_nodeFlags		( referenceTopology->get_nodeFlags() );
	set_connectivities 	( referenceTopology->get_connectivities() );
	m_flagValues 		= referenceTopology->get_flagValues();
	m_interiorFlagValues 	= referenceTopology->get_interiorFlagValues();
	m_boundaryFlagValues	= referenceTopology->get_boundaryFlagValues();
	m_flagToNodesMap	= referenceTopology->get_flagToNodesMap();
	m_boundaryNodes		= referenceTopology->get_boundaryNodes();
	m_interiorNodes		= referenceTopology->get_interiorNodes();
}


template<class T, size_t T2, class T3>
Topology<T, T2, T3>::Topology ( Topology<T, T2, T3>& referenceTopology ): Nodes<T>(referenceTopology), m_memAllocatedMembers(std::vector<bool>(3,1))
{
	m_connectivities 	= new connectivities_t( *(referenceTopology.get_connectivities()) );
	m_nodeFlags		= new nodeFlags_t ( *(referenceTopology.get_nodeFlags()) );
	m_flagValues 		= referenceTopology.get_flagValues();
	m_interiorFlagValues 	= referenceTopology.get_interiorFlagValues();
	m_boundaryFlagValues	= referenceTopology.get_boundaryFlagValues();
	m_flagToNodesMap	= referenceTopology.get_flagToNodesMap();
	m_boundaryNodes		= referenceTopology.get_boundaryNodes();
	m_interiorNodes		= referenceTopology.get_interiorNodes();
}


template<class T, size_t T2, class T3>
Topology<T, T2, T3>::~Topology()
{
	if ( m_memAllocatedMembers[0] ) delete m_nodeFlags;
	if ( m_memAllocatedMembers[1] ) delete m_connectivities;
}




template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::flagValues_t Topology<T, T2, T3>::find_flagValues( nodeFlags_t* nodeFlags  )
{
	flagValues_t flagValues;
	for ( T3 flag: *nodeFlags ) flagValues.insert ( flag );
	return flagValues;
}


template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::flagsMap_t Topology<T, T2, T3>::find_flagToNodesMap ( nodeFlags_t* nodeFlags, flagValues_t* flagValues )
{
	flagsMap_t mapValues; bool deleteFlag(0);
	if  ( !flagValues ) { flagValues = new flagValues_t( find_flagValues( nodeFlags ) ); deleteFlag = 1; }

	for ( auto& flagv : flagValues[0] ) mapValues[ flagv ] = std::vector<size_t>();

	for ( size_t index = 0; index<nodeFlags->size(); index++ ) mapValues[ nodeFlags[0][index] ].push_back(index);

	if ( deleteFlag ) delete flagValues;

	return mapValues;
}

template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::flagValues_t Topology<T, T2, T3>::find_complementarySet ( flagValues_t& fullSet, flagValues_t& subset )
{
	flagValues_t set_c ( fullSet );
	for ( auto& element: subset ) set_c.erase( element );

	return set_c;
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::find_complementarySet ( T3* set_begin, T3* set_end, T3* subset_begin, T3* subset_end, T3* complementary )
{
	flagValues_t set_c ( set_begin, set_end ), subset( subset_begin, subset_end );
	for ( auto& element: subset ) set_c.erase( element );

	std::copy( set_c.begin(), set_c.end(), complementary );
	return 0;
}


template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::find_complementaryNodes ( size_t* subset_begin, size_t* subset_end, size_t* complementary )
{
	std::set<size_t> subset ( subset_begin, subset_end );

	std::set<size_t> fullSet; std::set<size_t>::iterator it(fullSet.begin());
	for (size_t node = 0; node < this->get_numOfNodes(); node++) it = fullSet.insert( it, node );

	for ( auto& element: subset ) fullSet.erase( element );

	std::copy( fullSet.begin(), fullSet.end(), complementary );

	return 0;
}


template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::find_interiorNodes ( )
{
	m_interiorNodes.clear(); size_t index(0);
	for ( auto& flag: m_interiorFlagValues ) for ( auto& node: m_flagToNodesMap[flag] ) m_interiorNodes[ node ] = index++;
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::find_boundaryNodes ( )
{
	m_boundaryNodes.clear(); size_t index(0);
	for ( auto& flag: m_boundaryFlagValues ) for ( auto& node: m_flagToNodesMap[flag] ) m_boundaryNodes[ node ] = index++;
}



template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::nodeFlags_t* Topology<T, T2, T3>::get_nodeFlags () { return m_nodeFlags; }

template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::connectivities_t*	Topology<T, T2, T3>::get_connectivities	() { return m_connectivities; }

template<class T, size_t T2, class T3>
size_t Topology<T, T2, T3>::get_numOfConnectivities () { return m_connectivities->get_numOfConnectivities(); }

template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::flagValues_t&	Topology<T, T2, T3>::get_flagValues () { return m_flagValues; }

template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::flagValues_t&	Topology<T, T2, T3>::get_interiorFlagValues	() { return m_interiorFlagValues; }

template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::flagValues_t&	Topology<T, T2, T3>::get_boundaryFlagValues	() { return m_boundaryFlagValues; }

template<class T, size_t T2, class T3>
typename Topology<T, T2, T3>::flagsMap_t& Topology<T, T2, T3>::get_flagToNodesMap	() { return m_flagToNodesMap; }

template<class T, size_t T2, class T3>
size_t	Topology<T, T2, T3>::get_numOfNodesWithFlag ( flagValues_t& flags )
{
	size_t numOfNodes(0);
	for (auto& flag: flags ) numOfNodes += m_flagToNodesMap[ flag ].size();
	return numOfNodes;
}

template<class T, size_t T2, class T3>
size_t	Topology<T, T2, T3>::get_numOfNodesWithoutFlag ( flagValues_t& flags )
{
	flagValues_t flags_c = find_complementarySet ( m_flagValues,  flags);
	return get_numOfNodesWithFlag ( flags_c );
}


template<class T, size_t T2, class T3>
int	Topology<T, T2, T3>::get_nodesWithFlag ( flagValues_t& flags, std::vector<size_t>* containerToPopulate )
{
	for ( auto& flag: flags ) {
		containerToPopulate->insert( containerToPopulate->end(), m_flagToNodesMap[ flag ].begin(), m_flagToNodesMap[ flag ].end() );
	}
}

template<class T, size_t T2, class T3>
int	Topology<T, T2, T3>::get_nodesWithoutFlag 	( flagValues_t& flags, std::vector<size_t>* containerToPopulate )
{
	for ( auto& flag: find_complementarySet ( m_flagValues,  flags)  ) containerToPopulate->insert( containerToPopulate->end(), m_flagToNodesMap[ flag ].begin(), m_flagToNodesMap[ flag ].end() );
}

template<class T, size_t T2, class T3>
std::map<size_t, size_t>& Topology<T, T2, T3>::get_boundaryNodes ( )
{
	return m_boundaryNodes;
}

template<class T, size_t T2, class T3>
std::map<size_t, size_t>& Topology<T, T2, T3>::get_interiorNodes ( )
{
	return m_interiorNodes;
}

template<class T, size_t T2, class T3>
size_t Topology<T, T2, T3>::get_numOfBoundaryNodes ()
{
	return m_boundaryNodes.size();
}

template<class T, size_t T2, class T3>
size_t Topology<T, T2, T3>::get_numOfInteriorNodes	()
{
	return m_interiorNodes.size();
}




template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_numOfNodes(size_t newNumOfPoints) {
	Nodes<T>::set_numOfNodes(newNumOfPoints);
	m_nodeFlags->set_numOfNodes(newNumOfPoints);
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_numOfConnectivities(size_t numOfConnectivities) {
	m_connectivities->set_numOfConnectivities(numOfConnectivities); return 1; 
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_nodeFlags(std::istream& inStream) {
	inStream >> *(m_nodeFlags);
	m_flagValues = find_flagValues( m_nodeFlags );
	m_flagToNodesMap = find_flagToNodesMap ( m_nodeFlags, &m_flagValues );
	m_memAllocatedMembers[0] = 1; return 1;
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_nodeFlags( nodeFlags_t* data ) {
	m_nodeFlags = data;
	m_flagValues = find_flagValues( m_nodeFlags );
	m_flagToNodesMap = find_flagToNodesMap ( m_nodeFlags, &m_flagValues );
	m_memAllocatedMembers[0] = 0; return 1;
}


template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_nodeFlagWithCriterionOnComponentValue( size_t componentIndex, unary_criterion<T>& fun, const T3& newFlag )
{
	T* values = this->get_component( componentIndex );
	for ( size_t i = 0; i< this->get_numOfNodes(); i++ )
		if ( fun( values[i] ) ) m_nodeFlags[0][i] = newFlag;

	m_flagValues = find_flagValues(m_nodeFlags);
	m_flagToNodesMap = find_flagToNodesMap ( m_nodeFlags, &m_flagValues );
	return 1;
}


template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_connectivities(std::istream& inStream)
{ inStream >> *(m_connectivities); m_memAllocatedMembers[1] = 1; return 1; }

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_connectivities(std::istream& inStream, int nodeNumberingShift ) {
	inStream >> *(m_connectivities);
	if ( nodeNumberingShift ) for ( size_t &node : *(m_connectivities) ) node += nodeNumberingShift;
	m_memAllocatedMembers[1] = 1; return 1;
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_connectivities ( connectivities_t* conn  )
{ m_connectivities = conn; m_memAllocatedMembers[1] = 0; return 1; }

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_interiorFlagValues ( flagValues_t& newFlag ) 
{
	m_interiorFlagValues = newFlag;
	m_boundaryFlagValues = find_complementarySet( m_flagValues, newFlag );
	find_interiorNodes(); find_boundaryNodes();
	return 1;
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_interiorFlagValues ( T3* newFlag, int numOfFlags ) 
{
	m_interiorFlagValues.clear();
	m_interiorFlagValues.insert( newFlag, newFlag+numOfFlags);
	m_boundaryFlagValues = find_complementarySet( m_flagValues, m_interiorFlagValues );
	find_interiorNodes(); find_boundaryNodes();
	return 1;
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_boundaryFlagValues ( flagValues_t& newFlag ) 
{
	m_boundaryFlagValues = newFlag;
	m_interiorFlagValues = find_complementarySet( m_flagValues, m_boundaryFlagValues );
	find_interiorNodes(); find_boundaryNodes();
	return 1; 
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_boundaryFlagValues ( T3* newFlag, int numOfFlags ) 
{
	m_boundaryFlagValues.clear();
	m_boundaryFlagValues.insert( newFlag, newFlag+numOfFlags);
	m_interiorFlagValues = find_complementarySet( m_flagValues, m_boundaryFlagValues );
	find_interiorNodes(); find_boundaryNodes();
	return 1;
}

template<class T, size_t T2, class T3>
int Topology<T, T2, T3>::set_flagToNodesMap ( flagsMap_t& data) 
{ m_flagToNodesMap = data; }


#endif
