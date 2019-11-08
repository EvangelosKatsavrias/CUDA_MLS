//   CUDA_MLS Framework
//
//   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
//
//   This file is part of the CUDA_MLS Framework.
//
//   CUDA_MLS Framework is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License version 3 as published by
//   the Free Software Foundation.
//
//   CUDA_MLS Framework is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
//
//   Contact Info:
//   Evangelos D. Katsavrias
//   email/skype: vageng@gmail.com
// -----------------------------------------------------------------------

#ifndef TOPOLOGYADAPTERSHEADERFILE
#define TOPOLOGYADAPTERSHEADERFILE

#include<stdexcept>
#include<exception>
#include<iostream>
#include"Topology.h"


template<class nodalComponentType = float, size_t numOfConnectivityNodes = 3, class flagDataType = int>
class TopologyAdapter
{

public:
typedef Topology<nodalComponentType, numOfConnectivityNodes, flagDataType> 	topology_t;
typedef ContainerBase<size_t, std::vector<size_t> > 				perturbedNodes_t;
typedef ContainerBase<size_t, std::vector<size_t> > 				adaptedNodes_t;
typedef Nodes<nodalComponentType> 						perturbations_t;
typedef std::vector<size_t>							badConnectivities_t;
using flagValues_t = typename topology_t::flagValues_t;
typedef std::map< flagDataType, size_t >					perturbedFlag2NodesIndex_t;

protected:
	topology_t*			m_initialTopology = 0;
	topology_t*			m_adaptedTopology = 0;
	perturbedNodes_t*		m_perturbedNodes = 0;
	adaptedNodes_t*			m_adaptedNodes = 0;
	perturbations_t*		m_perturbations = 0;
	badConnectivities_t 		m_badConnectivities;
	float				m_adaptationTime;
	bool				m_adaptationCompletionFlag;


	flagValues_t			m_perturbedNodeFlags;
	flagValues_t			m_adaptedNodeFlags;
	perturbedFlag2NodesIndex_t	m_perturbedFlag2NodesIndex; // a map of perturbed node flags to the first perturbed node in the purturbedNodes container

	std::vector<bool>		m_memAllocatedMembers; // 0)perturbations, 1)adaptedNodes, 2) fixedNodes


public:
	TopologyAdapter( topology_t *initialTopology, topology_t *adaptedTopology = 0, perturbations_t* data=0, perturbedNodes_t* nodes=0, flagValues_t perturbedNodeFlags = flagValues_t() );
	~TopologyAdapter();


//	int			virtual find_adaptedNodes(); // perturbedNodes should be provided before or the node flags to be adapted
//	int			virtual find_adaptedNodeFlags(); // perturbedNodeFlags should be provided before
	int			virtual	adapt_nodes ();
	int			virtual	find_badConnectivities ();


	int			virtual	set_initialTopology 		( topology_t* topo) { m_initialTopology = topo; }
	int			virtual	set_adaptedTopology 		( topology_t* topo) { m_adaptedTopology = topo; }
	int 			virtual	set_perturbations		( perturbations_t* data );
	int 			virtual	set_perturbations		( );
	int			virtual	set_adaptedNodes		( adaptedNodes_t* adaptedNodes );
	int			virtual	set_allNodes			( );
	int			virtual	set_perturbedNodes		( perturbedNodes_t* perturbedNodes );
	int			virtual	set_adaptedNodeFlags		( flagValues_t& adaptedNodeFlags );
	int 			virtual set_adaptedNodeFlags 		( flagDataType* adaptedNodeFlags_begin, flagDataType* adaptedNodeFlags_end );
	int			virtual	set_perturbedNodeFlags		( flagValues_t& adaptedNodeFlags );
	int 			virtual set_perturbedNodeFlags 		( flagDataType* perturbedNodeFlags_begin, flagDataType* perturbedNodeFlags_end );
	int			virtual	set_adaptationTime 		( float time ) { m_adaptationTime = time; }
	int			virtual	set_adaptationCompletionFlag 	( bool flag ) { m_adaptationCompletionFlag=flag; }


	int			virtual	insert_perturbations		( nodalComponentType** perturbations, size_t numOfPerturbedNodes );
	int			virtual insert_perturbedNodesFlag	( flagDataType flag );
	int			virtual insert_adaptedNodesFlag		( flagDataType flag );
	int			virtual insert_perturbedNodes		( size_t* nodeIndex_begin, size_t* nodeIndex_end );
	int			virtual insert_adaptedNodes		( size_t* nodeIndex_begin, size_t* nodeIndex_end );

	
	int			virtual insert_perturbationData		( size_t* perturbedNodes, size_t numberOfPerturbedNodes, nodalComponentType** perturbations );
	int			virtual insert_perturbationData		( size_t* perturbedNodes, size_t numberOfPerturbedNodes, nodalComponentType* transformationMatrix );
	int			virtual insert_perturbationData		( flagValues_t flagsOfPerturbedNodes, nodalComponentType* transformationMatrix );


	int 			virtual	set_perturbationData		( size_t* perturbedNodes, size_t numberOfPerturbedNodes, nodalComponentType** perturbations, size_t numberOfComponents );
	int			virtual set_perturbationData		( size_t* perturbedNodes, size_t numberOfPerturbedNodes, nodalComponentType* transformationMatrix );
	int			virtual set_perturbationData		( flagValues_t flagsOfPerturbedNodes, nodalComponentType* transformationMatrix );



virtual	std::string 			get_typename 			() { return "Topology Adapter"; }
virtual	topology_t*			get_initialTopology 		() { return m_initialTopology; }
virtual	topology_t*			get_adaptedTopology 		() { return m_adaptedTopology; }
virtual perturbations_t* 		get_perturbations		() { return m_perturbations; }
virtual adaptedNodes_t*			get_adaptedNodes		() { return m_adaptedNodes; }
virtual perturbedNodes_t*		get_perturbedNodes		() { return m_perturbedNodes; }
	size_t				get_numOfPerturbedNodes		() { return m_perturbedNodes->size(); }
	size_t				get_numOfAdaptedNodes		() { return m_adaptedNodes->size(); }
	flagValues_t& 			get_adaptedNodeFlags		() { return m_adaptedNodeFlags; }
	flagValues_t&			get_perturbedNodeFlags		() { return m_perturbedNodeFlags; }
	size_t				get_numOfAdaptedNodeFlags	() { return m_adaptedNodeFlags.size(); }
	size_t				get_numOfPerturbedNodeFlags	() { return m_perturbedNodeFlags.size(); } 
	float				get_adaptationTime 		() { return m_adaptationTime; }
	bool				get_adaptationCompletionFlag 	() { return m_adaptationCompletionFlag; }
	std::vector<size_t>& 		get_badConnectivities 		() { return m_badConnectivities; }
	int 			virtual	get_adaptationParameters 	( std::ostream& out) { }
	int 			virtual	get_perturbations 		( std::string file ) {  }
	int 			virtual	get_perturbations 		( std::ostream& out ) {  }
	perturbedFlag2NodesIndex_t	get_perturbedFlag2NodesIndex	() { return m_perturbedFlag2NodesIndex; }

};



template<class T, size_t T2, class T3>
TopologyAdapter<T, T2, T3>::TopologyAdapter( 	topology_t*		initialTopology,
						topology_t*		adaptedTopology,
						perturbations_t* 	data,
						perturbedNodes_t* 	nodes,
						flagValues_t 		perturbedNodeFlags ):
	m_initialTopology(initialTopology), 
	m_adaptedTopology(adaptedTopology), 
	m_perturbations(data), 
	m_perturbedNodes(nodes), 
	m_adaptationTime(0), 
	m_adaptationCompletionFlag(0), 
	m_perturbedNodeFlags(perturbedNodeFlags), 
	m_memAllocatedMembers( 3, 0 ) 
{
	if ( !data ) 	set_perturbations();
	if ( !nodes ) 	set_allNodes();
	else set_perturbedNodes(nodes);
	if ( !perturbedNodeFlags.empty() && !nodes ) set_perturbedNodeFlags( perturbedNodeFlags );
}


template<class T, size_t T2, class T3>
TopologyAdapter<T, T2, T3>::~TopologyAdapter()
{
	if (m_memAllocatedMembers[0]) delete m_perturbations;
	if (m_memAllocatedMembers[1]) delete m_adaptedNodes;
	if (m_memAllocatedMembers[2]) delete m_perturbedNodes;
}

template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::adapt_nodes() { }


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::find_badConnectivities() { }


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbations( )
{
	if ( !m_memAllocatedMembers[0] ) {
		m_perturbations = new perturbations_t( m_initialTopology->get_numOfComponents(), size_t(0) );
		m_memAllocatedMembers[0] = 1; }
	m_perturbations->set_components();
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbations ( perturbations_t* data ) 
{ 
	if ( m_memAllocatedMembers[0] ) { delete m_perturbations; m_memAllocatedMembers[0] = 0; }
	m_perturbations = data;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_allNodes ( ) 
{
	if ( !m_memAllocatedMembers[2] ) {
		m_perturbedNodes = new perturbedNodes_t; 
		m_memAllocatedMembers[2] = 1; }
	m_perturbedNodes->clear();

	if ( !m_memAllocatedMembers[1] ) {
		m_adaptedNodes = new adaptedNodes_t; 
		m_memAllocatedMembers[1] = 1; }
	m_adaptedNodes->clear();
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbedNodes ( perturbedNodes_t* perturbedNodes ) 
{
	set_allNodes();
	if ( m_memAllocatedMembers[2] ) { delete m_perturbedNodes; m_memAllocatedMembers[2] = 0; }
	m_perturbedNodes = perturbedNodes; 

	m_adaptedNodes->resize( m_initialTopology->get_numOfNodes() -m_perturbedNodes->size() );
	m_initialTopology->find_complementaryNodes( m_perturbedNodes->begin(), m_perturbedNodes->end(), m_adaptedNodes->begin() );
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_adaptedNodes ( adaptedNodes_t* adaptedNodes ) 
{
	set_allNodes();
	if ( m_memAllocatedMembers[1] ) { delete m_adaptedNodes; m_memAllocatedMembers[1] = 0; }
	m_adaptedNodes = adaptedNodes; 

	m_perturbedNodes->resize( m_initialTopology->get_numOfNodes() -m_adaptedNodes->size() );
	m_initialTopology->find_complementaryNodes( m_adaptedNodes->begin(), m_adaptedNodes->end(), m_perturbedNodes->begin() );
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_adaptedNodeFlags ( T3* adaptedNodeFlags_begin, T3* adaptedNodeFlags_end )
{
	flagValues_t nodeFlags = m_initialTopology->get_flagValues();
	for ( T3* flag = adaptedNodeFlags_begin; flag != adaptedNodeFlags_end; flag++ ) if ( nodeFlags.find(*flag) == nodeFlags.end() ) std::cerr << "\nPlease, provide valid node flags.";
	m_adaptedNodeFlags.clear(); m_adaptedNodeFlags.insert( adaptedNodeFlags_begin, adaptedNodeFlags_end );

	m_perturbedNodeFlags = m_initialTopology->find_complementarySet( m_initialTopology->get_flagValues(), m_adaptedNodeFlags );

	m_perturbedFlag2NodesIndex.clear(); size_t counter(0);
	for ( auto& flag: m_perturbedNodeFlags ) {
		m_perturbedFlag2NodesIndex[ flag ] = counter;
		std::set<T3> flags = {flag};
		counter += m_initialTopology->get_numOfNodesWithFlag( flags );
	}

	set_allNodes(); std::vector<size_t> nodes;
	this->m_initialTopology->get_nodesWithFlag( m_adaptedNodeFlags, &nodes );
	m_adaptedNodes->set_content( nodes );

	nodes.clear();
	this->m_initialTopology->get_nodesWithFlag( m_perturbedNodeFlags, &nodes );
	m_perturbedNodes->set_content( nodes );

	return 1;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_adaptedNodeFlags ( flagValues_t& adaptedNodeFlags )
{
	flagValues_t nodeFlags = m_initialTopology->get_flagValues();
	for ( auto& flag: adaptedNodeFlags ) if ( nodeFlags.find(flag) == nodeFlags.end() ) std::cerr << "\nPlease, provide valid node flags.";
	m_adaptedNodeFlags = adaptedNodeFlags; 

	m_perturbedNodeFlags = m_initialTopology->find_complementarySet( m_initialTopology->get_flagValues(), m_adaptedNodeFlags );

	m_perturbedFlag2NodesIndex.clear(); size_t counter(0);
	for ( auto& flag: m_perturbedNodeFlags ) {
		m_perturbedFlag2NodesIndex[ flag ] = counter;
		std::set<T3> flags = {flag};
		counter += m_initialTopology->get_numOfNodesWithFlag( flags );
	}

	set_allNodes(); std::vector<size_t> nodes;
	this->m_initialTopology->get_nodesWithFlag( m_adaptedNodeFlags, &nodes );
	m_adaptedNodes->set_content( nodes );

	nodes.clear();
	this->m_initialTopology->get_nodesWithFlag( m_perturbedNodeFlags, &nodes );
	m_perturbedNodes->set_content( nodes );

	return 1;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbedNodeFlags ( T3* perturbedNodeFlags_begin, T3* perturbedNodeFlags_end )
{
	flagValues_t nodeFlags = m_initialTopology->get_flagValues();
	for ( T3* flag = perturbedNodeFlags_begin; flag != perturbedNodeFlags_end; flag++ ) if ( nodeFlags.find(*flag) == nodeFlags.end() ) std::cerr << "\nPlease, provide valid node flags.";
	m_perturbedNodeFlags.clear(); m_perturbedNodeFlags.insert( perturbedNodeFlags_begin, perturbedNodeFlags_end );

	m_perturbedFlag2NodesIndex.clear(); size_t counter(0);
	for ( auto& flag: m_perturbedNodeFlags ) {
		m_perturbedFlag2NodesIndex[ flag ] = counter;
		std::set<T3> flags = {flag};
		counter += m_initialTopology->get_numOfNodesWithFlag( flags );
	}


	m_adaptedNodeFlags = m_initialTopology->find_complementarySet( m_initialTopology->get_flagValues(), m_perturbedNodeFlags );

	set_allNodes(); std::vector<size_t> nodes;
	this->m_initialTopology->get_nodesWithFlag( m_adaptedNodeFlags, &nodes );
	m_adaptedNodes->set_content( nodes );

	nodes.clear();
	this->m_initialTopology->get_nodesWithFlag( m_perturbedNodeFlags, &nodes );
	m_perturbedNodes->set_content( nodes );

	return 1;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbedNodeFlags ( flagValues_t& fixedNodeFlags )
{
	flagValues_t nodeFlags = m_initialTopology->get_flagValues();
	for ( auto& flag: fixedNodeFlags ) if ( nodeFlags.find(flag) == nodeFlags.end() ) std::cerr << "\nPlease, provide valid node flags.";
	m_perturbedNodeFlags = fixedNodeFlags;

	m_perturbedFlag2NodesIndex.clear(); size_t counter(0);
	for ( auto& flag: m_perturbedNodeFlags ) {
		m_perturbedFlag2NodesIndex[ flag ] = counter;
//		counter += m_initialTopology->get_numOfNodesWithFlag( std::set<T3>(flag) );
	}

	m_adaptedNodeFlags = m_initialTopology->find_complementarySet( m_initialTopology->get_flagValues(), m_perturbedNodeFlags );

	set_allNodes(); std::vector<size_t> nodes;
	this->m_initialTopology->get_nodesWithFlag( m_adaptedNodeFlags, &nodes );
	m_adaptedNodes->set_content( nodes );

	nodes.clear();
	this->m_initialTopology->get_nodesWithFlag( m_perturbedNodeFlags, &nodes );
	m_perturbedNodes->set_content( nodes );

	return 1;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_perturbedNodes ( size_t* nodeIndex_begin, size_t* nodeIndex_end )
{
	int numOfExistingNodes = m_perturbedNodes->size();
	m_perturbedNodes->insert(m_perturbedNodes->end(), nodeIndex_begin, nodeIndex_end);

	return numOfExistingNodes;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_adaptedNodes ( size_t* nodeIndex_begin, size_t* nodeIndex_end )
{
	int numOfExistingNodes = m_adaptedNodes->size();
	m_adaptedNodes->insert(m_adaptedNodes->end(), nodeIndex_begin, nodeIndex_end);
	return numOfExistingNodes;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_perturbedNodesFlag ( T3 flag )
{
	int numOfExistingNodes = m_perturbedNodes->size();

	flagValues_t nodeFlags = m_initialTopology->get_flagValues();
	if ( nodeFlags.find(flag) == nodeFlags.end() ) std::cerr << "\nPlease, provide valid node flags.";
	m_perturbedNodeFlags.insert(flag);

	m_adaptedNodeFlags = m_initialTopology->find_complementarySet( m_initialTopology->get_flagValues(), m_perturbedNodeFlags );

	m_perturbedFlag2NodesIndex[ flag ] = m_perturbedNodes->size();

	std::vector<size_t> nodes; std::set<T3> flags{flag};
	this->m_initialTopology->get_nodesWithFlag( flags, &nodes );
	insert_perturbedNodes ( nodes.data(), nodes.data()+nodes.size() );


	nodes.clear();
	this->m_initialTopology->get_nodesWithFlag( m_adaptedNodeFlags, &nodes );
	m_adaptedNodes->set_content( nodes );

	return m_perturbedNodes->size() -numOfExistingNodes;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_adaptedNodesFlag ( T3 flag )
{
	int numOfExistingNodes = m_adaptedNodes->size();

	flagValues_t nodeFlags = m_initialTopology->get_flagValues();
	if ( nodeFlags.find(flag) == nodeFlags.end() ) std::cerr << "\nPlease, provide valid node flags.";
	m_adaptedNodeFlags.insert(flag);

	m_perturbedNodeFlags = m_initialTopology->find_complementarySet( m_initialTopology->get_flagValues(), m_adaptedNodeFlags );

	std::vector<size_t> nodes; std::set<T3> flags{flag};
	this->m_initialTopology->get_nodesWithFlag( flags, &nodes );
	insert_adaptedNodes ( nodes.data(), nodes.data()+nodes.size() );

	nodes.clear();
	this->m_initialTopology->get_nodesWithFlag( m_perturbedNodeFlags, &nodes );
	m_perturbedNodes->set_content( nodes );

	return m_adaptedNodes->size() -numOfExistingNodes;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_perturbations ( T** data, size_t numOfNodes )
{
	m_perturbations->insert_data( data, numOfNodes );
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_perturbationData ( size_t* perturbedNodes, size_t numberOfPerturbedNodes, T** perturbations )
{
	insert_perturbedNodes( perturbedNodes, perturbedNodes +numberOfPerturbedNodes );
	insert_perturbations( perturbations, numberOfPerturbedNodes );
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_perturbationData ( size_t* perturbedNodes, size_t numberOfPerturbedNodes, T* transformationMatrix )
{
	insert_perturbedNodes( perturbedNodes, perturbedNodes +numberOfPerturbedNodes );

	size_t numOfComponents = m_initialTopology->get_numOfComponents();
	T *initialNodeData[numOfComponents], *perturbationNodeData[numOfComponents];
	perturbations_t perturbations( numOfComponents, numberOfPerturbedNodes );
	perturbations.set_components();

	for ( size_t node = 0; node < numberOfPerturbedNodes; node++ ) {
		m_initialTopology->get_node( perturbedNodes[node], initialNodeData );
		perturbations.get_node( node, perturbationNodeData );
		for ( size_t comp = 0; comp < numOfComponents; comp++ ) {
			for ( size_t inProd = 0; inProd < numOfComponents; comp++ )
				perturbationNodeData[comp][0] += initialNodeData[inProd][0]*transformationMatrix[ comp*(numOfComponents+1) +inProd ];
			perturbationNodeData[comp][0] += transformationMatrix[ (comp+1)*(numOfComponents+1) ];
			perturbationNodeData[comp][0] -= initialNodeData[comp][0];
		}
	}

	insert_perturbations( perturbations.get_allComponents(), numberOfPerturbedNodes );
	return numberOfPerturbedNodes;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::insert_perturbationData ( flagValues_t flagsOfPerturbedNodes, T* transformationMatrix )
{

	size_t numberOfPerturbedNodes(0);
	for ( typename flagValues_t::iterator flag = flagsOfPerturbedNodes.begin(); flag != flagsOfPerturbedNodes.end(); flag++ )
		numberOfPerturbedNodes += insert_perturbedNodesFlag( *flag );
	size_t* perturbedNodes = m_perturbedNodes->end() -numberOfPerturbedNodes;

	size_t numOfComponents = m_initialTopology->get_numOfComponents();

	T *initialNodeData[numOfComponents], *perturbationNodeData[numOfComponents];
	perturbations_t perturbations( numOfComponents, numberOfPerturbedNodes );
	perturbations.set_components();
	
	for ( size_t node = 0; node < numberOfPerturbedNodes; node++ ) {
		m_initialTopology->get_node( perturbedNodes[node], initialNodeData );
		perturbations.get_node( node, perturbationNodeData );
		for ( size_t comp = 0; comp < numOfComponents; comp++ ) {
			for ( size_t inProd = 0; inProd < numOfComponents; inProd++ )
				perturbationNodeData[comp][0] += initialNodeData[inProd][0]*transformationMatrix[ comp*(numOfComponents+1) +inProd ];
			perturbationNodeData[comp][0] += transformationMatrix[ comp*(numOfComponents+1) +numOfComponents ];
			perturbationNodeData[comp][0] -= initialNodeData[comp][0];
		}
	}

	insert_perturbations( perturbations.get_allComponents(), numberOfPerturbedNodes );
	return numberOfPerturbedNodes;
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbationData ( size_t* perturbedNodes, size_t numberOfPerturbedNodes, T** perturbations, size_t numberOfComponents )
{
	set_perturbations(); set_allNodes();
	insert_perturbationData ( perturbedNodes, numberOfPerturbedNodes, perturbations );
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbationData ( size_t* perturbedNodes, size_t numberOfPerturbedNodes, T* transformationMatrix )
{
	set_perturbations(); set_allNodes();
	insert_perturbationData ( perturbedNodes, numberOfPerturbedNodes, transformationMatrix );
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::set_perturbationData ( flagValues_t flagsOfPerturbedNodes, T* transformationMatrix )
{
	set_perturbations(); set_allNodes();
	insert_perturbationData ( flagsOfPerturbedNodes, transformationMatrix );
}



/*
template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::find_adaptedNodes ( ) {
	set_allNodes();

	if ( m_perturbedNodes ) {
		m_adaptedNodes->resize( m_initialTopology->get_numOfNodes() -m_perturbedNodes->size() );
		m_initialTopology->find_complementaryNodes( m_perturbedNodes->begin(), m_perturbedNodes->end(), m_adaptedNodes->begin() );
		return 1;
	}

	if ( !m_adaptedNodeFlags.empty() ) {
		std::vector<size_t> nodes;
		this->m_initialTopology->get_nodesWithFlag( m_adaptedNodeFlags, &nodes );
		m_adaptedNodes->set_content( nodes );
	}
}


template<class T, size_t T2, class T3>
int TopologyAdapter<T, T2, T3>::find_adaptedNodeFlags()
{ 	m_adaptedNodeFlags = m_initialTopology->find_complementarySet( m_initialTopology->get_flagValues(), m_perturbedNodeFlags ); }
*/




#endif
