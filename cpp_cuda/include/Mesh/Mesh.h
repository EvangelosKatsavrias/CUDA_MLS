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

#ifndef MESHHEADERFILE
#define MESHHEADERFILE

#include<stdexcept>
#include<exception>
#include"ContainerBase.h"
#include"Topology.h"
#include"myBLAS.h"


template<class nodalComponentType = float, size_t numOfConnectivityNodes = 3, class flagType = int>
class Mesh: public Topology<nodalComponentType, numOfConnectivityNodes, flagType>
{

public:
typedef Topology<nodalComponentType, numOfConnectivityNodes, flagType>		topologyOfVertices_t;
typedef Mesh<nodalComponentType, numOfConnectivityNodes, flagType>		mesh_t;
typedef Connectivities<size_t, numOfConnectivityNodes>				vertexConnectivities_t;
typedef ContainerBase<nodalComponentType, std::vector<nodalComponentType> > 	meshMetric_t;
using	typename topologyOfVertices_t::nodeFlags_t;
using	typename topologyOfVertices_t::flagValues_t;
using	typename topologyOfVertices_t::flagsMap_t;


protected:
	std::vector<bool> m_memAllocatedMembers;


public:
	Mesh(	size_t 				numOfVertices = 4, 
		nodalComponentType** 		vertexData = 0, 
		vertexConnectivities_t* 	vertexConnectivities = 0, 
		nodeFlags_t* 			vertexNodeFlags = 0, 
		flagValues_t* 			interiorVertexFlagValues = 0, 
		flagsMap_t* 			vertexFlagToNodesMap = 0 );
	Mesh(mesh_t* referenceTopology);
	Mesh(mesh_t& referenceTopology);
	~Mesh();


virtual int			find_metrics 		 ( ) { }

virtual int 			set_vertexConnectivities ( vertexConnectivities_t* connectivities ) { return this->set_connectivities( connectivities ); }
virtual int			set_metrics 		 ( meshMetric_t* metrics ) {  }
	int			set_numOfVertices 	 ( size_t numOfNodes ) { return this->set_numOfNodes(numOfNodes); }


virtual	std::string		get_typename		( ) { return "Mesh"; }
virtual vertexConnectivities_t* get_vertexConnectivities( ) { return this->get_connectivities(); }
	size_t			get_numOfVertices 	( ) { return this->get_numOfNodes(); }
virtual size_t			get_numOfFaces		( ) { return this->get_numOfConnectivities(); }
virtual	meshMetric_t* 		get_metrics		( ) { }
virtual	meshMetric_t*		get_qualityMetrics	() {}
};


template<class T, size_t T2, class T3>
Mesh<T, T2, T3>::Mesh(	size_t 			numOfVertices, 
			T** 			vertexNodesData, 
			vertexConnectivities_t* vertexConnectivities, 
			nodeFlags_t* 		vertexNodeFlags, 
			flagValues_t* 		interiorVertexFlagValues, 
			flagsMap_t* 		vertexFlagToNodesMap ):
Topology< T, T2, T3> ( 3, numOfVertices, vertexNodesData, vertexConnectivities, vertexNodeFlags, interiorVertexFlagValues, vertexFlagToNodesMap ), m_memAllocatedMembers( std::vector<bool>(10, 0) )
{ }


template<class T, size_t T2, class T3>
Mesh<T, T2, T3>::Mesh(Mesh<T, T2, T3>* referenceTopology): Topology<T, T2, T3>(referenceTopology), m_memAllocatedMembers(std::vector<bool>(10, 0) )
{ }


template<class T, size_t T2, class T3>
Mesh<T, T2, T3>::Mesh(Mesh<T, T2, T3>& referenceTopology): Topology<T, T2, T3>(referenceTopology), m_memAllocatedMembers( std::vector<bool>(10, 1) )
{ }


template<class T, size_t T2, class T3>
Mesh<T, T2, T3>::~Mesh()
{ }


#endif
