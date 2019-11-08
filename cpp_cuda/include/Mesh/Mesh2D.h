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

#ifndef MESH2DHEADERFILE
#define MESH2DHEADERFILE

#include"Mesh.h"
#include"FacesArea.h"


template<class nodalComponentType = float, size_t numOfConnectivityNodes = 3, class flagType = int>
class Mesh2D: public Mesh<nodalComponentType, numOfConnectivityNodes, flagType>
{

public:
typedef Mesh<nodalComponentType, numOfConnectivityNodes, flagType>		mesh_t;
typedef Mesh2D<nodalComponentType, numOfConnectivityNodes, flagType> 		mesh2d_t;
typedef ElementsArea<nodalComponentType, size_t, numOfConnectivityNodes> 	elementsArea_t;
typedef ElementsQuality<nodalComponentType, size_t, numOfConnectivityNodes> 	qualityMetrics_t;

using	typename mesh_t::topologyOfVertices_t::nodeFlags_t;
using	typename mesh_t::topologyOfVertices_t::flagValues_t;
using	typename mesh_t::topologyOfVertices_t::flagsMap_t;
using	typename mesh_t::meshMetric_t;
using	typename mesh_t::vertexConnectivities_t;


protected:
	elementsArea_t*		m_elementsArea;
	qualityMetrics_t*	m_qualityMetrics;
	std::vector<bool> 	m_memAllocatedMembers;

public:
	Mesh2D(	size_t 			numOfNodes = 3, 
		nodalComponentType** 	nodesData=0, 
		vertexConnectivities_t* connectivities=0, 
		nodeFlags_t* 		nodeFlags=0, 
		flagValues_t* 		interiorFlagValues = 0, 
		flagsMap_t* 		flagToNodesMap = 0, 
		elementsArea_t* 	elementsArea=0,
	     	qualityMetrics_t*	qualityMetrics=0 );
	Mesh2D(mesh_t* referenceTopology);
	Mesh2D(mesh2d_t* referenceTopology);
	Mesh2D(mesh_t& referenceTopology);
	Mesh2D(mesh2d_t& referenceTopology);
	~Mesh2D();


	int 	virtual	find_elementsArea ();
	int 	virtual	find_qualityMetrics ();
	int		find_metrics()	{ find_elementsArea(); find_qualityMetrics(); }

	int 		set_elementsArea 	( elementsArea_t* elementsArea ){ m_elementsArea = elementsArea; m_memAllocatedMembers[0]=0; m_qualityMetrics->set_elementsArea( elementsArea ); }
	int 		set_qualityMetrics 	( qualityMetrics_t* qualityMetrics ){ m_qualityMetrics = qualityMetrics; m_memAllocatedMembers[1]=0; }
	int 		set_numOfElements 	( size_t numOfConnectivities) 	{ return set_numOfConnectivities( numOfConnectivities ); }
	int		set_numOfConnectivities	( size_t numOfConnectivities);


	std::string 	get_typename		( ) { return "Mesh2D";}
	size_t		get_numOfElements 	( ) { return this->get_numOfConnectivities(); }
	elementsArea_t* get_elementsArea 	( ) { return m_elementsArea; }
virtual	meshMetric_t* 	get_metrics 		( ) { return get_elementsArea(); }
virtual	meshMetric_t*	get_qualityMetrics	( ) { return m_qualityMetrics; }

};


template<class T, size_t T2, class T3>
Mesh2D<T, T2, T3>::Mesh2D (	size_t 			numOfNodes, 
				T** 			nodesData, 
				vertexConnectivities_t* connectivities, 
				nodeFlags_t* 		nodeFlags, 
				flagValues_t* 		interiorFlagValues, 
				flagsMap_t* 		flagToNodesMap, 
				elementsArea_t* 	elementsArea,
	       			qualityMetrics_t*	qualityMetrics	):
Mesh<T, T2, T3> ( numOfNodes, nodesData, connectivities, nodeFlags, interiorFlagValues, flagToNodesMap ), m_elementsArea(elementsArea), m_qualityMetrics(qualityMetrics), m_memAllocatedMembers(std::vector<bool>(10,0))
{
	if (!elementsArea) { m_elementsArea = new elementsArea_t(this, this->m_connectivities); m_memAllocatedMembers[0]=1; }
	if (!qualityMetrics) { m_qualityMetrics = new qualityMetrics_t(this, this->m_connectivities, m_elementsArea); m_memAllocatedMembers[1]=1; }
}

template<class T, size_t T2, class T3>
Mesh2D<T, T2, T3>::Mesh2D(Mesh<T, T2, T3>* referenceTopology): Mesh<T, T2, T3>(referenceTopology), m_memAllocatedMembers(std::vector<bool>(10,0))
{  m_elementsArea = new elementsArea_t(this, this->m_connectivities); m_memAllocatedMembers[0]=1; 
   m_qualityMetrics = new qualityMetrics_t(this, this->m_connectivities, m_elementsArea); m_memAllocatedMembers[1]=1;
}


template<class T, size_t T2, class T3>
Mesh2D<T, T2, T3>::Mesh2D(Mesh2D<T, T2, T3>* referenceTopology): Mesh<T, T2, T3>(referenceTopology) 
{ m_elementsArea = referenceTopology->get_elementsArea(); 
  m_qualityMetrics = referenceTopology->get_qualityMetrics();
  m_memAllocatedMembers = referenceTopology->m_memAllocatedMembers;
}

template<class T, size_t T2, class T3>
Mesh2D<T, T2, T3>::Mesh2D(Mesh<T, T2, T3>& referenceTopology): Mesh<T, T2, T3>(referenceTopology), m_memAllocatedMembers(std::vector<bool>(10,0))
{ m_elementsArea = new elementsArea_t(this, this->m_connectivities); m_memAllocatedMembers[0]=1;
  m_qualityMetrics = new qualityMetrics_t(this, this->m_connectivities, m_elementsArea); m_memAllocatedMembers[1]=1;}


template<class T, size_t T2, class T3>
Mesh2D<T, T2, T3>::Mesh2D(Mesh2D<T, T2, T3>& referenceTopology): Mesh<T, T2, T3>(referenceTopology), m_memAllocatedMembers(std::vector<bool>(10,0))
{ m_elementsArea = new elementsArea_t( *referenceTopology.get_elementsArea() ); m_memAllocatedMembers[0]=1;
 m_qualityMetrics = new qualityMetrics_t( *referenceTopology.get_qualityMetrics() ); m_memAllocatedMembers[1]=1;}


template<class T, size_t T2, class T3>
Mesh2D<T, T2, T3>::~Mesh2D()
{
	if (!m_memAllocatedMembers[0] ) { delete m_elementsArea; }
	if (!m_memAllocatedMembers[1] ) { delete m_qualityMetrics; }
}


template<class T, size_t T2, class T3>
int Mesh2D<T, T2, T3>::set_numOfConnectivities(size_t numOfConnectivities) {
	Topology<T, T2, T3>::set_numOfConnectivities( numOfConnectivities );
	m_elementsArea->resize( numOfConnectivities, 0);
}


template<class T, size_t T2, class T3>
int Mesh2D<T, T2, T3>::find_elementsArea () { m_elementsArea->find_areas(); }

template<class T, size_t T2, class T3>
int Mesh2D<T, T2, T3>::find_qualityMetrics () { m_qualityMetrics->find_qualities(); }


#endif
