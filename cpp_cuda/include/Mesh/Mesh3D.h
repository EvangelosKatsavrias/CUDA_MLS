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

#ifndef MESH3DHEADERFILE
#define MESH3DHEADERFILE

#include"Mesh.h"
#include"FacesArea.h"
#include"CellsVolume.h"


template<class nodalComponentType = float, size_t numOfVerticesPerFace = 3, class flagType = int, size_t numOfFacesPerCell = 4 >
class Mesh3D: public Mesh<nodalComponentType, numOfVerticesPerFace, flagType>
{

public:
typedef Mesh<nodalComponentType, numOfVerticesPerFace, flagType>				mesh_t;
typedef Mesh3D<nodalComponentType, numOfVerticesPerFace, flagType, numOfFacesPerCell> 		mesh3d_t;
typedef Topology<nodalComponentType, numOfFacesPerCell, flagType>				topologyOfFaces_t;
typedef Connectivities<size_t, numOfFacesPerCell>						faceConnectivities_t;
typedef CellsVolume<nodalComponentType, size_t, numOfVerticesPerFace, numOfFacesPerCell> 	cellsVolume_t;
typedef FacesArea< nodalComponentType, size_t, numOfVerticesPerFace > 				facesArea_t;
//typedef ConnectivitiesN										nodeFacesMap_t;
using	typename mesh_t::topologyOfVertices_t::nodeFlags_t;
using	typename mesh_t::topologyOfVertices_t::flagValues_t;
using	typename mesh_t::topologyOfVertices_t::flagsMap_t;
using	typename mesh_t::meshMetric_t;
using	typename mesh_t::vertexConnectivities_t;



protected:
	facesArea_t* 			m_facesArea;
	cellsVolume_t* 			m_cellsVolume;
//	nodeFacesMap_t			m_nodeFacesMap;
	faceConnectivities_t*		m_faceConnectivities;


public:
	Mesh3D(	size_t 			numOfVertices = 4, 
		nodalComponentType** 	vertexData = 0, 
		vertexConnectivities_t* vertexConnectivities = 0, 
		nodeFlags_t* 		vertexNodeFlags = 0, 
		flagValues_t* 		interiorVertexFlagValues = 0, 
		flagsMap_t* 		vertexFlagToNodesMap = 0, 
		facesArea_t* 		facesArea = 0, 
		faceConnectivities_t* 	faceConnectivities = 0, 
		cellsVolume_t* 		cellsVolume = 0 );
	Mesh3D(mesh_t* referenceTopology);
	Mesh3D(mesh3d_t* referenceTopology);
	Mesh3D(mesh_t& referenceTopology);
	Mesh3D(mesh3d_t& referenceTopology);
	~Mesh3D();


	int 		virtual	find_facesArea 		( );
	int 		virtual	find_cellsVolume 	( );
	int			find_metrics 		( ) 	{ find_facesArea(); /*find_cellsVolume();*/ }


	int 			set_numOfFaces 		( size_t numOfFaces );
	int 			set_numOfCells 		( size_t numOfCells );
//virtual int 			set_vertexConnectivities( vertexConnectivities_t* connectivities );// { return this->set_connectivities( connectivities ); }
	int 			set_faceConnectivities	( faceConnectivities_t* connectivities ) 	{ m_faceConnectivities = connectivities; }
	int 			set_facesArea 		( facesArea_t* facesArea ) 			{ m_facesArea = facesArea; }
	int 			set_cellsVolume 	( cellsVolume_t* cellsVolume ) 			{ m_cellsVolume = cellsVolume; }
	int			set_metrics 		( facesArea_t* metrics ) 			{ m_facesArea = metrics;  }


	std::string virtual	get_typename		( ) { return "Mesh3D"; }
	size_t 			get_numOfFaces 		( ) { return this->get_numOfConnectivities(); }
	size_t			get_numOfCells 		( ) { return m_faceConnectivities.size(); }
//virtual vertexConnectivities_t* get_vertexConnectivities( );// { return this->get_connectivities(); }
	faceConnectivities_t* 	get_faceConnectivities	( ) { return m_faceConnectivities; }
	facesArea_t* 		get_facesArea 		( ) { return m_facesArea; }
	cellsVolume_t* 		get_cellsVolume 	( ) { return m_cellsVolume; }
	meshMetric_t* 		get_metrics 		( ) { return get_facesArea(); }

};


template<class T, size_t T2, class T3, size_t T4>
Mesh3D<T, T2, T3, T4>::Mesh3D (	size_t 			numOfVertices, 
				T** 			vertexNodesData, 
				vertexConnectivities_t* vertexConnectivities, 
				nodeFlags_t* 		vertexNodeFlags, 
				flagValues_t* 		interiorVertexFlagValues, 
				flagsMap_t* 		vertexFlagToNodesMap, 
				facesArea_t*		facesArea,
				faceConnectivities_t*	faceConnectivities, 
				cellsVolume_t*		cellsVolume ):
Mesh< T, T2, T3> ( numOfVertices, vertexNodesData, vertexConnectivities, vertexNodeFlags, interiorVertexFlagValues, vertexFlagToNodesMap ), m_facesArea(facesArea), m_faceConnectivities( faceConnectivities ), m_cellsVolume(cellsVolume)
{
	if ( !facesArea ) 		{ m_facesArea 		= new facesArea_t(this, this->m_connectivities); 				this->m_memAllocatedMembers[0] = 1; }
	if ( !cellsVolume ) 		{ m_cellsVolume 	= new cellsVolume_t(this, this->m_connectivities, m_faceConnectivities);	this->m_memAllocatedMembers[1] = 1; }
	if ( !faceConnectivities ) 	{ m_faceConnectivities 	= new faceConnectivities_t; 							this->m_memAllocatedMembers[2] = 1; }
}


template<class T, size_t T2, class T3, size_t T4>
Mesh3D<T, T2, T3, T4>::Mesh3D(Mesh<T, T2, T3>* referenceTopology): Mesh<T, T2, T3>(referenceTopology)
{
	m_facesArea 		= new facesArea_t(this, this->m_connectivities); 				this->m_memAllocatedMembers[0] = 1;
	m_cellsVolume 		= new cellsVolume_t(this, this->m_connectivities, m_faceConnectivities);	this->m_memAllocatedMembers[1] = 1;
	m_faceConnectivities 	= new faceConnectivities_t; 							this->m_memAllocatedMembers[2] = 1;
}


template<class T, size_t T2, class T3, size_t T4>
Mesh3D<T, T2, T3, T4>::Mesh3D(Mesh3D<T, T2, T3, T4>* referenceTopology): Mesh<T, T2, T3>(referenceTopology)
{
	m_facesArea 		= referenceTopology->get_facesArea();
	m_cellsVolume 		= referenceTopology->get_cellsVolume();
	m_faceConnectivities 	= referenceTopology->get_faceConnectivities();
}


template<class T, size_t T2, class T3, size_t T4>
Mesh3D<T, T2, T3, T4>::Mesh3D(Mesh<T, T2, T3>& referenceTopology): Mesh<T, T2, T3>(referenceTopology)
{ 
	m_facesArea 		= new facesArea_t(this, this->m_connectivities); 				this->m_memAllocatedMembers[0] = 1;
	m_cellsVolume 		= new cellsVolume_t(this, this->m_connectivities, m_faceConnectivities);	this->m_memAllocatedMembers[1] = 1;
	m_faceConnectivities 	= new faceConnectivities_t; 							this->m_memAllocatedMembers[2] = 1;
}

	
template<class T, size_t T2, class T3, size_t T4>
Mesh3D<T, T2, T3, T4>::Mesh3D(Mesh3D<T, T2, T3, T4>& referenceTopology): Mesh<T, T2, T3>(referenceTopology)
{ 
	m_facesArea 		= new facesArea_t( *(referenceTopology.get_facesArea() ) );
	m_cellsVolume 		= new cellsVolume_t( *(referenceTopology.get_elementsVolume() ) );
	m_faceConnectivities 	= new faceConnectivities_t( *(referenceTopology.get_faceConnectivities() ) );
}


template<class T, size_t T2, class T3, size_t T4>
Mesh3D<T, T2, T3, T4>::~Mesh3D()
{
	if (!this->m_memAllocatedMembers[0] ) { delete m_facesArea; }
	if (!this->m_memAllocatedMembers[1] ) { delete m_cellsVolume; }
	if (!this->m_memAllocatedMembers[2] ) { delete m_faceConnectivities; }
}


template<class T, size_t T2, class T3, size_t T4>
int Mesh3D<T, T2, T3, T4>::set_numOfFaces(size_t numOfFaces) {
	Topology<T, T2, T3>::set_numOfConnectivities( numOfFaces );
	m_facesArea->resize( numOfFaces, 0);
}


template<class T, size_t T2, class T3, size_t T4>
int Mesh3D<T, T2, T3, T4>::set_numOfCells( size_t numOfCells ) {
	set_faceConnectivities ( numOfCells );
	m_cellsVolume->resize( numOfCells, 0);
}


template<class T, size_t T2, class T3, size_t T4>
int Mesh3D<T, T2, T3, T4>::find_cellsVolume () { m_cellsVolume->find_volumes(); }


template<class T, size_t T2, class T3, size_t T4>
int Mesh3D<T, T2, T3, T4>::find_facesArea () { m_facesArea->find_areas(); }


#endif
