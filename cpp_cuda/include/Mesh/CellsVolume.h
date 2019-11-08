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

#ifndef CELLSVOLUMEHEADERFILE
#define CELLSVOLUMEHEADERFILE

#include<stdexcept>
#include<exception>
#include"ContainerBase.h"
#include"Topology.h"
#include"myBLAS.h"
#include"CellsVolume.h"


template< typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace = 3, size_t numOfFaces = 4 >
class CellsVolume: public ContainerBase<nodalComponentType, std::vector<nodalComponentType> >
{

typedef ContainerBase<nodalComponentType, std::vector<nodalComponentType> > containerBase_t;

protected:
	Nodes<nodalComponentType>*						m_nodeCoords;
	Connectivities<connectivityType, numOfVerticesPerFace>*			m_vertexConnectivities;
	Connectivities<connectivityType, numOfFaces>*				m_faceConnectivities;
	PolyhedronVolume<nodalComponentType, numOfVerticesPerFace, numOfFaces>*	m_volumeEstimator;
	std::vector<bool>							m_allocatedMemberFlags;


public:
	CellsVolume ( Nodes<nodalComponentType>* nodes=0, Connectivities<connectivityType, numOfVerticesPerFace>* vertexConnectivities=0, Connectivities<connectivityType, numOfFaces>* faceConnectivities = 0, PolyhedronVolume<nodalComponentType, numOfVerticesPerFace, numOfFaces>* estimator=0 );
	CellsVolume ( CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>* referenceObject);
	CellsVolume ( CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>& referenceObject);
	~CellsVolume ();


	Nodes <nodalComponentType>*						get_nodeCoords() 		{ return m_nodeCoords; }
	Connectivities <connectivityType, numOfVerticesPerFace>*		get_vertexConnectivities()	{ return m_vertexConnectivities; }
	Connectivities <connectivityType, numOfFaces>*				get_faceConnectivities()	{ return m_faceConnectivities; }
	PolyhedronVolume <nodalComponentType, numOfVerticesPerFace, numOfFaces>*get_volumeEstimator()		{ return m_volumeEstimator; }


	int set_nodeCoords (Nodes<nodalComponentType>* nodeCoords) 							{ m_nodeCoords 		= nodeCoords; }
	int set_vertexConnectivities (Connectivities<connectivityType, numOfVerticesPerFace>* connectivities)		{ m_vertexConnectivities= connectivities; }
	int set_faceConnectivities (Connectivities<connectivityType, numOfFaces>* connectivities)			{ m_faceConnectivities 	= connectivities; }
	int set_volumeEstimator (PolyhedronVolume<nodalComponentType, numOfVerticesPerFace, numOfFaces>* volumeEstimator){ m_volumeEstimator 	= volumeEstimator; }

	int virtual find_volumes();
};


template<typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace, size_t numOfFaces>
CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>::CellsVolume(Nodes<nodalComponentType>* nodes, Connectivities<connectivityType, numOfVerticesPerFace>* vertexConnectivities, Connectivities<connectivityType, numOfFaces>* faceConnectivities, PolyhedronVolume<nodalComponentType, numOfVerticesPerFace, numOfFaces>* estimator ):
m_nodeCoords(nodes), m_vertexConnectivities(vertexConnectivities), m_faceConnectivities(faceConnectivities), m_volumeEstimator(estimator), m_allocatedMemberFlags(std::vector<bool>(4,0))
{
	if (!nodes) 			{ m_nodeCoords 		 = new Nodes<nodalComponentType>; 						m_allocatedMemberFlags[0]=1; }
	if (!vertexConnectivities) 	{ m_vertexConnectivities = new Connectivities<connectivityType, numOfVerticesPerFace>; 			m_allocatedMemberFlags[1]=1; }
	if (!faceConnectivities) 	{ m_faceConnectivities 	 = new Connectivities<connectivityType, numOfFaces>; 				m_allocatedMemberFlags[2]=1; }
	if (!estimator) 		{ m_volumeEstimator 	 = new PolyhedronVolume<nodalComponentType, numOfVerticesPerFace, numOfFaces>; 	m_allocatedMemberFlags[3]=1; }
}


template<typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace, size_t numOfFaces>
CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>::CellsVolume( CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>* referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType>>(referenceObject), m_allocatedMemberFlags(std::vector<bool>(4,0))
{
	m_nodeCoords 		= referenceObject->get_nodeCoords();
	m_vertexConnectivities	= referenceObject->get_vertexConnectivities();
	m_faceConnectivities	= referenceObject->get_faceConnectivities();
	m_volumeEstimator 	= referenceObject->get_volumeEstimator();
}


template<typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace, size_t numOfFaces>
CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>::CellsVolume( CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>& referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType> >(referenceObject), m_allocatedMemberFlags(std::vector<bool>(4,1))
{
	m_nodeCoords 		= new Nodes<nodalComponentType>( *(referenceObject.get_nodeCoords() ) );
	m_vertexConnectivities 	= new Connectivities<connectivityType, numOfVerticesPerFace>( *(referenceObject.get_vertexConnectivities() ) );
	m_faceConnectivities 	= new Connectivities<connectivityType, numOfFaces>( *(referenceObject.get_faceConnectivities() ) );
	m_volumeEstimator 	= new PolyhedronVolume<nodalComponentType, numOfVerticesPerFace, numOfFaces>;
}


template<typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace, size_t numOfFaces>
CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>::~CellsVolume()
{
	if (m_allocatedMemberFlags[0]) delete m_nodeCoords;
	if (m_allocatedMemberFlags[1]) delete m_vertexConnectivities;
	if (m_allocatedMemberFlags[2]) delete m_faceConnectivities;
	if (m_allocatedMemberFlags[3]) delete m_volumeEstimator;
}


template<typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace, size_t numOfFaces>
int CellsVolume<nodalComponentType, connectivityType, numOfVerticesPerFace, numOfFaces>::find_volumes()
{

	nodalComponentType *x[numOfVerticesPerFace], *y[numOfVerticesPerFace], *z[numOfVerticesPerFace];

	this->resize(m_faceConnectivities->get_numOfConnectivities());
	size_t vertexShift;
	connectivityType* cellVertices;
	connectivityType* cellFaces;

	for ( size_t cellIndex = 0; cellIndex < m_faceConnectivities->get_numOfConnectivities(); cellIndex++ ) {

		cellFaces = m_faceConnectivities->at(cellIndex);

		for ( size_t faceIndex = 0; faceIndex < numOfFaces; faceIndex++ ) {

			vertexShift 	= faceIndex*numOfVerticesPerFace;
			cellVertices = m_vertexConnectivities[0][ cellFaces[faceIndex] ];

			for (size_t vertexIndex = 0; vertexIndex < numOfVerticesPerFace; vertexIndex++ ) {
				x[vertexShift +vertexIndex] = m_nodeCoords[0][0] +cellVertices[vertexIndex];
				y[vertexShift +vertexIndex] = m_nodeCoords[0][1] +cellVertices[vertexIndex];
				z[vertexShift +vertexIndex] = m_nodeCoords[0][2] +cellVertices[vertexIndex];
			}
		}

		this->m_data[cellIndex] = (*m_volumeEstimator)(x, y, z);
	}
}


#endif
