#include "myQMesh.h"
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QAttribute>
#include <Qt3DCore/QEntity>


template<class T>
myQMesh<T>::myQMesh( Qt3DCore::QEntity *rootEntity, Mesh<T, 3, int>* mesh, QColor color ): m_rootEntity(rootEntity), m_mesh(mesh)
{

	size_t numOfNodes = m_mesh->get_numOfNodes();
	size_t numOfConnectivities = m_mesh->get_numOfConnectivities();
	T** nodes = m_mesh->get_allComponents();

	Connectivities<size_t, 3>& connect = m_mesh->get_connectivities()[0];
	size_t* elemVertices;
	for ( size_t element=0; element < numOfConnectivities; element++)
	{
		elemVertices = connect[element];
		m_uniqueEdges.insert( std::make_pair( elemVertices[0], elemVertices[1] ) );
		m_uniqueEdges.insert( std::make_pair( elemVertices[1], elemVertices[2] ) );
		m_uniqueEdges.insert( std::make_pair( elemVertices[2], elemVertices[0] ) );
	}

	numOfConnectivities = m_uniqueEdges.size();



	// =================== Create necessary objects and their relations 
	auto *lineEntity 	= new Qt3DCore::QEntity		( m_rootEntity );
	auto *line 		= new Qt3DRender::QGeometryRenderer( lineEntity );
	auto *geometry 		= new Qt3DRender::QGeometry	( line );

	auto *vertexBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::VertexBuffer, geometry);
	auto *colorBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::VertexBuffer, geometry);
	auto *indexBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::IndexBuffer, geometry);

	auto *positionAttribute = new Qt3DRender::QAttribute	( geometry );
	auto *colorAttribute 	= new Qt3DRender::QAttribute	( geometry );
	auto *indexAttribute 	= new Qt3DRender::QAttribute	( geometry );
	auto *material		= new Qt3DExtras::QPhongMaterial( m_rootEntity );


	// position vertices (start and end)
	m_vertexBufferData.resize( numOfNodes*3*sizeof(T)); 
	T *positions = reinterpret_cast<T*>(m_vertexBufferData.data());

	m_colorBufferData.resize( numOfNodes*3*sizeof(T)); 
	T *colors = reinterpret_cast<T*>(m_colorBufferData.data());

	// connectivity between vertices
	m_indexBufferData.resize( numOfConnectivities*2*sizeof(unsigned int));
	unsigned int *indices = reinterpret_cast<unsigned int*>(m_indexBufferData.data());


	vertexBuffer->setData( m_vertexBufferData);
	colorBuffer->setData( m_colorBufferData);
	indexBuffer->setData( m_indexBufferData);
	geometry->addAttribute(positionAttribute); // We add the vertices in the geometry
//	geometry->addAttribute(colorAttribute);
	geometry->addAttribute(indexAttribute); // We add the indices linking the points in the geometry
	line->setGeometry(geometry);
	line->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);
	lineEntity->addComponent(line);
	lineEntity->addComponent(material);


	// ==================== Define data and settings
	for ( size_t nodeI = 0; nodeI < numOfNodes; nodeI++ ) { 
		*positions++ = nodes[0][nodeI];
		*positions++ = nodes[1][nodeI];
		*positions++ = nodes[2][nodeI];
	}

	for ( auto& edgePoints: m_uniqueEdges ) {
		*indices++ = edgePoints.first;
		*indices++ = edgePoints.second;
	}


	positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
	positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
	positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
	positionAttribute->setVertexSize(3);
	positionAttribute->setBuffer(vertexBuffer);
	positionAttribute->setByteStride(3 * sizeof(float));
	positionAttribute->setCount( numOfNodes );

	/*
	colorAttribute->setName(Qt3DRender::QAttribute::defaultColorAttributeName());
	colorAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
	colorAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
	colorAttribute->setVertexSize(3);
	colorAttribute->setBuffer(vertexBuffer);
	colorAttribute->setByteStride(3 * sizeof(float));
	colorAttribute->setCount( numOfNodes );
*/

	indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
	indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedInt);
	indexAttribute->setBuffer( indexBuffer );
	indexAttribute->setCount( numOfConnectivities*2 );

	material->setAmbient( color );
}


