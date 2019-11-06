#include "myQMeshRendered.h"
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QAttribute>
#include <Qt3DCore/QEntity>


	
template<class T>
myQMeshRendered<T>::myQMeshRendered( Qt3DCore::QEntity *rootEntity, Mesh<T, 3, int>* mesh, QColor color ): m_rootEntity(rootEntity), m_mesh(mesh)
{

	size_t numOfNodes = m_mesh->get_numOfNodes();
	size_t numOfConnectivities = m_mesh->get_numOfConnectivities();
	T** nodes = m_mesh->get_allComponents();
	size_t* connectivities = m_mesh->get_connectivities()->data();



	// =================== Create necessary objects and their relations 
	auto *elementEntity 	= new Qt3DCore::QEntity		( m_rootEntity );
	auto *element 		= new Qt3DRender::QGeometryRenderer( elementEntity );
	auto *geometry 		= new Qt3DRender::QGeometry	( element );

	auto *vertexBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::VertexBuffer, geometry);
	auto *colorBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::VertexBuffer, geometry);
	auto *indexBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::IndexBuffer, geometry);

	auto *positionAttribute = new Qt3DRender::QAttribute	( geometry );
	auto *colorAttribute 	= new Qt3DRender::QAttribute	( geometry );
	auto *indexAttribute 	= new Qt3DRender::QAttribute	( geometry );
	auto *material		= new Qt3DExtras::QPerVertexColorMaterial( m_rootEntity );


	// position vertices (start and end)
	m_vertexBufferData.resize( numOfNodes*3*sizeof(T)); 
	T *positions = reinterpret_cast<T*>(m_vertexBufferData.data());

	m_colorBufferData.resize( numOfNodes*3*sizeof(T)); 
	T *colors = reinterpret_cast<T*>(m_colorBufferData.data());

	// connectivity between vertices
	m_indexBufferData.resize( numOfConnectivities*3*sizeof(unsigned int));
	unsigned int *indices = reinterpret_cast<unsigned int*>(m_indexBufferData.data());


	vertexBuffer->setData( m_vertexBufferData);
	colorBuffer->setData( m_colorBufferData);
	indexBuffer->setData( m_indexBufferData);
	geometry->addAttribute(positionAttribute); // We add the vertices in the geometry
	geometry->addAttribute(colorAttribute);
	geometry->addAttribute(indexAttribute); // We add the indices linking the points in the geometry
	element->setGeometry(geometry);
	element->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
	elementEntity->addComponent(element);
	elementEntity->addComponent(material);


	// ==================== Define data and settings
	for ( size_t nodeI = 0; nodeI < numOfNodes; nodeI++ ) { 
		*positions++ = nodes[0][nodeI];
		*positions++ = nodes[1][nodeI];
		*positions++ = nodes[2][nodeI];
		*colors++ = 1.f;
		*colors++ = 0.f;
		*colors++ = 0.f;
	}

	for ( size_t nodeI = 0; nodeI < numOfConnectivities*3; nodeI++ )
		*indices++ = connectivities[nodeI];


	positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
	positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
	positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
	positionAttribute->setVertexSize(3);
	positionAttribute->setBuffer(vertexBuffer);
	positionAttribute->setByteStride(3 * sizeof(float));
	positionAttribute->setCount( numOfNodes );


	colorAttribute->setName(Qt3DRender::QAttribute::defaultColorAttributeName());
	colorAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
	colorAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
	colorAttribute->setVertexSize(3);
	colorAttribute->setBuffer( colorBuffer);
	colorAttribute->setByteStride(3 * sizeof(float));
	colorAttribute->setCount( numOfNodes );


	indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
	indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedInt);
	indexAttribute->setBuffer( indexBuffer );
	indexAttribute->setCount( numOfConnectivities*3 );
	
}



