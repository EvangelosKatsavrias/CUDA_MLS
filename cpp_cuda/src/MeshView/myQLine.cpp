#include "myQLine.h"
#include <Qt3DRender/QGeometry>
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QMesh>



myQLine::myQLine(Qt3DCore::QEntity *rootEntity, QVector3D start, QVector3D end, QColor color): m_rootEntity(rootEntity), m_start(start), m_end(end), m_materialColor(color), m_vertexBufferData( QByteArray( 6*sizeof(float), 0 ) ), m_indexBufferData( QByteArray(2*sizeof(unsigned int), 0) )
{

	// =================== Create necessary objects and their relations 
	auto *m_lineEntity 	= new Qt3DCore::QEntity		( m_rootEntity );
	auto *line 		= new Qt3DRender::QGeometryRenderer( m_lineEntity );
	auto *geometry 		= new Qt3DRender::QGeometry	( line );

	auto *vertexBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::VertexBuffer, geometry);
//	auto *colorBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::VertexBuffer, geometry);
	auto *indexBuffer 	= new Qt3DRender::QBuffer	( Qt3DRender::QBuffer::IndexBuffer, geometry);

	auto *positionAttribute = new Qt3DRender::QAttribute	( geometry );
//	auto *colorAttribute 	= new Qt3DRender::QAttribute	( geometry );
	auto *indexAttribute 	= new Qt3DRender::QAttribute	( geometry );
	auto *m_material	= new Qt3DExtras::QPhongMaterial( m_rootEntity );


	// position vertices (start and end)
	m_positions = reinterpret_cast<float*>(m_vertexBufferData.data());
	// connectivity between vertices
	m_indices = reinterpret_cast<unsigned int*>(m_indexBufferData.data());
	*m_indices++ = 0; *m_indices++ = 1;
/*	QByteArray colorBufferData(6*sizeof(float), 0); // start.x, start.y, start.end + end.x, end.y, end.z
	float *m_colors = reinterpret_cast<float*>(colorBufferData.data());*/


	vertexBuffer->setData(m_vertexBufferData);
//	colorBuffer->setData(colorBufferData);
	indexBuffer->setData(m_indexBufferData);
	geometry->addAttribute(positionAttribute); // We add the vertices in the geometry
//	geometry->addAttribute(colorAttribute);
	geometry->addAttribute(indexAttribute); // We add the indices linking the points in the geometry
	line->setGeometry(geometry);
	line->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);	
	m_lineEntity->addComponent(line);
	m_lineEntity->addComponent(m_material);


	positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
	positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
	positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
	positionAttribute->setVertexSize( 3 );
	positionAttribute->setBuffer(vertexBuffer);
	positionAttribute->setByteStride( 3*sizeof(float));
	positionAttribute->setCount( 2 ); //number of vertices


	indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
	indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedInt);
	indexAttribute->setBuffer(indexBuffer);
	indexAttribute->setCount( 2 );

	setPositions( start, end );
	m_material->setAmbient( QColor( m_materialColor ) );

/*	colorAttribute->setName(Qt3DRender::QAttribute::defaultColorAttributeName());
	colorAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
	colorAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
	colorAttribute->setVertexSize( 3 );
	colorAttribute->setBuffer(vertexBuffer);
	colorAttribute->setByteStride( 3*sizeof(float));
	colorAttribute->setCount( 2 );*/
}


