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

#include"meshViewer.h"

#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>

#include <Qt3DRender/QMesh>
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QGeometry>
#include <Qt3DRender/QCamera>

#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QPerVertexColorMaterial>

#include <QOrbitCameraController>

#include "UCSAxesLines.h"


const double sizee(300);


template<typename T>
QRectF QGraphicsItem_Mesh2D_Edges<T>::boundingRect() const { return QRectF(-sizee, -sizee, 2*sizee, 2*sizee); }

template<typename T>
void QGraphicsItem_Mesh2D_Edges<T>::paint(QPainter* painter, const QStyleOptionGraphicsItem* styleOptions, QWidget* widget)
{
	T* x_coords = m_mesh->get_component(0);
	T* y_coords = m_mesh->get_component(1);

	QTransform tr = painter->transform(); double penThickness = 0.5/tr.m11();

	painter->setPen( QPen( Qt::red, penThickness ) );
	painter->drawLine( QPointF( 0.f, 0.f ), QPointF( 1.f, 0.f ) );
	painter->drawLine( QPointF( 1.f, 0.f ), QPointF( 0.7f, 0.2f ) );
	painter->drawLine( QPointF( 1.f, 0.f ), QPointF( 0.7f, -0.2f ) );
	painter->setPen( QPen( Qt::green, penThickness ) );
	painter->drawLine( QPointF( 0.f, 0.f ), QPointF( 0.f, 1.f ) );
	painter->drawLine( QPointF( 0.f, 1.f ), QPointF( 0.2f, 0.7f ) );
	painter->drawLine( QPointF( 0.f, 1.f ), QPointF( -0.2f, 0.7f ) );


	painter->setPen( QPen( Qt::black, penThickness ) );

	std::set< std::pair<size_t, size_t >, unorderLess<size_t> >& uniqueEdges = m_mesh->get_connectivities()->get_uniqueLinks();
	for ( auto& edgePoints: uniqueEdges )
		painter->drawLine( QPointF(x_coords[ edgePoints.first ], y_coords[ edgePoints.first ]), QPointF(x_coords[ edgePoints.second ], y_coords[ edgePoints.second ]) );
}



template<typename T>
QRectF QGraphicsItem_Mesh2D_Vertices<T>::boundingRect() const { return QRectF(-sizee, -sizee, 2*sizee, 2*sizee); }


template<typename T>
void QGraphicsItem_Mesh2D_Vertices<T>::paint(QPainter* painter, const QStyleOptionGraphicsItem* styleOptions, QWidget* widget)
{
	T* x_coords = m_mesh->get_component(0);
	T* y_coords = m_mesh->get_component(1);

	QTransform tr = painter->transform();

	double penThickness( 0.02 );

	painter->setPen( QPen( Qt::black, penThickness) );
	for (auto& i: m_mesh->get_interiorNodes() ) painter->drawPoint(QPointF(x_coords[i.second], y_coords[i.second]) );

	painter->setPen( QPen( Qt::blue, penThickness) );
	for (auto& i: m_mesh->get_boundaryNodes() ) painter->drawPoint(QPointF(x_coords[i.second], y_coords[i.second]) );


	painter->scale( 1/tr.m11(),  1/tr.m11() );

	QFont fontSet("Arial", 10); 
	painter->setFont(fontSet);
	qreal shift_x(tr.m11()/30), shift_y(tr.m11()/30);


//	painter->setPen( QPen( Qt::black, penThickness ) );
//	for (auto& i: m_mesh->get_interiorNodes() ) painter->drawText( x_coords[i.first]* tr.m11() +shift_x, y_coords[i.first]* tr.m11() +shift_y , QString::number(i.first) );


//	painter->setPen( QPen( Qt::blue, penThickness ) );
//	for (auto& i: m_mesh->get_boundaryNodes() ) painter->drawText( x_coords[i.first]* tr.m11() +shift_x, y_coords[i.first]* tr.m11() +shift_y , QString::number(i.first) );


}


template<typename T>
QRectF QGraphicsItem_Mesh2D_Vertex<T>::boundingRect() const { return QRectF(-sizee, -sizee, 2*sizee, 2*sizee); }


template<typename T>
void QGraphicsItem_Mesh2D_Vertex<T>::paint(QPainter* painter, const QStyleOptionGraphicsItem* styleOptions, QWidget* widget)
{

}


template<typename T>
QRectF QGraphicsItem_Mesh2D_VertexLabel<T>::boundingRect() const { return QRectF(-sizee, -sizee, 2*sizee, 2*sizee); }


template<typename T>
void QGraphicsItem_Mesh2D_VertexLabel<T>::paint(QPainter* painter, const QStyleOptionGraphicsItem* styleOptions, QWidget* widget)
{

}



template<typename T>
QRectF QGraphicsItem_Mesh2D_badElements<T>::boundingRect() const { return QRectF(-sizee, -sizee, 2*sizee, 2*sizee); }


template<typename T>
void QGraphicsItem_Mesh2D_badElements<T>::paint(QPainter* painter, const QStyleOptionGraphicsItem* styleOptions, QWidget* widget)
{
	T* x_coords = m_mesh->get_component(0);
	T* y_coords = m_mesh->get_component(1);

	painter->setPen( QPen( Qt::black, 0.01 ) );
	painter->setBrush( QBrush( Qt::red, Qt::SolidPattern ) );
	const size_t numOfVertices(3);
	Connectivities<size_t, numOfVertices>& conn = *m_mesh->get_connectivities();

	for ( size_t elementIndex = 0; elementIndex < m_numOfBadElements; elementIndex++ ) {
		size_t* elementVertices = conn[m_badElements[elementIndex]];
		QPointF points[numOfVertices];
		for (size_t vertexIndex = 0; vertexIndex < numOfVertices; vertexIndex++ ) {
			points[vertexIndex] = QPointF(x_coords[elementVertices[vertexIndex]], y_coords[elementVertices[vertexIndex]]);
		}
		painter->drawConvexPolygon(points, numOfVertices);
	}
}



template<typename T>
QRectF QGraphicsItem_Mesh2D_ElementsQuality<T>::boundingRect() const { return QRectF(-sizee, -sizee, 2*sizee, 2*sizee); }


template<typename T>
void QGraphicsItem_Mesh2D_ElementsQuality<T>::paint(QPainter* painter, const QStyleOptionGraphicsItem* styleOptions, QWidget* widget)
{
	T* x_coords = m_mesh->get_component(0);
	T* y_coords = m_mesh->get_component(1);

	painter->setPen( QPen( Qt::black, 0.01 ) );
	const size_t numOfVertices(3);
	Connectivities<size_t, numOfVertices>& conn = *m_mesh->get_connectivities();
	ContainerBase<T, std::vector< T > >& qualityMetrics = m_mesh->get_qualityMetrics()[0];

	for ( size_t elementIndex = 0; elementIndex < conn.get_numOfConnectivities(); elementIndex++ ) {
		size_t* elementVertices = conn[elementIndex];
		QPointF points[numOfVertices];
		T metricValue = qualityMetrics[ elementIndex ];
		if ( 0 <= metricValue &&  metricValue <= T(0.25) ) painter->setBrush( QBrush( Qt::green, Qt::SolidPattern ) );
		else if ( 0.25 <= metricValue &&  metricValue <= T(0.5) ) painter->setBrush( QBrush( Qt::cyan, Qt::SolidPattern ) );
		else if ( 0.5 < metricValue &&  metricValue <= T(0.8) ) painter->setBrush( QBrush( Qt::yellow, Qt::SolidPattern ) );
		else if ( 0.8 < metricValue &&  metricValue <= T(0.95) ) painter->setBrush( QBrush( Qt::magenta, Qt::SolidPattern ) );
		else if ( 0.95 < metricValue &&  metricValue <= T(0.99) ) painter->setBrush( QBrush( Qt::red, Qt::SolidPattern ) );
		else if ( 0.99 < metricValue &&  metricValue <= T(1) ) painter->setBrush( QBrush( Qt::darkRed, Qt::SolidPattern ) );
		for (size_t vertexIndex = 0; vertexIndex < numOfVertices; vertexIndex++ ) {
			points[vertexIndex] = QPointF(x_coords[elementVertices[vertexIndex]], y_coords[elementVertices[vertexIndex]]);
		}
		painter->drawConvexPolygon(points, numOfVertices);
	}
}



template<typename T>
QGraphicsItem_Mesh2D_fieldPlot<T>::QGraphicsItem_Mesh2D_fieldPlot( Qt3DExtras::Qt3DWindow *parent, Mesh<T, 3, int>* mesh, T* fieldValues ): m_mesh(mesh), m_fieldValues( fieldValues )
{
	T* x_coords = m_mesh->get_component(0);
	T* y_coords = m_mesh->get_component(1);
	size_t numOfNodes = m_mesh->get_numOfNodes();


	const size_t numOfVertices(3);
	size_t* connectivities = m_mesh->get_connectivities()->data();
	size_t numOfConnectivities = m_mesh->get_numOfConnectivities();


	Qt3DCore::QEntity	*rootEntity 	= new Qt3DCore::QEntity( );
	this->setRootEntity( rootEntity );


	// Plot UCS axes
	UCSAxesScene axes( rootEntity );



	// Set the data buffers
	QByteArray m_vertexBufferData;
	QByteArray m_normalsBufferData;
	QByteArray m_colorBufferData;
	QByteArray m_indexBufferData;


        m_vertexBufferData.resize( numOfNodes*3*sizeof(T));
        T *positions = reinterpret_cast<T*>(m_vertexBufferData.data());

	m_normalsBufferData.resize( numOfNodes*3*sizeof(T));
        T *normals = reinterpret_cast<T*>(m_normalsBufferData.data());

        m_colorBufferData.resize( numOfNodes*3*sizeof(T));
        float *colors = reinterpret_cast<float*>(m_colorBufferData.data());

        m_indexBufferData.resize( numOfConnectivities*3*sizeof(unsigned int));
        unsigned int *indices = reinterpret_cast<unsigned int*>(m_indexBufferData.data());


        auto *vertexBuffer      = new Qt3DRender::QBuffer       ( Qt3DRender::QBuffer::VertexBuffer);
        auto *normalsBuffer     = new Qt3DRender::QBuffer       ( Qt3DRender::QBuffer::VertexBuffer);
        auto *colorBuffer       = new Qt3DRender::QBuffer       ( Qt3DRender::QBuffer::VertexBuffer);
        auto *indexBuffer       = new Qt3DRender::QBuffer       ( Qt3DRender::QBuffer::IndexBuffer);
        vertexBuffer->setData( m_vertexBufferData);
        normalsBuffer->setData( m_normalsBufferData);
        colorBuffer->setData( m_colorBufferData);
        indexBuffer->setData( m_indexBufferData);

	
	// ==================== Define data and settings
        for ( size_t nodeI = 0; nodeI < numOfNodes; nodeI++ ) {
                *positions++ = x_coords[nodeI];
                *positions++ = y_coords[nodeI];
                *positions++ = T(0);
        }

	// FieldValues
	float val;
	auto maxValue = std::max_element( m_fieldValues, m_fieldValues +numOfNodes );
	auto minValue = std::min_element( m_fieldValues, m_fieldValues +numOfNodes );
	float map1 = float(1)/( float(*minValue) ); if ( *minValue >= 0 ) map1 = 0;
	float map2 = float(1)/( float(*maxValue) ); if ( *maxValue <= 0 ) map2 = 0;
	for ( size_t nodeI = 0; nodeI < numOfNodes; nodeI++ ) {
		if ( -10e-16 < float(m_fieldValues[nodeI]) && float(m_fieldValues[nodeI]) < 10e-16 ) { *colors++ = 1.f; *colors++ = 1.f; *colors++ = 0.f; continue; } 
		if ( float(m_fieldValues[nodeI]) < 0 ) { // zero: yellow, min: green, yellow rgb: 1 1 0
			val = map1*float(m_fieldValues[nodeI]);
			*colors++ = 1.f -val;
                	*colors++ = 1.f;
                	*colors++ = 0.f; }
		else {
			val = map2*float(m_fieldValues[nodeI]);
		 	*colors++ = 1.f; // zero: yellow, max: red
                	*colors++ = 1.f -val;
                	*colors++ = 0.f; }
        }

        for ( size_t nodeI = 0; nodeI < numOfConnectivities*3; nodeI++ )
                *indices++ = connectivities[nodeI];



	std::set< std::pair<size_t, size_t >, unorderLess<size_t> >& uniqueEdges = m_mesh->get_connectivities()->get_uniqueLinks();
	size_t numOfLineConnectivities = uniqueEdges.size();


	QByteArray		m_indexBufferData_l;
        m_indexBufferData_l.resize( numOfLineConnectivities*2*sizeof(unsigned int));
        unsigned int *indices_l	= reinterpret_cast<unsigned int*>( m_indexBufferData_l.data() );
        auto *indexBuffer_l     = new Qt3DRender::QBuffer       ( Qt3DRender::QBuffer::IndexBuffer);
        indexBuffer_l->setData( m_indexBufferData_l );

	for ( auto& edgePoints: uniqueEdges ) {
                *indices_l++ = edgePoints.first;
                *indices_l++ = edgePoints.second;
	}



	// Plot wireframe
/*	auto *lineEntity        = new Qt3DCore::QEntity         ( rootEntity );
        auto *line              = new Qt3DRender::QGeometryRenderer( lineEntity );
        auto *geometry_l        = new Qt3DRender::QGeometry     ( line );
        auto *positionAttribute_l = new Qt3DRender::QAttribute    ( geometry_l );
        auto *indexAttribute_l    = new Qt3DRender::QAttribute    ( geometry_l );
	Qt3DExtras::QPhongMaterial *material_l = new Qt3DExtras::QPhongMaterial( rootEntity);

        geometry_l->addAttribute(positionAttribute_l);
        geometry_l->addAttribute(indexAttribute_l);
        line->setGeometry( geometry_l );
        line->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);
        lineEntity->addComponent( line );
        lineEntity->addComponent( material_l );

        positionAttribute_l->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
        positionAttribute_l->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
        positionAttribute_l->setVertexBaseType(Qt3DRender::QAttribute::Float);
        positionAttribute_l->setVertexSize(3);
        positionAttribute_l->setBuffer(vertexBuffer);
        positionAttribute_l->setByteStride(3 * sizeof(float));
        positionAttribute_l->setCount( numOfNodes );

        indexAttribute_l->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
        indexAttribute_l->setVertexBaseType(Qt3DRender::QAttribute::UnsignedInt);
        indexAttribute_l->setBuffer( indexBuffer_l );
        indexAttribute_l->setCount( numOfLineConnectivities*2 );

	material_l->setAmbient( QColor(Qt::black) );
	material_l->setDiffuse( QColor(Qt::black) );
	material_l->setShininess( 0.f );
	material_l->setSpecular( QColor( Qt::black ) );
*/


	// Plot faces
	auto *elementEntity	= new Qt3DCore::QEntity         ( rootEntity );
	auto *element		= new Qt3DRender::QGeometryRenderer( elementEntity );
	auto *geometry		= new Qt3DRender::QGeometry     ( element );
	auto *positionAttribute = new Qt3DRender::QAttribute    ( geometry );
	auto *colorAttribute	= new Qt3DRender::QAttribute    ( geometry );
	auto *indexAttribute	= new Qt3DRender::QAttribute    ( geometry );
	Qt3DRender::QMaterial *material = new Qt3DExtras::QPerVertexColorMaterial( rootEntity);


        geometry->addAttribute(positionAttribute);
	geometry->addAttribute(colorAttribute);
        geometry->addAttribute(indexAttribute);
        element->setGeometry(geometry);
        element->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
        elementEntity->addComponent(element);
        elementEntity->addComponent(material);

        positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
        positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
        positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
        positionAttribute->setVertexSize(3);
        positionAttribute->setBuffer(vertexBuffer);
        positionAttribute->setByteStride(3 * sizeof(float));
        positionAttribute->setCount( numOfNodes );

        colorAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
        colorAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
        colorAttribute->setName(Qt3DRender::QAttribute::defaultColorAttributeName());
     	colorAttribute->setVertexSize(3);
        colorAttribute->setBuffer(colorBuffer);
        colorAttribute->setByteStride(3 * sizeof(float));
        colorAttribute->setCount( numOfNodes );

        indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
        indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedInt);
        indexAttribute->setBuffer( indexBuffer );
        indexAttribute->setCount( numOfConnectivities*3 );



	// Camera settings
      	Qt3DRender::QCamera *camera = this->camera();  
//    	camera->setPosition(QVector3D(0.f, 0.f, -50.f));
    	camera->setPosition(QVector3D(0.f, 0.f, 50.f));
//	camera->setUpVector(QVector3D(0.f, -1.f, 0.f));
	camera->setUpVector(QVector3D(0.f, 1.f, 0.f));
	camera->setViewCenter(QVector3D(0, 0, 0));
	Qt3DExtras::QOrbitCameraController *camController = new Qt3DExtras::QOrbitCameraController(rootEntity);
      	camController->setLinearSpeed( 200.0f );
       	camController->setLookSpeed( 500.0f );
       	camController->setCamera(camera);


}



void QGraphicsView_mod::wheelEvent ( QWheelEvent * event )
{
	double scaleFactor=1.15;
	if ( event->delta()<0) scaleFactor = 1/scaleFactor;

	scale(scaleFactor, scaleFactor);
}

