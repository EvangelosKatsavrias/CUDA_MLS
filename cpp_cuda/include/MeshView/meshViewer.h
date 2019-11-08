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

#ifndef MESHVIEWERHEADER
#define MESHVIEWERHEADER

#include <QtCore>
#include <QtGui>
#include <QtWidgets>
#include "App_MeshAdapter.h"
#include <utility>
#include <Qt3DExtras/Qt3DWindow>
#include"myAlgorithms.h"



template<typename T>
class QGraphicsItem_Mesh2D_Edges: public QGraphicsItem
{

protected:
	Mesh<T, 3, int>*						m_mesh;

	virtual QRectF boundingRect() const; virtual void paint(QPainter *,
			const QStyleOptionGraphicsItem*, QWidget*);

public: QGraphicsItem_Mesh2D_Edges(QGraphicsItem *parent = 0, Mesh<T, 3, int>* mesh=0): QGraphicsItem(parent), m_mesh(mesh) 
	{
	}

}; template class QGraphicsItem_Mesh2D_Edges<float>; template class
QGraphicsItem_Mesh2D_Edges<double>;


template<typename T> class QGraphicsItem_Mesh2D_Vertices: public QGraphicsItem
{

protected: Mesh<T, 3, int>*	m_mesh;

	virtual QRectF boundingRect() const; virtual void paint(QPainter *,
			const QStyleOptionGraphicsItem*, QWidget*);

public: QGraphicsItem_Mesh2D_Vertices(QGraphicsItem *parent = 0, Mesh<T, 3,
			int>* mesh=0): QGraphicsItem(parent), m_mesh(mesh) { }

}; template class QGraphicsItem_Mesh2D_Vertices<float>; template class
QGraphicsItem_Mesh2D_Vertices<double>;

template<typename T> class QGraphicsItem_Mesh2D_Vertex: public QGraphicsItem {

protected: Mesh<T, 3, int>*	m_mesh; int			m_pointIndex;

	virtual QRectF boundingRect() const; virtual void paint(QPainter *,
			const QStyleOptionGraphicsItem*, QWidget*);

public: QGraphicsItem_Mesh2D_Vertex(QGraphicsItem *parent = 0, Mesh<T, 3, int>*
			mesh=0): QGraphicsItem(parent), m_mesh(mesh) { }

}; template class QGraphicsItem_Mesh2D_Vertex<float>; template class
QGraphicsItem_Mesh2D_Vertex<double>;



template<typename T> class QGraphicsItem_Mesh2D_VertexLabel: public
							     QGraphicsItem {

protected: Mesh<T, 3, int>*	m_mesh;

	virtual QRectF boundingRect() const; virtual void paint(QPainter *,
			const QStyleOptionGraphicsItem*, QWidget*);

public: QGraphicsItem_Mesh2D_VertexLabel(QGraphicsItem *parent = 0, Mesh<T, 3,
			int>* mesh=0): QGraphicsItem(parent), m_mesh(mesh) { }

}; template class QGraphicsItem_Mesh2D_VertexLabel<float>; template class
QGraphicsItem_Mesh2D_VertexLabel<double>;



template<typename T> class QGraphicsItem_Mesh2D_badElements: public
							     QGraphicsItem {

protected: Mesh<T, 3, int>*	m_mesh; size_t*			m_badElements;
	   size_t			m_numOfBadElements;


	virtual QRectF boundingRect() const; virtual void paint(QPainter *,
			const QStyleOptionGraphicsItem*, QWidget*);

public: QGraphicsItem_Mesh2D_badElements(QGraphicsItem *parent = 0, Mesh<T, 3,
			int>* mesh=0, size_t numOfBadElements=0, size_t*
			badElements=0): QGraphicsItem(parent), m_mesh(mesh),
	m_numOfBadElements(numOfBadElements), m_badElements(badElements) { }

}; template class QGraphicsItem_Mesh2D_badElements<float>; template class
QGraphicsItem_Mesh2D_badElements<double>;



template<typename T> class QGraphicsItem_Mesh2D_ElementsQuality: public QGraphicsItem 
{

protected: Mesh<T, 3, int>*		m_mesh; 

	virtual QRectF boundingRect() const; virtual void paint(QPainter *,
			const QStyleOptionGraphicsItem*, QWidget*);

public: QGraphicsItem_Mesh2D_ElementsQuality(QGraphicsItem *parent = 0, Mesh<T, 3, int>* mesh=0): 
	QGraphicsItem(parent), m_mesh(mesh) { }

}; 
template class QGraphicsItem_Mesh2D_ElementsQuality<float>; 
template class QGraphicsItem_Mesh2D_ElementsQuality<double>;



class QGraphicsView_mod: public QGraphicsView
{


public:
	QGraphicsView_mod() 
	{
		setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
		setDragMode(QGraphicsView::ScrollHandDrag ); }
	virtual void wheelEvent(QWheelEvent *);

};


template<typename T>
class QGraphicsItem_Mesh2D_fieldPlot: public Qt3DExtras::Qt3DWindow
{

protected:
	Mesh<T, 3, int>*						m_mesh;
	T*								m_fieldValues;


public:
	QGraphicsItem_Mesh2D_fieldPlot(Qt3DExtras::Qt3DWindow *parent, Mesh<T, 3, int>* mesh, T* fieldValues );

};
template class QGraphicsItem_Mesh2D_fieldPlot<float>;
template class QGraphicsItem_Mesh2D_fieldPlot<double>;




//	void set_title(const std::string str) { m_title.insert( 0, QString::fromStdString(str.data() ) ); }

#endif
