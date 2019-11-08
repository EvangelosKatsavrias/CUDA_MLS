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

#ifndef MYQMESHRENDEREDHEADER
#define MYQMESHRENDEREDHEADER

#include <QApplication>
#include <Qt3DCore/QEntity>
#include <Qt3DRender/QMesh>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QPerVertexColorMaterial>
#include"myAlgorithms.h"
#include<set>
#include"Mesh.h"




template<class T>
class myQMeshRendered : public QObject
{

private:
	Qt3DCore::QEntity*	m_rootEntity;
	Qt3DCore::QEntity*	m_lineEntity;
	QByteArray		m_vertexBufferData;
	QByteArray		m_colorBufferData;
	T*			m_positions;
	QByteArray		m_indexBufferData;
	unsigned int*		m_indices;
	QColor			m_materialColor;
	Qt3DExtras::QPerVertexColorMaterial* m_material;
	Mesh<T, 3, int>*		m_mesh;


public:
	explicit myQMeshRendered(Qt3DCore::QEntity *rootEntity, Mesh<T, 3, int>* mesh, QColor color = QColor(Qt::red) ); 
	~myQMeshRendered() {}


public slots:
	void enableDraw(bool enabled) { m_lineEntity->setEnabled(enabled); }

};

template class myQMeshRendered<float>;
template class myQMeshRendered<double>;

#endif
