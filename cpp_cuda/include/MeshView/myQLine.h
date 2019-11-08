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

#ifndef MYQLINEHEADER
#define MYQLINEHEADER


#include <QApplication>
#include <QtWidgets/QApplication>
#include <QGuiApplication>
#include <Qt3DCore/QEntity>
#include <Qt3DExtras/QPhongMaterial>

#include <Qt3DCore/QTransform>

class myQLine : public QObject
{

private:
	Qt3DCore::QEntity*	m_rootEntity;
	Qt3DCore::QEntity*	m_lineEntity;
	QVector3D		m_start;
	QVector3D		m_end;
	QByteArray		m_vertexBufferData;
	float*			m_positions;
	QByteArray		m_indexBufferData;
	unsigned int*		m_indices;
	QColor			m_materialColor;
	Qt3DExtras::QPhongMaterial* m_material;	
	QVector3D		m_color_start;
	QVector3D		m_color_end;
	float*			m_colors;


public:
	explicit myQLine(Qt3DCore::QEntity *rootEntity, QVector3D start, QVector3D end, QColor color = QColor(Qt::red) ); 
	~myQLine() {}


public slots:
	void setPositions( QVector3D start, QVector3D end )
	{
		m_start = start; m_end = end;
		*m_positions++ = start.x(); 	*m_positions++ = start.y(); 	*m_positions++ = start.z();
		*m_positions++ = end.x(); 	*m_positions++ = end.y();	*m_positions++ = end.z();
	}

	void setColors( QVector3D start_c, QVector3D end_c )
	{
		m_color_start = start_c; m_color_end = end_c;
		*m_colors++ = start_c.x(); 	*m_colors++ = start_c.y(); 	*m_colors++ = start_c.z();
		*m_colors++ = end_c.x(); 	*m_colors++ = end_c.y();	*m_colors++ = end_c.z();
	}

	void setMaterialColor( QColor color )
	{
		m_materialColor = color;
		m_material->setAmbient( QColor( m_materialColor ) );
	}

	void enableDraw(bool enabled) { m_lineEntity->setEnabled(enabled); }

};



#endif
