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

#ifndef MYQARROWHEADER
#define MYQARROWHEADER

#include <QApplication>
#include <QtWidgets/QApplication>
#include <QGuiApplication>
#include <Qt3DCore/QEntity>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DCore/QTransform>

class myQArrow: public QObject
{
	private:
		Qt3DCore::QEntity*	m_rootEntity;
		Qt3DCore::QEntity*	m_coneEntity;
		Qt3DCore::QEntity*	m_cylinderEntity;
		QVector3D		m_start;
		QVector3D		m_end;
		float			m_scale;
		QColor			m_color;


	public:
		explicit myQArrow( Qt3DCore::QEntity *rootEntity, QVector3D start, QVector3D end, float scale = 1, QColor color = QColor(Qt::red) );
		~myQArrow() {}

};

#endif
