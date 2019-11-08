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

#ifndef UCSAXESLINESHEADER
#define UCSAXESLINESHEADER
#include"myQLine.h"
#include"myQArrow.h"


class UCSAxesLinesScene : public QObject
{

private:
	Qt3DCore::QEntity*	m_rootEntity;
	myQLine*		m_axis_x;
	myQLine*		m_axis_y;
	myQLine*		m_axis_z;


public:
	explicit UCSAxesLinesScene(Qt3DCore::QEntity *rootEntity): m_rootEntity(rootEntity)
	{
		m_axis_x = new myQLine( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(1.f, 0.f, 0.f), QColor(Qt::red) );
		m_axis_y = new myQLine( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 1.f, 0.f), QColor(Qt::green) );
		m_axis_z = new myQLine( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 0.f, 1.f), QColor(Qt::blue) );
	}
	~UCSAxesLinesScene() {}


public slots:
void enableDraw(bool enabled)
{
    //m_axis_x_->setEnabled(enabled);
}

};



class UCSAxesScene : public QObject
{

private:
	Qt3DCore::QEntity*	m_rootEntity;
	myQArrow*		m_axis_x;
	myQArrow*		m_axis_y;
	myQArrow*		m_axis_z;


public:
	explicit UCSAxesScene(Qt3DCore::QEntity *rootEntity): m_rootEntity(rootEntity)
	{
		m_axis_x = new myQArrow( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(1.f, 0.f, 0.f), 1.f, QColor(Qt::red) );
		m_axis_y = new myQArrow( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 1.f, 0.f), 1.f, QColor(Qt::green) );
		m_axis_z = new myQArrow( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 0.f, 1.f), 1.f, QColor(Qt::blue) );
	}
	~UCSAxesScene() {}


public slots:
void enableDraw(bool enabled)
{
    //m_axis_x_->setEnabled(enabled);
}

};



#endif
