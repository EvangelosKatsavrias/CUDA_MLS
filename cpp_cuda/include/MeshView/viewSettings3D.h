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

#ifndef VIEWSETTINGS3D
#define VIEWSETTINGS3D


#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DCore/QEntity>
#include <Qt3DRender/QCamera>
#include <iostream>


class view3DSettings: public QObject
{

public:
	explicit view3DSettings( Qt3DCore::QEntity* rootEntity, Qt3DExtras::Qt3DWindow* view );
	~view3DSettings() {}

	void set_sceneSize( float sceneSize ) { m_sceneSize = sceneSize; set_cameraDefaults(); }
	void set_backgroundColor( QColor color );
	void set_cameraDefaults ( );
	void set_lightingDefaults ( );


	void set_viewPlane ( int planeIndex ) {
		Qt3DRender::QCamera *camera = m_view->camera();
		switch (planeIndex) {
			case 0: 
				camera->setPosition(QVector3D(0.f, 0.f, m_sceneSize));
				camera->setUpVector(QVector3D(0, 0, 1));
				camera->setViewCenter(QVector3D(0, 0, 0));
				break;

			case 1: 
				camera->setPosition(QVector3D(0.f, m_sceneSize, 0.f));
				camera->setUpVector(QVector3D(0, 0, 1));
				camera->setViewCenter(QVector3D(0, 0, 0));
				break;

			case 2: 
				camera->setPosition(QVector3D(m_sceneSize, 0.f, 0.f));
				camera->setUpVector(QVector3D(0, 0, 1));
				camera->setViewCenter(QVector3D(0, 0, 0));
				break;
			case 3:
				camera->setPosition(QVector3D(m_sceneSize, m_sceneSize, m_sceneSize));
				camera->setUpVector(QVector3D(0, 0, 1));
				camera->setViewCenter(QVector3D(0, 0, 0));
				break;
			case 4:
				camera->setPosition(QVector3D(-m_sceneSize, m_sceneSize, m_sceneSize));
				camera->setUpVector(QVector3D(0, 0, 1));
				camera->setViewCenter(QVector3D(0, 0, 0));
				break;
			case 5:
				camera->setPosition(QVector3D(-m_sceneSize, -m_sceneSize, m_sceneSize));
				camera->setUpVector(QVector3D(0, 0, 1));
				camera->setViewCenter(QVector3D(0, 0, 0));
				break;
			case 6:
				camera->setPosition(QVector3D(m_sceneSize, -m_sceneSize, m_sceneSize));
				camera->setUpVector(QVector3D(0, 0, 1));
				camera->setViewCenter(QVector3D(0, 0, 0));
				break;
		}
	}

private:
	Qt3DExtras::Qt3DWindow*	m_view;
	Qt3DCore::QEntity*	m_rootEntity;
	float			m_sceneSize = 10;

};


#endif
