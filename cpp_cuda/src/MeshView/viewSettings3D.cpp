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

#include "viewSettings3D.h"
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DRender/QCamera>
#include <Qt3DRender/QCameraLens>
#include <QFirstPersonCameraController>
#include <QOrbitCameraController>
#include <Qt3DRender/QPointLight>
#include <Qt3DCore/QTransform>


view3DSettings::view3DSettings( Qt3DCore::QEntity* rootEntity, Qt3DExtras::Qt3DWindow* view ): m_rootEntity(rootEntity), m_view(view)
{
       	m_view->setRootEntity(m_rootEntity);
	set_backgroundColor( QColor(QRgb( 0x4d4d4f )) );
	set_cameraDefaults ( ); 
	set_lightingDefaults ( );


//	QWidget *container = QWidget::createWindowContainer(view);
//	QSize screenSize = m_view->screen()->size();
//	container->setMinimumSize(QSize(200, 100));
//	container->setMaximumSize(screenSize);

//	Qt3DInput::QInputAspect *input = new Qt3DInput::QInputAspect( m_rootEntity );
//	m_view->registerAspect(input);
//	input->setCamera( camera );

}


void view3DSettings::set_backgroundColor( QColor color ) 
{
	m_view->defaultFrameGraph()->setClearColor( color );
}


void view3DSettings::set_cameraDefaults ( ) 
{
      	Qt3DRender::QCamera *camera = m_view->camera();  
	camera->lens()->setPerspectiveProjection(45.0f, 16.0f/9.0f, 0.1f, 1000.0f);
    	camera->setPosition(QVector3D( m_sceneSize, 0, 0.0f));
	camera->setUpVector(QVector3D(0, 0, 1));
	camera->setViewCenter(QVector3D(0, 0, 0));

    	// For camera controls
       	Qt3DExtras::QOrbitCameraController *camController = new Qt3DExtras::QOrbitCameraController(m_rootEntity);
//	Qt3DExtras::QFirstPersonCameraController *camController = new Qt3DExtras::QFirstPersonCameraController(m_rootEntity);
      	camController->setLinearSpeed( 500.0f );
       	camController->setLookSpeed( 2000.0f );   
       	camController->setCamera(camera);
}


void view3DSettings::set_lightingDefaults ( ) 
{
      	Qt3DRender::QCamera *camera = m_view->camera();  
	// lighting
	Qt3DCore::QEntity *lightEntity = new Qt3DCore::QEntity(m_rootEntity);
	Qt3DCore::QEntity *lightEntity2 = new Qt3DCore::QEntity(m_rootEntity);
	Qt3DRender::QPointLight *light = new Qt3DRender::QPointLight(lightEntity);
	Qt3DRender::QPointLight *light2 = new Qt3DRender::QPointLight(lightEntity2);
	light->setColor("white");
	light->setIntensity(1.);
	Qt3DCore::QTransform *lightTransform = new Qt3DCore::QTransform(lightEntity);
	Qt3DCore::QTransform *lightTransform2 = new Qt3DCore::QTransform(lightEntity2);
	lightTransform->setTranslation(camera->position());
	lightTransform->setTranslation( QVector3D(0, 0, 20) );
	lightTransform2->setTranslation( QVector3D(0, 0, -20) );
	lightEntity->addComponent(light);
	lightEntity->addComponent(lightTransform);
	lightEntity2->addComponent(light2);
	lightEntity2->addComponent(lightTransform2);
}

