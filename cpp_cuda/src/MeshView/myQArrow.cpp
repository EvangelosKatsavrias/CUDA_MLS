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

#include "myQArrow.h"
#include <Qt3DRender/QGeometry>
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QMesh>
#include <Qt3DExtras/QCylinderMesh>
#include <Qt3DExtras/QConeMesh>

#include"myBLAS.h"
#include"math.h"
#include<iostream>


myQArrow::myQArrow( Qt3DCore::QEntity *rootEntity, QVector3D start, QVector3D end, float scale, QColor color): m_rootEntity(rootEntity),  m_start(start), m_end(end), m_scale(scale), m_color(color)
{

	QVector3D vec = end -start;
	float length = vec.length();

	auto vec_norm = vec/length;
	QQuaternion quat; quat = quat.fromAxisAndAngle( QVector3D::crossProduct( QVector3D( 0.f, 1.f, 0.f ), vec_norm ), acos(QVector3D::dotProduct( vec_norm, QVector3D( 0.f, 1.f, 0.f ) ) )*180/3.14159265 );

	Qt3DCore::QTransform cylinderTrans, cylinderTrans2;
	Qt3DCore::QTransform coneTrans, coneTrans2;

	cylinderTrans.setRotation(quat);
	coneTrans.setRotation(quat);

	cylinderTrans2.setTranslation( QVector3D( 0.f, 0.35f, 0.f) );
	cylinderTrans.setMatrix( cylinderTrans.matrix() *cylinderTrans2.matrix() );

	coneTrans2.setTranslation( QVector3D(0.f, 0.85f, 0.f) );
	coneTrans.setMatrix( coneTrans.matrix() *coneTrans2.matrix() );


	// Cylinder shape data
	m_cylinderEntity = new Qt3DCore::QEntity(m_rootEntity);
	Qt3DExtras::QPhongMaterial *cylinderMaterial = new Qt3DExtras::QPhongMaterial();

	Qt3DExtras::QCylinderMesh *cylinder = new Qt3DExtras::QCylinderMesh();
	cylinder->setRadius( 0.05 );
	cylinder->setLength( 0.7 );
	cylinder->setRings(100);
	cylinder->setSlices(20);
	
	// CylinderMesh Transform
	Qt3DCore::QTransform *cylinderTransform = new Qt3DCore::QTransform();
	cylinderTransform->setMatrix( cylinderTrans.matrix() );
	
	cylinderMaterial->setDiffuse( m_color );
	
	// Cylinder
	m_cylinderEntity->addComponent(cylinder);
	m_cylinderEntity->addComponent(cylinderMaterial);
	m_cylinderEntity->addComponent(cylinderTransform);



	// Cone shape data
	m_coneEntity = new Qt3DCore::QEntity(m_rootEntity);
	Qt3DExtras::QPhongMaterial *coneMaterial = new Qt3DExtras::QPhongMaterial();


	Qt3DExtras::QConeMesh *cone = new Qt3DExtras::QConeMesh( m_coneEntity );
	cone->setTopRadius(0);
	cone->setBottomRadius( 0.1 );
	cone->setLength( 0.3 );
	cone->setRings(50);
	cone->setSlices(20);
	
	// ConeMesh Transform
	Qt3DCore::QTransform *coneTransform = new Qt3DCore::QTransform();
	coneTransform->setMatrix( coneTrans.matrix() );

	coneMaterial->setDiffuse( m_color );
	
	// Cone
	m_coneEntity->addComponent(cone);
	m_coneEntity->addComponent(coneMaterial);
	m_coneEntity->addComponent(coneTransform);

}



