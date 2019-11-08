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

#include <QApplication>
#include <QGuiApplication>

#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QCommandLinkButton>
#include <QtGui/QScreen>


#include <Qt3DCore/QEntity>
#include <Qt3DExtras/Qt3DWindow>


#include"App_MeshAdapter.h"
#include"viewSettings3D.h"
#include"myAlgorithms.h"
#include"UCSAxesLines.h"
#include"myQMesh.h"
#include"myQMeshRendered.h"

#include"myFileMenu_gui.h"
#include <iostream>


int app3Dtest2(  int argc, char *argv[]  )
{
	int argc1(1); char *argv1[1]; char title[16]={"Mesh adaptation"}; argv1[0] = title;
	QApplication app(argc1, argv1);


	std::string fileName = argv[1];

	app_meshAdapter<float> app1;
	app1.readProjectSettings( fileName );
	app1.set_initialMesh(  );
	app1.set_adaptedMesh(  );
	app1.set_meshAdapter();
	app1.write_nodesInfo ( std::cout );
	app1.adaptMesh();
	app1.write_adaptationProcessInfo( std::cout );



/*	
	app_meshAdapter<float> app1;
	app1.readProjectSettings( "Alfa147/alfa147.msh" );
	app1.set_initialMesh(  );
//	Isless_equal<float> fun( -80.f );
//	IsInRange<float> fun2( 2.f, -2.f );
//	app1.get_initialMesh()->set_nodeFlagWithCriterionOnComponentValue( 1, fun, 1 );
//	app1.get_initialMesh()->set_nodeFlagWithCriterionOnComponentValue( 1, fun2, 2 );
	std::set<int> flags = app1.get_initialMesh()->get_flagValues();
	std::cout << " =============================================================> flag values: "; for ( auto& el: flags ) std::cout << el<< std::endl;
	app1.set_adaptedMesh(  );
	app1.set_meshAdapter();
	app1.write_nodesInfo ( std::cout );
	app1.adaptMesh();
	app1.write_adaptationProcessInfo( std::cout );
*/




	size_t numOfNodes = app1.get_initialMesh()->get_numOfNodes();
	float* xcoords = app1.get_initialMesh()->get_component(0);
	float* ycoords = app1.get_initialMesh()->get_component(1);
	float* zcoords = app1.get_initialMesh()->get_component(2);
	float* maxX = std::max_element( xcoords, xcoords+numOfNodes );
	float* maxY = std::max_element( ycoords, ycoords+numOfNodes );
	float* maxZ = std::max_element( zcoords, zcoords+numOfNodes );

	float maxSize = std::max(std::max( *maxX, *maxY ), *maxZ );

	// Setting 3d visualization
	Qt3DExtras::Qt3DWindow* view 		= new Qt3DExtras::Qt3DWindow();
	Qt3DCore::QEntity*	rootEntity 	= new Qt3DCore::QEntity();
	view3DSettings* 	viewSetup 	= new view3DSettings(rootEntity, view);
	viewSetup->set_sceneSize( 10.f*maxSize );
	UCSAxesScene*		axesLines	= new UCSAxesScene(rootEntity);	
	myQMesh<float>* qmesh_l = new myQMesh<float>( rootEntity, app1.get_initialMesh(), QColor(Qt::white) );
	myQMeshRendered<float>* qmesh = new myQMeshRendered<float>( rootEntity, app1.get_initialMesh(), QColor(Qt::white) );


	Qt3DExtras::Qt3DWindow* view2 		= new Qt3DExtras::Qt3DWindow();
	Qt3DCore::QEntity*	rootEntity2 	= new Qt3DCore::QEntity();
	view3DSettings* 	viewSetup2 	= new view3DSettings(rootEntity2, view2);
	viewSetup2->set_sceneSize( 10.f*maxSize );
	UCSAxesScene*		axesLines2	= new UCSAxesScene(rootEntity2);	
	myQMesh<float>* qmesh2_l = new myQMesh<float>( rootEntity2, app1.get_adaptedMesh(), QColor(Qt::white) );
	myQMeshRendered<float>* qmesh2 = new myQMeshRendered<float>( rootEntity2, app1.get_adaptedMesh(), QColor(Qt::white) );



	QWidget *container = QWidget::createWindowContainer(view);
	QWidget *container2 = QWidget::createWindowContainer(view2);
	
	QWidget *widget = new QWidget;
	QHBoxLayout *hLayout = new QHBoxLayout(widget);
	QVBoxLayout *vLayout = new QVBoxLayout();
	vLayout->setAlignment(Qt::AlignTop);
	hLayout->addWidget(container, 1);
	hLayout->addWidget(container2, 1);
	hLayout->addLayout(vLayout);
	widget->setWindowTitle(QStringLiteral("Mesh adaptation"));


	// Create control widgets
	QCommandLinkButton *info = new QCommandLinkButton();
	info->setText(QStringLiteral("Mesh adaptation Visualization"));
	info->setDescription(QString::fromLatin1("Initial mesh on the left side, adapted mesh on the right side."));
	info->setIconSize(QSize(0,0));


	QButtonGroup* viewButtons_ini = new QButtonGroup(widget);
	QRadioButton* viewXY_ini = new QRadioButton(widget);
	QRadioButton* viewXZ_ini = new QRadioButton(widget);
	QRadioButton* viewYZ_ini = new QRadioButton(widget);
	QRadioButton* viewIso_ini = new QRadioButton(widget);
	QRadioButton* viewIso1_ini = new QRadioButton(widget);
	QRadioButton* viewIso2_ini = new QRadioButton(widget);
	QRadioButton* viewIso3_ini = new QRadioButton(widget);
	viewButtons_ini->addButton( viewXY_ini, 0 );
	viewButtons_ini->addButton( viewXZ_ini, 1 );
	viewButtons_ini->addButton( viewYZ_ini, 2 );
	viewButtons_ini->addButton( viewIso_ini, 3 );
	viewButtons_ini->addButton( viewIso1_ini, 4 );
	viewButtons_ini->addButton( viewIso2_ini, 5 );
	viewButtons_ini->addButton( viewIso3_ini, 6 );
	viewYZ_ini->setChecked(true);

	viewXY_ini->setText( QStringLiteral("XY plane view") );
	viewXZ_ini->setText( QStringLiteral("XZ plane view") );
	viewYZ_ini->setText( QStringLiteral("YZ plane view") );
	viewIso_ini->setText( QStringLiteral("Isometric view 1") );
	viewIso1_ini->setText( QStringLiteral("Isometric view 2") );
	viewIso2_ini->setText( QStringLiteral("Isometric view 3") );
	viewIso3_ini->setText( QStringLiteral("Isometric view 4") );

	QObject::connect( viewButtons_ini, static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked), viewSetup, static_cast<void (view3DSettings::*)(int)>(&view3DSettings::set_viewPlane ) );


	QButtonGroup* viewButtons_ada = new QButtonGroup(widget);
	QRadioButton* viewXY_ada = new QRadioButton(widget);
	QRadioButton* viewXZ_ada = new QRadioButton(widget);
	QRadioButton* viewYZ_ada = new QRadioButton(widget);
	QRadioButton* viewIso_ada = new QRadioButton(widget);
	QRadioButton* viewIso1_ada = new QRadioButton(widget);
	QRadioButton* viewIso2_ada = new QRadioButton(widget);
	QRadioButton* viewIso3_ada = new QRadioButton(widget);
	viewButtons_ada->addButton( viewXY_ada, 0 );
	viewButtons_ada->addButton( viewXZ_ada, 1 );
	viewButtons_ada->addButton( viewYZ_ada, 2 );
	viewButtons_ada->addButton( viewIso_ada, 3 );
	viewButtons_ada->addButton( viewIso1_ada, 4 );
	viewButtons_ada->addButton( viewIso2_ada, 5 );
	viewButtons_ada->addButton( viewIso3_ada, 6 );
	viewYZ_ada->setChecked(true);

	viewXY_ada->setText( QStringLiteral("XY plane view") );
	viewXZ_ada->setText( QStringLiteral("XZ plane view") );
	viewYZ_ada->setText( QStringLiteral("YZ plane view") );
	viewIso_ada->setText( QStringLiteral("Isometric view 1") );
	viewIso1_ada->setText( QStringLiteral("Isometric view 2") );
	viewIso2_ada->setText( QStringLiteral("Isometric view 3") );
	viewIso3_ada->setText( QStringLiteral("Isometric view 4") );

	QObject::connect( viewButtons_ada, static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked), viewSetup2, static_cast<void (view3DSettings::*)(int)>(&view3DSettings::set_viewPlane ) );




	QLabel* viewSettings_ini = new QLabel( "Initial mesh view options", widget );
	QLabel* viewSettings_ada = new QLabel( "Adapted mesh view options", widget );
	QLabel* blank = new QLabel( "", widget );

	vLayout->addWidget(info);
	vLayout->addWidget( viewSettings_ini );
	vLayout->addWidget( viewXY_ini );
	vLayout->addWidget( viewXZ_ini );
	vLayout->addWidget( viewYZ_ini );
	vLayout->addWidget( viewIso_ini );
	vLayout->addWidget( viewIso1_ini );
	vLayout->addWidget( viewIso2_ini );
	vLayout->addWidget( viewIso3_ini );
	vLayout->addWidget( blank );
	vLayout->addWidget( viewSettings_ada );
	vLayout->addWidget( viewXY_ada );
	vLayout->addWidget( viewXZ_ada );
	vLayout->addWidget( viewYZ_ada );
	vLayout->addWidget( viewIso_ada );
	vLayout->addWidget( viewIso1_ada );
	vLayout->addWidget( viewIso2_ada );
	vLayout->addWidget( viewIso3_ada );
	



	QMenu 		*fileMenu 	= new QMenu("&File");
	OpenFileAction 	*openAction 	= new OpenFileAction("&Open", fileMenu);
	SaveFileAction 	*saveAction 	= new SaveFileAction("&Save", fileMenu);
	QAction 	*exitAction 	= new QAction("E&xit", fileMenu);
	fileMenu->addAction( openAction ); fileMenu->addAction( saveAction ); fileMenu->addAction( exitAction );

	QMenuBar* menu = new QMenuBar;
	menu->addMenu(fileMenu);

	QObject::connect(exitAction, SIGNAL(triggered()), qApp, SLOT(quit()));
	QObject::connect(openAction, &OpenFileAction::triggered, openAction, &OpenFileAction::openFile);
	QObject::connect(saveAction, &SaveFileAction::triggered, saveAction, &SaveFileAction::saveFile);



	QMainWindow mainWindow;
	mainWindow.setMenuBar(menu);
	mainWindow.setCentralWidget(widget);
	mainWindow.resize(1750, 850);
	mainWindow.show();


       	return app.exec();
}

