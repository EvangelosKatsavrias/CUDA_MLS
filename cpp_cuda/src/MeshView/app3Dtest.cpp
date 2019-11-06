#include <QApplication>

#include <Qt3DCore/QEntity>
#include <Qt3DExtras/Qt3DWindow>

#include"App_MeshAdapter.h"
#include"viewSettings3D.h"
#include"myAlgorithms.h"
#include"UCSAxesLines.h"
#include"myQMesh.h"
#include"myQMeshRendered.h"


#include <iostream>


int app3Dtest( )
{
	int argc(1); char *argv[1]; char title[16]={"Mesh adaptation"}; argv[0] = title;
	QApplication app(argc, argv);


/*
	app_meshAdapter<float> app1;
	app1.readProjectSettings( "Cube/cube.msh" );
	app1.set_initialMesh(  );
	app1.set_adaptedMesh(  );
	app1.set_meshAdapter();
	app1.write_nodesInfo ( std::cout );
	app1.adaptMesh();
	app1.write_adaptationProcessInfo( std::cout );
	myQMesh<float>* qmesh = new myQMesh<float>( rootEntity, app1.get_initialMesh() );
	myQMesh<float>* qmesh2 = new myQMesh<float>( rootEntity, app1.get_adaptedMesh(), QColor(Qt::blue) );
*/


	app_meshAdapter<float> app1;
	app1.readProjectSettings( "Alfa147/alfa147.msh" );
	app1.set_initialMesh(  );
	Isless_equal<float> fun( -90.f );
	IsInRange<float> fun2( 2.f, -2.f );
	app1.get_initialMesh()->set_nodeFlagWithCriterionOnComponentValue( 1, fun, 1 );
	app1.get_initialMesh()->set_nodeFlagWithCriterionOnComponentValue( 1, fun2, 2 );
	std::set<int> flags = app1.get_initialMesh()->get_flagValues();
	std::cout << " =============================================================> flag values: "; for ( auto& el: flags ) std::cout << el<< std::endl;
	app1.set_adaptedMesh(  );
	app1.set_meshAdapter();
	app1.write_nodesInfo ( std::cout );
	app1.adaptMesh();
	app1.write_adaptationProcessInfo( std::cout );


	Qt3DExtras::Qt3DWindow* view 		= new Qt3DExtras::Qt3DWindow();
	Qt3DCore::QEntity*	rootEntity 	= new Qt3DCore::QEntity();
	view3DSettings* 	viewSetup 	= new view3DSettings(rootEntity, view);
	UCSAxesLinesScene*	axesLines	= new UCSAxesLinesScene(rootEntity);	
	myQMeshRendered<float>* qmesh = new myQMeshRendered<float>( rootEntity, app1.get_initialMesh(), QColor(Qt::white) );
	myQMesh<float>* qmesh_l = new myQMesh<float>( rootEntity, app1.get_initialMesh(), QColor(Qt::white) );
	view->show();
	view->resize(1200, 800);


	Qt3DExtras::Qt3DWindow* view2 		= new Qt3DExtras::Qt3DWindow();
	Qt3DCore::QEntity*	rootEntity2 	= new Qt3DCore::QEntity();
	view3DSettings* 	viewSetup2 	= new view3DSettings(rootEntity2, view2);
	UCSAxesLinesScene*	axesLines2	= new UCSAxesLinesScene(rootEntity2);	
	myQMeshRendered<float>* qmesh2 = new myQMeshRendered<float>( rootEntity2, app1.get_adaptedMesh(), QColor(Qt::white) );
	myQMesh<float>* qmesh2_l = new myQMesh<float>( rootEntity2, app1.get_adaptedMesh(), QColor(Qt::white) );
	view2->show();
	view2->resize(1200, 800);

       	return app.exec();
}
