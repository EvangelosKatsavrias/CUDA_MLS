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

#include<QApplication>
#include"meshViewer.h"
#include"App_MeshAdapter.h"
#include"myFileMenu_gui.h"
#include"Colorbar.h"

#include<numeric>
#include<sstream>


int meshViewerTest(  int argc, char *argv[]  )
{
	int argc1(1); char *argv1[1]; char title[16]={"Mesh adaptation"}; argv1[0] = title;
	QApplication app(argc1, argv1 );


	typedef double data_t;

	std::string fileName = argv[1];

	app_meshAdapter<data_t> app1;
	app1.readProjectSettings( fileName );
	app1.set_initialMesh(  );
	app1.set_adaptedMesh(  );
	app1.set_meshAdapter();
	app1.adaptMesh();

	app1.write_adaptationProcessInfo( std::cout );
	app1.write_adaptedNodesFile();
	app1.write_invalidElementsFile();

	//app1.saveProjectSettings( fileName );
//	app1.assure_results("../temp/plegmata/res.nod");


	// ====================== Evaluate statistics on elements quality ====================

	Connectivities<size_t, 3>* conn = app1.get_initialMesh()->get_connectivities();
	ContainerBase<data_t, std::vector< data_t > >* qualityMetrics = app1.get_initialMesh()->get_qualityMetrics();
	data_t meanQualityValue(0);
	for ( auto& value: qualityMetrics[0] ) {
		meanQualityValue += value;
	}
	auto minmaxQuality = std::minmax_element( qualityMetrics->begin(), qualityMetrics->end() );
	meanQualityValue /= conn->get_numOfConnectivities();

	data_t qualityStdDeviation(0);
	for ( auto& value: qualityMetrics[0] ) {
		qualityStdDeviation += (value-meanQualityValue)*(value-meanQualityValue);
	}
	qualityStdDeviation /= conn->get_numOfConnectivities()-1;
	qualityStdDeviation = sqrt(qualityStdDeviation);


	ContainerBase<data_t, std::vector< data_t > >* qualityMetrics_ad = app1.get_adaptedMesh()->get_qualityMetrics();
	data_t meanQualityValue_ad(0);
	for ( auto& value: qualityMetrics_ad[0] ) {
		meanQualityValue_ad += value;
	}
	auto minmaxQuality_ad = std::minmax_element( qualityMetrics_ad->begin(), qualityMetrics_ad->end() );
	meanQualityValue_ad /= conn->get_numOfConnectivities();

	data_t qualityStdDeviation_ad(0);
	for ( auto& value: qualityMetrics_ad[0] ) {
		qualityStdDeviation_ad += (value-meanQualityValue_ad)*(value-meanQualityValue_ad);
	}

	qualityStdDeviation_ad /= conn->get_numOfConnectivities()-1;
	qualityStdDeviation_ad = sqrt(qualityStdDeviation_ad);




	// ==================== Visualize Initial and Adapted meshes ====================

	QGraphicsItem_Mesh2D_Edges<data_t>* initialMesh_edges = new QGraphicsItem_Mesh2D_Edges<data_t>(0, app1.get_initialMesh() );
	QGraphicsItem_Mesh2D_Edges<data_t>* adaptedMesh_edges = new QGraphicsItem_Mesh2D_Edges<data_t>(0, app1.get_adaptedMesh() );
	QGraphicsItem_Mesh2D_Vertices<data_t>* initialMesh_vertices = new QGraphicsItem_Mesh2D_Vertices<data_t>(0, app1.get_initialMesh() );
	QGraphicsItem_Mesh2D_Vertices<data_t>* adaptedMesh_vertices = new QGraphicsItem_Mesh2D_Vertices<data_t>(0, app1.get_adaptedMesh() );
	QGraphicsItem_Mesh2D_badElements<data_t>* adaptedMesh_badElements = new QGraphicsItem_Mesh2D_badElements<data_t>(0, app1.get_adaptedMesh(), app1.get_meshAdapter()->get_badElements().size(), app1.get_meshAdapter()->get_badElements().data() );
	QGraphicsItem_Mesh2D_ElementsQuality<data_t>* initialMesh_ElementsQuality = new QGraphicsItem_Mesh2D_ElementsQuality<data_t>(0, app1.get_initialMesh() );
	QGraphicsItem_Mesh2D_ElementsQuality<data_t>* adaptedMesh_ElementsQuality = new QGraphicsItem_Mesh2D_ElementsQuality<data_t>(0, app1.get_adaptedMesh() );

	const data_t sizee(15);

	// Initial mesh
	QGraphicsScene* sceneInitial = new QGraphicsScene; sceneInitial->setSceneRect(-sizee, -sizee, 2*sizee, 2*sizee);
	sceneInitial->addItem(initialMesh_edges);
	sceneInitial->addItem(initialMesh_vertices);
	sceneInitial->addItem(initialMesh_ElementsQuality);
	QGraphicsTextItem * sceneTitle_initial = sceneInitial->addText("Initial Mesh");
	sceneTitle_initial->setPos(QPoint(-0.8*sizee,-0.8*sizee)); sceneTitle_initial->setFlag( QGraphicsItem::ItemIgnoresTransformations );

	QGraphicsView* viewInitial = new QGraphicsView_mod;
	viewInitial->setScene(sceneInitial); viewInitial->scale(30,30);
	viewInitial->centerOn(initialMesh_edges);
	viewInitial->scale(1, -1);


//	QWidget *initialMeshWidget = QWidget::createWindowContainer( meshgraphicsitem );
	QLabel *initialMeshLabel = new QLabel( QString("Initial Mesh\n")
						+QString("Mesh quality: Min value ") + QString::number( *minmaxQuality.first )
						+QString(", Max value ") + QString::number( *minmaxQuality.second )
						+QString(", Mean value ") + QString::number( meanQualityValue )
						+QString(", Std deviation ") + QString::number( qualityStdDeviation )
						, viewInitial );

	QFont fontMesh = initialMeshLabel->font();
	fontMesh = initialMeshLabel->font();
	fontMesh.setPointSize(10); fontMesh.setBold(true);
	initialMeshLabel->setFont(fontMesh);
	initialMeshLabel->setFixedHeight( 30 );
	QualityColorbar* initialMeshColorbar = new QualityColorbar;

	QLabel *QualityColorbarLabel1 = new QLabel( QString("0") );
	QLabel *QualityColorbarLabel2 = new QLabel( QString("0.5") );
	QLabel *QualityColorbarLabel3 = new QLabel( QString("1") );

	QHBoxLayout* labelsStack = new QHBoxLayout;
	labelsStack->addWidget(QualityColorbarLabel1);
	labelsStack->addStretch();
	labelsStack->addWidget(QualityColorbarLabel2);
	labelsStack->addStretch();
	labelsStack->addWidget(QualityColorbarLabel3);

	QVBoxLayout* initialMeshStack = new QVBoxLayout;
	initialMeshStack->addWidget( initialMeshLabel );
	initialMeshStack->addWidget( viewInitial );
	initialMeshStack->addLayout( labelsStack );
	initialMeshStack->addWidget( initialMeshColorbar );



	// Adapted mesh
	QGraphicsScene* sceneAdapted = new QGraphicsScene; sceneAdapted->setSceneRect(-sizee,-sizee,2*sizee,2*sizee);
	sceneAdapted->addItem(adaptedMesh_edges);
//	sceneAdapted->addItem(adaptedMesh_vertices);
//	sceneAdapted->addItem(adaptedMesh_badElements);
	sceneAdapted->addItem(adaptedMesh_ElementsQuality);
	QGraphicsTextItem * sceneTitle_adapted = sceneAdapted->addText("Adapted Mesh");
	sceneTitle_adapted->setPos(QPoint(-0.8*sizee,-0.8*sizee)); sceneTitle_adapted->setFlag( QGraphicsItem::ItemIgnoresTransformations );

	QGraphicsView* viewAdapted = new QGraphicsView_mod;
	viewAdapted->setScene(sceneAdapted); viewAdapted->scale(30,30);
	viewAdapted->centerOn( adaptedMesh_edges );
	viewAdapted->scale(1, -1);

	QLabel *adaptedMeshLabel = new QLabel( QString("Adapted Mesh\n")
						+QString("Mesh quality: Min value ") + QString::number( *minmaxQuality_ad.first )
						+QString(", Max value ") + QString::number( *minmaxQuality_ad.second )
						+QString(", Mean value ") + QString::number( meanQualityValue_ad )
						+QString(", L2 norm ") + QString::number( qualityStdDeviation_ad )
						, viewAdapted);
	adaptedMeshLabel->setFont(fontMesh);
	adaptedMeshLabel->setFixedHeight( 30 );
	QualityColorbar* adaptedMeshColorbar = new QualityColorbar;

	QLabel *QualityColorbarLabel21 = new QLabel( QString("0") );
	QLabel *QualityColorbarLabel22 = new QLabel( QString("0.5") );
	QLabel *QualityColorbarLabel23 = new QLabel( QString("1") );

	QHBoxLayout* labelsStack2 = new QHBoxLayout;
	labelsStack2->addWidget(QualityColorbarLabel21);
	labelsStack2->addStretch();
	labelsStack2->addWidget(QualityColorbarLabel22);
	labelsStack2->addStretch();
	labelsStack2->addWidget(QualityColorbarLabel23);

	QVBoxLayout* adaptedMeshStack = new QVBoxLayout;
	adaptedMeshStack->addWidget( adaptedMeshLabel );
	adaptedMeshStack->addWidget( viewAdapted );
	adaptedMeshStack->addLayout( labelsStack2 );
	adaptedMeshStack->addWidget( adaptedMeshColorbar );


	QHBoxLayout* stack = new QHBoxLayout;
//	stack->addWidget(viewInitial); stack->addWidget(viewAdapted);
	stack->addLayout( initialMeshStack ); stack->addLayout( adaptedMeshStack );

	QWidget* centralWidget = new QWidget;
	centralWidget->setLayout(stack);





	// ====================== Evaluate displacement and strain fields ==========================
	size_t numOfNodes = app1.get_initialMesh()->get_numOfNodes();
	ContainerBase<data_t, std::vector< data_t > >* elementsArea = app1.get_initialMesh()->get_metrics();

	// displacement field x
	data_t* init_x = app1.get_initialMesh()->get_component(0);
	data_t* adapt_x = app1.get_adaptedMesh()->get_component(0);
	data_t disp_x[app1.get_adaptedMesh()->get_numOfNodes()];
	data_t meanDisp_x(0);
	for ( size_t i = 0; i < numOfNodes; i++ ) {
		disp_x[i] = adapt_x[i] -init_x[i];
		meanDisp_x += disp_x[i];
	}
	auto maxDisp_x = std::max_element( disp_x, disp_x+numOfNodes );
	auto minDisp_x = std::min_element( disp_x, disp_x+numOfNodes );
	meanDisp_x /= numOfNodes;

	data_t L2norm_dispx(0); size_t* vert;
	for ( size_t elem = 0; elem < conn->get_numOfConnectivities(); elem++ ) {
		vert = conn[0][elem];
		data_t totDisp = ( disp_x[vert[0]] +disp_x[vert[1]] +disp_x[vert[2]] )*0.3333333;
		L2norm_dispx += totDisp*totDisp*elementsArea[0][elem];
	}
	L2norm_dispx = sqrt( L2norm_dispx );

	// displacement field y
	data_t* init_y = app1.get_initialMesh()->get_component(1);
	data_t* adapt_y = app1.get_adaptedMesh()->get_component(1);
	data_t disp_y[app1.get_adaptedMesh()->get_numOfNodes()];
	data_t meanDisp_y(0);
	for ( size_t i = 0; i < numOfNodes; i++ ) {
		disp_y[i] = adapt_y[i] -init_y[i];
		meanDisp_y += disp_y[i];
	}
	auto maxDisp_y = std::max_element( disp_y, disp_y+numOfNodes );
	auto minDisp_y = std::min_element( disp_y, disp_y+numOfNodes );
	meanDisp_y /= numOfNodes;

	data_t L2norm_dispy(0);
	for ( size_t elem = 0; elem < conn->get_numOfConnectivities(); elem++ ) {
		vert = conn[0][elem];
		data_t totDisp = ( disp_y[vert[0]] +disp_y[vert[1]] +disp_y[vert[2]] )*0.3333333;
		L2norm_dispy += totDisp*totDisp*elementsArea[0][elem];
	}
	L2norm_dispy = sqrt(L2norm_dispy);


	std::ofstream fileStrOut; fileStrOut.open("./strainsData.dat");

	// strain field dux/dx
	std::vector< std::set<size_t> >& vertexVertex_conn = app1.get_initialMesh()->get_connectivities()->get_nodeLinks();
	EvaluationData<data_t>* derivatives = app1.get_meshAdapter()->get_derivativesData();
	data_t* gradux_dx = derivatives->get_fieldComponent( 0 );


//	data_t gradux_dx[ app1.get_adaptedMesh()->get_numOfNodes() ];
	
//	for ( size_t i = 0; i < numOfNodes; i++ ) std::cout << gradux_dx[i] << std::endl;

	size_t node(0);
	data_t nodeDisp, nodePos_x; 
	data_t meandux(0);
	for ( auto& vertex: vertexVertex_conn ) {
//		nodeDisp = disp_x[node];
//		nodePos_x = init_x[node];
//		gradux_dx[node] = 0;
//		size_t counter(0); data_t maxDist(0); size_t maxIndex(0);
//		for ( auto& vv: vertex ) {
//			data_t dx = std::abs( init_x[vv] -nodePos_x ); if ( dx > maxDist ) { maxDist = dx; maxIndex = vv; }
//			counter++;
//		}
//		gradux_dx[node] = std::abs(disp_x[ maxIndex ] -nodeDisp) / maxDist;

		meandux += gradux_dx[node];
		node++;
	}
	auto maxdux = std::max_element( gradux_dx, gradux_dx+numOfNodes );
	auto mindux = std::min_element( gradux_dx, gradux_dx+numOfNodes );
	meandux /= numOfNodes;

	data_t L2norm_dudx(0);
	for ( size_t elem = 0; elem < conn->get_numOfConnectivities(); elem++ ) {
		vert = conn[0][elem];
		data_t totGrad = ( gradux_dx[vert[0]] +gradux_dx[vert[1]] +gradux_dx[vert[2]] )*0.333333;
		L2norm_dudx += totGrad*totGrad*elementsArea[0][elem];
		if ( totGrad*totGrad < 1e-8 ) { gradux_dx[vert[0]] = 0; gradux_dx[vert[1]] = 0; gradux_dx[vert[2]] = 0; L2norm_dudx = 0;  }
	}
	L2norm_dudx = sqrt(L2norm_dudx);


	// strain field duy/dy
	data_t* graduy_dy = derivatives->get_fieldComponent( 1 );

//	data_t graduy_dy[ app1.get_adaptedMesh()->get_numOfNodes() ];
	node = 0;
	data_t nodePos_y;
	data_t meanduy(0);
	for ( auto& vertex: vertexVertex_conn ) {
//		nodeDisp = disp_y[node];
//		nodePos_y = init_y[node];
		//graduy_dy[node] = 0;
//		size_t counter(0); data_t maxDist(0); size_t maxIndex(0);
//		for ( auto& vv: vertex ) {
//			data_t dy = std::abs( init_y[vv] -nodePos_y ); if ( dy > maxDist ) { maxDist = dy; maxIndex = vv; }
//			counter++;
//		}
//		graduy_dy[node] = std::abs(disp_y[ maxIndex ] -nodeDisp) / maxDist;

		meanduy += graduy_dy[node];
		node++;
	}
	auto maxduy = std::max_element( graduy_dy, graduy_dy+numOfNodes );
	auto minduy = std::min_element( graduy_dy, graduy_dy+numOfNodes );
	meanduy /= numOfNodes;

	data_t L2norm_dvdy(0);
	for ( size_t elem = 0; elem < conn->get_numOfConnectivities(); elem++ ) {
		vert = conn[0][elem];
		data_t totGrad = ( graduy_dy[vert[0]] +graduy_dy[vert[1]] +graduy_dy[vert[2]] )*0.333333;
		L2norm_dvdy += totGrad*totGrad*elementsArea[0][elem];
		if ( totGrad*totGrad < 1e-6 ) { graduy_dy[vert[0]] = 0; graduy_dy[vert[1]] = 0; graduy_dy[vert[2]] = 0; L2norm_dvdy = 0; }
	}
	L2norm_dvdy = sqrt(L2norm_dvdy);


	// shear strain field

	data_t* shear_xy = derivatives->get_fieldComponent( 2 );
//	data_t shear_xy[ app1.get_adaptedMesh()->get_numOfNodes() ];
	node = 0;
	data_t nodeDisp_y;
	data_t meanshear(0);
	for ( auto& vertex: vertexVertex_conn ) {
//		nodeDisp = disp_x[node];
//		nodeDisp_y = disp_y[node];
//		nodePos_x = init_x[node];
//		nodePos_y = init_y[node];
//		shear_xy[node] = 0;
//		data_t tempShear;
//		size_t counter(0); data_t maxDist(0); size_t maxIndex(0);
//		for ( auto& vv: vertex ) {
//			data_t dx = std::abs( init_x[vv] -nodePos_x );
//			data_t dy = std::abs( init_y[vv] -nodePos_y ); 
//			tempShear = ( disp_x[ vv ] -nodeDisp) / dy + ( disp_y[ vv ] -nodeDisp_y ) / dx;
//			if ( (counter == 0 || tempShear < shear_xy[node]) && std::isfinite(tempShear) ) shear_xy[node] = tempShear;
//			counter++;
//		}
		//shear_xy[node] = ( disp_x[ maxIndex ] -nodeDisp) / ( init_y[ maxIndex ] -nodePos_y ) + ( disp_y[ maxIndex ] -nodeDisp_y ) / ( init_x[ maxIndex ] -nodePos_x );

		meanshear += shear_xy[node];
		node++;
	}
	auto maxshear = std::max_element( shear_xy, shear_xy+numOfNodes );
	auto minshear= std::min_element( shear_xy, shear_xy+numOfNodes );
	meanshear /= numOfNodes;

	data_t L2norm_shear(0);
	for ( size_t elem = 0; elem < conn->get_numOfConnectivities(); elem++ ) {
		vert = conn[0][elem];
		data_t totGrad = ( shear_xy[vert[0]] +shear_xy[vert[1]] +shear_xy[vert[2]] )*0.333333;
		L2norm_shear += totGrad*totGrad*elementsArea[0][elem];
		if ( totGrad*totGrad < 1e-6 ) { shear_xy[vert[0]] = 0; shear_xy[vert[1]] = 0; shear_xy[vert[2]] = 0; L2norm_shear = 0; }
	}
	L2norm_shear = sqrt(L2norm_shear );



	// =========================== Visualize displacement and strain fields ============
	// x disp field
	QGraphicsItem_Mesh2D_fieldPlot<data_t>* fieldPlot_dispx = new QGraphicsItem_Mesh2D_fieldPlot<data_t>(0, app1.get_adaptedMesh(), disp_x );
	QWidget *xDispWidget = QWidget::createWindowContainer( fieldPlot_dispx );
	QLabel *xDispLabel = new QLabel( QString("Displacements x\n")
					+QString("Min value ") + QString::number( *minDisp_x )
					+QString(", Max value ") + QString::number( *maxDisp_x )
					+QString(", Mean value ") + QString::number( meanDisp_x )
					+QString(", L2 norm ") + QString::number( L2norm_dispx )
					, xDispWidget );
	QFont font = xDispLabel->font();
	font.setPointSize(10); font.setBold(true);
	xDispLabel->setFont(font);
	xDispLabel->setFixedHeight( 30 );
	QVBoxLayout* displacementsStack_x = new QVBoxLayout;


	Colorbar* colorbar1 = new Colorbar;
	displacementsStack_x->addWidget( xDispLabel );
	displacementsStack_x->addWidget( xDispWidget );
	displacementsStack_x->addWidget( colorbar1 );


	// y disp field
	QGraphicsItem_Mesh2D_fieldPlot<data_t>* fieldPlot_dispy = new QGraphicsItem_Mesh2D_fieldPlot<data_t>(0, app1.get_adaptedMesh(), disp_y );
	QWidget *yDispWidget = QWidget::createWindowContainer( fieldPlot_dispy );
	QLabel *yDispLabel = new QLabel( QString("Displacements y\n")
					+QString("Min value ") + QString::number( *minDisp_y )
					+QString(", Max value ") + QString::number( *maxDisp_y )
					+QString(", Mean value ") + QString::number( meanDisp_y )
					+QString(", L2 norm ") + QString::number( L2norm_dispy )
					, yDispWidget );
	QFont fonty = yDispLabel->font();
	fonty.setPointSize(10); fonty.setBold(true);
	yDispLabel->setFont(fonty);
	yDispLabel->setFixedHeight( 30 );
	QVBoxLayout* displacementsStack_y = new QVBoxLayout;

	Colorbar* colorbar2 = new Colorbar;
	displacementsStack_y->addWidget( yDispLabel );
	displacementsStack_y->addWidget( yDispWidget );
	displacementsStack_y->addWidget( colorbar2 );


	// disp window
	QHBoxLayout* displacementsStack = new QHBoxLayout;
	displacementsStack->addLayout( displacementsStack_x );
	displacementsStack->addLayout( displacementsStack_y );
	QWidget* displacementsWidget = new QWidget;
	displacementsWidget->setLayout(displacementsStack);



	// du/dx strain field
	QGraphicsItem_Mesh2D_fieldPlot<data_t>* fieldPlot_gradux_dx = new QGraphicsItem_Mesh2D_fieldPlot<data_t>(0, app1.get_adaptedMesh(), gradux_dx );
	QWidget *dudxStrainWidget = QWidget::createWindowContainer( fieldPlot_gradux_dx );
	QLabel *dudxStrainLabel = new QLabel( QString("Grad field du/dx\n")
						+QString("Min value ") + QString::number( *mindux )
						+QString(", Max value ") + QString::number( *maxdux )
						+QString(", Mean value ") + QString::number( meandux )
						+QString(", L2 norm ") + QString::number( L2norm_dudx )
						, dudxStrainWidget );
	font = dudxStrainLabel->font();
	font.setPointSize(10); font.setBold(true);
	dudxStrainLabel->setFont(font);
	dudxStrainLabel->setFixedHeight( 30 );
	QVBoxLayout* strainsStack_dudx = new QVBoxLayout;
	Colorbar* colorbar3 = new Colorbar;
	strainsStack_dudx->addWidget( dudxStrainLabel );
	strainsStack_dudx->addWidget( dudxStrainWidget );
	strainsStack_dudx->addWidget( colorbar3 );


	// dv/dy strain field
	QGraphicsItem_Mesh2D_fieldPlot<data_t>* fieldPlot_graduy_dy = new QGraphicsItem_Mesh2D_fieldPlot<data_t>(0, app1.get_adaptedMesh(), graduy_dy );
	QWidget *dvdyStrainWidget = QWidget::createWindowContainer( fieldPlot_graduy_dy );
	QLabel *dvdyStrainLabel = new QLabel( QString("Grad field dv/dy\n")
						+QString("Min value ") + QString::number( *minduy )
						+QString(", Max value ") + QString::number( *maxduy )
						+QString(", Mean value ") + QString::number( meanduy )
						+QString(", L2 norm ") + QString::number( L2norm_dvdy )
						, dvdyStrainWidget );
	font = dvdyStrainLabel->font();
	font.setPointSize(10); font.setBold(true);
	dvdyStrainLabel->setFont(font);
	dvdyStrainLabel->setFixedHeight( 30 );
	QVBoxLayout* strainsStack_dvdy = new QVBoxLayout;
	Colorbar* colorbar4 = new Colorbar;
	strainsStack_dvdy->addWidget( dvdyStrainLabel );
	strainsStack_dvdy->addWidget( dvdyStrainWidget );
	strainsStack_dvdy->addWidget( colorbar4 );



	// shear strain field
	QGraphicsItem_Mesh2D_fieldPlot<data_t>* fieldPlot_shear_xy = new QGraphicsItem_Mesh2D_fieldPlot<data_t>(0, app1.get_adaptedMesh(), shear_xy );
	QWidget *shearStrainWidget = QWidget::createWindowContainer( fieldPlot_shear_xy );
	QLabel *shearStrainLabel = new QLabel( QString("Grad field dv/dx + du/dy\n")
						+QString("Min value ") + QString::number( *minshear )
						+QString(", Max value ") + QString::number( *maxshear )
						+QString(", Mean value ") + QString::number( meanshear )
						+QString(", L2 norm ") + QString::number( L2norm_shear )
						, shearStrainWidget );
	font = shearStrainLabel->font();
	font.setPointSize(10); font.setBold(true);
	shearStrainLabel->setFont(font);
	shearStrainLabel->setFixedHeight( 30 );
	QVBoxLayout* strainsStack_shear = new QVBoxLayout;
	Colorbar* colorbar5 = new Colorbar;
	strainsStack_shear->addWidget( shearStrainLabel );
	strainsStack_shear->addWidget( shearStrainWidget );
	strainsStack_shear->addWidget( colorbar5 );



	// strains window
	QHBoxLayout* strainsStack = new QHBoxLayout;
	strainsStack->addLayout( strainsStack_dudx );
	strainsStack->addLayout( strainsStack_dvdy );
	strainsStack->addLayout( strainsStack_shear );
	QWidget* strainsWidget = new QWidget;
	strainsWidget->setLayout( strainsStack );


	// tabs widget
	QTabWidget* tabs = new QTabWidget;
	tabs->insertTab( 0, centralWidget, QString( "&Meshes" ) );
	tabs->insertTab( 1, displacementsWidget, QString( "&Displacement fields" ) );
	tabs->insertTab( 2, strainsWidget, QString("&Strain fields") );






	// ================== Visualize info in side panel =================
	computingMode computeMode = app1.get_option_computingMode();
	QString computeModeStr( get_computingMode_name( computeMode ).data() );
	DistributionType distType = app1.get_option_weightingFunctionType();
	QString distTypeStr( get_DistributionType_name( distType ).data() );
	solutionMethod solMethod = app1.get_option_solutionMethod();
	QString solutionMethodStr( get_solutionMethod_name( solMethod ).data() );
	LinearSolverType solverType = app1.get_option_linearSolverType	( );
	QString solverTypeStr( get_LinearSolverType_name( solverType ).data() );
	std::map<size_t, std::vector<displacementsSource> > dataSources = app1.get_dataSources( );
	size_t		polynomialDegree 	= app1.get_polynomialDegree		( );
	data_t		weightingFunSpan 	= app1.get_weightingFunctionSpan	( );
	data_t		cutOffMult 		= app1.get_cutOffMultiplier		( );
	data_t		interpolationFactor	= app1.get_interpolationFactor		( );
	size_t		domainCardinality 	= app1.get_domainCardinality		( );
	float 		adaptationTime = app1.get_meshAdapter()->get_adaptationTime( );


	QCommandLinkButton *info = new QCommandLinkButton();
	info->setText(QStringLiteral("Mesh adaptation Visualization"));
	info->setDescription(
			QString("Initial mesh on the left side, adapted mesh on the right side. \nOn the adapted mesh, the degenerated elements are colored as red faces.\nNumber of degenerated elements: ") 				+QString::number( app1.get_meshAdapter()->get_badElements().size() )
			+QString("\n\nAdaptation Info:")
			+QString("\nDomain cardinality: ") +QString::number( domainCardinality )
			+QString("\nSolution method: ") +solutionMethodStr
			+QString("\nPolynomial degree: ") +QString::number( polynomialDegree )
			+QString("\nComputing mode: ") +computeModeStr
			+QString("\nWeighting function type: ") +distTypeStr
			+QString("\nWeighting function span: ") +QString::number( weightingFunSpan )
			+QString("\nCut off multiplier: ") +QString::number( cutOffMult )
			+QString("\nInterpolation factor: ") +QString::number( interpolationFactor )
			+QString("\nSolver type: ") +solverTypeStr
			+QString("\nAdaptation time: ") +QString::number( adaptationTime*1000 )+QString("ms")
			+QString("\n\n\nDisplacement and strain metrics:")
			+QString("\nDisplacements total L2 norm ") + QString::number( sqrt(L2norm_dispx*L2norm_dispx +L2norm_dispy*L2norm_dispy) )
			+QString("\nDisplacements in x L2 norm ") + QString::number( L2norm_dispx )
			+QString("\nDisplacements in y L2 norm ") + QString::number( L2norm_dispy )
			+QString("\nTotal energy seminorm ") + QString::number( L2norm_dudx*L2norm_dudx +L2norm_dvdy*L2norm_dvdy +L2norm_shear*L2norm_shear )
			+QString("\nGrad du/dx L2 norm ") + QString::number( L2norm_dudx )
			+QString("\nGrad dv/dy L2 norm ") + QString::number( L2norm_dvdy )
			+QString("\nGrad du/dy+dv/dx L2 norm ") + QString::number( L2norm_shear )
			+QString("\n\n\nMesh quality info: ")
			+QString("\nInitial Mesh, mesh quality statistics: ")
			+QString("\n  min: ") + QString::number( *minmaxQuality.first )
			+QString("\n  max: ") + QString::number( *minmaxQuality.second )
			+QString("\n  mean: ") + QString::number( meanQualityValue )
			+QString("\n  std.dev.: ") + QString::number( qualityStdDeviation )
			+QString("\nAdapted Mesh, mesh quality statistics: ")
			+QString("\n  min: ") + QString::number( *minmaxQuality_ad.first )
			+QString("\n  max: ") + QString::number( *minmaxQuality_ad.second )
			+QString("\n  mean: ") + QString::number( meanQualityValue_ad )
			+QString("\n  std.dev.: ") + QString::number( qualityStdDeviation_ad )
			);
	info->setIconSize(QSize(0,0));
	stack->addWidget(info);



	// =================== Create an app menu bar =====================
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
	mainWindow.setCentralWidget(tabs);
	mainWindow.resize(1700, 1000);
	mainWindow.show();

	return app.exec();
}
