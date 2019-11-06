#include"App_MeshAdapter.h"
#include"Mesh3D.h"
#include"fileFunctions.h"


typedef float data_t;

void gridDataTest()
{

// ==================== 2D example =============================

app_meshAdapter<data_t> app1;
//Mesh2DDisplacedBoundaryAdapter<double, int, int, 3>* name = new Mesh2DDisplacedBoundaryAdapter<double, int, int, 3>;

app1.readProjectSettings( "Leaf/Leaf.msh" );
//app1.createProject_console();

app1.set_initialMesh( );
app1.set_adaptedMesh( );
app1.set_meshAdapter( );
app1.write_nodesInfo ( std::cout );
app1.adaptMesh();

app1.write_adaptationProcessInfo( std::cout );


//app1.write_adaptedNodesFile();
//app1.write_invalidElementsFile();
//app1.assure_results("../temp/plegmata/res.nod");

//app1.saveProjectSettings( "Leaf/Leaf_test.msh" );


// ================== 3D example =======================
//convert_objFile<float>( "Leaf/teapot.obj" );
//convert_objFile<float>( "Leaf/alfa147.obj" );
//convert_objFile<float>( "Leaf/cube.obj" );

//Mesh3D<float, 3, int, 4>* mesh3d = new Mesh3D<float, 3, int, 4>; delete mesh3d;
/*
app_meshAdapter<float> app1;
app1.readProjectSettings( "Cube/cube.msh" );
app1.set_initialMesh(  );
app1.set_adaptedMesh(  );
app1.set_meshAdapter();
app1.write_nodesInfo ( std::cout );
app1.adaptMesh();
app1.write_adaptationProcessInfo( std::cout );
*/
//std::cout << *app1.get_initialMesh();
//std::cout << *app1.get_adaptedMesh();




// ================== plot data
//
/*
std::cout << "project name " << app1.get_projectName() << std::endl;
std::cout << "nodes file " << app1.get_nodesFile() << std::endl;
std::cout << "conn file " << app1.get_connectivitiesFile() << std::endl;
std::cout << "solution method " << app1.get_option_solutionMethod() << std::endl;
std::cout << "weighting function type " << app1.get_option_weightingFunctionType() << std::endl;
std::cout << "computing mode " << app1.get_option_computingMode() << std::endl;
std::cout << "solver type " << app1.get_option_linearSolverType() << std::endl;
std::cout << "polynomial degree " << app1.get_polynomialDegree() << std::endl;
std::cout << "weighting function span " << app1.get_weightingFunctionSpan() << std::endl;
std::cout << "cut off multiplier " << app1.get_cutOffMultiplier() << std::endl;
std::cout << "domain cardinality " << app1.get_domainCardinality() << std::endl;
std::cout << "number of perturbed node sets " << app1.get_perturbedNodeSets().size() << std::endl;
std::cout << "perturbed node sets \n"; for ( auto& el: app1.get_perturbedNodeSets() ) std::cout << el << std::endl;
std::cout << "data sources \n"; for ( auto& set: app1.get_dataSources() ) for ( auto& el: set.second ) std::cout << el<< std::endl;
std::cout << "transformations \n"; for ( auto& set: app1.get_transformations() ) for ( auto& el: set.second ) std::cout << el<< std::endl;
std::cout << "translations \n"; for ( auto& set: app1.get_translations() ) for ( auto& el: set.second ) std::cout << el<< std::endl;
std::cout << "rotations euler angles\n"; for ( auto& set: app1.get_rotationsEulerAngles() ) for ( auto& el: set.second ) std::cout << el<< std::endl;
std::cout << "rotations axis angle\n"; for ( auto& set: app1.get_rotationsAxisAngle() ) for ( auto& el: set.second ) std::cout << el<< std::endl;
std::cout << "rotations quaternion\n"; for ( auto& set: app1.get_rotationsQuaternion() ) for ( auto& el: set.second ) std::cout << el<< std::endl;
*/


}

