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

#include"App_MeshAdapter.h"
#include"fileFunctions.h"


template<class T>
int app_meshAdapter<T>::read_perturbedNodesData( std::istream& streamIn )
{

std::string line; std::istringstream iss;
size_t nodeSetTemp, nodeSetsCounter(0);
std::string stringTemp; 
T tempValue;
displacementsSource dataSourceTemp;

size_t numOfComponents = m_domainCardinality; //get_initialMesh()->get_numOfComponents();
streamIn.ignore();
T dx, dy, dz, x_0, y_0, z_0, theta_x, theta_y, theta_z, u_x, u_y, u_z, theta, w, x, y, z;
TransformationMatrix<T> trMatrixTemp( numOfComponents );

int flag = get_dataLine( streamIn, line );


while( flag ) {

	iss.clear(); iss.str(line); 
	iss >> nodeSetTemp; insert_perturbedNodeSet( nodeSetTemp );
	iss >> dataSourceTemp; insert_dataSource( nodeSetTemp, dataSourceTemp );

	switch ( dataSourceTemp ) {	
		case ProvidedInFile:
			stringTemp.clear();
			iss >> stringTemp; insert_dataFile(nodeSetTemp, stringTemp);
			break;

		case AffineTransformation:
			trMatrixTemp = TransformationMatrix<T>( iss, numOfComponents ); insert_transformation( nodeSetTemp, trMatrixTemp );
			break;

		case Translation:
		       	iss >> dx >> dy >> dz; 
			insert_translations( nodeSetTemp, std::vector<T>{dx, dy, dz} );
			break;

		case Rotation_EulerAngles:
			iss >> x_0 >> y_0 >> z_0 >> theta_x >> theta_y >> theta_z;
			insert_rotationsEulerAngles( nodeSetTemp, std::vector<T>{x_0, y_0, z_0, theta_x, theta_y, theta_z} );
			break;

		case Rotation_AxisAngle:
			iss >> u_x >> u_y >> u_z >> theta;
			insert_rotationsAxisAngle( nodeSetTemp, std::vector<T>{u_x, u_y, u_z, theta} );
			break;

		case Rotation_Quaternion:
			iss >> w >> x >> y >> z;
			insert_rotationsQuaternion( nodeSetTemp, std::vector<T>{w, x, y, z} );
			break;

		default: std::cerr << "\nInvalid option selected.\n"; break;
	}

	nodeSetsCounter++;

	flag = get_consecutiveData( streamIn, line );
}

return nodeSetsCounter;

}



template<class T>
int app_meshAdapter<T>::read_weightingData( std::istream& streamIn )
{


std::string line; std::istringstream iss;
size_t nodeSetTemp, nodeSetsCounter(0);
T weightTemp, cutOffTemp, interpFactorTemp;

streamIn.ignore();

std::set<size_t> nodesets = get_perturbedNodeSets();

int flag = get_dataLine( streamIn, line );

while( flag ) {

	iss.clear(); iss.str(line); 
	iss >> nodeSetTemp;
	iss >> weightTemp;
	iss >> cutOffTemp;
	iss >> interpFactorTemp;

	//std::cout << nodeSetTemp << std::endl;
	if ( nodesets.find( nodeSetTemp ) == nodesets.end() ) { std::cout << "Node sets specified in weighting data are not compatible with the node sets in the adaptation data.\n"; throw; return 0; }

	std::vector<T> newData = { weightTemp, cutOffTemp, interpFactorTemp };

	insert_weightingData( nodeSetTemp, newData );

	nodeSetsCounter++;

	flag = get_consecutiveData( streamIn, line );
}

if ( nodesets.size( ) > nodeSetsCounter ) { std::cout << "Weighting data has not been given for every adapted node set.\n"; throw; return 0; }

return nodeSetsCounter;

}



template<class T>
int app_meshAdapter<T>::createProject_console ()
{

std::cout << "\n\nSpecify the project name: ";
std::cin >> m_projectName; m_projectName += ".msh";

std::cout << "\n\nSpecify the domain cardinality: ";
std::cin >> m_domainCardinality;


std::cout << "\n\nData of perturbed boundary nodes: (node set number) (fixation type) (parameters for current node set and fixation type)\ne.g.(node set 3 is first translated and then rotated) \n3 2 0.5 0.5 1\n3 3 1 1 0.5 0.2 0.3 0.4\n";
std::cout << "Node set number: 0-n\n";
std::cout << "Displacement type:\n";
std::cout << "\t0)\t From file (under the name 'dxdy.dat').\n";
std::cout << "\t1)\t Affine transformation of boundary node set.\n";
std::cout << "\t2)\t Linear translation of boundary node set.\n";
std::cout << "\t3)\t Rotation of boundary node set, with Euler angles around a specified point.\n";
std::cout << "\t4)\t Rotation of boundary node set, with an axis and an angle.\n";
std::cout << "\t5)\t Rotation of boundary node set, with a quaternion.\n";
std::cout << "Parameters of the fixation type of the current node set:\n";
std::cout << "\tNode set perturbed from file (under the name 'dx.dat').\n";
std::cout << "\tNode set perturbed with an affine transformation, transformation parameters: a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44.\n";
std::cout << "\tNode set perturbed with a translation, translation parameters: dx, dy, dz.\n";
std::cout << "\tNode set perturbed with a rotation, Euler angles as parameters: x0, y0, z0, theta_x, theta_y, theta_z.\n";
std::cout << "\tNode set perturbed with a rotation, a rotation axis and an angle as parameters: u_x, u_y, u_z, theta.\n";
std::cout << "\tNode set perturbed with a rotation, a quaternion as parameters: w, x, y, z.\n";


read_perturbedNodesData( std::cin );


std::cout << "\n\nOptions for the solution method:\n";
std::cout << "0)\tMoving least squares in LCS (MLS_LCS).\n";
std::cout << "Please, enter an option number: ";
std::cin >> m_option_solutionMethod;

std::cout << "\n\nPlease, enter the desired polynomial degree: ";
std::cin >> m_polynomialDegree;

std::cout << "\n\nOptions for the weighting function type:\n";
std::cout << "0)\tWendland distribution\n";
std::cout << "1)\tGaussian distribution\n";
std::cout << "2)\tInverse distance distribution\n";
std::cout << "Please, enter an option number: ";
std::cin >> m_option_weightingFunctionType;


std::cout << "\n\nPlease, enter the desired span, cutoff multiplier, and interpolation factor of the weighting function (in one line per each node set): ";
//std::cin >> m_weightingFunctionSpan >> m_cutOffMultiplier;
//if ( m_weightingFunctionSpan < 1e-10 ) std::cerr << "Please, provide a non zero radius.\n";
read_weightingData( std::cin );


std::cout << "\n\nOptions for the computing mode:\n";
std::cout << "0)\tCPU only\n";
std::cout << "1)\tCPU multithreading\n";
std::cout << "2)\tCUDA GPU only\n";
std::cout << "3)\tCPU + CUDA GPU\n";
std::cout << "4)\tCUDA MultiGPU\n";
std::cout << "Please, enter an option number: ";
std::cin >> m_option_computingMode;


std::cout << "\n\nOptions for the linear solver:\n";
std::cout << "1)\tGauss elimination with back substitution\n";
std::cout << "2)\tConjugate Gradient\n";
std::cout << "3)\tLU Doolittle\n";
std::cout << "4)\tLU Crout\n";
std::cout << "5)\tLDU\n";
std::cout << "6)\tLDL\n";
std::cout << "7)\tCholesky bandedn\n";
std::cout << "8)\tCholesky skyline\n";
std::cout << "Please, enter an option number: ";
std::cin >> m_option_linearSolverType;



std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
std::cout << "\n\nPlease, specify the node data file (press enter for default projectName+\".nod\"): ";
if (std::cin.peek() != '\n') std::cin >> m_nodesFile;
else m_nodesFile = m_projectName.substr(0, m_projectName.rfind('.')) + ".nod";
std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
std::cout << m_nodesFile << std::endl;

std::cout << "\n\nPlease, specify the connectivities data file (press enter for default projectName+\".ele\"): ";
if (std::cin.peek() != '\n') std::cin >> m_connectivitiesFile;
else m_connectivitiesFile = m_projectName.substr(0, m_projectName.rfind('.'))+ ".ele";
std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
std::cout << m_connectivitiesFile << std::endl;

std::cout << "\n\nPlease, specify a file to store the adapted node data (press enter for default projectName+\".res\"): ";
if (std::cin.peek() != '\n') std::cin >> m_fileOfAdaptedMeshNodes;
else m_fileOfAdaptedMeshNodes = m_projectName.substr(0, m_projectName.rfind('.'))+ ".res";
std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
std::cout << m_fileOfAdaptedMeshNodes << std::endl;


}




template<class T>
int app_meshAdapter<T>::readProjectSettings ( std::string fileName ) {

std::ifstream newFileStream; newFileStream.exceptions(std::ios_base::badbit); newFileStream.open( fileName );
std::string line; std::istringstream iss;

std::cerr << "File open status: " << strerror(errno) << ": " << fileName << std::endl;


size_t lineOfInputIndex(0);

while( lineOfInputIndex < 4 )
{
	if ( !get_dataLine( newFileStream, line) ) break;

	iss.clear(); iss.str(line); lineOfInputIndex++;

	switch (lineOfInputIndex) {
		case 1:
			iss >> m_projectName; continue;

		case 2:
			iss >> m_nodesFile; continue;

		case 3:
			iss >> m_connectivitiesFile; continue;

		case 4:
			iss >> m_domainCardinality; continue;
	}
}


size_t numOfLines = read_perturbedNodesData( newFileStream );
lineOfInputIndex = 0;


while( lineOfInputIndex < 3 )
{

	if ( !get_dataLine( newFileStream, line) ) break;

	iss.clear(); iss.str(line); lineOfInputIndex++;

	switch (lineOfInputIndex) {
	
		case 1:
			iss >> m_option_solutionMethod; continue;

		case 2:
			iss >> m_polynomialDegree; continue;

		case 3:
			iss >> m_option_weightingFunctionType; continue;
	}
}

read_weightingData( newFileStream );

lineOfInputIndex = 0;

while( get_dataLine( newFileStream, line) )
{
	iss.clear(); iss.str(line); lineOfInputIndex++;

	switch (lineOfInputIndex) {
	
		case 1:
			iss >> m_option_computingMode; continue;

		case 2:
			iss >> m_option_linearSolverType; continue;

		case 3:
			iss >> m_fileOfAdaptedMeshNodes; continue;
	}
}


if (m_nodesFile == "-") m_nodesFile = m_projectName.substr(0, m_projectName.rfind('.')) + ".nod";
if (m_connectivitiesFile == "-") m_connectivitiesFile = m_projectName.substr(0, m_projectName.rfind('.'))+ ".ele";
if (m_fileOfAdaptedMeshNodes == "-") m_fileOfAdaptedMeshNodes = m_projectName.substr(0, m_projectName.rfind('.'))+ ".res";

m_weightingFunctionSpan = (*m_weightingData.begin()).second[0];
m_cutOffMultiplier = (*m_weightingData.begin()).second[1];
m_interpolationFactor = (*m_weightingData.begin()).second[2];

}



template<class T>
int app_meshAdapter<T>::saveProjectSettings (std::string file) {

std::ofstream newFileStream; newFileStream.exceptions(std::ios_base::badbit); newFileStream.open(file);
newFileStream << "! Project file\n";
newFileStream << m_projectName << "\n\n";
newFileStream << "! Node data file (use '-' for a file path similar to the project file, with .nod extension) \n";
newFileStream << m_nodesFile << "\n\n";
newFileStream << "! Connectivity data file (use '-' for a file path similar to project file, with .ele extension) \n";
newFileStream << m_connectivitiesFile << "\n\n";
newFileStream << "! Domain cardinality \n";
newFileStream << m_domainCardinality << "\n\n";
newFileStream << "! Perturbed nodes data: (node set) (data source type){0:file, 1:affine transformation, 2:Linear displacement, 3:Rotation_EulerAngles, 4:Rotation_AxisAngle, 5:Rotation_Quaternion} (parameters for current node set and fixation type)\n";

for ( auto& i : m_perturbedNodeSets ) {

	auto sources 		= m_dataSources[i].begin();
	auto transformations 	= m_transformations[i].begin();
	auto translations 	= m_translations[i].begin();
	auto rotationsEulerAngles= m_rotationsEulerAngles[i].begin();
	auto rotationsAxisAngle	= m_rotationsAxisAngle[i].begin();
	auto rotationsQuaternion= m_rotationsQuaternion[i].begin();
	auto dataFile		= m_dataFiles[i].begin();

	for ( std::vector<displacementsSource>::iterator j = sources; j != m_dataSources[i].end(); j++ ) {
		newFileStream << i;
		switch ( *j ) {
			case 0:
				newFileStream << " 0";
				newFileStream << " " << *(dataFile++);
				continue;
			case 1:
				newFileStream << " 1";
				for ( auto& el: *(transformations++) ) newFileStream << " " << el;
				continue;
			case 2:
				newFileStream << " 2";
				for ( auto& el: *(translations++) ) newFileStream << " " << el;
				continue;
			case 3:
				newFileStream << " 3";
				for ( auto& el: *(rotationsEulerAngles++) ) newFileStream << " " << el;
				continue;
			case 4:
				newFileStream << " 4";
				for ( auto& el: *(rotationsAxisAngle++) ) newFileStream << " " << el;
				continue;
			case 5:
				newFileStream << " 5";
				for ( auto& el: *(rotationsQuaternion++) ) newFileStream << " " << el;
				continue;
		}
	}

	newFileStream << std::endl;
}


newFileStream << std::endl;
newFileStream << "! Solution method, option number (0:MLS_LCS)\n";
newFileStream << size_t(m_option_solutionMethod) << "\n\n";
newFileStream << "! Polynomial degree\n";
newFileStream << m_polynomialDegree << "\n\n";
newFileStream << "! Weighting function type, option number (0:Wendland, 1:Gaussian, 2:Inverse Distance)\n";
newFileStream << size_t(m_option_weightingFunctionType) << "\n\n";


newFileStream << "! Weighting function span and cutoff multiplier (in one line per each node set)\n";
for ( auto& i : m_perturbedNodeSets ) {

	auto nodesetData = m_weightingData[i];
	newFileStream << i << "\t";
	newFileStream << nodesetData[0] << "\t";
	newFileStream << nodesetData[1];
	newFileStream << std::endl;
}

newFileStream << "! Computing mode, option number (0:CPU, 1:CPU_MT, 2:CUDA, 3:CPU+CUDA, 4:CUDA MultiGPU)\n";
newFileStream << size_t(m_option_computingMode) << "\n\n";
newFileStream << "! Linear solver type (1:GaussEliminationWithBackSubst, 2:Conjugate Gradient, 3:LU Doolittle, 4:LU Crout, 5:LDU, 6:LDL, 7:Cholesky, 8:Cholesky banded, 9:Cholesky skyline)\n";
newFileStream << size_t(m_option_linearSolverType) << "\n\n";
newFileStream << "! File to store the adapted node data (use '-' for a file path similar to project file, with .res extension) \n";
newFileStream << m_fileOfAdaptedMeshNodes << "\n";

}




template<class T>
int app_meshAdapter<T>::saveProjectSettings ( ) {
	return saveProjectSettings( m_projectName );
}




template<class T>
int app_meshAdapter<T>::write_nodesInfo (std::ostream& out)
{
	out << "\nTotal number of nodes:\t\t " 		<< m_initialMesh->get_numOfNodes() << std::endl;
	out << "Number of adapted nodes:\t " 		<< m_meshBoundaryAdapter->get_numOfAdaptedNodes() << std::endl;
	out << "Number of perturbed nodes:\t " 		<< m_meshBoundaryAdapter->get_numOfPerturbedNodes () << std::endl;
	out << "Total number of boundary nodes:\t " 	<< m_initialMesh->get_numOfBoundaryNodes() << std::endl;
}



template<class T>
int app_meshAdapter<T>::write_adaptationProcessInfo (std::ostream& out)
{
	out << "Adaptation process completed " 	<< ( (m_meshBoundaryAdapter->get_adaptationCompletionFlag())?"successfully":"unsuccessfully") << std::endl;
//	out << "Adaptation type: " 		<< m_meshBoundaryAdapter->get_typename() << std::endl;
//	out << "Adaptation parameters: "; 	m_meshBoundaryAdapter->write_adaptationParameters(out);
//	out << "Polynomial degree: " 		<< m_polynomialDegree << std::endl;
//	out << "Weighting function span: " 	<< m_weightingFunctionSpan << std::endl;
//	out << "Total number of polynomial terms: " << m_meshBoundaryAdapter->get_numOfPolynomialTerms() << std::endl;
	out << "Execution time: " 		<< m_meshBoundaryAdapter->get_adaptationTime() << "sec" << std::endl;
	out << "Number of bad elements: " 	<< m_meshBoundaryAdapter->get_badElements().size() << std::endl;

}
