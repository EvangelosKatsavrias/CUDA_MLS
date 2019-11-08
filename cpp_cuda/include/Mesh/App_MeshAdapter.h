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

#ifndef APPMESHADAPTER
#define APPMESHADAPTER


#include<stdexcept>
#include<exception>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include"Mesh2D.h"
#include"Mesh3D.h"
#include"MeshAdapters.h"
#include"TransformationMatrix.h"
#include"vectorPlotOutStreamer.h"



template<class T>
class app_meshAdapter
{


public:
	typedef std::pair<size_t, displacementsSource > 		dispSrcPair_t;
	typedef std::pair<size_t, TransformationMatrix<T> > 		trPair_t;
	typedef std::pair<size_t, std::string > 			filePair_t;
	typedef std::pair<size_t, std::vector<T> > 			parameters_t;

	typedef std::set<size_t> 					perturbedNodeSets_t;
	typedef std::map<size_t, std::vector<displacementsSource> > 	dataSources_t;
	typedef std::map<size_t, std::vector<TransformationMatrix<T> > > transformations_t;
	typedef std::map<size_t, std::vector<std::vector<T> > >		translations_t;
	typedef std::map<size_t, std::vector<std::vector<T> > >		rotationsEulerAngles_t;
	typedef std::map<size_t, std::vector<std::vector<T> > >		rotationsAxisAngle_t;
	typedef std::map<size_t, std::vector<std::vector<T> > >		rotationsQuaternion_t;
	typedef std::map<size_t, std::vector<std::string> >		dataFiles_t;
	typedef std::map<size_t, std::vector<T> >			weightingData_t;
	typedef MeshDisplacedBoundaryAdapter<T, 3, int>			meshAdapter_t;
	typedef Mesh<T, 3, int>						mesh_t;


protected:
	std::string 		m_projectName;
	std::string		m_nodesFile;
	std::string		m_connectivitiesFile;
	std::string		m_fileOfAdaptedMeshNodes;


	solutionMethod		m_option_solutionMethod;
	DistributionType	m_option_weightingFunctionType;
	computingMode		m_option_computingMode;
	LinearSolverType	m_option_linearSolverType;
	size_t			m_polynomialDegree;
	T			m_weightingFunctionSpan;
	T			m_cutOffMultiplier;
	T			m_interpolationFactor;
	size_t			m_domainCardinality;


	perturbedNodeSets_t 	m_perturbedNodeSets;
	dataSources_t 		m_dataSources;
	transformations_t 	m_transformations;
	translations_t		m_translations;
	rotationsEulerAngles_t	m_rotationsEulerAngles;
	rotationsAxisAngle_t	m_rotationsAxisAngle;
	rotationsQuaternion_t	m_rotationsQuaternion;
	dataFiles_t		m_dataFiles;
	weightingData_t		m_weightingData;


	mesh_t*			m_initialMesh;
	mesh_t*			m_adaptedMesh;
	meshAdapter_t*		m_meshBoundaryAdapter;

	std::vector<bool>	m_memAllocatedMembers;

public:
	app_meshAdapter ( mesh_t* initialMesh=0, meshAdapter_t* boundaryAdapter=0 );
	~app_meshAdapter (  );


	int 	createProject_console		( );
	int 	readProjectSettings 		( std::string );
	int 	saveProjectSettings 		( );
	int 	saveProjectSettings 		( std::string );
	int	set_nodesFile			( std::string file) { m_nodesFile = file; }
	int	set_connectivitiesFile		( std::string file) { m_connectivitiesFile = file; }
	int	set_fileOfAdaptedMeshNodes	( std::string file ) { m_fileOfAdaptedMeshNodes = file; }
	int	set_projectName			( std::string name ) { m_projectName = name; }


	int	adaptMesh			( );
	int	assure_results			( std::string file );


	int	set_option_solutionMethod 	( solutionMethod solMethod = MLS_LCS ) 			{ m_option_solutionMethod = solMethod; }
	int	set_option_weightingFunctionType( DistributionType weightingFunction = Wendland ) 	{ m_option_weightingFunctionType = weightingFunction; }
	int	set_option_computingMode 	( computingMode computingMode = CPU ) 			{ m_option_computingMode = computingMode; }
	int	set_option_linearSolverType	( LinearSolverType solverType = ConjugateGradient )	{ m_option_linearSolverType = solverType; }
	int	set_polynomialDegree		( size_t newDegree )					{ m_polynomialDegree = newDegree; }
	int	set_weightingFunctionSpan	( T newSpan )						{ m_weightingFunctionSpan = newSpan; }
	int	set_cutOffMultiplier		( T cutOffMultiplier )					{ m_cutOffMultiplier = cutOffMultiplier; }
	int	set_interpolationFactor		( T factor )						{ m_interpolationFactor = factor; }
	int	set_domainCardinality		( size_t cardinality )					{ m_domainCardinality = cardinality; }


	int	insert_perturbedNodeSet		( size_t nodeSet )					{ m_perturbedNodeSets.insert( nodeSet ); }
	int	insert_dataSource		( size_t nodeSet, displacementsSource source )		{ m_dataSources[nodeSet].push_back(source); }
	int	insert_transformation		( size_t nodeSet, TransformationMatrix<T> transformationMatrix ) { m_transformations[nodeSet].push_back(transformationMatrix); }
	int	insert_translations		( size_t nodeSet, std::vector<T> values ) 		{ m_translations[nodeSet].push_back(values); }
	int	insert_rotationsEulerAngles	( size_t nodeSet, std::vector<T> values ) 		{ m_rotationsEulerAngles[nodeSet].push_back(values); }
	int	insert_rotationsAxisAngle	( size_t nodeSet, std::vector<T> values ) 		{ m_rotationsAxisAngle[nodeSet].push_back(values); }
	int	insert_rotationsQuaternion	( size_t nodeSet, std::vector<T> values ) 		{ m_rotationsQuaternion[nodeSet].push_back(values); }
	int	insert_dataFile			( size_t nodeSet, std::string dataFile ) 		{ m_dataFiles[nodeSet].push_back(dataFile); }
	int	insert_weightingData		( size_t nodeSet, std::vector<T> values ) 		{ m_weightingData[nodeSet] = values; }

	
	int	set_perturbedNodeSets 		( std::set<size_t> perturbedNodeSets ) 			{ m_perturbedNodeSets = perturbedNodeSets; }
	int	set_dataSources 		( dataSources_t dataSources ) 				{ m_dataSources = dataSources; }
	int	set_transformations 		( transformations_t transformations ) 			{ m_transformations = transformations; }
	int	set_dataFiles			( dataFiles_t dataFiles ) 				{ m_dataFiles = dataFiles; }
	int	set_initialMesh			( mesh_t* );
	int	set_initialMesh			( );
	int	set_adaptedMesh			( mesh_t* );
	int	set_adaptedMesh			( );
	int	set_meshAdapter 		( );
	int	set_meshAdapter 		( meshAdapter_t* adapter );


	std::string 		get_projectName			() { return m_projectName;}
	std::string		get_nodesFile			() { return m_nodesFile; }
	std::string		get_connectivitiesFile		() { return m_connectivitiesFile; }
	std::string		get_fileOfAdaptedMeshNodes	() { return m_fileOfAdaptedMeshNodes;}

	solutionMethod 		get_option_solutionMethod 	( ) { return m_option_solutionMethod; }
	DistributionType	get_option_weightingFunctionType( ) { return m_option_weightingFunctionType; }
	computingMode		get_option_computingMode 	( ) { return m_option_computingMode; }
	LinearSolverType 	get_option_linearSolverType	( ) { return m_option_linearSolverType; }
	size_t			get_polynomialDegree		( ) { return m_polynomialDegree; }
	T			get_weightingFunctionSpan	( ) { return m_weightingFunctionSpan; }
	T			get_cutOffMultiplier		( ) { return m_cutOffMultiplier; }
	T			get_interpolationFactor		( ) { return m_interpolationFactor; }
	size_t			get_domainCardinality		( ) { return m_domainCardinality; }


	std::set<size_t>	get_perturbedNodeSets 		( ) { return m_perturbedNodeSets; }
	dataSources_t 		get_dataSources 		( ) { return m_dataSources; }
	transformations_t 	get_transformations 		( ) { return m_transformations; }
	translations_t 		get_translations 		( ) { return m_translations; }
	rotationsEulerAngles_t	get_rotationsEulerAngles	( ) { return m_rotationsEulerAngles; }
	rotationsAxisAngle_t	get_rotationsAxisAngle		( ) { return m_rotationsAxisAngle; }
	rotationsQuaternion_t	get_rotationsQuaternion		( ) { return m_rotationsQuaternion; }
	dataFiles_t		get_dataFiles			( ) { return m_dataFiles; }
	weightingData_t		get_weightingData		( ) { return m_weightingData; }


	mesh_t* 		get_initialMesh 		( ) { return m_initialMesh; }
	mesh_t* 		get_adaptedMesh 		( ) { return m_adaptedMesh; }
	meshAdapter_t* 		get_meshAdapter 		( ) { return m_meshBoundaryAdapter; }


	int			write_nodesInfo 		( std::ostream& out );
	int			write_boundaryNodesDisplacements( std::ostream& out );
	int			write_adaptedNodesFile 		( );
	int			write_invalidElementsFile 	( );
	int			write_adaptationProcessInfo	( std::ostream& out );
	int			read_perturbedNodesData		( std::istream& streamIn );
	int			read_weightingData		( std::istream& streamIn );

};


template class app_meshAdapter<double>;
template class app_meshAdapter<float>;

#endif
