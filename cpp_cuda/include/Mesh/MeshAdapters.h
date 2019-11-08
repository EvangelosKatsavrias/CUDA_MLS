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

#ifndef MESHADAPTERSHEADERFILE
#define MESHADAPTERSHEADERFILE


#include<math.h>
#include<cmath>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<time.h>
#include<chrono>
#include<memory>
#include"myBLAS.h"
#include"Mesh.h"
#include"TopologyAdapters.h"
#include"LeastSquares.h"


enum displacementsSource { ProvidedInFile, AffineTransformation, Translation, Rotation_EulerAngles, Rotation_AxisAngle, Rotation_Quaternion };
enum solutionMethod { MLS_LCS };
std::istream& operator >> (std::istream& stream, displacementsSource& var );
std::istream& operator >> (std::istream& stream, solutionMethod& var );
std::ostream& operator << (std::ostream& stream, displacementsSource& var );
std::ostream& operator << (std::ostream& stream, solutionMethod& var );
std::string get_displacementsSource_name( displacementsSource& source );
std::string get_solutionMethod_name( solutionMethod& method );


template<class nodalComponentType = float, size_t numOfConnectivityNodes = 3, class flagDataType = int >
class MeshDisplacedBoundaryAdapter: public TopologyAdapter<nodalComponentType, numOfConnectivityNodes, flagDataType>
{

public:
	typedef Mesh<nodalComponentType, numOfConnectivityNodes, flagDataType> 			mesh_t;
	typedef TopologyAdapter<nodalComponentType, numOfConnectivityNodes, flagDataType> 	topologyAdapter_t;
using	typename topologyAdapter_t::perturbations_t;
using	typename topologyAdapter_t::perturbedNodes_t;
using	typename topologyAdapter_t::adaptedNodes_t;
using	typename topologyAdapter_t::badConnectivities_t;
using	typename topologyAdapter_t::flagValues_t;
using	meshMetric_t = typename mesh_t::meshMetric_t;


protected:
	mesh_t*					m_initialMesh		= this->m_initialTopology;
	mesh_t*					m_adaptedMesh 		= this->m_adaptedTopology;
	perturbations_t*& 			m_displacements		= this->m_perturbations;
	badConnectivities_t& 			m_badElements 		= this->m_badConnectivities;
	std::vector<nodalComponentType>		m_weightingFunctionSpans;
	std::vector<nodalComponentType>		m_weightingFunctionCutOffLimits;
	std::vector<nodalComponentType>		m_weightingFunctionInterpolationFactors;
	EvaluationData<nodalComponentType>*	m_derivativesData 	= 0;


	int			m_polynomialDegree 	= 0;
	nodalComponentType	m_weightingFunctionSpan = 10;
	nodalComponentType	m_cutOffMultiplier 	= 10;
	nodalComponentType	m_interpolationFactor 	= 0;
	solutionMethod		m_solutionMethod	= MLS_LCS;
	DistributionType	m_weightingFunctionType = Wendland;
	computingMode		m_computingMode 	= CPU;
	LinearSolverType	m_linearSolverType 	= ConjugateGradient;


private:
	MeshDisplacedBoundaryAdapter ( 	mesh_t* 		initialMesh,
					mesh_t* 		adaptedMesh,
					mesh_t* 		ptrInitial,
					mesh_t* 		ptrAdapted,
		       			perturbations_t*	displacements,
					perturbedNodes_t*	perturbedNodes,
					flagValues_t		perturbedNodeFlags):
	topologyAdapter_t( initialMesh, adaptedMesh, displacements, perturbedNodes, perturbedNodeFlags ), m_initialMesh(ptrInitial), m_adaptedMesh(ptrAdapted){ }


public:
	MeshDisplacedBoundaryAdapter ( mesh_t* initialMesh, mesh_t* adaptedMesh = 0, perturbations_t* displacements = 0, perturbedNodes_t* perturbedNodes = 0, flagValues_t perturbedBoundaryFlags = flagValues_t(), solutionMethod method=MLS_LCS ):
	MeshDisplacedBoundaryAdapter ( initialMesh, adaptedMesh, initialMesh, adaptedMesh, displacements, perturbedNodes, perturbedBoundaryFlags )
	{ m_solutionMethod = method; }
	~MeshDisplacedBoundaryAdapter ( ) { if ( m_derivativesData ) delete m_derivativesData;  }


	int			adapt_nodes 			( );
	int			find_badElements 		( );


	int			set_polynomialDegree 		( int num) 					{ m_polynomialDegree = num; }
	int			set_weightingFunctionSpan	( nodalComponentType span ) 			{ m_weightingFunctionSpan = span; }
	int			set_cutOffMultiplier 		( nodalComponentType cutOff ) 			{ m_cutOffMultiplier = cutOff; }
	int			set_interpolationFactor		( nodalComponentType factor ) 			{ m_interpolationFactor = factor; }
	int			set_solutionMethod 		( solutionMethod solMethod = MLS_LCS ) 		{ m_solutionMethod = solMethod; }
	int			set_weightingFunctionType	( DistributionType weightingFunction = Wendland ) { m_weightingFunctionType = weightingFunction; }
	int			set_computingMode 		( computingMode computingMode = CPU ) 		{ m_computingMode = computingMode; }
	int			set_linearSolverType 		( LinearSolverType solverType = ConjugateGradient ) { m_linearSolverType = solverType; }

	int			set_weightingFunctionSpans	( flagDataType nodeSet, nodalComponentType value );
	int			set_weightingFunctionCutOffLimits( flagDataType nodeSet, nodalComponentType value );
	int			set_weightingFunctionInterpolationFactors( flagDataType nodeSet, nodalComponentType value );
	int			set_derivativesData( EvaluationData<nodalComponentType>* data ) {  m_derivativesData = data; return 1; }


	std::string 		get_typename			( ) { return "Mesh Displaced Boundary Adapter\n";}
	int			get_polynomialDegree		( ) { return m_polynomialDegree;}
	nodalComponentType	get_weightingFunctionSpan	( ) { return m_weightingFunctionSpan; }
	nodalComponentType	get_cutOffMultiplier		( ) { return m_cutOffMultiplier; }
	nodalComponentType	get_interpolationFactor		( ) { return m_interpolationFactor; }
	solutionMethod 		get_solutionMethod 		( ) { return m_solutionMethod; }
	DistributionType	get_weightingFunctionType	( ) { return m_weightingFunctionType; }
	computingMode		get_computingMode 		( ) { return m_computingMode; }
	LinearSolverType	get_linearSolverType 		( ) { return m_linearSolverType; }
	badConnectivities_t& 	get_badElements 		( ) { return this->get_badConnectivities(); }
	std::vector<nodalComponentType>& get_weightingFunctionSpans		( ) { return m_weightingFunctionSpans; }
	std::vector<nodalComponentType>& get_weightingFunctionCutOffLimits	( ) { return m_weightingFunctionCutOffLimits; }
	std::vector<nodalComponentType>& get_weightingFunctionInterpolationFactors( ) { return m_weightingFunctionInterpolationFactors; }
	EvaluationData<nodalComponentType>* get_derivativesData(  ) {  return m_derivativesData; }

	
	int 		virtual	write_adaptationParameters 	( std::ostream& out) { out << "Defined vector of displacements\n"; }
	int 		virtual	write_nodeDisplacements 	( std::string file);
	int 		virtual	write_nodeDisplacements		( std::ostream& out);
};






// ================================= adaptation algorithms
template<class T, size_t T2, class T3>
int MeshDisplacedBoundaryAdapter<T, T2, T3>::adapt_nodes ( ) {


typename topologyAdapter_t::adaptedNodes_t& 	adaptedNodes 		= this->get_adaptedNodes()[0];
typename topologyAdapter_t::perturbedNodes_t& 	perturbedNodes 		= this->get_perturbedNodes()[0];
T** 						domainComponents 	= this->m_initialTopology->get_allComponents();
T** 						adaptedDomainComponents = this->m_adaptedTopology->get_allComponents();
T** 						displacementComponents	= this->m_displacements->get_allComponents();
ScatteredData<T> 				samplePoints		(this->m_initialTopology->get_numOfComponents(), m_displacements->get_numOfComponents(), perturbedNodes.size() );
EvaluationData<T> 				evaluationData		(this->m_initialTopology->get_numOfComponents(), m_displacements->get_numOfComponents(), adaptedNodes.size() );
						

m_derivativesData = new EvaluationData<T>(this->m_initialTopology->get_numOfComponents(), 3, this->m_initialTopology->get_numOfNodes() );


clock_t clocks; size_t counter(0);


// set the field data
samplePoints.set_AllFieldComponents( displacementComponents );
evaluationData.set_AllFieldComponents();
m_derivativesData->set_AllFieldComponents();


// set the domain data
samplePoints.set_AllDomainComponents();
evaluationData.set_AllDomainComponents();
T	**samplePointNodes 	= samplePoints.get_AllDomainComponents(),
	**evaluationDataNodes 	= evaluationData.get_AllDomainComponents();

for ( size_t component = 0; component < this->m_initialTopology->get_numOfComponents(); component++ ) {
	counter = 0;
	for ( auto& node: perturbedNodes ) samplePointNodes[component][counter++] = domainComponents[component][node];
	counter = 0;
	for ( auto& node: adaptedNodes ) evaluationDataNodes[component][counter++] = domainComponents[component][node];
}


m_derivativesData->set_AllDomainComponents( this->m_initialTopology->get_allComponents() );


switch ( m_solutionMethod ) {

	case MLS_LCS: {

		BasisFunctions<T> 		basisFunctions		( this->m_initialTopology->get_numOfComponents(), m_polynomialDegree ); 
		WeightingFunction<T> 		weightingFunction	( m_weightingFunctionType, m_weightingFunctionSpan, m_cutOffMultiplier, m_interpolationFactor ); 
		LinearSolver<T> 		solver			( m_linearSolverType, basisFunctions.get_numOfMonomials() );

		T* spans = m_weightingFunctionSpans.data();
		T* limits = m_weightingFunctionCutOffLimits.data();
		T* interpolationFactors = m_weightingFunctionInterpolationFactors.data();

		weightingFunction.set_variableSpan( spans, this->get_numOfPerturbedNodes() );
		weightingFunction.set_variableCutOff( limits, this->get_numOfPerturbedNodes() );
		weightingFunction.set_variableInterpolationFactor( interpolationFactors, this->get_numOfPerturbedNodes() );

		MovingLeastSquares_LCS<T> mls( &weightingFunction, &samplePoints, &evaluationData, &basisFunctions, &solver, m_derivativesData );
		mls.set_computingMode( m_computingMode );
		mls.evaluatePoints();
//		mls.evaluateDerivatives();


		// set the adaptations to the adapted mesh
		T** adaptations = evaluationData.get_AllFieldComponents();
		for ( size_t component = 0; component < this->m_initialTopology->get_numOfComponents(); component++ ) {
			counter = 0;
			for ( auto& node: perturbedNodes ) 	adaptedDomainComponents[component][node] += displacementComponents[component][ counter++ ];
			counter = 0;
			for ( auto& node: adaptedNodes ) 	{ adaptedDomainComponents[component][node] += adaptations[component][ counter++ ]; 
				//std::cout << adaptations[component][counter-1] << std::endl; 
				}
		}

		clocks = mls.get_calculationClocks() + clock();
		this->set_adaptationCompletionFlag(1);
		this->set_adaptationTime( (float)clocks/CLOCKS_PER_SEC );
		mls.streamAllProperties(std::cout); std::cout << std::endl;

		break;
	}
}


}



// ===========================================


template<class T1, size_t T2, class T3>
int MeshDisplacedBoundaryAdapter<T1, T2, T3>::set_weightingFunctionSpans( T3 nodeSet, T1 value )
{
	ContainerBase<size_t, std::vector<size_t> >* 	perturbedNodes = this->get_perturbedNodes();
	std::set<T3>					perturbedFlagValues = this->get_perturbedNodeFlags();
	std::map<T3, size_t> 				flagsMap = this->get_perturbedFlag2NodesIndex();

	if ( perturbedFlagValues.find( nodeSet ) == perturbedFlagValues.end() ) { std::cout << "The nodes under the flag " << nodeSet << ", are not perturbed." ; throw; }

	std::set<T3> flags = {nodeSet}; size_t numOfNodes = m_initialMesh->get_numOfNodesWithFlag( flags );

	m_weightingFunctionSpans.resize( this->get_numOfPerturbedNodes() );

	size_t nodesIndex = flagsMap[nodeSet];

	for ( size_t node = nodesIndex; node < nodesIndex+numOfNodes; node++ ) 
		m_weightingFunctionSpans[node] = value;

}


template<class T1, size_t T2, class T3>
int MeshDisplacedBoundaryAdapter<T1, T2, T3>::set_weightingFunctionCutOffLimits( T3 nodeSet, T1 value )
{
	ContainerBase<size_t, std::vector<size_t> >* 	perturbedNodes = this->get_perturbedNodes();
	std::set<T3>					perturbedFlagValues = this->get_perturbedNodeFlags();
	std::map<T3, size_t> 				flagsMap = this->get_perturbedFlag2NodesIndex();

	if ( perturbedFlagValues.find( nodeSet ) == perturbedFlagValues.end() ) { std::cout << "The nodes under the flag " << nodeSet << ", are not perturbed." ; throw; }

	std::set<T3> flags = {nodeSet}; size_t numOfNodes = m_initialMesh->get_numOfNodesWithFlag( flags );

	m_weightingFunctionCutOffLimits.resize( this->get_numOfPerturbedNodes() );

	size_t nodesIndex = flagsMap[nodeSet];

	for ( size_t node = nodesIndex; node < nodesIndex+numOfNodes; node++ )
		m_weightingFunctionCutOffLimits[node] = value;

}


template<class T1, size_t T2, class T3>
int MeshDisplacedBoundaryAdapter<T1, T2, T3>::set_weightingFunctionInterpolationFactors( T3 nodeSet, T1 value )
{
	ContainerBase<size_t, std::vector<size_t> >* 	perturbedNodes = this->get_perturbedNodes();
	std::set<T3>					perturbedFlagValues = this->get_perturbedNodeFlags();
	std::map<T3, size_t> 				flagsMap = this->get_perturbedFlag2NodesIndex();

	if ( perturbedFlagValues.find( nodeSet ) == perturbedFlagValues.end() ) { std::cout << "The nodes under the flag " << nodeSet << ", are not perturbed." ; throw; }

	std::set<T3> flags = {nodeSet}; size_t numOfNodes = m_initialMesh->get_numOfNodesWithFlag( flags );

	m_weightingFunctionInterpolationFactors.resize( this->get_numOfPerturbedNodes() );

	size_t nodesIndex = flagsMap[nodeSet];

	for ( size_t node = nodesIndex; node < nodesIndex+numOfNodes; node++ )
		m_weightingFunctionInterpolationFactors[node] = value;

}




template<class T1, size_t T2, class T3>
int MeshDisplacedBoundaryAdapter<T1, T2, T3>::find_badElements() {

	m_badElements.clear();
	meshMetric_t* metrics_initial = m_initialMesh->get_metrics();
	meshMetric_t* metrics_adapted = m_adaptedMesh->get_metrics();


for (size_t elemIndex = 0; elemIndex < m_initialMesh->get_numOfConnectivities(); elemIndex++ )
	if ( metrics_initial[0][elemIndex]*metrics_adapted[0][elemIndex] <= T1(0) )  m_badElements.push_back(elemIndex);

}



template<class T1, size_t T2, class T3>
int MeshDisplacedBoundaryAdapter<T1, T2, T3>::write_nodeDisplacements ( std::string file) {
	std::ofstream out; out.exceptions(std::ios_base::badbit);
	out.open(file); write_nodeDisplacements (out); out.close(); }


template<class T1, size_t T2, class T3>
int MeshDisplacedBoundaryAdapter<T1, T2, T3>::write_nodeDisplacements (std::ostream& out)
{

out << this->m_initialTopology->get_numOfBoundaryNodes() << std::endl;

for ( auto& node: this->m_initialTopology->get_boundaryNodes() ) {
	out << node.first;
	for (size_t component = 0; component < m_displacements->get_numOfComponents(); component++)
		out << "\t\t" << m_displacements[component][node.second];
	out << std::endl;
}

}


#endif
