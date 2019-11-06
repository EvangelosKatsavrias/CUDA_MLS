#include"LeastSquares.h"
#include"weightedLeastSquares_cuda.h"
#include"myBLAS.h"
#include<cmath>


template<class T>
std::string WeightedLeastSquares<T>::getLSType() { return "Weighted least squares"; }


template<class T>
void WeightedLeastSquares<T>::copy_properties(WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject)
{
	LeastSquares<T>::copy_properties(referenceWeightedLeastSquaresObject);
	m_weightingFunction 		= referenceWeightedLeastSquaresObject->m_weightingFunction;
	m_stationaryPoints 		= referenceWeightedLeastSquaresObject->m_stationaryPoints;
	m_computingMode			= referenceWeightedLeastSquaresObject->m_computingMode;
	std::fill(m_flags_createdData_w.begin(), m_flags_createdData_w.end(), 0); 
}


template <typename T>
void WeightedLeastSquares<T>::streamAllProperties(std::ostream& out)
{
	LeastSquares<T>::streamAllProperties(out);
	out << "Number of stationary points: " << get_numOfStationaryPoints() << std::endl;
	out << "Weighting function type: " << m_weightingFunction->get_distributionTypeName();
	out << ", Weighting function span: " << m_weightingFunction->get_span();
	out << ", Weighting function cutoff multiplier: " << m_weightingFunction->get_cutOffMultiplier() << std::endl;
	out << ", Weighting function interpolation factor: " << m_weightingFunction->get_interpolationFactor() << std::endl;
	out << "Computing mode: " << m_computingMode << std::endl;
}


template<class T>
int WeightedLeastSquares<T>::set_weightingFunction(WeightingFunction<T>* newWeightingFunction)
{ 
	if (m_flags_createdData_w[0]) { delete(m_weightingFunction); m_flags_createdData_w[0] = 0;}
	m_weightingFunction = newWeightingFunction;
	return 1;
}


template<class T>
int WeightedLeastSquares<T>::set_stationaryPoints(StationaryPoints<T>* newStationaryPoints)
{ 
	if (m_flags_createdData_w[1]) { delete(m_stationaryPoints); m_flags_createdData_w[1] = 0;}
	m_stationaryPoints = newStationaryPoints;
	return 1;
}

template<class T>
int WeightedLeastSquares<T>::set_computingMode(computingMode mode)
{	m_computingMode = mode; }


template<class T>
void buildSystemAtSamplePoint(T samplePoint, T sampleValue, T* Sum_bbT, T* Sum_bf, T stationaryPoint, int numOfMonomials, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction)
{
	T bbT[numOfMonomials*numOfMonomials];
	T b[numOfMonomials];
	T theta;
	T theta_b[numOfMonomials];

	theta = weightingFunction(fabs(samplePoint - stationaryPoint));
	basisFunctions(b, samplePoint);
	multiplicationOfVectorbyScalar(b, b+numOfMonomials, theta_b, theta);
	tensorProductOfVectors(theta_b, theta_b+numOfMonomials, b, b+numOfMonomials, bbT);
	sumOfVectors(bbT, bbT+numOfMonomials*numOfMonomials, Sum_bbT, Sum_bbT);
	multiplicationOfVectorbyScalar(theta_b, theta_b+numOfMonomials, theta_b, sampleValue);
	sumOfVectors(theta_b, theta_b+numOfMonomials, Sum_bf, Sum_bf);
}


template<class T>
void findCoefficientsAtStationaryPoint(T stationaryPoint, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, LinearSolver<T> linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, T* Sum_bf)
{
	T Sum_bbT[numOfMonomials*numOfMonomials]; std::fill(Sum_bbT, Sum_bbT+numOfMonomials*numOfMonomials, T(0));
	
	for (int samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++)
		buildSystemAtSamplePoint<T>(samplePoints[samplePointIndex], sampleValues[samplePointIndex], Sum_bbT, Sum_bf, stationaryPoint, numOfMonomials, basisFunctions, weightingFunction);

	linearSolver(Sum_bbT, Sum_bf);
}

template<class T>
void WeightedLeastSquares<T>::find_coefficients()
{
	int numOfSamplePoints 			= this->m_scatteredData->get_numOfNodes();
	T* samplePoints 			= this->m_scatteredData->get_domainComponent(0);
	T* sampleValues 			= this->m_scatteredData->get_fieldComponent(0);

	T* stationaryPoints			= this->m_stationaryPoints->get_component(0);
	int numOfStationaryPoints		= this->get_numOfStationaryPoints();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	LinearSolver<T>& linearSolver		= *(this->m_linearSolver);
	linearSolver.set_numberOfColumns(numOfMonomials);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);


	if (m_computingMode)
	{
		this->m_coefficients.resize(numOfStationaryPoints*numOfMonomials);
		cudaVersion_findCoefficients<T>(weightingFunction.get_span(), numOfStationaryPoints, stationaryPoints, numOfSamplePoints, samplePoints, sampleValues, numOfMonomials, this->m_coefficients.data());
	}
	else
	{
		this->m_coefficients.reserve(numOfStationaryPoints*numOfMonomials);
		T Sum_bf[numOfMonomials];
		for (int stationaryPointIndex = 0; stationaryPointIndex < get_numOfStationaryPoints(); stationaryPointIndex++)
		{
			std::fill(Sum_bf, Sum_bf+numOfMonomials, T(0));
			findCoefficientsAtStationaryPoint<T>(stationaryPoints[stationaryPointIndex], numOfSamplePoints, samplePoints, sampleValues, numOfMonomials, linearSolver, basisFunctions, weightingFunction, Sum_bf);

			for (int index=0; index<numOfMonomials;index++) this->m_coefficients.push_back(Sum_bf[index]);
		}
	}
}


template<class T>
T evaluateSinglePoint(T* stationaryPoints, T evaluationPoint, T* coefficients_local, int numOfStationaryPoints, int numOfMonomials, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction)
{
	T weightPerStationaryPoint[numOfStationaryPoints], weight_sum(0);
	for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPoints; stationaryPointIndex++)
	{
		weightPerStationaryPoint[stationaryPointIndex] = weightingFunction(fabs(evaluationPoint - stationaryPoints[stationaryPointIndex]));
		weight_sum += weightPerStationaryPoint[stationaryPointIndex];
	}

	T coefficients[numOfMonomials], temp[numOfMonomials];  std::fill(coefficients, coefficients+numOfMonomials, 0);
	for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPoints; stationaryPointIndex++)
	{
		multiplicationOfVectorbyScalar(coefficients_local+stationaryPointIndex*numOfMonomials, coefficients_local+(stationaryPointIndex+1)*numOfMonomials, temp, weightPerStationaryPoint[stationaryPointIndex]/weight_sum);
		sumOfVectors(coefficients, coefficients+numOfMonomials, temp, coefficients);
	}

	T b[numOfMonomials];
	basisFunctions(b, evaluationPoint);
	return  dotProductOfVectors(b, b+numOfMonomials, coefficients);
}


template<class T>
void WeightedLeastSquares<T>::evaluationMethod()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	T* stationaryPoints			= this->m_stationaryPoints->get_component(0);
	T* coefficients_local			= this->m_coefficients.data();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);
	int numOfStationaryPoints		= get_numOfStationaryPoints();

	if (m_computingMode)
		cudaVersion_evaluationMethod(numOfEvaluationPoints, evaluationPoints, weightingFunction.get_span(), numOfStationaryPoints, stationaryPoints, numOfMonomials, coefficients_local, evaluations);
	else
	{
		for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
			evaluations[evaluationPointIndex] = evaluateSinglePoint<T>(stationaryPoints, evaluationPoints[evaluationPointIndex], coefficients_local, numOfStationaryPoints, numOfMonomials, basisFunctions, weightingFunction);
	}
}



template<class T>
void WeightedLeastSquares<T>::derivativesEvaluationMethod()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	T* stationaryPoints			= this->m_stationaryPoints->get_component(0);
	T* coefficients_local			= this->m_coefficients.data();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);
	int numOfStationaryPoints		= get_numOfStationaryPoints();

//	if (m_computingMode)
//		cudaVersion_evaluationMethod(numOfEvaluationPoints, evaluationPoints, weightingFunction.get_span(), numOfStationaryPoints, stationaryPoints, numOfMonomials, coefficients_local, evaluations);
//	else
//	{
//		for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
//			evaluations[evaluationPointIndex] = evaluateSinglePoint<T>(stationaryPoints, evaluationPoints[evaluationPointIndex], coefficients_local, numOfStationaryPoints, numOfMonomials, basisFunctions, weightingFunction);
//	}
}
