#include"LeastSquares.h"
#include"weightedLeastSquares_cuda.h"
#include"myBLAS.h"
#include<cmath>

template<class T>
std::string MovingLeastSquares<T>::getLSType() { return "Moving least squares"; }


template <typename T>
void MovingLeastSquares<T>::streamAllProperties(std::ostream& out)
{
	LeastSquares<T>::streamAllProperties(out);
	out << "Weighting function type: " << m_weightingFunction->get_distributionTypeName();
	out << ", Weighting function span: " << m_weightingFunction->get_span();
	out << ", Weighting function cutoff multiplier: " << m_weightingFunction->get_cutOffMultiplier();
	out << ", Weighting function interpolation factor: " << m_weightingFunction->get_interpolationFactor() << std::endl;
	out << "Computing mode: " << m_computingMode << std::endl;
}


template<class T>
MovingLeastSquares<T>::MovingLeastSquares(WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData ): WeightedLeastSquares<T>(weightingFunction, scatteredData, evaluationData, nullptr, basisFunctions, linearSolver, derivativesData)
{
	delete(this->m_stationaryPoints);
	this->m_stationaryPoints = new StationaryPoints<T>(this->m_evaluationData->get_domain());
}

template<class T>
MovingLeastSquares<T>::MovingLeastSquares(LeastSquares<T>* referenceLeastSquaresObject, WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData  ): WeightedLeastSquares<T>(referenceLeastSquaresObject, weightingFunction, scatteredData, evaluationData, nullptr, basisFunctions, linearSolver, derivativesData )
{
	delete(this->m_stationaryPoints);
	this->m_stationaryPoints = new StationaryPoints<T>(this->m_evaluationData->get_domain());
}


template<class T>
MovingLeastSquares<T>::MovingLeastSquares(WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject, WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData  ): WeightedLeastSquares<T>(referenceWeightedLeastSquaresObject, weightingFunction, scatteredData, evaluationData, nullptr, basisFunctions, linearSolver, derivativesData )
{
	this->m_stationaryPoints = new StationaryPoints<T>(this->m_evaluationData->get_domain());
	this->m_flags_createdData_w[1] = 1;
}


template<class T>
void MovingLeastSquares<T>::evaluationMethod()
{

	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	T* coefficients_local			= this->m_coefficients.data();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);


	if (m_computingMode)
		cudaVersion_evaluationMethodMovingLS(numOfEvaluationPoints, evaluationPoints, numOfMonomials, coefficients_local, evaluations);
	else
	{
		T b[numOfMonomials];
		for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
		{
			basisFunctions(b, evaluationPoints[evaluationPointIndex]);
			evaluations[evaluationPointIndex] = dotProductOfVectors(b, b+numOfMonomials, coefficients_local+evaluationPointIndex*numOfMonomials);
		}
	}

}


template<class T>
void MovingLeastSquares<T>::derivativesEvaluationMethod()
{

	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	T* coefficients_local			= this->m_coefficients.data();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);


//	if (m_computingMode)
//		cudaVersion_evaluationMethodMovingLS(numOfEvaluationPoints, evaluationPoints, numOfMonomials, coefficients_local, evaluations);
//	else
//	{
//		T b[numOfMonomials];
//		for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
//		{
//			basisFunctions(b, evaluationPoints[evaluationPointIndex]);
//			evaluations[evaluationPointIndex] = dotProductOfVectors(b, b+numOfMonomials, coefficients_local+evaluationPointIndex*numOfMonomials);
//		}
//	}

}
