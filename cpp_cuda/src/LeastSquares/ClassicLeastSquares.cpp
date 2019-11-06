#include"LeastSquares.h"
#include"myBLAS.h"
#include<cmath>


template<class T>
std::string ClassicLeastSquares<T>::getLSType() { return "Classic least squares"; }


template<class T>
void ClassicLeastSquares<T>::find_coefficients()
{
	int numOfSamplePoints 			= this->m_scatteredData->get_numOfNodes();
	T* samplePoints 			= this->m_scatteredData->get_domainComponent(0);
	T* sampleValues 			= this->m_scatteredData->get_fieldComponent(0);

	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	LinearSolver<T>& linearSolver		= *(this->m_linearSolver);
	linearSolver.set_numberOfColumns(numOfMonomials);


	T b[numOfMonomials], bbT[numOfMonomials*numOfMonomials], Sum_bbT[numOfMonomials*numOfMonomials], Sum_bf[numOfMonomials];
	std::fill(Sum_bbT, Sum_bbT+numOfMonomials*numOfMonomials, T(0));
	std::fill(Sum_bf, Sum_bf+numOfMonomials, T(0));


	for (int samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++)
	{
		basisFunctions(b, samplePoints[samplePointIndex]);
		tensorProductOfVectors(b, b+numOfMonomials, b, b+numOfMonomials, bbT);
		sumOfVectors(bbT, bbT+numOfMonomials*numOfMonomials, Sum_bbT, Sum_bbT);
		multiplicationOfVectorbyScalar(b, b+numOfMonomials, b, sampleValues[samplePointIndex]);
		sumOfVectors(b, b+numOfMonomials, Sum_bf, Sum_bf);
	}

	linearSolver(Sum_bbT, Sum_bf); this->m_coefficients.assign(Sum_bf, Sum_bf+numOfMonomials);
}


template<class T>
void ClassicLeastSquares<T>::evaluationMethod()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	T* coefficients				= this->m_coefficients.data();

	T b[numOfMonomials];

	for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
	{
		basisFunctions(b, evaluationPoints[evaluationPointIndex]);
		evaluations[evaluationPointIndex] = dotProductOfVectors(b, b+numOfMonomials, coefficients);
	}
}

	
template<class T>
void ClassicLeastSquares<T>::derivativesEvaluationMethod()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	T* coefficients				= this->m_coefficients.data();

	T b[numOfMonomials];

//	for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
//	{
//		basisFunctions(b, evaluationPoints[evaluationPointIndex]);
//		evaluations[evaluationPointIndex] = dotProductOfVectors(b, b+numOfMonomials, coefficients);
//	}
}

