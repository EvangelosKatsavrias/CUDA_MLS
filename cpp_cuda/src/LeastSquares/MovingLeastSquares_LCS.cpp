#include"LeastSquares.h"
#include"myBLAS.h"
#include<cmath>
#include"vectorPlotOutStreamer.h"

template<class T>
MovingLeastSquares_LCS<T>::MovingLeastSquares_LCS(WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData ):
WeightedLeastSquares<T>(weightingFunction, scatteredData, evaluationData, nullptr, basisFunctions, linearSolver, derivativesData), MovingLeastSquares<T>(weightingFunction, scatteredData, evaluationData, basisFunctions, linearSolver, derivativesData)
{}



template<class T>
MovingLeastSquares_LCS<T>::MovingLeastSquares_LCS(LeastSquares<T>* referenceLeastSquaresObject, WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData ):
WeightedLeastSquares<T>(referenceLeastSquaresObject, weightingFunction, scatteredData, evaluationData, nullptr, basisFunctions, linearSolver, derivativesData), MovingLeastSquares<T>(referenceLeastSquaresObject, weightingFunction, scatteredData, evaluationData, basisFunctions, linearSolver, derivativesData)
{}


template<class T>
MovingLeastSquares_LCS<T>::MovingLeastSquares_LCS(WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject, WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData ):
WeightedLeastSquares<T>(referenceWeightedLeastSquaresObject, weightingFunction, scatteredData, evaluationData, nullptr, basisFunctions, linearSolver, derivativesData), MovingLeastSquares<T>(referenceWeightedLeastSquaresObject, weightingFunction, scatteredData, evaluationData, basisFunctions, linearSolver, derivativesData)
{}


template<class T>
std::string MovingLeastSquares_LCS<T>::getLSType()
{ return "Moving least squares in LCS"; }


template<class T>
void MovingLeastSquares_LCS<T>::evaluationMethod_1Field()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_numOfMonomials();
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* coefficients_local			= this->m_coefficients.data();

	for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
		evaluations[evaluationPointIndex] = *(coefficients_local+evaluationPointIndex*numOfMonomials);

}


template<class T>
void MovingLeastSquares_LCS<T>::evaluationMethod_2Field()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_numOfMonomials();
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations_1			= this->m_evaluationData->get_fieldComponent(0);
	T* evaluations_2			= this->m_evaluationData->get_fieldComponent(1);
	T* coefficients_local			= this->m_coefficients.data();


	for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
	{
		evaluations_1[evaluationPointIndex] = *(coefficients_local+evaluationPointIndex*2*numOfMonomials);
		evaluations_2[evaluationPointIndex] = *(coefficients_local+(evaluationPointIndex*2+1)*numOfMonomials);
	}
}


template<class T>
void MovingLeastSquares_LCS<T>::evaluationMethod_3Field()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_numOfMonomials();
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations_1			= this->m_evaluationData->get_fieldComponent(0);
	T* evaluations_2			= this->m_evaluationData->get_fieldComponent(1);
	T* evaluations_3			= this->m_evaluationData->get_fieldComponent(2);
	T* coefficients_local			= this->m_coefficients.data();

	for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++)
	{
		evaluations_1[evaluationPointIndex] = *(coefficients_local+evaluationPointIndex*3*numOfMonomials);
		evaluations_2[evaluationPointIndex] = *(coefficients_local+(evaluationPointIndex*3+1)*numOfMonomials);
		evaluations_3[evaluationPointIndex] = *(coefficients_local+(evaluationPointIndex*3+2)*numOfMonomials);
	}
}


template<class T>
void MovingLeastSquares_LCS<T>::evaluationMethod()
{
	switch ( this->m_evaluationData->get_numOfFieldsPerPoint() )
	{
		case 1: { evaluationMethod_1Field(); break; }
		case 2: { evaluationMethod_2Field(); break; }
		case 3: { evaluationMethod_3Field(); break; }
	}
}








template<class T>
T weightDerivativeFunction(T distance, T reciprocal_span)
{
//	if ( distance*reciprocal_span > 1 ) return T(0);
	T r = distance/reciprocal_span;
	r *= r;
	return -T(4)/( r*r + 1e-12 );
}


template<class T>
void findDerivativesAtStationaryPoint_LCS2D(T stationaryPoint_x, T stationaryPoint_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, LinearSolver<T>& linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, std::vector<T>& samplePointsDistance, T* dudx, T* dvdy, T* shear )
{
	if ( numOfMonomials == 1) {

		T sumTheta(0), F1x(0), F1y(0), F1dx(0), F1dy(0), F1xdy(0), F1ydx(0), sumdThetax(0), sumdThetay(0);
		for ( int samplePointIndex=0; samplePointIndex < numOfSamplePoints; samplePointIndex++ )
		{
			T dx = (samplePoints_x[samplePointIndex]-stationaryPoint_x);
			T dy = (samplePoints_y[samplePointIndex]-stationaryPoint_y);
			T distance = sqrt( dx*dx +dy*dy );

			if ( distance > 1e-6 ) {
				T theta = weightingFunction( distance );
				T theta_dx = theta*dx/( distance*distance );
				T theta_dy = theta*dy/( distance*distance );
	
				sumTheta += theta;
				sumdThetax += theta_dx;
				sumdThetay += theta_dy;
	
	
				F1x += theta*sampleValues_x[samplePointIndex];
				F1y += theta*sampleValues_y[samplePointIndex];
	
				F1dx += theta_dx*sampleValues_x[samplePointIndex];
				F1dy += theta_dy*sampleValues_y[samplePointIndex];
				
				F1xdy += theta_dy*sampleValues_x[samplePointIndex];
				F1ydx += theta_dx*sampleValues_y[samplePointIndex];
			}
			else { 
				F1x = 0; F1y = 0; F1dx = 0; F1dy = 0; F1xdy = 0; F1ydx = 0; 
				sumTheta = 1; sumdThetax = 1; sumdThetay = 1;
				break; 
			}

		}

		*dudx = F1dx/sumTheta - F1x/(sumTheta*sumTheta)*sumdThetax;
		*dvdy = F1dy/sumTheta - F1y/(sumTheta*sumTheta)*sumdThetay;
		*shear = (F1ydx/sumTheta - F1y/(sumTheta*sumTheta)*sumdThetax)
			+ (F1xdy/sumTheta - F1x/(sumTheta*sumTheta)*sumdThetay);
	}


	if ( numOfMonomials > 1 ) {
		T A11(0), A12(0), A13(0), A22(0), A23(0), A33(0), F1x(0), F2x(0), F3x(0), F1y(0), F2y(0), F3y(0);
		for ( int samplePointIndex=0; samplePointIndex < numOfSamplePoints; samplePointIndex++ )
		{
			T dx = (samplePoints_x[samplePointIndex]-stationaryPoint_x);
			T dy = (samplePoints_y[samplePointIndex]-stationaryPoint_y);
			T distance = sqrt( dx*dx +dy*dy );
	
			T theta = weightingFunction( distance );
			A11 += theta;
			A12 += theta*dx;
			A13 += theta*dy;
			A23 += theta*dx*dy;
			A22 += theta*dx*dx;
			A33 += theta*dy*dy;
	
			F1x += theta*sampleValues_x[samplePointIndex];
			F2x += theta*sampleValues_x[samplePointIndex]*dx;
			F3x += theta*sampleValues_x[samplePointIndex]*dy;
	
			F1y += theta*sampleValues_y[samplePointIndex];
			F2y += theta*sampleValues_y[samplePointIndex]*dx;
			F3y += theta*sampleValues_y[samplePointIndex]*dy;
		}
	
	
		T P11 = A12/A11;
		T P12 = A13/A11;
		T P21 = ( A23 -P12*A12 )/( A22 -P11*A12 );
	
		T C3x = ( F3x -P12*F1x -P21*(F2x -P11*F1x) ) / ( A33 -P12*A13 -P21*( A23-P11*A13 ) );
		T C2x = ( F2x -P11*F1x -C3x*(A23 -P11*A13) ) / ( A22 -P11*A12 );
		T C1x = ( F1x -C2x*A12 -C3x*A13 )/A11;
	
		T C3y = ( F3y -P12*F1y -P21*(F2y -P11*F1y) ) / ( A33 -P12*A13 -P21*( A23-P11*A13 ) );
		T C2y = ( F2y -P11*F1y -C3y*(A23 -P11*A13) ) / ( A22 -P11*A12 );
		T C1y = ( F1y -C2y*A12 -C3y*A13 )/A11;
	
	
		*dudx = C2x;
		*dvdy = C3y;
		*shear = C3x + C2y;

	}

//	std::cout << *shear << std::endl;

/*  ======= 1D case
	T A11(0), A12(0), A22(0), F1(0), F2(0); 
	for ( int samplePointIndex=0; samplePointIndex < numOfSamplePoints; samplePointIndex++ )
	{
		T dx = (samplePoints_x[samplePointIndex]-stationaryPoint_x);
		T distance = std::abs(dx);

		T theta = weightingFunction( distance );
		A11 += theta;
		A12 += theta*dx;
		A22 += theta*dx*dx;
		F1 += theta*sampleValues_x[samplePointIndex];
		F2 += theta*sampleValues_x[samplePointIndex]*dx;
	}

	T P1 = A12/A11;

	*dudx = (F2x -P1*F1x)/(A22 -P1*A12) ;
	*/
}


template<class T>
void cpuVersion_findDerivatives_LCS2D(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, T* dudx, T* dvdy, T* shear, LinearSolver<T>& linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, bool progressIndexFlag )
{
	T distance;
	std::vector<T> samplePointsDistance; samplePointsDistance.resize( numOfSamplePoints );
	for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPoints; stationaryPointIndex++)
	{
		int counter(0); std::fill(samplePointsDistance.begin(), samplePointsDistance.end(), T(-1) );
		for ( size_t samplePointIndex = 0; samplePointIndex<numOfSamplePoints; samplePointIndex++ ) {
			distance = distance2D(stationaryPoints_x[stationaryPointIndex], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex], samplePoints_y[samplePointIndex]);
			if ( distance < weightingFunction.get_span( samplePointIndex ) ) { samplePointsDistance[samplePointIndex]=distance; counter++; }
		}
		findDerivativesAtStationaryPoint_LCS2D ( stationaryPoints_x[stationaryPointIndex], stationaryPoints_y[stationaryPointIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, linearSolver, basisFunctions, weightingFunction, samplePointsDistance, dudx+stationaryPointIndex, dvdy+stationaryPointIndex, shear+stationaryPointIndex );
	}

}


template<class T>
void MovingLeastSquares_LCS<T>::derivativesEvaluationMethod_2Field()
{
	int numOfSamplePoints 			= this->m_scatteredData->get_numOfNodes();
	T* samplePoints_x 			= this->m_scatteredData->get_domainComponent(0);
	T* samplePoints_y 			= this->m_scatteredData->get_domainComponent(1);
	T* sampleValues_x 			= this->m_scatteredData->get_fieldComponent(0);
	T* sampleValues_y 			= this->m_scatteredData->get_fieldComponent(1);

	T* stationaryPoints_x			= this->m_derivativesData->get_domainComponent(0);
	T* stationaryPoints_y			= this->m_derivativesData->get_domainComponent(1);
	int numOfStationaryPoints		= this->m_derivativesData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_numOfMonomials();
	int polDegree				= this->m_basisFunctions->get_degree();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	LinearSolver<T>& linearSolver		= *(this->m_linearSolver);
	linearSolver.set_numberOfColumns(numOfMonomials);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);

	this->m_derivativesData->set_fieldComponent(0);
	T* dudx	= this->m_derivativesData->get_fieldComponent(0);
	T* dvdy	= this->m_derivativesData->get_fieldComponent(1);
	T* shear= this->m_derivativesData->get_fieldComponent(2);

	cpuVersion_findDerivatives_LCS2D<T>(numOfStationaryPoints, stationaryPoints_x, stationaryPoints_y, numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, dudx, dvdy, shear, linearSolver, basisFunctions, weightingFunction, 1 );

}


template<class T>
void MovingLeastSquares_LCS<T>::derivativesEvaluationMethod()
{
	switch ( this->m_evaluationData->get_numOfFieldsPerPoint() )
	{
		case 1: { break; }
		case 2: { derivativesEvaluationMethod_2Field(); break; }
		case 3: { break; }
	}
}
