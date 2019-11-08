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

#include"LeastSquares.h"

template<class T>
WeightedLeastSquares<T>::WeightedLeastSquares(WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, StationaryPoints<T>* stationaryPoints, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData ): LeastSquares<T>(scatteredData, evaluationData, basisFunctions, linearSolver, derivativesData), m_weightingFunction(weightingFunction), m_stationaryPoints(stationaryPoints)
{
	m_flags_createdData_w = std::vector<bool>(5, 0);
	if (!weightingFunction) { m_weightingFunction = new WeightingFunction<T>; m_flags_createdData_w[0]=1;}
	if (!stationaryPoints) { m_stationaryPoints = new StationaryPoints<T>; m_flags_createdData_w[1]=1;}

	m_computingMode = CPU;
}


template<class T>
WeightedLeastSquares<T>::WeightedLeastSquares(LeastSquares<T>* referenceLeastSquaresObject, WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, StationaryPoints<T>* stationaryPoints, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData ): LeastSquares<T>(referenceLeastSquaresObject, scatteredData, evaluationData, basisFunctions, linearSolver, derivativesData), m_weightingFunction(weightingFunction), m_stationaryPoints(stationaryPoints)
{
	m_flags_createdData_w = std::vector<bool>(5, 0);

	if (!weightingFunction) { m_weightingFunction = new WeightingFunction<T>; m_flags_createdData_w[0]=1;}
	if (!stationaryPoints) { m_stationaryPoints = new StationaryPoints<T>; m_flags_createdData_w[1]=1;}

	m_computingMode = CPU;
}


template<class T>
WeightedLeastSquares<T>::WeightedLeastSquares(WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject, WeightingFunction<T>* weightingFunction, ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, StationaryPoints<T>* stationaryPoints, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData ):
LeastSquares<T>(referenceWeightedLeastSquaresObject, scatteredData, evaluationData, basisFunctions, linearSolver, derivativesData )
{
	m_flags_createdData_w = std::vector<bool>(5, 0);

	if (weightingFunction) 	m_weightingFunction = weightingFunction;
	else			m_weightingFunction = referenceWeightedLeastSquaresObject->m_weightingFunction;

	if (stationaryPoints) 	 m_stationaryPoints = stationaryPoints;
	else			 m_stationaryPoints = referenceWeightedLeastSquaresObject->m_stationaryPoints;

	m_computingMode = referenceWeightedLeastSquaresObject->m_computingMode;
}


template<class T>
WeightedLeastSquares<T>::~WeightedLeastSquares()
{
	if (m_flags_createdData_w[0]) delete(this->m_weightingFunction);
	if (m_flags_createdData_w[1]) delete(this->m_stationaryPoints);
}
