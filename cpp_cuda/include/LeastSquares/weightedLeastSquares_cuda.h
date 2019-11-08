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

#include<cuda.h>
#include"LinearSolver.h"
#include"WeightingFunction.h"
//#include"BasisFunctions.h"


// 3D
template<class T>
void cudaVersion_findCoefficientsLCS3D(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& solver, WeightingFunction<T>& wFun, size_t cudaDeviceID );

template<class T>
void cudaVersion_findCoefficientsLCS3D2(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID );


// 2D
template<class T>
void cudaVersion_findCoefficientsLCS2D(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& solver, WeightingFunction<T>& wFun, size_t cudaDeviceID );

template<class T>
void cudaVersion_findCoefficientsLCS2D2(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& solver, WeightingFunction<T>& wFun, size_t cudaDeviceID );

template<class T>
void cudaVersion_findCoefficientsLCS2DConstantMemory2(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& solver, WeightingFunction<T>& wFun, size_t cudaDeviceID );


template<class T>
void cudaVersion_findCoefficientsLCS2D3(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& solver, WeightingFunction<T>& wFun, size_t cudaDeviceID );

template<class T>
void cudaVersion_findCoefficientsLCS2D4(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& solver, WeightingFunction<T>& wFun, size_t cudaDeviceID );



// 1D
template<class T>
void cudaVersion_findCoefficients(T span, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* coefficients);

template<class T>
void cudaVersion_evaluationMethod(int numOfEvaluationPoints, T* evaluationPoints, T span, int numOfStationaryPoints, T* stationaryPoints, int numOfMonomials, T* host_coefficients, T* host_evaluations);

template<class T>
void cudaVersion_findCoefficientsLCS(T span, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* coefficients);

template<class T>
void cudaVersion_evaluationMethodLCS(int numOfEvaluationPoints, T* evaluationPoints, T span, int numOfStationaryPoints, T* stationaryPoints, int numOfMonomials, T* host_coefficients, T* host_evaluations);

template<class T>
void cudaVersion_evaluationMethodMovingLS(int numOfEvaluationPoints, T* evaluationPoints, int numOfMonomials, T* host_coefficients, T* host_evaluations);
