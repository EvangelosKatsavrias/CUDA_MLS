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

#ifndef LEASTSQUARESHEADER
#define LEASTSQUARESHEADER


#include<iostream>
#include<vector>
#include<functional>
#include<algorithm>
#include"BasisFunctions.h"
#include"ScatteredData.h"
#include"EvaluationData.h"
#include"StationaryPoints.h"
#include"LinearSolver.h"
#include"WeightingFunction.h"


enum computingMode { CPU, CPU_MT, CUDA, CPU_CUDA, MultiGPU_CUDA, CPU_MultiGPU_CUDA, OPENCL };
std::ostream& operator << ( std::ostream& out, computingMode& mode );
std::istream& operator >> (std::istream& stream, computingMode& var );
std::string get_computingMode_name( computingMode& mode );


template<class T=float>
class LeastSquares 
{
protected:
	ScatteredData<T>*	m_scatteredData;
	EvaluationData<T>*	m_evaluationData;
	BasisFunctions<T>*	m_basisFunctions;
	LinearSolver<T>*	m_linearSolver;
	std::vector<T>		m_coefficients;
	EvaluationData<T>*	m_derivativesData;


	bool			m_coefficientsEval_CompletionFlag = 0;
	bool			m_pointEvaluation_CompletionFlag = 0;
	bool			m_evaluationNeedFlag = 1;
	bool			m_updateNeedFlag = 0;
	std::vector<bool>	m_flags_createdData;
	float			m_evalTime_coefficients =0;
	float			m_evalTime_evalPoints =0;
	clock_t			m_calculationClocks = 0;

public:
	LeastSquares (ScatteredData<T>* scatteredData=0,  EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );
	LeastSquares (LeastSquares<T>* referenceLeastSquaresObject, ScatteredData<T>* scatteredData=0,  EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData = 0 );
	~LeastSquares();

	int			set_scatteredData (ScatteredData<T>*);
	int			set_evaluationData (EvaluationData<T>*);
	int			set_derivativesData (EvaluationData<T>*);
	int			set_basisFunctions (BasisFunctions<T>* newBasisFunctions);
	int			set_linearSolver (LinearSolver<T>*);
	std::vector<T>& 	get_coefficients ();
	ScatteredData<T>* 	get_scatteredData () {return m_scatteredData;}
	BasisFunctions<T>* 	get_basisFunctions () {return m_basisFunctions;}
	EvaluationData<T>* 	get_evaluationData () {return m_evaluationData;}
	EvaluationData<T>* 	get_derivativesData () {return m_derivativesData;}
	LinearSolver<T>* 	get_linearSolver () {return m_linearSolver;}
	float			get_coeffEvalTime () {return m_evalTime_coefficients;}
	float			get_pointsEvalTime () {return m_evalTime_evalPoints;}
	clock_t			get_calculationClocks () {return m_calculationClocks;}

	void 			evaluatePoints ();
	void			evaluateDerivatives ();
	void 		virtual find_coefficients () =0;
	void 		virtual evaluationMethod () =0;
	void		virtual	derivativesEvaluationMethod () = 0;
	std::string 	virtual getLSType () =0;
	void 		virtual streamAllProperties (std::ostream& out);
	void 		virtual copy_properties (LeastSquares<T>* referenceLeastSquaresObject);
};

template class LeastSquares<float>;
template class LeastSquares<double>;


// =====================================================

template<class T=float>
class ClassicLeastSquares: public LeastSquares<T>
{

public:
	using LeastSquares<T>::LeastSquares;
	using LeastSquares<T>::copy_properties;

	void find_coefficients ();
	void evaluationMethod ();
	void derivativesEvaluationMethod ();
	std::string getLSType ();
};

template class ClassicLeastSquares<float>;
template class ClassicLeastSquares<double>;


// =====================================================

template<class T=float>
class WeightedLeastSquares: public LeastSquares<T>
{

protected:
	WeightingFunction<T>*	m_weightingFunction;
	StationaryPoints<T>*	m_stationaryPoints;
	std::vector<bool>	m_flags_createdData_w;
	computingMode		m_computingMode;

public:
	WeightedLeastSquares (WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, StationaryPoints<T>* stationaryPoints=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData = 0 );
	WeightedLeastSquares (LeastSquares<T>* referenceLeastSquaresObject, WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, StationaryPoints<T>* stationaryPoints=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );
	WeightedLeastSquares (WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject, WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, StationaryPoints<T>* stationaryPoints=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData = 0 );
	~WeightedLeastSquares ();

	int			set_weightingFunction (WeightingFunction<T>*);
	int			set_stationaryPoints (StationaryPoints<T>* newStationaryPoints);
	int			set_computingMode (computingMode);
	StationaryPoints<T>*  	get_stationaryPoints () {return m_stationaryPoints;}
	WeightingFunction<T>* 	get_weightingFunction () {return m_weightingFunction;}
	computingMode		get_computingMode () {return m_computingMode;}
	int 			get_numOfStationaryPoints () { return m_stationaryPoints->get_numOfNodes(); }

	void 			streamAllProperties (std::ostream& out);
	void 			copy_properties (WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject);
	void 			find_coefficients ();
	void 			evaluationMethod ();
	void 			derivativesEvaluationMethod ();
	std::string 		getLSType ();
};

template class WeightedLeastSquares<float>;
template class WeightedLeastSquares<double>;



// ==========================================================================
template<class T=float>
class WeightedLeastSquares_LCS: virtual public WeightedLeastSquares<T>
{
protected:
	using WeightedLeastSquares<T>::m_weightingFunction;
	using WeightedLeastSquares<T>::m_stationaryPoints;
	using WeightedLeastSquares<T>::m_computingMode;
	void  findCoefficients_LCS1D();
	void  findCoefficients_LCS2D();
	void  findCoefficients_LCS3D();

public:
	using WeightedLeastSquares<T>::WeightedLeastSquares;
	using WeightedLeastSquares<T>::get_numOfStationaryPoints;

	void 		find_coefficients ();
	void 		evaluationMethod ();
	void 		derivativesEvaluationMethod ();
	std::string 	getLSType ();
};

template class WeightedLeastSquares_LCS<float>;
template class WeightedLeastSquares_LCS<double>;


// ==========================================================================
template<class T=float>
class MovingLeastSquares: virtual public WeightedLeastSquares<T>
{
protected:
	using WeightedLeastSquares<T>::m_weightingFunction;
	using WeightedLeastSquares<T>::m_stationaryPoints;
	using WeightedLeastSquares<T>::m_computingMode;

public:
	MovingLeastSquares (WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );
	MovingLeastSquares (LeastSquares<T>* referenceLeastSquaresObject, WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );
	MovingLeastSquares (WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject, WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );

	using		WeightedLeastSquares<T>::find_coefficients;
	void 		streamAllProperties(std::ostream&);
	void 		evaluationMethod ();
	void 		derivativesEvaluationMethod ();
	std::string 	getLSType ();
};

template class MovingLeastSquares<float>;
template class MovingLeastSquares<double>;


// ==========================================================================
template<class T=float>
class MovingLeastSquares_LCS: public MovingLeastSquares<T>, public WeightedLeastSquares_LCS<T>
{

protected:
	using WeightedLeastSquares<T>::m_weightingFunction;
	using WeightedLeastSquares<T>::m_stationaryPoints;
	using WeightedLeastSquares<T>::m_computingMode;

public:
	MovingLeastSquares_LCS(WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );
	MovingLeastSquares_LCS(LeastSquares<T>* referenceLeastSquaresObject, WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );
	MovingLeastSquares_LCS(WeightedLeastSquares<T>* referenceWeightedLeastSquaresObject, WeightingFunction<T>* weightingFunction=0, ScatteredData<T>* scatteredData=0, EvaluationData<T>* evaluationData=0, BasisFunctions<T>* basisFunctions=0, LinearSolver<T>* linearSolver=0, EvaluationData<T>* derivativesData=0 );

	using MovingLeastSquares<T>::streamAllProperties;
	using WeightedLeastSquares_LCS<T>::find_coefficients;
	void 		evaluationMethod();
	void 		evaluationMethod_1Field();
	void 		evaluationMethod_2Field();
	void 		evaluationMethod_3Field();
	void 		derivativesEvaluationMethod();
	void 		derivativesEvaluationMethod_2Field();
	std::string 	getLSType();
};

template class MovingLeastSquares_LCS<float>;
template class MovingLeastSquares_LCS<double>;


// ========================================================================
template <typename T>
std::ostream& operator<< (std::ostream& out, LeastSquares<T>& givenLS)
{
	givenLS.streamAllProperties(out);
	return out << *(givenLS.get_evaluationData());
}

#endif
