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

#ifndef WEIGHTINGFUNCTIONHEADER
#define WEIGHTINGFUNCTIONHEADER

#include<cuda_runtime.h>
#include<cmath>
#include<string>
#include<iostream>

template<class T>
__host__ __device__
T wendlandDistribution(T reciprocal_span, T distance, T cutOff = 1, T interpolationFactor = 0)
{
	T r = distance*reciprocal_span;
	if ( r > cutOff ) return T(0);
	T prop = (T(1)-r); prop *= prop;
	return ( prop*prop ) *( T(4)*r +T(1) ) / ( std::pow(r, interpolationFactor) +1e-12);
}


template<class T>
__host__ __device__
T gaussianDistribution(T reciprocal_span, T distance, T cutOffMult=1, T interpolationFactor = 0)
{
	T r = distance*reciprocal_span;
//	if ( r > cutOffMult ) return T(0);
//	r *= r;
	return exp( -r ) /( std::pow(r, interpolationFactor) +1e-12 );
}


template<class T>
__host__ __device__
T inverseDistanceDistribution(T reciprocal_span, T distance, T cutOffMult = 1, T interpolationFactor = 0)
{
	T r = distance*reciprocal_span;
//	if ( r > cutOffMult ) return T(0);
//	r *= r;
//	return T(1)/( r*r + 1e-12 );
	return T(1)/( std::pow(r, interpolationFactor) + 1e-12 );

}


enum DistributionType{ Wendland, Gaussian, InverseDistance, NonPredefined };
std::istream& operator >> (std::istream& stream, DistributionType& var );
std::ostream& operator << (std::ostream& stream, DistributionType& var );
std::string get_DistributionType_name( DistributionType& type );


template<class T=float>
class WeightingFunction
{
protected:
	T 			m_span;
	T			m_reciprocalSpan;
	DistributionType 	m_distributionType;
	T 			m_cutOffMultiplier=1;
	T			m_interpolationFactor = 0;
	T*			m_variableSpanData = 0;
	T*			m_variableCutOffData = 0;
	T*			m_variableInterpolationData = 0;
	size_t			m_variableCounter = 0;
	size_t			m_numOfVariableData = 0;
	T(*m_evalFunc)(T, T, T, T);

public:
	__device__ __host__ WeightingFunction(DistributionType type =Wendland, T span=1, T cutOffMultiplier = 1, T interpolationFactor = 0)
	{ set_span(span); set_distributionType(type); set_cutOffMultiplier(cutOffMultiplier); set_interpolationFactor(interpolationFactor); }

	__device__ __host__ WeightingFunction( WeightingFunction& wfun )
	{ set_span( wfun.m_span ); set_distributionType( wfun.m_distributionType ); set_cutOffMultiplier( wfun.m_cutOffMultiplier ); set_interpolationFactor( wfun.m_interpolationFactor ); 
       		m_variableSpanData = wfun.m_variableSpanData; m_variableCutOffData = wfun.m_variableCutOffData; m_variableInterpolationData = wfun.m_variableInterpolationData; m_numOfVariableData = wfun.m_numOfVariableData; }


	__device__ __host__ ~WeightingFunction() {}


	__device__ __host__ void set_span(T newSpan)
	{ m_span = newSpan; m_reciprocalSpan = T(1)/newSpan; }

	__device__ __host__ void set_variableSpan(T* newSpans, size_t numOfSpans )
	{ m_variableSpanData = newSpans; m_numOfVariableData = numOfSpans; m_variableCounter = 0; }

	__device__ __host__ void set_variableCutOff(T* newCutOff, size_t numOfData )
	{ m_variableCutOffData = newCutOff; m_numOfVariableData = numOfData; m_variableCounter = 0; }

	__device__ __host__ void set_variableInterpolationFactor(T* newFactor, size_t numOfData )
	{ m_variableInterpolationData = newFactor; m_numOfVariableData = numOfData; m_variableCounter = 0; }

	__device__ __host__ void set_numOfVariableData(size_t numOfData )
	{ m_numOfVariableData = numOfData; m_variableCounter = 0; }

	__device__ __host__ void set_interpolationFactor(T newFactor)
	{ m_interpolationFactor = newFactor; }


	__device__ __host__ void set_cutOffMultiplier(T newCutOffMultiplier)
	{
		m_cutOffMultiplier = newCutOffMultiplier;
		if (newCutOffMultiplier==-1)
			switch (m_distributionType)
			{
				case Wendland:		{m_cutOffMultiplier =  1; break;}
				case Gaussian:		{m_cutOffMultiplier =  1; break;}
				case InverseDistance:	{m_cutOffMultiplier =  1; break;}
			};
	}


	__device__ __host__ void set_distributionType(DistributionType type)
	{
		switch (type)
		{
			case Wendland:		{ m_evalFunc = wendlandDistribution; break; }
			case Gaussian:		{ m_evalFunc = gaussianDistribution; break; }
			case InverseDistance:	{ m_evalFunc = inverseDistanceDistribution; break; }
		};
		m_distributionType = type;
	}


	__host__ __device__ void set_distributionFunction( T(*distributionFunction)(T, T, T, T) )
	{ m_evalFunc = distributionFunction; m_distributionType = NonPredefined; }


	__device__ __host__ T (*get_evaluationFunction())(T, T, T, T) { return m_evalFunc; }


	__host__ __device__ T operator() (T radialDistance)
	{
		//return m_evalFunc( m_reciprocalSpan, radialDistance, m_cutOffMultiplier);
		if ( m_numOfVariableData ) {
			if (m_variableCounter == m_numOfVariableData) m_variableCounter = 0;

			set_span( m_variableSpanData[ m_variableCounter ] );
			set_cutOffMultiplier( m_variableCutOffData[ m_variableCounter ] );
			set_interpolationFactor( m_variableInterpolationData[ m_variableCounter ] );
			m_variableCounter++;

		}
		switch ( m_distributionType ) {
			case Wendland: 		return wendlandDistribution( m_reciprocalSpan, radialDistance, m_cutOffMultiplier, m_interpolationFactor);
			case Gaussian: 		return gaussianDistribution( m_reciprocalSpan, radialDistance, m_cutOffMultiplier, m_interpolationFactor);
			case InverseDistance: 	return inverseDistanceDistribution( m_reciprocalSpan, radialDistance, m_cutOffMultiplier, m_interpolationFactor); }
	}


	T 			get_span		()		{return m_span;}
	T 			get_span		( size_t index ){return m_variableSpanData[ index ]; }
	T 			get_cutOffMultiplier	()		{return m_cutOffMultiplier;}
	T 			get_cutOffMultiplier	( size_t index ){return m_variableCutOffData[index]; }
	T 			get_interpolationFactor ()		{return m_interpolationFactor;}
	T 			get_interpolationFactor	( size_t index ){return m_variableInterpolationData[index]; }
	DistributionType 	get_distributionType	()		{return m_distributionType;}
	std::string		get_distributionTypeName();
	size_t 			get_counterValue	( )		{return m_variableCounter; }

};

template class WeightingFunction<float>;
template class WeightingFunction<double>;

#endif
