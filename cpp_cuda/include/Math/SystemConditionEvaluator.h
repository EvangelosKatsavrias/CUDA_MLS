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

#ifndef SYSTEMCONDITIONEVALUATORHEADER
#define SYSTEMCONDITIONEVALUATORHEADER 

template<class T> 
class SystemConditionEvaluator
{

	T*	m_data;
	int	m_numOfColumns;
	T	m_determinant_preestimate;
	T	m_HadamardCondNumber;
	T	m_spectralCondNumber;

public:
	SystemConditionEvaluator(T* data, int numOfColumns): m_data(data), m_numOfColumns(numOfColumns) {  }
	~SystemConditionEvaluator() {}


	static T preestimate_spectralCondNumber(T* data, int numOfColumns) { 
		T maxVal(0);
		for (int i=0; i<numOfColumns*numOfColumns; i+numOfColumns)
			if ( maxVal < data[i] ) maxVal = data[i];

		
	}

};


template class SystemConditionEvaluator<float>;
template class SystemConditionEvaluator<double>;

#endif
