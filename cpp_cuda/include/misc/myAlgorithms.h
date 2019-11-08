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

#ifndef MYALGORITHMSHEADER
#define MYALGORITHMSHEADER

#include<algorithm>
#include<functional>


template <typename T>
struct unorderLess
{
    bool operator () (const std::pair<T, T>& lhs, const std::pair<T, T>& rhs) const
    {
        const auto lhs_order = lhs.first < lhs.second ? lhs : std::tie(lhs.second, lhs.first);
        const auto rhs_order = rhs.first < rhs.second ? rhs : std::tie(rhs.second, rhs.first);

        return lhs_order < rhs_order;
    }
};


template<class T>
class unary_criterion
{


public:
	virtual bool operator() ( T values ) { return 1; }

	T	m_value;

};



template<class T>
class Isless_equal : public unary_criterion<T> 
{
public:
	T m_limit;

	Isless_equal(T lim = 0): m_limit(lim) { this->m_value = m_limit; }
	Isless_equal( Isless_equal& ref) { m_limit = ref.m_limit; }

	bool operator() (T value ) { return value <= m_limit; }
};


template<class T>
class IsInRange : public unary_criterion<T> 
{
public:
	T m_lowerBound;
	T m_upperBound;


	IsInRange(T lowerBound = 0, T upperBound = 0): m_lowerBound( lowerBound ), m_upperBound(upperBound) {  }
	IsInRange( IsInRange & ref ) { 
		m_lowerBound = ref.m_lowerBound; 
		m_upperBound = ref.m_upperBound; 
	}

	bool operator() (T value ) { return m_lowerBound <= value <= m_upperBound; }
};


#endif
