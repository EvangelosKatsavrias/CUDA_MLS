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

#ifndef VECTORPLOTERHEADER
#define VECTORPLOTERHEADER

#include<vector>
#include<set>
#include<iostream>

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& givenVector)
{
	for (auto& element: givenVector) out << element << std::endl;

	return out;
}


//#include<iterator>
// if ( !v.empty() ) std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));


template <typename T>
std::istream& operator>> (std::istream& streamIn, std::vector<T>& givenVector)
{
	for (auto& element: givenVector) streamIn >> element;

	return streamIn;
}


template <typename T>
std::ostream& operator<< (std::ostream& out, const std::set<T>& givenVector)
{
	for (auto& element: givenVector) out << element << std::endl;

	return out;
}


template <typename T>
std::istream& operator>> (std::istream& streamIn, std::set<T>& givenVector)
{
	for (auto& element: givenVector) streamIn >> element;

	return streamIn;
}



#endif
