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

#include"BasisFunctions.h"

template<class T>
BasisFunctions<T>::BasisFunctions(int dimensions, int degree)
{ set_dimensions(dimensions); set_degree(degree); }


template<class T>
int BasisFunctions<T>::get_degree() {return m_degree;}
template<class T>
int BasisFunctions<T>::get_dimensions() {return m_dimensions;}
template<class T>
void BasisFunctions<T>::set_degree(int newDegree) {m_degree = newDegree;}
template<class T>
void BasisFunctions<T>::set_dimensions(int newDimensions) {m_dimensions = newDimensions;}

template<class T>
void BasisFunctions<T>::operator () (T* monomials, T evalPoint)
{ univariateMonomials(monomials, evalPoint, m_degree); }

template<class T>
void BasisFunctions<T>::operator () (T* monomials, T evalPoint_x, T evalPoint_y)
{ bivariateMonomials(monomials, evalPoint_x, evalPoint_y, m_degree );  }

template<class T>
void BasisFunctions<T>::operator () (T* monomials, T evalPoint_x, T evalPoint_y, T evalPoint_z)
{ trivariateMonomials(monomials, evalPoint_x, evalPoint_y, evalPoint_z, m_degree );  }
