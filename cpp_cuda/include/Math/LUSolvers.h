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

#include"LUDecomposers.h"

template<typename T>
void LU_Crout_solver(T* L, T* U, T* b, T* x, int n);

template<typename T>
void LU_Doolittle_solver(T* L, T* U, T* b, T* x, int n);

template<typename T>
void LDU_solver(T* L, T* U, T* D, T* b, T* x, int n);

template<typename T>
void LDL_solver(T* L, T* D, T* b, T* x, int n);

template<typename T>
void Cholesky_solver(T* L, T* b, T* x, int n);

template<typename T>
void Cholesky_solver_banded(T* L, T* b, T* x, int n, int n_b);

template<typename T>
void Cholesky_solver_skyline(T* L, T* b, T* x, int n, int* n_max);

