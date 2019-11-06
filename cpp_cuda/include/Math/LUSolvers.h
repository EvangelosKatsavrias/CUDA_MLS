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

