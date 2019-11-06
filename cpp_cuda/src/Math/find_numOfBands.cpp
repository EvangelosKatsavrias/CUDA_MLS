#include"matrix_memoryStorage.h"


template<class T>
int find_numOfBands(T* A, int n)
{ for (int i = n-1; i >= 0; i-- ) if ( A[i] !=0 ) return i+1; }


template int find_numOfBands(float* A, int n);
template int find_numOfBands(double* A, int n);
