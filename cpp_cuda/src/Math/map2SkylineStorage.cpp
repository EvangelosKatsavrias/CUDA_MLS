#include<limits>
#include<numeric>

template<typename T>
void find_skyline_mapVector( T* A, int* maxa, int n)
{
	maxa[0] = 1; maxa[1] = 2;
	for ( int j=2; j<n; j++ ) for ( int i=1; i<j; i++ ) if ( A[i*n+j] != 0 ) { maxa[j] = maxa[j-1] +j-i+1; break; }
}

template void find_skyline_mapVector( double* A, int* maxa, int n);
template void find_skyline_mapVector( float* A, int* maxa, int n);


int map2SkylineStorage(int i, int j, int* maxa)
{
	if (i > j) { int temp(i); i=j; j=temp; }; // if the queried element belongs to the lower triangular part of the matrix, consider the corresponding reciprocal element
	int firstElement = j -( maxa[j+1] -maxa[j] ); // the row of the most top non zero element of the queried matrix column
	if ( firstElement > i ) return std::numeric_limits<int>::quiet_NaN();
	return maxa[j] +(j-i);
}
