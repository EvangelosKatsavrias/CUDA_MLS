#include"linearizers.h"

template<class T>
T linearize_rowWiseZeroBased_2d(T i, T j, T size_x)
{
	return i + j*size_x;
}

template<class T>
T linearize_rowWiseZeroBased_3d(T i, T j, T k, T size_x, T size_y)
{
	return i + j*size_x + k*size_x*size_y;
}

#include"linearizers.ipp"
