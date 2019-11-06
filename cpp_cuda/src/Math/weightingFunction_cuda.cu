#include"weightingFunction_cuda.h"


template<class T>
__device__
T weightingFunction3(T reciprocal_span, T distance)
{
	T r = distance*reciprocal_span;
	if ( r > 1 ) return 0;
	T prop = T(1)-r; prop *= prop;
	return (prop*prop) *( T(4)*r +T(1) );
}

template __device__ float weightingFunction3(float reciprocal_span, float distance);
template __device__ double weightingFunction3(double reciprocal_span, double distance);
