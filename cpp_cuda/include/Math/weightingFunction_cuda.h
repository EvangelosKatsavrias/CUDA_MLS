#include<cuda.h>

template<class T>
__device__
T weightingFunction3(T reciprocal_span, T distance);
