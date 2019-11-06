#include"immintrin.h"

void SSESqrt(const double * pOut, const double * pIn)
{
	
	__m256d in = _mm256_loadu_pd(pIn);
//	__m256d in2 = _mm256_invsqrt_pd(in);
	
}
