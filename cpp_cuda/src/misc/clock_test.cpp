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

#include<memory>
#include<chrono>
#include<time.h>


#include"pmmintrin.h"
inline void SSESqrt_Recip_Times_X( float * pOut, float * pIn )
{
   __m128 in = _mm_load_ss( pIn );
   _mm_store_ss( pOut, _mm_mul_ss( in, _mm_rsqrt_ss( in ) ) );
   // compiles to movss, movaps, rsqrtss, mulss, movss
}


#include<x86intrin.h>
#include<cpuid.h>
long long ReadTSC()
{
	int dummy[4];
	volatile int DontSkip;
	long long clock;
	__cpuid(0,dummy[0], dummy[1],dummy[2],dummy[3]);
	DontSkip = dummy[0];
	clock = __rdtsc();
	return clock;
}


void clock_test()
{

/*
	std::chrono::system_clock sysClock;
	std::chrono::steady_clock stdyClock;
	std::chrono::high_resolution_clock hgResClock;

	std::chrono::time_point<std::chrono::system_clock,std::chrono::nanoseconds> curP = sysClock.now();
	auto point1 = hgResClock.now();
	auto point1_1 = stdyClock.now();


	float in[100], out[100];
	for (int index=0; index<100; index++) in[index]=4;


	clock_t clocks = -clock();
	SSESqrt_Recip_Times_X(out, in);
	for (int index=0; index<100000000;index++) 3*4;
	clocks += clock();
	auto point2 = hgResClock.now();
	auto point2_1 = stdyClock.now();


	auto duration = point2-point1;
	auto duration_1 = point2_1-point1_1;
	time_t tt1 = sysClock.to_time_t(curP);
	std::cout << clocks << ",  " << CLOCKS_PER_SEC << ",  " << (double)clocks/CLOCKS_PER_SEC << std::endl << curP.time_since_epoch().count() << ",  " << ctime(&tt1) << std::endl << duration.count() << std::endl << duration_1.count();
*/

}
