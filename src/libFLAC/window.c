/* libFLAC - Free Lossless Audio Codec library
 * Copyright (C) 2006-2009  Josh Coalson
 * Copyright (C) 2011-2016  Xiph.Org Foundation
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * - Neither the name of the Xiph.org Foundation nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <math.h>
#include "share/compat.h"
#include "FLAC/assert.h"
#include "FLAC/format.h"
#include "private/window.h"

#ifndef FLAC__INTEGER_ONLY_LIBRARY

float __M_PI = (float) M_PI;
float __M_2PI = 2.0f * (float) M_PI;

void FLAC__window_bartlett(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	FLAC__int32 n;
	const float k1 = 2.0f/(float)N;

	if (L & 1) {
		for (n = 0 ; n <= N/2; n++)
			window[n] = k1 * n;
		
		for (; n <= N; n++)
			window[n] = 2.0f - k1 * n;
	}
	else {
		for (n = 0; n <= L/2-1; n++)
			window[n] = k1 * n;
		for (; n <= N; n++)
			window[n] = k1 * n;
	}
}

void FLAC__window_bartlett_hann(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	const float k1 = __M_2PI/(float)N;
	const float k2 = 1.0f/((float)N - 0.5f);
	FLAC__int32 n;

	for (n = 0; n < L; n++)
		window[n] = (FLAC__real)(0.62f - 0.48f * fabsf((float)n * k2 ) - 0.38f * cosf(k1 * ((float)n)));
}

void FLAC__window_blackman(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	FLAC__int32 n;
	const float k1 = 2.0f * __M_PI/(float)N;
	//const float k2 = 2.0f * k1;
	float fsin = 0.0f;
	float fcos = 0.0f;

	for (n = 0; n < L; n++)
		sincosf(k1 * n, &fsin, &fcos);
//		window[n] = (FLAC__real)(0.42f - 0.5f * cosf( k1 * n) + 0.08f * cosf( k2 * n ));
		window[n] = (FLAC__real)(0.42f - 0.5f * fcos + 0.08f * fcos * fcos - fsin * fsin);
}

/* 4-term -92dB side-lobe */
void FLAC__window_blackman_harris_4term_92db_sidelobe(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	FLAC__int32 n;
	const float k1 = __M_2PI/(float) N;

	for (n = 0; n <= N; n++){
		float n_i = k1 * n;
		window[n] = (FLAC__real)(0.35875f - 0.48829f * cosf(n_i ) + 0.14128f * cosf(2.0f * n_i ) - 0.01168f * cosf(3.0f * n_i ));
	}
}

void FLAC__window_connes(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	const float N2 = (float)(N / 2);
	FLAC__int32 n;

	for (n = 0; n <= N; n++) {
		float k1 = 1.0 / N2;
		float k = ((float)n - N2) *  k1;
		k = 1.0f - k * k;
		window[n] = (FLAC__real)(k * k);
	}
}

void FLAC__window_flattop(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	FLAC__int32 n;
	const float k1 = __M_2PI/(float) N;
	for (n = 0; n < L; n++)
		window[n] = (FLAC__real)(0.21557895f - 0.41663158f * cosf(k1 * n ) + 0.277263158f * cosf(2.0f * k1 * n ) - 0.083578947f * cosf(3.0f * k1 * n ) + 0.006947368f * cosf(4.0f * k1 * n ));
}

void FLAC__window_gauss(FLAC__real *window, const FLAC__int32 L, const FLAC__real stddev)
{
	const FLAC__int32 N = L - 1;
	const float N2 = (float)N / 2.0f;
	const float k1 = 1.0f/(stddev * N2);
	FLAC__int32 n;

	for (n = 0; n <= N; n++) {
		const float k = ((float)n - N2) * k1;
		window[n] = (FLAC__real)expf(-0.5f * k * k);
	}
}

void FLAC__window_hamming(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	const float k1 =__M_2PI/ (float)N ;
	FLAC__int32 n;

	for (n = 0; n < L; n++)
		window[n] = (FLAC__real)(0.54f - 0.46f * cosf( k1 * n));
}

void FLAC__window_hann(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	const float k1 =__M_2PI/ (float)N ;
	FLAC__int32 n;

	for (n = 0; n < L; n++)
		window[n] = (FLAC__real)(0.5f - 0.5f * cosf(k1 * n));
}

void FLAC__window_kaiser_bessel(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	FLAC__int32 n;
	const float k1 =__M_2PI/ (float)N ;

	for (n = 0; n < L; n++){
		float n_i = k1 * (float)n;
		window[n] = (FLAC__real)(0.402f - 0.498f * cosf( n_i) + 0.098f * cosf(2.0f * n_i ) - 0.001f * cosf(3.0f * n_i));
	}
}

void FLAC__window_nuttall(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	const float k1 =__M_2PI/ (float)N ;
	FLAC__int32 n;

	for (n = 0; n < L; n++){
		float n_i = k1 * (float)n;
		window[n] = (FLAC__real)(0.3635819f - 0.4891775f*cosf(n_i) + 0.1365995f*cosf(2.0f*n_i) - 0.0106411f*cosf(3.0f*n_i));
	}
}

void FLAC__window_rectangle(FLAC__real *window, const FLAC__int32 L)
{
	FLAC__int32 n;

	for (n = 0; n < L; n++)
		window[n] = 1.0f;
}

void FLAC__window_triangle(FLAC__real *window, const FLAC__int32 L)
{
	FLAC__int32 n;
	const float k1 = 1.0f /((float)L + 1.0f);
	const float k2 = 2.0f * k1;
	const float k3 = 2.0f * (float)L + 2.0f;
//	const float k4 = 2.0f * k1;
	// 2L - 2n +2 -> 2L+2 - 2n
	if (L & 1) {
		for (n = 1; n <= (L+1)/2; n++)
			window[n-1] = n * k2;
		for (; n <= L; n++)
//			window[n-1] = (float)(2 * (L - n + 1)) * k1;
			window[n-1] = k3 -  n * k1;
	}
	else {
		for (n = 1; n <= L/2; n++)
			window[n-1] =  n * k2;
		for (; n <= L; n++)
			window[n-1] = k3 -  n * k1;
//			window[n-1] = (float)(2 * (L - n + 1)) * k1;
	}
}

void FLAC__window_tukey(FLAC__real *window, const FLAC__int32 L, const FLAC__real p)
{
	if (p <= 0.0)
		FLAC__window_rectangle(window, L);
	else if (p >= 1.0)
		FLAC__window_hann(window, L);
	else {
		const FLAC__int32 Np = (FLAC__int32)(p / 2 * L) - 1;
		FLAC__int32 n;
		/* start with rectangle... */
		FLAC__window_rectangle(window, L);
		/* ...replace ends with hann */
		if (Np > 0) {
			const float k1 = __M_PI / (float) Np;
			for (n = 0; n <= Np; n++) {
				window[n] = (FLAC__real)(0.5f - 0.5f * cosf(k1 * n));
				window[L-Np-1+n] = (FLAC__real)(0.5f - 0.5f * cosf( k1 * (n+Np) ));
			}
		}
	}
}

void FLAC__window_partial_tukey(FLAC__real *window, const FLAC__int32 L, const FLAC__real p, const FLAC__real start, const FLAC__real end)
{
	const FLAC__int32 start_n = (FLAC__int32)(start * L);
	const FLAC__int32 end_n = (FLAC__int32)(end * L);
	const FLAC__int32 N = end_n - start_n;
	FLAC__int32 Np, n, i;

	if (p <= 0.0f)
		FLAC__window_partial_tukey(window, L, 0.05f, start, end);
	else if (p >= 1.0f)
		FLAC__window_partial_tukey(window, L, 0.95f, start, end);
	else {

		const float k1 = (__M_2PI * (float)N / p); 
		Np = ((FLAC__int32)p / 2 * N);

		for (n = 0; n < start_n && n < L; n++)
			window[n] = 0.0f;
		for (i = 1; n < (start_n+Np) && n < L; n++, i++)
			window[n] = (FLAC__real)(0.5f - 0.5f * cosf(k1 * i));
		for (; n < (end_n-Np) && n < L; n++)
			window[n] = 1.0f;
		for (i = Np; n < end_n && n < L; n++, i--)
			window[n] = (FLAC__real)(0.5f - 0.5f * cosf(k1 * i));
		for (; n < L; n++)
			window[n] = 0.0f;
	}
}

void FLAC__window_punchout_tukey(FLAC__real *window, const FLAC__int32 L, const FLAC__real p, const FLAC__real start, const FLAC__real end)
{
	const FLAC__int32 start_n = (FLAC__int32)(start * L);
	const FLAC__int32 end_n = (FLAC__int32)(end * L);
	FLAC__int32 Ns, Ne, n, i;

	if (p <= 0.0f)
		FLAC__window_punchout_tukey(window, L, 0.05f, start, end);
	else if (p >= 1.0f)
		FLAC__window_punchout_tukey(window, L, 0.95f, start, end);
	else {

		Ns = (FLAC__int32)(p / 2.0f * start_n);
		Ne = (FLAC__int32)(p / 2.0f * (L - end_n));

		for (n = 0, i = 1; n < Ns && n < L; n++, i++)
			window[n] = (FLAC__real)(0.5f - 0.5f * cos(M_PI * i / Ns));
		for (; n < start_n-Ns && n < L; n++)
			window[n] = 1.0f;
		for (i = Ns; n < start_n && n < L; n++, i--)
			window[n] = (FLAC__real)(0.5f - 0.5f * cos(M_PI * i / Ns));
		for (; n < end_n && n < L; n++)
			window[n] = 0.0f;
		for (i = 1; n < end_n+Ne && n < L; n++, i++)
			window[n] = (FLAC__real)(0.5f - 0.5f * cos(M_PI * i / Ne));
		for (; n < L - (Ne) && n < L; n++)
			window[n] = 1.0f;
		for (i = Ne; n < L; n++, i--)
			window[n] = (FLAC__real)(0.5f - 0.5f * cos(M_PI * i / Ne));
	}
}

void FLAC__window_welch(FLAC__real *window, const FLAC__int32 L)
{
	const FLAC__int32 N = L - 1;
	const float N2 = (float)(N / 2);
	const float k1 = 1.0f/N2;
	FLAC__int32 n;

	for (n = 0; n <= N; n++) {
		const float k = ((float)n - N2) * k1;
		window[n] = (FLAC__real)(1.0f - k * k);
	}
}

#endif /* !defined FLAC__INTEGER_ONLY_LIBRARY */
