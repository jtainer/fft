// 
// Fast Fourier Transform using Cooley-Tukey algorithm
//
// 2022, Jonathan Tainer
//

#include "fft.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Calculate corresponding index for reorder
static unsigned int reverse_bits(unsigned int orig, unsigned int bits) {
	unsigned int flip = 0;
	for (unsigned int i = 0; i < bits; i++) {
		flip <<= 1;
		flip |= orig & 1;
		orig >>= 1;
	}
	return flip;
}

// Reorder data so FFT can be computer in-place
static void fft_reorder(float complex* buf, unsigned int n) {
	unsigned int bits = 0;
	while ((n >> bits) > 1) bits++;
	for (unsigned int i = 0; i < n; i++) {
		unsigned int i0 = i;
		unsigned int i1 = reverse_bits(i0, bits);
		if (i1 <= i0) {
			float complex tmp = buf[i0];
			buf[i0] = buf[i1];
			buf[i1] = tmp;
		}
	}
}

// Input time domain signal in buf
// Convert to freq domain and store result in buf
// n (num of samples) must be a power of 2
// s (stride length) must be 1 for initial function call
static void fft_recurse(float complex* buf, unsigned int n, unsigned int s) {
	if (n == 1) return;

	// Recursively split the DFT in halves
	fft_recurse(buf, n/2, s*2);
	fft_recurse(buf+n/2, n/2, s*2);

	// Combine halves into full DFT
	for (unsigned int k = 0; k < n/2; k++) {
		float complex p = buf[k];
		float complex u = cexpf(I*-2*M_PI*k/n);
		float complex q = u * buf[k+n/2];
		buf[k] = p + q;
		buf[k+n/2] = p - q;
	}
}

// sig_td = input time domain signal (not modified)
// sig_fd = buffer to store resulting freq domain signal
// n = sample count (MUST be a power of 2)
void fft(float* sig_td, float complex* sig_fd, unsigned int n) {
	// Copy time domain input signal into freq domain output signal
	// so the FFT can be calculated in-place (better cache coherency)
	for (unsigned int i = 0; i < n; i++) {
		sig_fd[i] = sig_td[i] + 0*I;
	}

	fft_reorder(sig_fd, n);
	fft_recurse(sig_fd, n, 1);
}

// Input freq domain signal in buf
// Convert to time domain and store result in buf
// n (number of samples) must be a power of 2
// s (stride length) must be 1 for initial function call
static void ifft_recurse(float complex* buf, unsigned int n, unsigned int s) {
	if (n == 1) return;

	// Split DFT in half
	for (unsigned int k = 0; k < n/2; k++) {
		float complex u = cexpf(I*-2.f*M_PI*k/n);
		float complex q = (buf[k] - buf[k+n/2]) / 2.f;
		float complex p = buf[k] - q;
		buf[k] = p;
		buf[k+n/2] = q / u;
	}

	// Recursively split DFT in half until n == 1
	ifft_recurse(buf, n/2, s*2);
	ifft_recurse(buf+n/2, n/2, s*2);
}

// sig_td = buffer to store output time domain signal
// sig_fd = input freq domain signal (modified due to in-place calculation)
// n = sample count (MUST be a power of 2)
void ifft(float* sig_td, float complex* sig_fd, unsigned int n) {
	// Calculate FFT in-place in the freq domain signal buffer
	ifft_recurse(sig_fd, n, 1);
	fft_reorder(sig_fd, n);

	// Copy real component of freq domain buffer to time domain buffer
	for (unsigned int i = 0; i < n; i++) {
		sig_td[i] = crealf(sig_fd[i]);
	}
}

// Complex time domain signal in sig is overwritten
// with freq domain signal
// n = sample count
void fft_inpl(float complex* sig, unsigned int n) {
	fft_reorder(sig, n);
	fft_recurse(sig, n, 1);
}

// Complex freq domain signal in sig is overwritten
// with time domain signal
// n = sample count
void ifft_inpl(float complex* sig, unsigned int n) {
	ifft_recurse(sig, n, 1);
	fft_reorder(sig, n);
}

fft_param_cache fft_create_param_cache(unsigned int n) {
	fft_param_cache param;
	param.n = n;
	param.lut = malloc(sizeof(float complex)*n/2);
	fft_init_param_cache(param);
	return param;
}

void fft_init_param_cache(fft_param_cache param) {
	for (unsigned int k = 0; k < param.n/2; k++) {
		param.lut[k] = cexpf(I*-2*M_PI*k/param.n);
	}
}

void fft_destroy_param_cache(fft_param_cache param) {
	free(param.lut);
}

static void fft_recurse_cached(float complex* buf, fft_param_cache param, unsigned int s) {
	if (param.n == 1) return;

	// Recursively split the DFT in halves
	fft_param_cache param_next = { param.n/2, param.lut };
	fft_recurse_cached(buf, param_next, s*2);
	fft_recurse_cached(buf+param.n/2, param_next, s*2);

	// Combine halves into full DFT
	for (unsigned int k = 0; k < param.n/2; k++) {
		float complex p = buf[k];
		float complex u = param.lut[k*s];
		float complex q = u * buf[k+param.n/2];
		buf[k] = p + q;
		buf[k+param.n/2] = p - q;
	}

}

static void ifft_recurse_cached(float complex* buf, fft_param_cache param, unsigned int s) {
	if (param.n == 1) return;

	// Split DFT in half
	for (unsigned int k = 0; k < param.n/2; k++) {
		float complex u = param.lut[k*s];
		float complex q = (buf[k] - buf[k+param.n/2]) / 2.f;
		float complex p = buf[k] - q;
		buf[k] = p;
		buf[k+param.n/2] = q / u;
	}
	
	// Recursively split DFT in half until n == 1
	fft_param_cache param_next = { param.n/2, param.lut };
	ifft_recurse_cached(buf, param_next, s*2);
	ifft_recurse_cached(buf+param.n/2, param_next, s*2);
}

void fft_inpl_cached(float complex* sig, fft_param_cache param) {
	fft_reorder(sig, param.n);
	fft_recurse_cached(sig, param, 1);
}

void ifft_inpl_cached(float complex* sig, fft_param_cache param) {
	ifft_recurse_cached(sig, param, 1);
	fft_reorder(sig, param.n);
}
