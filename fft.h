// 
// Fast Fourier Transform using Cooley-Tukey algorithm
//
// 2022, Jonathan Tainer
//

#ifndef FFT_H
#define FFT_H

#include <complex.h>

//
// Simple FFT implementation (no parameter caching)
//

// Reads time domain signal and writes to freq domain signal
// n (number of samples) must be a power of 2
void fft(float* sig_td, float complex* sig_fd, unsigned int n);

// Reads freq domain signal and writes to time domain signal
// n (number of samples) must be a power of 2
void ifft(float* sig_td, float complex* sig_fd, unsigned int n);

// In-place implementations of the transform functions
// Input signal stored in sig is overwritten
// n must be a power of 2
void fft_inpl(float complex* sig, unsigned int n);
void ifft_inpl(float complex* sig, unsigned int n);

// 
// FFT Implementation using cached parameters
//
// Faster, but uses more memory
// 27% speed increase on i7-10700k
// Increase will probably be greater on systems with no trig instructions
//

typedef struct fft_param_cache {
	unsigned int n;
	float complex* lut;
} fft_param_cache;

// Allocates memory and calculates lut values
// n must be a power of 2
fft_param_cache fft_create_param_cache(unsigned int n);

// Only calculates lut values
// param.n must be initialized to a power of 2
// param.lut must point to a sufficient amount of memory
void fft_init_param_cache(fft_param_cache param);

// Only deallocates memory, doesnt set param.lut to NULL
void fft_destroy_param_cache(fft_param_cache param);

// In-place implementation of transform using cached parameters
void fft_inpl_cached(float complex* sig, fft_param_cache param);

// In-place implementation of inverse transform using cached parameters
void ifft_inpl_cached(float complex* sig, fft_param_cache param);

#endif
