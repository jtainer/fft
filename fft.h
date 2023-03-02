// 
// Fast Fourier Transform using Cooley-Tukey algorithm
//
// 2022, Jonathan Tainer
//

#ifndef FFT_H
#define FFT_H

#include <complex.h>

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

#endif
