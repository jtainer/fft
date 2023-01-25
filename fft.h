// 
// Fast Fourier Transform using Cooley-Tukey algorithm
//
// 2022, Jonathan Tainer
//

#ifndef FFT_H
#define FFT_H

#include <complex.h>

// Reads time domain signal from input buffer x
// Writes freq domain signal to output buffer X
// n (number of samples) must be a power of 2
void fft(float* x, float complex* X, unsigned int n);

// Reads freq. domain signal from X
// Writes time domain signal to x
// n (number of samples) must be a power of 2
void ifft(float* c, float complex* X, unsigned int n);

#endif
