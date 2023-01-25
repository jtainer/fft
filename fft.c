// 
// Fast Fourier Transform using Cooley-Tukey algorithm
//
// 2022, Jonathan Tainer
//

#include "fft.h"
#include <math.h>

// Reads time domain signal from input buffer x
// Writes freq domain signal to output buffer X
// n (num of samples) must be a power of 2
// s (stride length) must be 1 for initial function call
static void fft_recurse(float* x, float complex* X, unsigned int n, unsigned int s) {
	if (n == 1) {
		X[0] = x[0] + 0 * I;
	}
	else {
		// Recursively split the DFT in halves
		fft_recurse(x, X, n/2, s*2);
		fft_recurse(x+s, X+n/2, n/2, 2*s);

		// Combine halves into full DFT
		for (unsigned int k = 0; k < n/2; k++) {
			float complex p = X[k];
			float complex u = cexpf(I*-2*M_PI*k/n);
			float complex q = u * X[k+n/2];
			X[k] = p + q;
			X[k+n/2] = p - q;
		}
	}
}

// Wrapper for recursive function call
void fft(float* x, float complex* X, unsigned int n) {
	// TODO: check params for validity here (not inside fft_recurse)
	fft_recurse(x, X, n, 1);
}

// Reads freq domain signal from X
// Writes time domain signal to x
// n (number of samples) must be a power of 2
// s (stride length) must be 1 for initial function call
static void ifft_recurse(float* x, float complex* x, unsigned int n, unsigned int s) {
	if (n == 1) {
		x[0] = crealf(X[0]);
	}
	else {
		for (unsigned int k = 0; k < n/2; k++) {
			float complex u = cexpf(I*-2.f*M_PI*k/n);
			float complex q = (X[k] - X[k+n/2]) / 2.f;
			float complex p = X[k] - q;
			X[k] = p;
			X[k+n/2] = q / u;
		}
		ifft_recurse(x, X, n/2, s*2);
		ifft_recurse(x+s, X+n/2, n/2, s*2);
	}
}

// Wrapper for recursive function call
void ifft(float* x, float complex* X, unsigned int n) {
	// TODO: check params for validity
	ifft_recurse(x, X, n, 1);
}
