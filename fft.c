// 
// Fast Fourier Transform using Cooley-Tukey algorithm
//
// 2022, Jonathan Tainer
//

#include "fft.h"
#include <math.h>

// Reads signal from input buffer x
// Writes result to output buffer X
// n (num of samples) must be a power of 2
// s (stride length) must be 1
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

