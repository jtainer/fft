// 
// Fast Fourier Transform using Cooley-Tukey algorithm
//
// 2022, Jonathan Tainer
//

#include "fft.h"
#include <math.h>

// Input time domain signal in buf
// Convert to freq domain and store result in buf
// n (num of samples) must be a power of 2
// s (stride length) must be 1 for initial function call
static void fft_recurse(float complex* buf, unsigned int n, unsigned int s) {
	if (n == 1) return;

	// Recursively split the DFT in halves
	fft_recurse(buf, n/2, s*2);
	fft_recurse(buf+n/2, n/2, 2*s);

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

	fft_recurse(sig_fd, n, 1);
}

// Input freq domain signal in buf
// Convert to time domain and store result in buf
// n (number of samples) must be a power of 2
// s (stride length) must be 1 for initial function call
static void ifft_recurse(float complex* buf, unsigned int n, unsigned int s) {
	if (n == 1) return;

	for (unsigned int k = 0; k < n/2; k++) {
		float complex u = cexpf(I*-2.f*M_PI*k/n);
		float complex q = (buf[k] - buf[k+n/2]) / 2.f;
		float complex p = buf[k] - q;
		buf[k] = p;
		buf[k+n/2] = q / u;
	}
	ifft_recurse(buf, n/2, s*2);
	ifft_recurse(buf+n/2, n/2, s*2);
}

// sig_td = buffer to store output time domain signal
// sig_fd = input freq domain signal (modified due to in-place calculation)
// n = sample count (MUST be a power of 2)
void ifft(float* sig_td, float complex* sig_fd, unsigned int n) {
	// Calculate FFT in-place in the freq domain signal buffer
	ifft_recurse(sig_fd, n, 1);

	// Copy real component of freq domain buffer to time domain buffer
	for (unsigned int i = 0; i < n; i++) {
		sig_td[i] = crealf(sig_fd[i]);
	}
}

// Complex time domain signal in sig is overwritten
// with freq domain signal
// n = sample count
void fft_inpl(float complex* sig, unsigned int n) {
	fft_recurse(sig, n, 1);
}

// Complex freq domain signal in sig is overwritten
// with time domain signal
// n = sample count
void ifft_inpl(float complex* sig, unsigned int n) {
	ifft_recurse(sig, n, 1);
}
