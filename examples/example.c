// 
// Basic FFT example
//
// 2022, Jonathan Tainer
//

#include "../fft.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
	// buffer size must be a power of 2
	const unsigned int n = 8;
	float* input = malloc(sizeof(float) * n);
	float complex* output = malloc(sizeof(float complex) * n);

	// Store some data in the input buffer
	for (unsigned int i = 0; i < n; i++) {
		input[i] = (float) rand() / RAND_MAX;
	}

	fft(input, output, n);

	for (unsigned int i = 0; i < n; i++) {
		printf("%f\t%f\n", crealf(output[i]), cimagf(output[i]));
	}

	free(input);
	free(output);

	return 0;
}
