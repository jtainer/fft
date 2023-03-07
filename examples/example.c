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
	const unsigned int n = 1024*1024*1024;
	float complex* signal = malloc(sizeof(float complex) * n);
	fft_param_cache param = fft_create_param_cache(n);

	// Test using cached parameters
	fft_inpl_cached(signal, param);

	// Test without cached paramterers
//	fft_inpl(signal, n);

	free(signal);

	return 0;
}
