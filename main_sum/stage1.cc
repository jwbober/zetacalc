/*
 * Function for computing the "chunk"
 *
 *      sum_{n=v}^{v + N - 1} exp(i*y*log(n)) / sqrt(n)
 *
 * for y = t, t+delta, t+2*delta, ... , t+(M-1)*delta, using "stage 1", which 
 * computes this sum directly, with each summand computed to (nearly) full double 
 * precision, regardless of the size of n and t.
 * 
 * These functions are relatively slow compared to the stage 2 and stage 3 functions.
 *
 */


#include <queue>
#include <fstream>
#include "theta_sums.h"
#include "main_sum.h"
#include "log.h"

using namespace std;

Complex zeta_block_stage1(mpz_t v, unsigned int N, mpfr_t t, Double delta, int M, Complex * S) {
    //
    // Compute the "chunk"
    //
    //      sum_{n=v}^{v + N - 1} exp(i*y*log(n)) / sqrt(n)
    //
    // for y = t, t+delta, t+2*delta, ... , t+(M-1)*delta, using "stage 1", and place
    // the result into S[0], S[1], ..., S[M-1]
    //
    // We also return the value of S[0], which is the value of the sum at y = t.
    // Currently if M = 1 there is some slight inefficiency, since we care more about 
    // the case where M is larger.
    //
   

    // Special case for dealing with the empty sum:
    if(N == 0) return 0.0;
    
    mpz_t n;
    mpz_init(n);
    mpz_set(n, v);

    for(unsigned int k = 0; k < N; k++) {
        Complex current_term = exp_itlogn(n)/sqrt(mpz_get_d(n));
        Complex multiplier = exp(I * delta * log(mpz_get_d(n)));
        
        for(int l = 0; l < M; l++) {
            S[l] += current_term;
            current_term *= multiplier;
        }

        mpz_add_ui(n, n, 1u); // n = n + 1
    }

    mpz_clear(n);

    return S[0];
}
