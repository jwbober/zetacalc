/*
 * Functions for computing the "chunk"
 *
 *          sum_{n = v}^{v + N - 1} exp(i*t*log(n)) / sqrt(n)
 * 
 * for y = t, t+delta, t+2*delta, ... , t+(M-1)*delta.
 * These "stage 3" functions use the "theta algorithm", which 
 * is very fast when v is v >> t^(1/3).
 *
 */

#include <queue>
#include <fstream>
#include "theta_sums.h"
#include "main_sum.h"
#include "log.h"

using namespace std;

unsigned int stage_3_block_size(Double v, Double t) {
    //
    // The block size is calculated according to the formula:
    //
    //      K = min( .9 * v / t^(1/3), v / 500^2 )
    //
    // The factor 0.9 is to better control the number of theta-sum 
    // derivatives that are needed to compute the block
    //
    
    unsigned int block_size = (unsigned int)(  min(.9 * v * pow(t, -.3333333333333333333333), v/500.0 * 500.0) );
    return block_size;
}


Complex zeta_block_stage3_basic(mpz_t v, unsigned int *K, mpfr_t t, Complex ZZ[30], Double epsilon, int Kmin) {
    //
    // function to compute the block
    //
    //      sum_{n = v}^{v + K - 1} exp(i*t*log(n)) / sqrt(n)
    //
    // using the theta algorithm. Notice that K is re-calculated inside,
    // so its value may differ from the one passed.
    //


    // if sum is empty, return 0
    if(*K == 0) return 0.0;
   
    // if sum consists of 1 term only, then return the value of the term directly
    if(*K == 1) return exp_itlogn(v)/sqrt(mpz_get_d(v));
    
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    
    // calculate an appropriate block size
    unsigned int block_size = min(stage_3_block_size(vv, tt), *K);

    if(block_size < 50 && *K >= 50) {
        cout << "Error: in stage 3, computed a block size of " << block_size;
        cout << ", which is smaller than 50. Shouldn't reach stage two until we can use a larger block size." << endl;
        cout << "Refusing to continue, and returning NAN, without setting K." << endl;
        return 0.0/0.0;
    }

    *K = block_size;
    Double w = (block_size-1)/vv;
    Double w_power = 1;

    Complex Z[30];

    // notice that j = 18 is hard-coded here; eventually, j should be passed instead
    int j = 18;

    for(int l = 0; l <= j; l++) {
        Z[l] = ZZ[l] * w_power;
        w_power *= w;
    }

    // Compute the linea and quadratic arguments a & b in the theta sum
    int precision = mpfr_get_prec(t);

    MPFR_DECL_INIT(a, precision);
    MPFR_DECL_INIT(b, precision);
    MPFR_DECL_INIT(x, precision);

    mpfr_const_pi(x, GMP_RNDN); // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN); // x = 2 pi
    mpfr_mul_z(x, x, v, GMP_RNDN);  // x = 2 pi v
    mpfr_div(a, t, x, GMP_RNDN); // a = t / (2 pi v)

    mpfr_div_z(b, a, v, GMP_RNDN); // b = a / v
    mpfr_div_si(b, b, -2, GMP_RNDN); // b = - a / (2*v)

    // call the theta algorithm to compute the block
    // Notice that we pass "block_size - 1" to the theta algorithm since it takes
    // the end-point of the sum, rather than its length 
    Complex S = compute_exponential_sums(a, b, j, block_size-1, Z, epsilon, Kmin, 0);

    // on multiplying S by the leading term, we obtain the value of the block
    Complex z = exp_itlogn(v);
    z = z / sqrt(mpz_get_d(v));
    S = S * z;

    return S;
}


Complex zeta_block_stage3(mpz_t v0, unsigned int N, mpfr_t t, Complex Z[30], Double delta, int M, Complex * S, int Kmin) {
    //
    // function to compute the "chunk"
    //      
    //      sum_{n = v}^{v + N - 1} exp(i*t*log(n)) / sqrt(n)
    //
    // in blocks of length K, where each block is computed using the
    // theta algorithm
    //

    // if sum is empty, return 0
    if(N == 0) return 0;

    mpz_t v;
    mpz_init(v);

    mpz_set(v, v0);
    unsigned int K = N;

    while(N > 0) {
        // NOTICE: epsilon = exp(-20) has been hard-coded here; it should be passed eventually.
        Complex current_term = zeta_block_stage3_basic(v, &K, t, Z, exp(-20), Kmin);
        Complex multiplier = exp(I * delta * log(mpz_get_d(v)));

        for(int l = 0; l < M; l++) {
            S[l] += current_term;
            current_term *= multiplier;
        }

        N = N - K;
        mpz_add_ui(v, v, K); // v = v + K
        K = N;
    }

    mpz_clear(v);

    return S[0];
}


void compute_taylor_coefficients(mpfr_t t, Complex Z[30]) {
    //
    // Auxilliary coeffients for the taylor expansion of each
    // block into a lnear combination of theta sums. Notice only
    // the first two terms are used from the taylor expansion of
    // 1 / sqrt(1 + k / v)
    //

    Double tt1 = mpfr_get_d(t, GMP_RNDN);
    Double tt2 = tt1 * tt1;
    Double tt3 = tt2 * tt1;
    Double tt4 = tt3 * tt1;
    Double tt5 = tt4 * tt1;
    Double tt6 = tt5 * tt1;
    Double tt7 = tt6 * tt1;

    Z[0] = 1.0;
    Z[1] = -.5;
    Z[2] =  .375;
    Z[3] = I * tt1 / 3.0;
    Z[4] = -I * tt1 * 5.0/12.0;
    Z[5] = I * tt1 * 9.0/20.0;
    Z[6] = - tt2 / 18.0;
    Z[7] =  tt2 / 9.0;
    Z[8] = - tt2 * 77.0/480.0;
    Z[9] = -I * tt3 / 162.0;
    Z[10] = I  * tt3 * 11.0/648.0;
    Z[11] = -I * tt3 * 133.0/4320.0;
    Z[12] = tt4 / 1944.0;
    Z[13] = 0.0;
    Z[14] = 0.0;
    Z[15] = tt5 / 29160.0;
    Z[16] = 0.0;
    Z[17] = 0.0;
    Z[18] = tt6 / 524880.0;
    Z[19] = 0;
    Z[20] = 0;
    Z[21] = tt7 / 11022480.0;
    Z[22] = 0;

}
