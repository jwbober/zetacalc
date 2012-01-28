#include <vector>
#include <string>

#include "theta_sums.h"
#include "main_sum.h"

using namespace std;

complex<double> zeta_block_mpfr(mpz_t v, unsigned int K, mpfr_t t) {
    //
    // Compute the sum_{k=v}^{v + K - 1) n^{-.5 + it}
    //
    // We use machine doubles for the terms in the sum, but we use
    // mpfr to accurately calculate the quantity t log n mod 2 pi
    // for each term in the sum.

    // This function computes the sum accurately, and it is relatively
    // simple, but it is quite slow. It is written for testing purposes,
    // so that other code can be compared to it.
 
    //
    // First we figure out how much precision we will need from mpfr.
    //
    // We want to accurately calculate t log n mod 2pi to 53 bits, which
    // means that we need to compute t log n to 53 + log_2(t log n) bits.
    //
    // For safety we add an extra 2 bits.
    //

    mpfr_t x, y, n, twopi;
    mpfr_init2(x, 53);
    mpfr_init2(y, 53);

    mpfr_log2(x, t, GMP_RNDN);                          // Now x = log2(t)
    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log(y, y, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // And y = log2(log n) (We calculate these quantities
                                                        // to low precision because we are only interested
                                                        // in their size, really. There is probably 
                                                        // a clever way to do faster using less precision.

    mpfr_add(x, x, y, GMP_RNDN);
    int mod_precision = mpfr_get_ui(x, GMP_RNDN) + 55;  // This is the precision that we need when
                                                        // need to calculate the quantity t log n
                                                        // when we mod by 2 pi

    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // y = log2(v + K) now.
    int n_precision = mpfr_get_ui(y, GMP_RNDN) + 2;     // This is the precision that we need to exactly
                                                        // represent the largest integer that will occur
                                                        // in the summation.
 
    mpfr_init2(n, n_precision);

    mpfr_init2(twopi, mod_precision);
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_clear(x);
    mpfr_init2(x, mod_precision);

    complex<double> S = 0.0;
    
    mpfr_set_z(n, v, GMP_RNDN);             // The summation starts with n = v
    for(unsigned int k = 0; k <= K-1; k++) {
        mpfr_log(x, n, GMP_RNDN);           // x = log(n)
        mpfr_mul(x, x, t, GMP_RNDN);        // a = t log n
        mpfr_fmod(x, x, twopi, GMP_RNDN);
        complex<double> z = exp(I * mpfr_get_d(x, GMP_RNDN));
        z = z/sqrt(mpfr_get_d(n, GMP_RNDN));
        S = S + z;
    
        mpfr_add_ui(n, n, 1, GMP_RNDN);     // now on the next iteration, n will be v + k + 1
    }

    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(n);
    mpfr_clear(twopi);

    return S;
}

complex<double> zeta_block_mpfr(mpfr_t v, unsigned int K, mpfr_t t) {
    //
    // This is a function written for convenience that does exactly
    // the same thing as the above function but accepts an mpfr_t
    // for v.
    //
    mpz_t vv;
    mpz_init(vv);

    mpfr_get_z(vv, v, GMP_RNDN);

    complex<double> S = zeta_block_mpfr(vv, K, t);

    mpz_clear(vv);
    return S;
}

complex<double> zeta_block_mpfr(unsigned int v, unsigned int K, mpfr_t t) {
    //
    // This is a function written for convenience that does exactly
    // the same thing as the above function but accepts an unsigned int
    // for v.

    mpz_t vv;
    mpz_init(vv);

    mpz_set_ui(vv, v);
    complex<double> S = zeta_block_mpfr(vv, K, t);

    mpz_clear(vv);
    return S;
}

double compare_zeta_sums(mpz_t v, mpz_t length, mpfr_t t) {
    complex<double> answer1;
    complex<double> answer2;
    
    answer1 = zeta_block_mpfr(v, mpz_get_ui(length), t);
    partial_zeta_sum(v, length, t, 1.0, 1, &answer2, 800, 1, 1, 0);
    
    return abs(answer1 - answer2);
}

void test() {
    mpfr_t t;
    mpz_t v;
    mpz_t length;

    mpfr_init2(t, 300);
    mpz_init(v);
    mpz_init(length);

    mpfr_set_str(t, "1e30", 10, GMP_RNDN);

    complex<double> answer1;
    complex<double> answer2[100];

    mpz_set_ui(length, 1000000);

    mpz_set_ui(v, 1);
    cout << compare_zeta_sums(v, length, t) << endl;
    mpz_set_str(v, "94568330", 10);
    cout << compare_zeta_sums(v, length, t) << endl;
    mpz_set_str(v, "11999999500000", 10);
    cout << compare_zeta_sums(v, length, t) << endl;
    mpz_set_str(v, "1061566261010081", 10);
    cout << compare_zeta_sums(v, length, t) << endl;

    mpfr_clear(t);
    mpz_clear(v);
    mpz_clear(length);
}

int main(int argc, char ** argv_) {
    vector<string> argv(argv_, argv_ + argc);

    unsigned int seed = time(NULL);
    
    cout << "Seeding rand() and gmp with " << seed << "." << endl;
    srand(seed);
    
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, seed);
    test();

    gmp_randclear(rand_state);
}
