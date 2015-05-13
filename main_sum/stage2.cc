/*
 * Functions for computing the "chunk"
 *
 *      sum_{n=v}^{v + N - 1} exp(i*y*log(n)) / sqrt(n)
 * 
 * for y = t, t+delta, t+2*delta, ... , t+(M-1)*delta.
 * These "stage 2" functions use taylor expansions on the summands
 * to achieve significant speedup as compared to individually
 * computing each summand to full precision. This method is only
 * efficient when n is >> t^(1/4). 
 *
 */


#include <queue>
#include <fstream>
#include "theta_sums.h"
#include "main_sum.h"
#include "log.h"

using namespace std;

unsigned int stage_2_block_size(Double v, Double t) {
    //
    // For a given v and t, we determine a good block size to use
    // in stage 2 via the formula
    //
    //      block_size = min{ 5000, min( v / t^0.25 , v / 500^2) }
    //
    // The numbers 5000 and 500 in the formula are somwhat arbitrary and
    // can be adjusted.

    unsigned int block_size = (unsigned int)( min(v * pow(t, -.25), v/(500.0 * 500.0)  ) );
    if(block_size > 5000) return 5000;

    return block_size;
}


Complex zeta_block_stage2_basic(mpz_t v, unsigned int *K, mpfr_t t, Double epsilon) {
    //
    // A function to compute the block
    //
    //      sum_{n = v}^{v + K - 1} exp(i*t*log n) / sqrt(n)
    //
    // where K is calculated by calling stage_2_block_size.
    //

    Complex S = 0;

    // if sum is emptry, return 0
    if(*K == 0) return 0.0;

    // if sum consists of 1 term, return value directly
    if(*K == 1) return exp_itlogn(v)/sqrt(mpz_get_d(v));
    
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    
    // compute an appropriate block size 
    unsigned int block_size = min(stage_2_block_size(vv, tt), *K);

    if(block_size < 2) {
        cout << "Error: in stage 2, computed a block size that is too small " << endl;
        cout << "Refusing to continue, and returning NAN, without setting K." << endl;
        return 0.0/0.0;
    }

    *K = block_size;

    Double x = block_size/vv;

    // The following code to estimate the number of terms that we need in the taylor expansion 
    // is useful (and fast) because it doesn't take any logarithms (instead we just grab the 
    // exponent of the relevant numbers to estimate the log base 2). 
    int logepsilon = fastlog2(epsilon);
    int logtt = fastlog2(tt);
    int logx = fastlog2(x);

    // estimate the number of terms needed in the taylor expansions of 
    // t * log(1 + k / v) and -0.5 * log(1 + k / v)
    int number_of_log_terms = (logepsilon - logtt)/logx;
    int number_of_sqrt_terms = logepsilon/logx;

    // arrays that will hold auxilliary coefficents in the Taylor expansions
    // of t * log(1 + k / v) and -0.5 * log(1 + k / v)
    Double a[number_of_log_terms + 1];
    Double b[number_of_sqrt_terms + 1];

    // We want to accurately compute the quantities t/v mod pi. To do so
    // it should be enough to use log2(t) - log2(v) + 53 bits.
    int vsize = mpz_sizeinbase(v, 2);
    int precision = mpfr_get_exp(t) - vsize + 53;       

    mpfr_t z, mp_v_power, twopi_l;
    mpfr_init2(z, precision);               // We are going to change the
    mpfr_init2(mp_v_power, precision);      // precision on these variables,
    mpfr_init2(twopi_l, precision);         // so we can't use MPFR_DECL_INIT.

    MPFR_DECL_INIT(one_over_v, precision);
    MPFR_DECL_INIT(twopi, precision);
    MPFR_DECL_INIT(z1, 53);


    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_set_z(one_over_v, v, GMP_RNDN);
    mpfr_ui_div(one_over_v, 1, one_over_v, GMP_RNDN);

    int sign = 1;
    Double one_over_vv = 1.0/vv;
    Double v_power = one_over_vv;
    mpfr_set(mp_v_power, one_over_v, GMP_RNDN); // mp_v_power = 1 / v
    mpfr_set(twopi_l, twopi, GMP_RNDN); // two_pi_l = two_pi

    // Code to compute auxilliary coefficients in the taylor expansion
    // of t * log(1 + k / v) mod 2_pi. The mpfr precision is adjusted on 
    // each iteration to improve efficiency.
    //
    for(int l = 1; l <= number_of_log_terms; l++) {
        if(precision >= 53) {
            // if required precision is >= 53 bits, then use mpfr 
            mpfr_set_prec(z, precision);
            mpfr_prec_round(mp_v_power, precision, GMP_RNDN);
            mpfr_prec_round(twopi_l, precision, GMP_RNDN);
            mpfr_mul(z, t, mp_v_power, GMP_RNDN); // z = t / v^l
            
            mpfr_fmod(z1, z, twopi_l, GMP_RNDN); // z1 = t / v^l mod 2_pi_l
            a[l] = sign * mpfr_get_d(z1, GMP_RNDN)/l; // a[l] = (-1)^(l+1) * z1 / l
            mpfr_mul(mp_v_power, mp_v_power, one_over_v, GMP_RNDN); // mp_v_power = mp_v_power / v
            v_power = v_power * one_over_vv; 
            mpfr_add(twopi_l, twopi_l, twopi, GMP_RNDN); // twopi_l = twopi_l + twopi 
            precision = precision - vsize;
        }
        else {
            // if required precision falls below 53 bits, then use double arithmetic directly.
            a[l] = tt * sign * v_power/l;
            v_power = v_power * one_over_vv;
        }

        sign = -sign;
    }

    // code to compute auxilliary coefficients in the taylor expansion of -0.5 * log(1 + k / v)
    // there is no need for mpfr here (because we're not doing and mod-ing by 2_pi)
    //
    Double s = 1;
    for(int l = 1; l <= number_of_sqrt_terms; l++) {
        s = s * (-1.0/vv); // s = (-1)^l / v^l 
        b[l] = .5 * s/l; // b[l] = 0.5 * (-1)^l / (l * v^l)
    }

    for(unsigned int k = 0; k < block_size; k++) {
        Double k_power = 1;
        Double x = 0;
        Double y = 0;

        for(int l = 1; l <= number_of_log_terms; l++) {
            k_power = k_power * k;
            y = y + a[l] * k_power;

            if(l <= number_of_sqrt_terms)
                x = x + b[l] * k_power;
        }

        // x + I*y is the value of (I*t - 0.5) * log(1 + k / v) calculated via 
        // taylor expansions
        S = S + exp(x + I * y);
    }

    S = S / sqrt(vv);
    S = S * exp_itlogn(v);

    mpfr_clear(mp_v_power);
    mpfr_clear(z);
    //mpfr_clear(twopi);
    mpfr_clear(twopi_l);
    //mpfr_clear(one_over_v);
    //mpfr_clear(z1);

    return S;
}

Complex zeta_block_stage2(mpz_t v0, unsigned int N, mpfr_t t, Double delta, int M, Complex * S) {
    //
    // A function to compute the the "chunk"
    //
    //      sum_{n=v0}^{v0 + N - 1} exp(it log n) / sqrt(n)
    //
    // in blocks of length K (this is done by calling zeta_block_stage2_basic repeatedly, which 
    // computes an appropriate K on each call)
    //

    // if the sum is empty return 0 
    if(N == 0) return 0;

    mpz_t v;
    mpz_init(v);

    mpz_set(v, v0);
    unsigned int K = N;

    while(N > 0) {
        Complex current_term = zeta_block_stage2_basic(v, &K, t, exp(-20));
        Complex multiplier = exp(I * delta * log(mpz_get_d(v)));

        for(int l = 0; l < M; l++) {
            S[l] += current_term;
            current_term *= multiplier;
        }
        
        N = N - K;
        mpz_add_ui(v, v, K);
        K = N;
    }

    mpz_clear(v);

    return S[0];
}
