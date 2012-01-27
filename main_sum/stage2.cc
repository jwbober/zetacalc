/*
 * Functions for computing
 *   SUM(v, K, t) = sum_{n=v}^{v + K - 1} n^{-(1/2 + it)}
 *
 * These "stage 2" functions use taylor expansions on the summands
 * to achieve significant speedup as compared to individually
 * computing each summand to full precision. This method is only
 * effective once n is >> t^(1/4). (Our current choice of parameters
 * cause stage 2 to be used once n > 3 t^(1/4), so that we can
 * compute the sum at least 3 terms at a time.)
 *
 * For more information on the implementation, see comments below.
 * 
 */



#include <queue>
#include <fstream>

#include "theta_sums.h"
#include "main_sum.h"

#include "log.h"

using namespace std;

inline unsigned int stage_2_block_size(Double v, Double t) {
    //
    // For a given v and t determine a good block size to use
    // in stage 2.
    //
    
    // For now we set the block size to be v/t^{1/4}

    //unsigned int block_size = min((unsigned int)( v * pow(t, -.25) ), (unsigned int)(pow(t, 1.0/12.0)));
    //unsigned int block_size = (unsigned int)( v * pow(t, -.25) );
    //unsigned int block_size = (unsigned int)( min(v * pow(t, -.25), sqrt(v)/500.0  ) );
    //unsigned int block_size = (unsigned int)( min(v * pow(t, -.25), pow(v, .75)/sqrt(300000.0)  ) );
    unsigned int block_size = (unsigned int)( min(v * pow(t, -.25), v/(500.0 * 500.0)  ) );

    if(block_size > 5000)
        return 5000;

    return block_size;
}

Complex zeta_block_stage2_basic(mpz_t v, unsigned int *K, mpfr_t t, Double epsilon) {
    //
    // This routine calculates the sum
    //
    // sum_{n=v}^{v + K - 1} n^{.5 + it) = sum_{n=v}^{v + K} exp(it log n)/sqrt(n)
    //
    // to a nominal precision of epsilon*K/sqrt(v)
    //
    // Unless *K = 0 or 1, we always compute a good value of *K, and store it
    // in the passed variable. The *K that is passed is the maximum block
    // size that we will use, however.
    //
    // To deal with precision issues in the calculation of the exponential
    // and to get significant speedup, we do a change of variables to
    // write this sum as 
    //
    // exp(it log v)/sqrt(v) sum_{k=0}^{K-1} exp(it log(1 + k/v))/sqrt(1 + k/v).
    //
    // Then we write 1/sqrt(v + k) as exp(log -.5(1 + k/v)),
    // so the sum is
    //
    // exp(i t log(1 + k/v) + .5 log(1 + k/v))
    //
    // We Taylor expand the logarithms with a few terms, and instead
    // of computing each term using mpfr, we just calculate the first
    // terms mod 2pi using mpfr, and then multiply by powers of k
    // in the inner sum using double arithmetic.
    //
    // Let x = K/v.
    // The number of terms we compute in the taylor expansion of the log is
    //
    // -(log(t) - log(epsilon))/log(x)
    //
    // We compute the initial terms using mpfr up to -log(t)/log(x).
    //
    // The number of terms of the in the taylor expansion for the square root
    // term is log(epsilon)/log(x).

    Complex S = 0;

    if(*K == 0) {
        return 0.0;
    }
    if(*K == 1) {
        if(0)
            cout << v << "  " << v << endl;
        return exp_itlogn(v)/sqrt(mpz_get_d(v));
    }  
    
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    
    unsigned int block_size = min(stage_2_block_size(vv, tt), *K);
    if(block_size < 2) {
        cout << "Error: in stage 2, computed a block size of " << block_size << ", which is smaller than 2. Shouldn't reach stage two until we can use a larger block size." << endl;
        cout << "Refusing to continue, and returning NAN, without setting K." << endl;
        return 0.0/0.0;
    }

    *K = block_size;

    if(0) {
        mpz_t n;
        mpz_init(n);
        mpz_add_ui(n, v, block_size - 1);
        cout << v << " " << n << endl;
        mpz_clear(n);
    }

    Double x = block_size/vv;

    // It seems like the following code to estimate the number
    // of terms that we need in the taylor expansion might be
    // useful (and fast) because it doesn't take any logarithms.
    // Instead we just grab the exponent of the relevant numbers
    // to estimate the log base 2. In practice it seems like
    // this slows things down a little bit, though, perhaps
    // because we sometimes (often?) take one extra term
    // in the taylor expansion, so we don't use it right now.

    // UPDATE: After writing the above comment, I changed the
    // way that we decide whether or not to use mpfr to calculate
    // the terms, and now this this method seems about as good,
    // and maybe better, so I am using it again.
    //

    int logepsilon = fastlog2(epsilon);
    int logtt = fastlog2(tt);
    int logx = fastlog2(x);

    int number_of_log_terms = (logepsilon - logtt)/logx;
    int number_of_sqrt_terms = logepsilon/logx;

    Double a[number_of_log_terms + 1];
    Double b[number_of_sqrt_terms + 1];

    int vsize = mpz_sizeinbase(v, 2);
    int precision = mpfr_get_exp(t) - vsize + 53;       
                                                // We want to accurately compute
                                                // the quantities t/v mod pi. To do so
                                                // it should be enough
                                                // to use log2(t) - log2(v) + 53 bits.
    

    mpfr_t mp_v_power, z, twopi, one_over_v, twopi_l, z1;
    
    mpfr_init2(mp_v_power, precision);
    mpfr_init2(one_over_v, precision);
    mpfr_init2(z, precision);
    mpfr_init2(twopi, precision);
    mpfr_init2(twopi_l, precision);
    mpfr_init2(z1, 53);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_set_z(one_over_v, v, GMP_RNDN);
    mpfr_ui_div(one_over_v, 1, one_over_v, GMP_RNDN);

    int sign = 1;
    Double one_over_vv = 1.0/vv;
    Double v_power = one_over_vv;
    mpfr_set(mp_v_power, one_over_v, GMP_RNDN);
    mpfr_set(twopi_l, twopi, GMP_RNDN);
    for(int l = 1; l <= number_of_log_terms; l++) {
        //if(l <= number_of_log_terms_mpfr) {
        if(precision >= 53) {
            //
            // The following calls to reduce the precision on each iteration
            // seem to slow things down just a very little bit in the tests I have run, and maybe
            // they shouldn't be there. But I suspect that in cases
            // where v is very large there should be some gain
            // in using them, so I am leaving them for now.
            //
            mpfr_set_prec(z, precision);
            mpfr_prec_round(mp_v_power, precision, GMP_RNDN);
            mpfr_prec_round(twopi_l, precision, GMP_RNDN);
            mpfr_mul(z, t, mp_v_power, GMP_RNDN);
            
            //mpfr_div_si(z, z, l, GMP_RNDN);
            //mpfr_frac(z, z, GMP_RNDN);
            mpfr_fmod(z1, z, twopi_l, GMP_RNDN);
            a[l] = sign * mpfr_get_d(z1, GMP_RNDN)/l;
            mpfr_mul(mp_v_power, mp_v_power, one_over_v, GMP_RNDN);
            v_power = v_power * one_over_vv;
            mpfr_add(twopi_l, twopi_l, twopi, GMP_RNDN);
            precision = precision - vsize;
        }
        else {
            a[l] = tt * sign * v_power/l;
            v_power = v_power * one_over_vv;
        }
        sign = -sign;
    }


    Double s = 1;
    for(int l = 1; l <= number_of_sqrt_terms; l++) {
        s = s * (-1.0/vv);       // RECENTLY FIXED. VERIFY.
        b[l] = .5 * s/l;
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
        S = S + exp(x + I * y);
    }
    S = S / sqrt(vv);

    S = S * exp_itlogn(v);

    mpfr_clear(mp_v_power);
    mpfr_clear(z);
    mpfr_clear(twopi);
    mpfr_clear(twopi_l);
    mpfr_clear(one_over_v);
    mpfr_clear(z1);

    return S;
}

Complex zeta_block_stage2(mpz_t n, unsigned int N, mpfr_t t, Double delta, int M, Complex * S) {
    //for(int l = 0; l < M; l++) {
    //    S[l] = 0.0;
    //}
    
    if(N == 0) {
        return 0;
    }


    mpz_t v;
    mpz_init(v);

    mpz_set(v, n);
    unsigned int K = N;
    while(N > 0) {
        Complex current_term = zeta_block_stage2_basic(v, &K, t, exp(-20));
        Complex multiplier = exp(I * delta * log(mpz_get_d(v)));
        for(int l = 0; l < M; l++) {
            S[l] += current_term;
            current_term *= multiplier;
        }
        //cout << v << "  " << K << endl;
        N = N - K;
        mpz_add_ui(v, v, K);
        K = N;
    }

    mpz_clear(v);

    return S[0];
}

Complex zeta_sum_stage2(mpz_t n, mpz_t N, mpfr_t t, Double delta, int M, Complex * S, int verbose) {
    //
    // Compute and return the sum
    //
    // \sum_{k=n}^{n + N - 1} n^{.5 + it}
    //
    // We do this by repeatedly calling zeta_block_stage2 with a block size of 1000000
    // and add up all the terms. We could just call it directly, but we anticipate
    // shortly changing this routine to use multiple threads.
    //

    time_t start_wall_time = time(NULL);
    clock_t last_cpu_time = clock();

    for(int l = 0; l < M; l++) {
        S[l] = 0.0;
    }


    const unsigned int block_size = 10000000;

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, N, block_size);

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_set(v, n);
    mpz_sub_ui(v, v, block_size);

    Complex * S2 = new Complex[M];
    for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
        mpz_add_ui(v, v, block_size);
        zeta_block_stage2(v, block_size, t, delta, M, S2);
        for(int l = 0; l < M; l++) {
            S[l] += S2[l];
        }
        //S = S + zeta_block_mpfr(v, block_size, t);
        if(mpz_divisible_ui_p(k, 20u)) {
            time_t current_wall_time = time(NULL);
            clock_t current_cpu_time = clock();
            time_t elapsed_wall_time = current_wall_time = start_wall_time;
            double elapsed_cpu_time = ((double)current_cpu_time - (double)last_cpu_time)/CLOCKS_PER_SEC;
            if(verbose)
                cout << "In stage2, completed " << k << " large blocks out of " << number_of_blocks << ". Spent " << elapsed_wall_time << " seconds so far and " << elapsed_cpu_time << " cpu seconds this block. " << endl;
            last_cpu_time = current_cpu_time;
        }
    }
    mpz_add_ui(v, v, block_size);

    zeta_block_stage2(v, remainder, t, delta, M, S2);
    for(int l = 0; l < M; l++) {
        S[l] += S2[l];
    }
    delete [] S2;

    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);

    time_t end_wall_time = time(NULL);
    time_t elapsed_wall_time = end_wall_time - start_wall_time;
    if(verbose)
        cout << "Spent " << elapsed_wall_time << " seconds in stage 2." << endl;

    return S[0];

}
