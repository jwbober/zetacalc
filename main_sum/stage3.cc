/*
 * Functions for computing
 *   SUM(v, K, t) = sum_{n=v}^{v + K - 1} n^{-(1/2 + it)}
 *
 * These "stage 3" functions use Ghaith Hiary's "theta algorithm" and
 * are very fast when v is large. They should only be used when
 * v >> t^(1/3) (our current choice is v > 1100 t^(1/3), but the specific
 * choice of 1100 depends on the quality of the theta algorithm
 * implementation and may also be machine dependent, so it is subject
 * to change.)
 *
 */

#include <queue>
#include <fstream>

#include "theta_sums.h"
#include "main_sum.h"
#include "log.h"

using namespace std;

unsigned int stage_3_block_size(Double v, Double t) {
    // With this choice, the block size will hit N when
    // v = N/.9 * t^(1/3)
    //
    // So if we want to start stage 3 with a block size of 100, then we should
    // set stage_2_bound to 112 * t^(1/3)
    //unsigned int block_size = (unsigned int)(  min(.9 * v * pow(t, -.3333333333333333333333), pow(v, .75)/sqrt(500)) );
    unsigned int block_size = (unsigned int)(  min(.9 * v * pow(t, -.3333333333333333333333), v/500.0 * 500.0) );
    //unsigned int block_size = (unsigned int)(2 * v * pow(t, -.3333333333333333333333));
    return block_size;
}

Complex zeta_block_stage3_basic(mpz_t v, unsigned int *K, mpfr_t t, Complex ZZ[30], Double epsilon, int Kmin) {
    if(*K == 0) {
        return 0.0;
    }
    if(*K == 1) {
        if(0) {
            cout << v << " " << v << endl;
        }
        return exp_itlogn(v)/sqrt(mpz_get_d(v));
    }
    
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    
    unsigned int block_size = min(stage_3_block_size(vv, tt), *K);
    if(block_size < 50 && *K >= 50) {
        cout << "Error: in stage 3, computed a block size of " << block_size << ", which is smaller than 50. Shouldn't reach stage two until we can use a larger block size." << endl;
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




    Double w = (block_size-1)/vv;
    Double w_power = 1;

    Complex Z[30];

    for(int l = 0; l < 19; l++) {
        Z[l] = ZZ[l] * w_power;
        w_power *= w;
    }

    int j = 18;

    // Compute Z[l]
 
    mpfr_t a, b, x;

    int precision = mpfr_get_prec(t);

    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(x, precision);

    mpfr_const_pi(x, GMP_RNDN);             // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);         // x = 2 pi
    mpfr_mul_z(x, x, v, GMP_RNDN);            // x = 2 pi v
    mpfr_div(a, t, x, GMP_RNDN);            // a = t / (2 pi v)

//    mpfr_mul_si(x, x, -2, GMP_RNDN);        // x = -4 pi v
//    mpfr_mul(x, x, v, GMP_RNDN);            // x = -4 pi v^2
//    mpfr_div(b, t, x, GMP_RNDN);            // b = -t/ (4 pi v^2)

//    mpfr_mul_si(b, v, -2, GMP_RNDN);
    mpfr_div_z(b, a, v, GMP_RNDN);
    mpfr_div_si(b, b, -2, GMP_RNDN);

//    cout << mpfr_get_d(b, GMP_RNDN) << endl;

    Complex S = compute_exponential_sums(a, b, j, block_size-1, Z, epsilon, Kmin, 0);

    // we don't need the previous values of a and b anymore, so
    // we can erase them.

    Complex z = exp_itlogn(v);
    z = z / sqrt(mpz_get_d(v));
    S = S * z;

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(x);


//    cout << *K << endl;

    return S;
}

Complex zeta_block_stage3(mpz_t n, unsigned int N, mpfr_t t, Complex Z[30], Double delta, int M, Complex * S, int Kmin) {
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
        Complex current_term = zeta_block_stage3_basic(v, &K, t, Z, exp(-20), Kmin);
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

Complex zeta_sum_stage3(mpz_t n, mpz_t N, mpfr_t t, Double delta, int M, Complex * S, int Kmin, int verbose) {
    //
    // Compute and return the sum
    //
    // \sum_{k=n}^{n + N - 1} n^{.5 + it}
    //
    // We do this by repeatedly calling zeta_block_stage3 with a block size of 1000000
    // and add up all the terms. We could just call it directly, but we anticipate
    // shortly changing this routine to use multiple threads.
    //

    time_t start_wall_time = time(NULL);
    clock_t last_cpu_time = clock();
    double total_cpu_time = 0;

    for(int l = 0; l < M; l++) {
        S[l] = 0.0;
    }


    Complex * S2 = new Complex[M];
    for(int l = 0; l < M; l++) {
        S[l] = 0.0;
        S2[l] = 0.0;
    }


    const unsigned int block_size = 100000000;

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, N, block_size);

    Complex Z[30];
    compute_taylor_coefficients(t, Z);

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_set(v, n);
    mpz_sub_ui(v, v, block_size);

    for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
        mpz_add_ui(v, v, block_size);
        zeta_block_stage3(v, block_size, t, Z, delta, M, S2, Kmin);
        for(int l = 0; l < M; l++) {
            S[l] += S2[l];
        }
        //S = S + zeta_block_mpfr(v, block_size, t);
        if(mpz_divisible_ui_p(k, 20u)) {
            time_t current_wall_time = time(NULL);
            clock_t current_cpu_time = clock();
            time_t elapsed_wall_time = current_wall_time - start_wall_time;
            double elapsed_cpu_time = ((double)current_cpu_time - (double)last_cpu_time)/CLOCKS_PER_SEC;
            total_cpu_time += elapsed_cpu_time;
            if(verbose) {
                cout << "In stage3, completed " << k << " large blocks out of " << number_of_blocks << "." << endl;
                cout << "        In stage3 thus far: " << elapsed_wall_time << " real seconds; " << total_cpu_time << " cpu seconds; " << elapsed_cpu_time << "cpu seconds this chunk. " << endl;
            }
            last_cpu_time = current_cpu_time;
            int current_blocksize = stage_3_block_size(mpz_get_d(v), mpfr_get_d(t, GMP_RNDN));
            if(verbose)
                cout << "        Current blocksize ~= " << current_blocksize << endl;
        }
    }
    mpz_add_ui(v, v, block_size);

    zeta_block_stage3(v, remainder, t, Z, delta, M, S2, Kmin);
    for(int l = 0; l < M; l++) {
        S[l] += S2[l];
    }
    //S = S + zeta_block_mpfr(v, remainder, t);

    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);

    time_t end_wall_time = time(NULL);
    time_t elapsed_wall_time = end_wall_time - start_wall_time;
    if(verbose)
        cout << "Spent " << elapsed_wall_time << " seconds in stage 3." << endl;

    delete [] S2;

    return S[0];

}


void compute_taylor_coefficients(mpfr_t t, Complex Z[30]) {
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
