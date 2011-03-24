/*
 * Functions for computing
 *   SUM(v, K, t) = sum_{n=v}^{v + K - 1} n^{-(1/2 + it)}
 *
 * These "stage 1" functions compute this sum directly, with each summand
 * computed to (nearly) full double precision, regardless of the size of n and t.
 * 
 * These functions are relatively slow compared to the stage 2 and stage 3 functions,
 * but for small n there isn't much of an alternative.
 *
 */


#include <queue>
#include <fstream>

#include "theta_sums.h"
#include "rs_sum.h"

#include "log.h"

using namespace std;

const int MAX_THREADS = 30;
extern string NUM_THREADS_FILE;

const int stage1_block_size = 10000;

Complex zeta_sum_stage1_version2(mpz_t start, mpz_t length, mpfr_t t, Double delta, int M, Complex * S) {
    //
    // Compute the sums S(start, length, t + delta m) for m = 0, 1, ..., M - 1
    // and put the results in S[0], S[1], ..., S[M - 1]
    //
    // Also, return S[0].
    //
    // Before calling this function, the initialization for the exp_itlogn() function
    // has to be done by calling create_exp_itlogn_table(t).
    //
    // We repeatedly call zeta_block_stage1 with a block size of stage1_block_size
    // and add up all of the terms. We could just call zeta_block_stage1
    // directly, but we anticipate shortly changing this routine
    // to use multiple threads.
    

    time_t start_time = time(NULL);

    for(int l = 0; l < M; l++) {
        S[l] = 0;
    }

    Complex S2[M];

    const unsigned int block_size = 10000;

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, length, block_size);

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_sub_ui(v, start, block_size);
    for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
        mpz_add_ui(v, v, block_size);
        zeta_block_stage1(v, block_size, t, delta, M, S2);
        for(int l = 0; l < M; l++) {
            S[l] += S2[l];
        }
    }
    mpz_add_ui(v, v, block_size);

    zeta_block_stage1(v, remainder, t, delta, M, S2);
    for(int l = 0; l < M; l++) {
        S[l] += S2[l];
    }

    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);

    time_t end_time = time(NULL);

    time_t elapsed_time = end_time - start_time;
    cout << "Spent " << elapsed_time << " seconds in stage 1." << endl;

    return S[0];
}


Complex zeta_block_stage1(mpz_t v, unsigned int K, mpfr_t t, Double delta, int M, Complex * S) {
    //
    // compute SUM(v, K, t + m delta) for m = 0, 1, 2, ... M - 1 and place
    // the result into S[0], S[1], ..., S[M-1]
    //
    // We also return the value of S[0]. Currently if M = 1 there is some slight
    // inefficiency, since we care more about the case where M is larger.
    //
    // S should already be initialized, but we will explicitly set it to 0 here.

    
    // Start by zeroing the return array.
    for(int l = 0; l < M; l++) {
        S[l] = 0;
    }
    
    // Special case for dealing with the empty sum:
    if(K == 0) {
        return 0.0;
    }
    
    mpz_t n;
    mpz_init(n);
    mpz_set(n, v);

    for(unsigned int k = 0; k < K; k++) {
        Complex current_term = exp_itlogn(n)/sqrt(mpz_get_d(n));
        Complex multiplier = exp(I * delta * log(mpz_get_d(n)));
        
        for(int l = 0; l < M; l++) {
            S[l] += current_term;
            current_term *= multiplier;
        }

        mpz_add_ui(n, n, 1u);
    }

    mpz_clear(n);

    return S[0];
}

Complex zeta_sum_stage1(mpz_t N, mpfr_t t, Double delta, int M, Complex * S) {
    // THIS FUNCTION SHOULD NO LONGER BE USED, AND SHOULD BE REMOVED.
    //
    // Compute and return the sum
    //
    // \sum_{n=1}^{N - 1} n^{-.5 + it)
    //
    // NOTE: Unlike the other zeta_sum_* functions, N is NOT the
    // number of terms that will be computed. This will compute N - 1
    // terms. (Or, it can be thought of as computing N terms,
    // but with the 0th term being 0).
    //
    // We repeatedly call zeta_block_stage1 with a block size of 10000
    // and add up all of the terms. We could just call zeta_block_stage1
    // directly, but we anticipate shortly changing this routine
    // to use multiple threads.

    time_t start_time = time(NULL);


    for(int l = 0; l < M; l++) {
        S[l] = 0;
    }

    Complex S2[M];

    const unsigned int block_size = 10000;

    mpz_t number_of_terms;
    mpz_init(number_of_terms);
    mpz_sub_ui(number_of_terms, N, 1u);

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, number_of_terms, block_size);

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_set_si(v, 1 - (int)block_size);
    for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
        mpz_add_ui(v, v, block_size);
        zeta_block_stage1(v, block_size, t, delta, M, S2);
        for(int l = 0; l < M; l++) {
            S[l] += S2[l];
        }
        //S = S + zeta_block_mpfr(v, block_size, t);
    }
    mpz_add_ui(v, v, block_size);

    zeta_block_stage1(v, remainder, t, delta, M, S2);
    for(int l = 0; l < M; l++) {
        S[l] += S2[l];
    }

    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);
    mpz_clear(number_of_terms);

    time_t end_time = time(NULL);

    time_t elapsed_time = end_time - start_time;
    cout << "Spent " << elapsed_time << " seconds in stage 1." << endl;

    return S[0];
}
