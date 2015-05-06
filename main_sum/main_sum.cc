/*
 * 
 * A description of main_sum.cc: 
 *
 * 1. zetacalc is first asked to compute a partial zeta sum. This is done by feeding zetacalc an input file, 
 * or from the commandline. 
 *
 * 2. The function partial_zeta_sum is then called. It decides which stage (i.e. method) to use for 
 * which piece of the partial sum; here:
 * 
 *   stage 1 is direct evaluation using exp_itlogn and mpfr
 *   stage 2 is direct evaluation using taylor expansions and mpfr
 *   stage 3 uses the theta algorithm
 *      
 * 3. partial_zeta_sum_stage (templated function) handles the pieces corresponding to stages 1, 2, and 3
 * (notice, almost always, only one piece is not empty). This is done by
 *
 *      a. calling zeta_sum_thread to create threads to compute "chunks" 
 *      b. creating sum_data_t structs (templated), which keep track of the next available "chunk", 
 *         and handle thread locking/unlocking & "reporting".
 *
 * 4. In turn, zeta_sum_thread calls zeta_block_stage1, or zeta_block_stage2, or zeta_block_stage3. These functions 
 * grab the next available "chunk" in the partial sum (by calling sum_data_t->next_block), they call
 * the function compute_taylor_coefficients, then compute the "chunk" in smaller blocks using the appropriate 
 * stage (i.e. method).
 * 
 */

#include <iostream>
#include <random>
#include "theta_sums.h"
#include "main_sum.h"
#include "log.h"

using namespace std;

typedef std::mt19937 RNG;

    // This function is declared here and defined at the end of the file
template<int stage> Complex partial_zeta_sum_stage(mpz_t start, mpz_t length, mpfr_t t, double delta, int M, Complex * S, int number_of_threads, double epsilon, int fraction, int verbose);


void stage_2_start(mpz_t v, mpfr_t t) {
    //
    // Computes the starting point in the main sum for stage 2 and 
    // stores it in v. The starting point is computed according to 
    // the formula:
    //
    //      v =  max{ floor( 3 * t^{1/4} ) , 1100000 } 
    //
    // With this choice, the block size in stage 2 will be at least 3.
    // The appearance of 1100000 in the formula is to ensure that the
    // Taylor expansions used in stage 2 don't get too long (which they
    // can if v is too small).
    //

    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));

    mpfr_root(x, t, 4, GMP_RNDN); // x = t^(1/4)
    mpfr_mul_ui(x, x, 3u, GMP_RNDN); // x = 3 * t^(1/4)
    mpfr_get_z(v, x, GMP_RNDN); // v = floor(3 * t^(1/4)

    // if v < 1100000, then set v = 1100000
    if(mpz_cmp_ui(v, 1100000u) < 0) mpz_set_ui(v, 1100000u);

    mpfr_clear(x);
}


void stage_3_start(mpz_t v, mpfr_t t) {
    //
    // A function to compute the starting point in the main sum for the stage 3 and 
    // store it in v. The starting point is computed according to the formula:
    //
    //      v = 1200 * t^(1/3)
    //
    // With this choice, the block size in stage 3 will be at least 1200.
    //
    
    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));

    mpfr_cbrt(x, t, GMP_RNDN); // x = t^(1/3)
    mpfr_mul_ui(x, x, 1200u, GMP_RNDN);  // x = stage_3_start * t^(1/3)
    mpfr_get_z(v, x, GMP_RNDN); // v = floor( stage_3_start * t^(1/3) )

    mpfr_clear(x);
}


Complex partial_zeta_sum(mpz_t start, mpz_t length, mpfr_t t, Double delta, int M, Complex * S, int Kmin, int number_of_threads, int fraction, int verbose) {
    //
    // This function computes the partial sum:
    //
    //          \sum_{n=start}^{start + length - 1} exp(i*y*log n)/sqrt(n)
    //
    // for y = t,t+delta,t+2*delta,...,t+(M-1)*delta, then stores the resulting 
    // values in S. It also returns the value of the partial sum at y = t
    //

    mpz_t n2, n3, end, N1, N2, N3;
    mpz_init(n2);
    mpz_init(n3);
    mpz_init(end);
    mpz_init(N1);
    mpz_init(N2);
    mpz_init(N3);

    stage_2_start(n2, t); // n2 is the starting point of stage_2 
    stage_3_start(n3, t); // n3 is the starting point of stage_3

    mpz_add(end, start, length); // end = start+length, is the endpoint of the partial sum 
                                 // (so the last term in the sum is n = end - 1)

    // S1, S2, and S3 will store the contributions of stages 1, 2, and 3 (if any) to the 
    // partial zeta sum.
    Complex * S1 = new Complex[M]; 
    Complex * S2 = new Complex[M]; 
    Complex * S3 = new Complex[M]; 

    // stage 1 contribution to the partial sum comes from start <= n < n2,
    // stage 2 from n2 <= n < n3, and stage 3 from n3 <= n < end. The "if 
    // statements" below reset n2 & n3 appropriately to mark the start & 
    // end points of stages 2 & 3.
    //
    if(mpz_cmp(end,n3) < 0) 
        mpz_set(n3, end); // if end < n3, set n3 = end
    if(mpz_cmp(start,n3) > 0) 
        mpz_set(n3, start); // if start > n3, set n3 = start
    if(mpz_cmp(end,n2) < 0) 
        mpz_set(n2, end); // if end < n2, set n2 = end
    if(mpz_cmp(start,n2) > 0) 
        mpz_set(n2, start); // if start > n2, set n2 = start

    if(verbose) {
        cout << "For t = " << mpfr_get_d(t, GMP_RNDN) << ": " << endl;
        //cout << "   Using stage1 up to n = " << n1 << endl;
        cout << "       Starting sum at n = " << start << endl;
        cout << "    Starting stage2 at n = " << n2 << endl;
        cout << "    Starting stage3 at n = " << n3 << endl;
        cout << "         Ending sum at n = " << end << endl;
    }

    mpz_sub(N1, n2, start); // N1 = n2 - start, is the length of stage 1 sum
    mpz_sub(N2, n3, n2); // N2 = n3 - n2, is the length of stage 2 sum
    mpz_sub(N3, end, n3); // N3 = end - n3, is the length of stage 3 sum

    // we carry out a precomputation, which is needed in conjunction with a 
    // modification of Fenynman's algorithm for computing exp(i*t*log(n)); see 
    // log.cc
    create_exp_itlogn_table(t);

    // compute the stages 1, 2, & 3 sums
    partial_zeta_sum_stage<1>(start, N1, t, delta, M, S1, number_of_threads, exp(-20.0), fraction, verbose);
    if(verbose) cout << "Done with stage 1. Sum was: " << S1[0] << endl;

    partial_zeta_sum_stage<2>(n2, N2, t, delta, M, S2, number_of_threads, exp(-20.0), fraction, verbose);
    if(verbose) cout << "Done with stage 2. Sum was: " << S2[0] << endl;
    
    partial_zeta_sum_stage<3>(n3, N3, t, delta, M, S3, number_of_threads, exp(-20.0), fraction, verbose);
    if(verbose) cout << "Done with stage 3. Sum was: " << S3[0] << endl;

    // Store the contributions of stages 1, 2 & 3 in S, which is the array passed to partial_zeta_sum.
    // We first initialize S to zero, then add up the contributions of S1, S2, & S3
    for(int l = 0; l < M; l++) {
        S[l] = 0.0;
    }
    for(int l = 0; l < M; l++) {
        S[l] = S1[l] + S2[l] + S3[l];
    }

    mpz_clear(n2);
    mpz_clear(n3);
    mpz_clear(end);
    mpz_clear(N1);
    mpz_clear(N2);
    mpz_clear(N3);

    delete [] S1;
    delete [] S2;
    delete [] S3;

    return S[0]; // return the value of the partial sum at t
}


template<int stage> struct sum_data_t {
    // 
    // A struct to help keep track of threads. It contains the function next_block, 
    // which increments the starting point of the next available "chunk", and the function 
    // report (to extract the result of the "chunk" computation)
    //
    
    int M; // # of gird points where partial sum will be computed
    complex<double> * S; // storage space for the values of the partial sum at t, t + delta, ..., t + (M-1)*delta
    pthread_mutex_t * report_mutex; // to handle thread locking/unlocking in preparation for a block computation
    pthread_mutex_t * next_mutex; // to manage threads 
    
    mpz_t next; // the starting point of the next available block
    mpz_t end; // the end point of the partial sum (which consists of several blocks)
    double length; // the length (= # of terms) of the next available block
    mpfr_t t; // the height t where the block is computed
    double tt; // a double version of t for double-precision arithmetic
    double delta; // to compute the block at t, t + delta, t + 2*delta, ... , t + (M-1) * delta
    double epsilon; // a precision parameter that is passed to the theta algorithm 
    
    // to help with debugging 
    int verbose; 
    int percent_finished;
    int fraction;

    RNG rng;
    std::uniform_int_distribution<unsigned long> randint;

    // a struct constructor: sets variable values and initializes a thread
    sum_data_t(mpz_t start, mpz_t _length, mpfr_t _t, double _delta , int _M, complex<double> * _S, double _epsilon, int _fraction, int _verbose, unsigned int seed) {
        mpz_init(next);
        mpz_init(end);
        mpfr_init2(t, mpfr_get_prec(_t));

        mpz_set(next, start); // next = start
        mpz_add(end, start, _length); // end = start + _length
        mpfr_set(t, _t, GMP_RNDN);  
        tt = mpfr_get_d(t, GMP_RNDN); 
        length = mpz_get_d(_length);
        fraction = _fraction;
        rng = RNG(seed);
        randint = std::uniform_int_distribution<unsigned long>(0, 2 * fraction);

        delta = _delta;
        M = _M;
        S = _S;
        epsilon = _epsilon;
        verbose = _verbose;
        report_mutex = new pthread_mutex_t;
        next_mutex = new pthread_mutex_t;
        pthread_mutex_init(report_mutex, NULL);
        pthread_mutex_init(next_mutex, NULL);
        percent_finished = -1;
    }

    // a struct destructor
    ~sum_data_t() {
        pthread_mutex_destroy(report_mutex);
        pthread_mutex_destroy(next_mutex);
        delete report_mutex;
        delete next_mutex;

        mpz_clear(next);
        mpz_clear(end);
        mpfr_clear(t);
    }

    // locks thread, computes and returns block_size and stores the starting point of next 
    // available "chunk", then unlocks thread 
    unsigned long next_block(mpz_t start) {
        pthread_mutex_lock(next_mutex);
        unsigned long max_block_size;
        
        if(stage == 1) { max_block_size = 10000;}
        else if(stage == 2) { max_block_size = 1000000;}
        else if(stage == 3) { max_block_size = 10000000;}
        else {
            cout << "this segment of the code should never be reached, exiting" << endl;
            exit(-1);
        }
        
        unsigned int block_size;
        mpz_sub(start, end, next); // start = end - next (notice this is
                                   // an abuse of notation since now start
                                   // denotes the length of the remaining
                                   // partial sum)

        if(mpz_cmp_ui(start, max_block_size) < 0) { // if start < max_block_size,
            block_size = mpz_get_ui(start);         // set block_size = start
        }
        else {
            block_size = max_block_size;
        }

        double remainder = mpz_get_d(start);
        int current_percent_finished = 1000 * (1.0 - remainder/length);
        if(percent_finished != current_percent_finished && verbose) {
            cout << "stage" << stage << " percent complete: "
                 << current_percent_finished/10.0 << endl;
        }
        percent_finished = current_percent_finished;
        
        mpz_set(start, next);               // start = next (start is now the
                                            // beginning point of the remainder
                                            // of the partial sum) 
        
        mpz_add_ui(next, next, block_size); // next = next + block_size
        
        if(fraction > 1) {
            mpz_add_ui(next, next, block_size*randint(rng));
            if(mpz_cmp(next, end) > 0) {
                mpz_set(next, end);
            }
        }

        pthread_mutex_unlock(next_mutex);
        return block_size;
    }

    // once a thread is done, it "reports" its output to be added to the array S 
    // that was created with the struct 
    void report(complex<double> * S2) {
        pthread_mutex_lock(report_mutex);
        
        for(int m = 0; m < M; m++) {
            S[m] += S2[m];
        }

        pthread_mutex_unlock(report_mutex);
    }
};


template<int stage> void * zeta_sum_thread(void * data) {
    //
    // creates a sum_data_t struct, computes auxilliary coefficients in the linear 
    // combination of quadratic sums that results from taylor-expanding the terms in 
    // the main sum, and starts a computation of the partial sum (in "chunks")
    //
    
    sum_data_t<stage> * sum_data = (sum_data_t<stage>*)data;

    mpz_t v;
    mpz_init(v);

    // returns size of currently available block, stores its starting point in v, and 
    // increments "next" in sum_data to "next + length"
    unsigned long length = sum_data->next_block(v);

    // array where thread output will be stored
    complex<double> * S = new complex<double>[sum_data->M];

    // initializes S[0], S[1], ... , S[M-1] to zero
    for(int l = 0; l < sum_data->M; l++) S[l] = 0;

    // array where the auxilliary coefficients in the linear combination of quadratic 
    // sums is stored (however, it doesn't seem this needs to be called each time since
    // it only depends on t ...)
    Complex Z[30];
    if(stage == 3)
        compute_taylor_coefficients(sum_data->t, Z);

    // compute the partial sum in "chunks" 
    while(length != 0) {
            if(stage==1) 
                zeta_block_stage1(v, length, sum_data->t, sum_data->delta, sum_data->M, S);
            if(stage==2) 
                zeta_block_stage2(v, length, sum_data->t, sum_data->delta, sum_data->M, S);
            if(stage==3) 
                zeta_block_stage3(v, length, sum_data->t, Z, sum_data->delta, sum_data->M, S);
        length = sum_data->next_block(v);
    }

    // add the computation output stored in S to the internal array in sum_data struct 
    // (which is the same as the array originally passed to partial_zeta_sum_stage)
    sum_data->report(S);

    mpz_clear(v);
    pthread_exit(NULL);
}

template<int stage> Complex partial_zeta_sum_stage(mpz_t start, mpz_t length, mpfr_t t, double delta, int M, Complex * S, int number_of_threads, double epsilon, int fraction, int verbose) {
    // 
    // Computes the partial sum:
    //
    //      \sum_{n = start}^{start + length - 1} exp(i*t*log(n)) / sqrt(n)
    //
    // at y = t, t+delta, t+2*delta,...,t+(M-1)*delta, and stores the result 
    // in S[0], S[1], ... , S[M-1]. The stage number determines which method 
    // to use for computing the partial main sum.
    //
    // The function attempts to compute the sum to within epsilon
    // using the given number of threads. It also returns S[0], which
    // is the value of the partial sum at t.
    //

    for(int l = 0; l < M; l++) S[l] = 0;

    sum_data_t<stage> sum(start, length, t, delta, M, S, epsilon, fraction, verbose, 0);
    pthread_t threads[number_of_threads];

    for(int n = 0; n < number_of_threads; ++n) 
        pthread_create(&threads[n], NULL, zeta_sum_thread<stage>, (void *)(&sum));

    for(int n = 0; n < number_of_threads; ++n) 
        pthread_join(threads[n], NULL);

    return S[0];
}
