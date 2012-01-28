#include <queue>
#include <fstream>

#include "theta_sums.h"
#include "main_sum.h"

#include "log.h"

using namespace std;

const int MAX_THREADS = 30;
int default_number_of_threads = 4;

unsigned int stage3_start = 1200;

template<int stage> Complex zeta_sum(mpz_t start, mpz_t length, mpfr_t t, double delta, int M, Complex * S, int number_of_threads, double epsilon, int fraction, int verbose);

void stage_1_bound(mpz_t v, mpfr_t t) {
    //
    // Compute the endpoint for the stage 1 sum for computing zeta(.5 + it)
    // and put it into v.
    //

    // Right now this is 3 * t^{1/4}. With this choice, stage 2 will
    // start up when it can use a block size of about 3.

    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));

    mpfr_root(x, t, 4, GMP_RNDN);
    mpfr_mul_ui(x, x, 3u, GMP_RNDN);
    mpfr_get_z(v, x, GMP_RNDN);

    if(mpz_cmp_ui(v, 1100000u) < 0) {
        mpz_set_ui(v, 1100000u);
    }

    mpfr_clear(x);
}


void stage_2_bound(mpz_t v, mpfr_t t) {
    //
    // Compute the endpoint for the stage 2 sum for computing zeta(.5 + it)
    // and put it into v.
    //

    // Right now this is 1120 t^{1/3}. With this choice, we won't
    // start stage 3 until it can use a block size of approximately 1000.
    // (It is possible that we want to make this even larger right now.)

    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));

    mpfr_cbrt(x, t, GMP_RNDN);
    //mpfr_mul_ui(x, x, 1120u, GMP_RNDN);
    //mpfr_mul_ui(x, x, 51u, GMP_RNDN);
    //mpfr_mul_ui(x, x, 890u, GMP_RNDN);
    mpfr_mul_ui(x, x, stage3_start, GMP_RNDN);
    
    // temporary multiplication to avoid entering stage3
    //mpfr_mul_ui(x, x, 100000u, GMP_RNDN);
    
    mpfr_get_z(v, x, GMP_RNDN);

    mpfr_clear(x);
}

void stage_3_bound(mpz_t v, mpfr_t t) {
    //
    // Compute the endpoint for stage 3.
    //
    // This is just floor( sqrt(t/2pi) ) + 1;
    // The "+ 1" is because we stop BEFORE we reach stage_3_bound.
    
    mpfr_t twopi, x;
    mpfr_init2(twopi, mpfr_get_prec(t));
    mpfr_init2(x, mpfr_get_prec(t));
   
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2u, GMP_RNDN);
    
    mpfr_div(x, t, twopi, GMP_RNDN);
    mpfr_sqrt(x, x, GMP_RNDN);
    mpfr_floor(x, x);

    mpfr_get_z(v, x, GMP_RNDN);
    mpz_add_ui(v, v, 1u);

    mpfr_clear(twopi);
    mpfr_clear(x);
}

Complex partial_zeta_sum(mpz_t start, mpz_t length, mpfr_t t, Double delta, int N, Complex * S, int Kmin, int number_of_threads, int fraction, int verbose) {
    //
    // Evaluate the sum \sum_{n=start}^{start + length - 1} exp(i(t + k delta)log n)/sqrt(n)
    // for 0 <= k < N in whatever way seems best possible. The answers will be put into S,
    // which should have space for N entries, and the return value will be the
    // value of the sum for k = 0
    //

    for(int l = 0; l < N; l++) {
        S[l] = 0.0;
    }

    mpz_t n1, n2, n3, N1, N2, N3;

    mpz_init(n1);
    mpz_init(n2);
    mpz_init(n3);
    mpz_init(N1);
    mpz_init(N2);
    mpz_init(N3);

    stage_1_bound(n1, t);
    stage_2_bound(n2, t);
    mpz_add(n3, start, length);

    Complex * S1 = new Complex[N];
    Complex * S2 = new Complex[N];
    Complex * S3 = new Complex[N];

    if(mpz_cmp(n2, n3) > 0) {      // if n2 > n3, set n2 = n3
        mpz_set(n2, n3);
    }
    if(mpz_cmp(n2, start) < 0) {
        mpz_set(n2, start);
    }
    if(mpz_cmp(n1, n3) > 0) {      // if n1 > n3, set n1 = n3
        mpz_set(n1, n3);
    }
    if(mpz_cmp(n1, start) < 0) {
        mpz_set(n1, start);
    }

    if(verbose) {
        cout << "For t = " << mpfr_get_d(t, GMP_RNDN) << ": " << endl;
        cout << "   Using stage1 up to n = " << n1 << endl;
        cout << "   Using stage2 up to n = " << n2 << endl;
        cout << "   Using stage3 up to n = " << n3 << endl;
    }

    mpz_sub(N1, n1, start);
    mpz_sub(N2, n2, n1);        // N2 and N3 hold the lengths of the sums for stage 2 and stage 3
    mpz_sub(N3, n3, n2);

    create_exp_itlogn_table(t);

    //zeta_sum_stage1(start, N1, t, delta, N, S1, verbose);
    zeta_sum<1>(start, N1, t, delta, N, S1, number_of_threads, exp(-20.0), fraction, verbose);
    if(verbose)
        cout << "Done with stage 1. Sum was: " << S1[0] << endl;

    //zeta_sum_stage2(n1, N2, t, delta, N, S2, verbose);
    zeta_sum<2>(n1, N2, t, delta, N, S2, number_of_threads, exp(-20.0), fraction, verbose);
    if(verbose)
        cout << "Done with stage 2. Sum was: " << S2[0] << endl;
    
    //create_exp_itlogn_table(t);

    //zeta_sum_stage3(n2, N3, t, delta, N, S3, verbose, Kmin);
    zeta_sum<3>(n2, N3, t, delta, N, S3, number_of_threads, exp(-20.0), fraction, verbose);
    if(verbose)
        cout << "Done with stage 3. Sum was: " << S3[0] << endl;


    mpz_clear(n1);
    mpz_clear(n2);
    mpz_clear(n3);
    mpz_clear(N1);
    mpz_clear(N2);
    mpz_clear(N3);

    for(int l = 0; l < N; l++) {
        S[l] = S1[l] + S2[l] + S3[l];
    }

    delete [] S1;
    delete [] S2;
    delete [] S3;

    return S[0];
}


template<int stage> struct sum_data_t {
    int M;
    complex<double> * S;
    pthread_mutex_t * report_mutex;
    pthread_mutex_t * next_mutex;
    
    mpz_t next;
    mpz_t end;
    double length;
    mpfr_t t;
    double tt;
    double delta;
    double epsilon;
    int verbose;
    int percent_finished;
    int fraction;

    sum_data_t(mpz_t start, mpz_t _length, mpfr_t _t, double _delta , int _M, complex<double> * _S, double _epsilon, int _fraction, int _verbose) {
        mpz_init(next);
        mpz_init(end);
        mpfr_init2(t, mpfr_get_prec(_t));

        mpz_set(next, start);
        mpz_add(end, start, _length);
        mpfr_set(t, _t, GMP_RNDN);
        tt = mpfr_get_d(t, GMP_RNDN);
        length = mpz_get_d(_length);
        fraction = _fraction;

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

    ~sum_data_t() {
        pthread_mutex_destroy(report_mutex);
        pthread_mutex_destroy(next_mutex);
        delete report_mutex;
        delete next_mutex;

        mpz_clear(next);
        mpz_clear(end);
        mpfr_clear(t);
    }

    unsigned long next_block(mpz_t start) {
        pthread_mutex_lock(next_mutex);
        unsigned int max_block_size;
        if(stage == 1) {
            max_block_size = 10000;
        }
        else if(stage == 2) {
            max_block_size = 1000000;
        }
        else if(stage == 3) {
            max_block_size = 10000000;
        }
        else {
            cout << "this code should never be reached" << endl;
            exit(-1);
        }
        
        unsigned int block_size;
        mpz_sub(start, end, next);
        if(mpz_cmp_ui(start, max_block_size) < 0) {
            block_size = mpz_get_ui(start);
        }
        else {
            block_size = max_block_size;
        }
        double remainder = mpz_get_d(start);
        int current_percent_finished = 1000 * (1.0 - remainder/length);
        if(percent_finished != current_percent_finished) {
            cout << "stage" << stage << " percent complete: " << current_percent_finished/10.0 << endl;
        }
        percent_finished = current_percent_finished;
        mpz_set(start, next);
        mpz_add_ui(next, next, block_size);

        pthread_mutex_unlock(next_mutex);
        return block_size;
    }

    void report(complex<double> * S2) {
        pthread_mutex_lock(report_mutex);
            

        for(int m = 0; m < M; m++) {
            S[m] += S2[m];
        }

        pthread_mutex_unlock(report_mutex);
    }
};

template<int stage> void * zeta_sum_thread(void * data) {
    sum_data_t<stage> * sum_data = (sum_data_t<stage>*)data;

    mpz_t v;
    mpz_init(v);
    unsigned long length = sum_data->next_block(v);

    unsigned int seed = mpz_fdiv_ui(v, 123456789) + time(NULL);

    complex<double> * S = new complex<double>[sum_data->M];
    
    Complex Z[30];
    if(stage == 3)
        compute_taylor_coefficients(sum_data->t, Z);

    while(length != 0) {
        if(sum_data->fraction > 0) {
            int n = rand_r(&seed);
            if(n % sum_data->fraction == 0) {
                if(stage==1)
                    zeta_block_stage1(v, length, sum_data->t, sum_data->delta, sum_data->M, S);
                if(stage==2)
                    zeta_block_stage2(v, length, sum_data->t, sum_data->delta, sum_data->M, S);
                if(stage==3)
                    zeta_block_stage3(v, length, sum_data->t, Z, sum_data->delta, sum_data->M, S);
            }
        }
        else {
            if(stage==1)
                zeta_block_stage1(v, length, sum_data->t, sum_data->delta, sum_data->M, S);
            if(stage==2)
                zeta_block_stage2(v, length, sum_data->t, sum_data->delta, sum_data->M, S);
            if(stage==3)
                zeta_block_stage3(v, length, sum_data->t, Z, sum_data->delta, sum_data->M, S);
        }

        length = sum_data->next_block(v);
    }
    sum_data->report(S);

    mpz_clear(v);
    pthread_exit(NULL);
}

template<int stage> Complex zeta_sum(mpz_t start, mpz_t length, mpfr_t t, double delta, int M, Complex * S, int number_of_threads, double epsilon, int fraction, int verbose) {
    //
    // Compute the sums S(start, length, t + delta m) for m = 0, 1, ..., M - 1
    // and put the results in S[0], S[1], ..., S[M - 1]
    //
    // Also, return S[0].


    for(int l = 0; l < M; l++) {
        S[l] = 0;
    }

    sum_data_t<stage> sum(start, length, t, delta, M, S, epsilon, fraction, verbose);

    pthread_t threads[number_of_threads];

    for(int n = 0; n < number_of_threads; ++n) {
        pthread_create(&threads[n], NULL, zeta_sum_thread<stage>, (void *)(&sum));
    }
    for(int n = 0; n < number_of_threads; ++n) {
        pthread_join(threads[n], NULL);
    }

    return S[0];

}
