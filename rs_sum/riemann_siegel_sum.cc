#include <queue>
#include <fstream>

#include "theta_sums.h"
#include "rs_sum.h"

#include "log.h"

using namespace std;

const int MAX_THREADS = 30;
//const char * NUM_THREADS_FILE = "/home/bober/math/experiments/theta_sums/number_of_threads";
string NUM_THREADS_FILE;
bool use_num_threads_file;
int default_number_of_threads = 4;

unsigned int stage3_start = 1200;

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

Complex partial_zeta_sum(mpz_t start, mpz_t length, mpfr_t t, Double delta, int N, Complex * S, string number_of_threads_filename, int Kmin) {
    //
    // Evaluate the sum \sum_{n=start}^{start + length - 1} exp(i(t + k delta)log n)/sqrt(n)
    // for 0 <= k < N in whatever way seems best possible. The answers will be put into S,
    // which should have space for N entries, and the return value will be the
    // value of the sum for k = 0
    //
    
    if(number_of_threads_filename != "") {
        NUM_THREADS_FILE = number_of_threads_filename;
        use_num_threads_file = true;
    }
    else {
        use_num_threads_file = false;
    }

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

    Complex S1[N];
    Complex S2[N];
    Complex S3[N];

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

    cout << "For t = " << mpfr_get_d(t, GMP_RNDN) << ": " << endl;
    cout << "   Using stage1 up to n = " << n1 << endl;
    cout << "   Using stage2 up to n = " << n2 << endl;
    cout << "   Using stage3 up to n = " << n3 << endl;

    mpz_sub(N1, n1, start);
    mpz_sub(N2, n2, n1);        // N2 and N3 hold the lengths of the sums for stage 2 and stage 3
    mpz_sub(N3, n3, n2);

    create_exp_itlogn_table(t);

    zeta_sum_stage1_version2(start, N1, t, delta, N, S1);
    cout << "Done with stage 1. Sum was: " << S1[0] << endl;

    zeta_sum_stage2(n1, N2, t, delta, N, S2);
    cout << "Done with stage 2. Sum was: " << S2[0] << endl;
    
    //create_exp_itlogn_table(t);

    zeta_sum_stage3(n2, N3, t, delta, N, S3, Kmin);
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

    return S[0];
}



Complex zeta_sum2(mpfr_t t, Double delta, int N, Complex * S) {
    mpz_t length;
    mpz_t one;
    mpz_init(length);
    mpz_init(one);
    create_exp_itlogn_table(t);

    stage_3_bound(length, t);
    mpz_set_ui(one, 1u);

//    Complex Z = partial_zeta_sum(one, length, t, delta, N, S);
    
    Complex Z = 0;

    mpz_clear(length);
    mpz_clear(one);
    return Z;
}


Complex zeta_sum(mpfr_t t, Double delta, int N, Complex * S) {
    for(int l = 0; l < N; l++) {
        S[l] = 0.0;
    }


    create_exp_itlogn_table(t);

    mpz_t n1, n2, n3, N2, N3;

    mpz_init(n1);
    mpz_init(n2);
    mpz_init(n3);
    mpz_init(N2);
    mpz_init(N3);

    stage_1_bound(n1, t);
    stage_2_bound(n2, t);
    stage_3_bound(n3, t);

    Complex S1[N];
    Complex S2[N];
    Complex S3[N];

    if(mpz_cmp(n2, n3) > 0) {      // if n2 > n3, set n2 = n3
        mpz_set(n2, n3);
    }
    if(mpz_cmp(n1, n3) > 0) {      // if n1 > n3, set n1 = n3
        mpz_set(n1, n3);
    }

    cout << "For t = " << mpfr_get_d(t, GMP_RNDN) << ": " << endl;
    cout << "   Using stage1 up to n = " << n1 << endl;
    cout << "   Using stage2 up to n = " << n2 << endl;
    cout << "   Using stage3 up to n = " << n3 << endl;

    mpz_sub(N2, n2, n1);        // N2 and N3 hold the lengths of the sums for stage 2 and stage 3
    mpz_sub(N3, n3, n2);

    zeta_sum_stage1(n1, t, delta, N, S1);
    cout << "Done with stage 1. Sum was: " << S1[0] << endl;

    zeta_sum_stage2(n1, N2, t, delta, N, S2);
    cout << "Done with stage 2. Sum was: " << S2[0] << endl;
    
    //create_exp_itlogn_table(t);

    zeta_sum_stage3(n2, N3, t, delta, N, S3, 800);
    cout << "Done with stage 3. Sum was: " << S3[0] << endl;

    mpz_clear(n1);
    mpz_clear(n2);
    mpz_clear(n3);
    mpz_clear(N2);
    mpz_clear(N3);

    for(int l = 0; l < N; l++) {
        S[l] = S1[l] + S2[l] + S3[l];
    }

    return S[0];
}
