#include <iostream>
#include <cmath>

#include "mpfr.h"
#include "gmp.h"


typedef double Double;

using namespace std;



Double * tlog_table;
int tlog_table_size = 0;

mpfr_t twopi_table[20];

void make_tlog_table(mpfr_t t, int number_of_entries) {
    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));
//    mpfr_root(x, t, 5, GMP_RNDN);
//    tlog_table_size = mpfr_get_si(x, GMP_RNDD);


    tlog_table_size = number_of_entries;
    tlog_table = new Double [tlog_table_size + 1];

    mpfr_t twopi;
    mpfr_init2(twopi, mpfr_get_prec(t));

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);

    for(int n = 1; n <= tlog_table_size; n++) {
        mpfr_set_si(x, n, GMP_RNDN);
        mpfr_log(x, x, GMP_RNDN);
        mpfr_mul(x, x, t, GMP_RNDN);
        mpfr_fmod(x, x, twopi, GMP_RNDN);
        tlog_table[n] = mpfr_get_d(x, GMP_RNDN);
    }

    mpfr_init2(twopi_table[0], mpfr_get_prec(t));
    mpfr_set_ui(twopi_table[0], 0, GMP_RNDN);
    for(int n = 1; n <= 20; n++) {
        mpfr_init2(twopi_table[n], mpfr_get_prec(t));
        mpfr_add(twopi_table[n], twopi_table[n-1], twopi, GMP_RNDN);
    }
    
//    mpfr_out_str(0, 10, 0, twopi_table[0], GMP_RNDN);

    mpfr_clear(x);
    mpfr_clear(twopi);
}

Double tlog(mpfr_t t, mpz_t v) {
    //
    // Return t log v mod 2pi.
    //
    // Uses a precomputed table of values of t log v mod pi, which
    // has to be filled first with a call to make_tlog_table
    //

    mpz_t N;
    mpz_init(N);
    mpz_set(N, v);

    int l = mpz_scan1(N, 0);
    if(l > 0) {
        mpz_fdiv_q_2exp(N, N, l);
    }
    //while(mpz_divisible_2exp_p(N, 1) != 0) {
    //    l++;
    //    mpz_fdiv_q_2exp(N, N, 1);
    //}

    // now v = N * 2^l, with N odd.
    
    if(mpz_fits_sint_p(N) != 0) {
        int z = mpz_get_si(N);
        if(z <= tlog_table_size) {
            mpz_clear(N);
            return l * tlog_table[2] + tlog_table[z];
        }
    }

    mpz_sub_ui(N, N, 1);

    Double x = l * tlog_table[2] + tlog(t, N);

    Double tt = mpfr_get_d(t, GMP_RNDN);
    Double NN = mpz_get_d(N);

    int number_of_terms = 1; // (int)( -37 - log(tt))/(-log(NN));
    int number_of_terms_mpfr = 1; //(int)ceil(-log(tt)/(-log(NN)));

    NN = 1.0/NN;

    Double a[number_of_terms + 1];
    mpfr_t mp_Npower;
    Double Npower;
    mpfr_init2(mp_Npower, mpfr_get_prec(t));

    mpfr_set_si(mp_Npower, 1, GMP_RNDN);
    Npower = 1;

    mpfr_t z, one_over_N;
    mpfr_init2(z, mpfr_get_prec(t));
    mpfr_init2(one_over_N, mpfr_get_prec(t));
    int sign = 1;

    mpfr_set_ui(one_over_N, 1, GMP_RNDN);
    mpfr_div_z(one_over_N, one_over_N, N, GMP_RNDN);

    for(int l = 1; l <= number_of_terms; l++) {
        if(l <= number_of_terms_mpfr) {
            mpfr_mul(mp_Npower, mp_Npower, one_over_N, GMP_RNDN);
            Npower = Npower * NN;
            mpfr_mul(z, t, mp_Npower, GMP_RNDN);
            mpfr_fmod(z, z, twopi_table[l], GMP_RNDN);
            a[l] = sign * mpfr_get_d(z, GMP_RNDN)/l;
        }
        else {
            Npower = Npower * NN;
            a[l] = tt * sign * Npower/l;
        }
        sign = -sign;
    }

    for(int l = 1; l <= number_of_terms; l++) {
        x = x + a[l];
    }

    mpfr_clear(mp_Npower);
    mpfr_clear(z);
    mpz_clear(N);

    return x;
}
