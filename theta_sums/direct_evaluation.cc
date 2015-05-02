#include "theta_sums.h"

#include <cmath>
#include <iostream>

using namespace std;

Complex compute_exponential_sums_directly(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int method) {
    //
    // method == 1: automatically choose whether or not to use mpfr
    // method == 4: use mpfr
    // method == 5: don't use mpfr
    //
    // (these numbers are nonconsecutive so that they can be directly passed
    // from the function compute_exponential_sums();
    //
    Complex S = 0;

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    if ( (method == 1 && K > mpfr_Kmin) || method == 4) {
        for(int l = 0; l <= j; l++) {
            Complex x = 0.0;
            if(v[l] != 0.0) {
                x = v[l] * direct_exponential_sum_evaluation2(mp_a, mp_b, l, 0, K);
            }
            S = S + x;
        }
    }
    else {
        Double one_over_K = 1. / (Double)K;
        for(int n = 0; n <= K; n++) {
            Complex x = 0;
            Double n_over_K_powers = 1.;
            //Complex common_exp = EXP( two_pi_i * (Double)n * (a + b * (Double)n) );
            Double z = 2 * PI * n * ( a + b * n);
            Complex common_exp = Complex( cos(z), sin(z) );
            Complex S2 = 0.0;
            for(int l = 0; l <= j; l++) {
                Complex y = 0.0;
                if(v[l] != 0.) {
                    y = v[l] * n_over_K_powers;
                }
                S2 = S2 + y;

                n_over_K_powers *= (Double)n * one_over_K;
            }
            x = S2 * common_exp;
            S += x;
        }
    }
 
    return S;

}


Complex direct_exponential_sum_evaluation2(Double a, Double b, int j, int m, int M, int working_precision) {
    // return the sum 
    //
    // Sum_{n=0}^N exp(2 pi i alpha t + 2 pi i beta t^2)
    //
    // computed by direct evaluation
    
    Complex S = (Complex)0.0;

    /*
    if(working_precision > 53) {
        mpfr_t a, b;
        mpfr_init2(a, working_precision);
        mpfr_init2(b, working_precision);
        mpfr_set_d(a, to_double(alpha), GMP_RNDN);
        mpfr_set_d(b, to_double(beta), GMP_RNDN);
        S = direct_exponential_sum_evaluation(a, b, m, M);
        mpfr_clear(a);
        mpfr_clear(b);
        return S;
    }
    */

    for(int n = m; n <= M; n++) {
        S = S + pow(n, j) * EXP( (Complex)2.0 * PI * I * (Double)n * (a + b * (Double)n) );
    }

    S = S/pow(M, j);

    return S;

}

Complex direct_exponential_sum_evaluation2(mpfr_t a, mpfr_t b, int j, int m, int M) {
    mpfr_t real_part, imaginary_part;
    mpfr_t t;
    mpfr_t t2;
    mpfr_t t3;
    mpfr_init2(real_part, mpfr_get_prec(a));
    mpfr_init2(imaginary_part, mpfr_get_prec(a));
    mpfr_set_ui(real_part, 0, GMP_RNDN);
    mpfr_set_ui(imaginary_part, 0, GMP_RNDN);

    mpfr_init2(t, mpfr_get_prec(a));
    mpfr_init2(t2, mpfr_get_prec(a));
    mpfr_init2(t3, mpfr_get_prec(a));

    for(int k = m; k <= M; k++) {
        mpfr_mul_si(t, a, k, GMP_RNDN);         // t = ak
        mpfr_mul_si(t2, b, k, GMP_RNDN);        // t2 = bk
        mpfr_mul_si(t2, t2, k, GMP_RNDN);       // now t2 = bk^2
        mpfr_add(t2, t, t2, GMP_RNDN);          // now t2 = ak + bk^2
        mpfr_const_pi(t, GMP_RNDN);             // t = pi
        mpfr_mul_2ui(t, t, 1, GMP_RNDN);        // now t = 2pi
        mpfr_mul(t, t, t2, GMP_RNDN);           // now t = 2pi(ak + bk^2)
        mpfr_sin_cos(t, t2, t, GMP_RNDN);       // t = sin(2 pi(ak + bk^2)), t2 = cos(2 pi (ak + bk^2))

        mpfr_set_si(t3, k, GMP_RNDN);
        mpfr_pow_ui(t3, t3, j, GMP_RNDN);
        mpfr_mul(t, t, t3, GMP_RNDN);
        mpfr_mul(t2, t2, t3, GMP_RNDN);

        mpfr_add(real_part, real_part, t2, GMP_RNDN);
        mpfr_add(imaginary_part, imaginary_part, t, GMP_RNDN);
    }

    mpfr_set_si(t3, M, GMP_RNDN);
    mpfr_pow_ui(t3, t3, j, GMP_RNDN);
    mpfr_div(real_part, real_part, t3, GMP_RNDN);
    mpfr_div(imaginary_part, imaginary_part, t3, GMP_RNDN);

    Complex S(mpfr_get_d(real_part, GMP_RNDN), mpfr_get_d(imaginary_part, GMP_RNDN));
    mpfr_clear(real_part);
    mpfr_clear(imaginary_part);
    mpfr_clear(t);
    mpfr_clear(t2);
    mpfr_clear(t3);
    return S;
}
