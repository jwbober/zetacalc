#include "theta_sums.h"
#include "precomputed_tables.h"
#include "log.h"

#include <cmath>
#include <iostream>

using namespace std;

complex<double> compute_exponential_sums_using_theta_algorithm(
        mpfr_t mp_a,
        mpfr_t mp_b,
        int J,
        int K,
        complex<double> * v,
        double epsilon,
        int _Kmin) {
    // Compute the linear combination of exponential sums
    //
    // sum_{j=0}^J v[j] * 1/K^j sum_{k=0}^K k^j exp(2 pi i a k + 2 pi i b k^2)
    //
    // using Hiary's "theta sum algorithm" to an "attempted" precision of
    // epsilon.
    //
    // "attempted" precision is defined here as: if all of our code were
    // correct and every arithmetic operation were exact, then the answer
    // returned would be accurate to within an absolute precision of +/-
    // epsilon. The second assumption is, of course, false, and it is not hard
    // to come up with examples (where the value of the sum is large, e.g.)
    // where we don't evaulate the sum that accurately, but we should still
    // get pretty good relative precision in such cases.

    theta_cache * cache = build_theta_cache(mp_a, mp_b, J, K);

    double a = cache->a;
    double b = cache->b;

    int q = cache->q;
    double w = a + 2 * b * K - (double)q;
    if(w < 0) {
        // It is "impossible" for this to happen, because q = floor(2 + abK).
        // But because of unpredictable rounding with the FPU, this can
        // actually happen.
        w = w + 1;
        q = q - 1;
        cache->q = q;
    }

    int p = int(ceil(a));
    double w1 = ceil(a) - a;
    int p1 = q - p;

    complex<double> z[J+1];
    complex<double> z2[J+1];
    double z_epsilon[J+1];
    double v_epsilon[J+1];

    for(int j = 0; j <= J; j++) {
        z[j] = 0.0;
        for(int s = j; s <= J; s++) {
            z[j] += v[s] * binomial_coefficient(s, j);
        }
        z[j] *= I_power(j);
        z_epsilon[j] = (1.0)/(12 * (J + 1.0)) * epsilon/abs(z[j]);
        v_epsilon[j] = (1.0)/(12 * (J + 1.0)) * epsilon/abs(v[j]);
        if(K < z_epsilon[j]) z[j] = 0.0;
        if(K < v_epsilon[j]) v[j] = 0.0;
    }

    complex<double> jbulk1[J+1];
    complex<double> jbulk2[J+1];
    complex<double> jboundary1[J+1];
    complex<double> jboundary2[J+1];

    JBulk(jbulk1, w, b, J, p1, K, cache, z_epsilon);
    JBulk(jbulk2, w1, b, J, p1, K, cache, v_epsilon);
    JBoundary(jboundary1, 2*b*K - w1, 1 - w, b, J, K, cache, z_epsilon);
    JBoundary(jboundary2, 2*b*K - w, 1 - w1, b, J, K, cache, v_epsilon);

    complex<double> S = 0.0;

    complex<double> c = I*cache->ExpABK;
    complex<double> c2 = cache->ExpAK * exp(-2*PI*w*K);
    complex<double> c8 = -I * cache->ExpBK_inverse;

    bool smallw = false;
    if(2.0 * PI * w * K <= -fastlog(epsilon) + (J + 1) * log(2.0)) {
        smallw = true;
    }

    if(smallw) {
        for(int j = 0; j <= J; j++) {
            z2[j] = 0.0;
            for(int s = j; s <= J; s++) {
                complex<double> zz = v[s] * binomial_coefficient(s, j);
                z2[j] += zz * pow(2.0, (s + 1.0)/2.0) * exp_i_pi4(s + 1);
            }
        }
    }
    for(int j = 0; j <= J; j++) {
        if(z[j] != 0.0)
            S += z[j] * c * (   - IC7(K, j, w, b, cache, z_epsilon[j])
                                - jbulk1[j]
                                + minus_one_power(j) * 1.0 * IC9H(K, j, 1 - w, b, cache, z_epsilon[j])
                                - jboundary1[j]    );
        if(v[j] != 0.0)
            S += v[j] * I_power(j + 1) * (   minus_one_power(j + 1) * 1.0 * jbulk2[j]
                            + minus_one_power(j + 1) * 1.0 * IC7(K, j, w1, b, cache, v_epsilon[j])
                            + IC9H(K, j, 1 - w1, b, cache, v_epsilon[j])
                            + minus_one_power(j+1) * 1.0 * jboundary2[j] );
        if(smallw) {
            if(z[j] != 0.0) {
                S += z[j] * c * IC1c(K, j, w, b, c8, cache, z_epsilon[j]);
            }
            if(z2[j] != 0.0) {
                S -= c2 * z2[j] * IC9E(K, j, w, b, cache, 1.0/(12 * (J + 1.0) ) * epsilon/abs(c2*z2[j]));
            }
        }
    }

    c = .5 * cache->ExpABK;
    for(int j = 0; j <= J; j++) {
        S += c * v[j];
    }

    S += .5 * v[0];


    MPFR_DECL_INIT(a1, mpfr_get_prec(mp_a));
    MPFR_DECL_INIT(b1, mpfr_get_prec(mp_b));
    mpfr_ui_div(b1, 1, mp_b, GMP_RNDN);
    mpfr_div_2ui(b1, b1, 1, GMP_RNDN);
    mpfr_mul(a1, mp_a, b1, GMP_RNDN);

    mpfr_div_2ui(b1, b1, 1, GMP_RNDN);
    mpfr_neg(b1, b1, GMP_RNDN);

    complex<double> v2[J+1];
    compute_subsum_coefficients(v2, v, cache);
    if(p == 1) S -= v2[0];
    free_theta_cache(cache);

    complex<double> subsum = compute_exponential_sums(a1, b1, J, q, v2, epsilon, _Kmin);

    S += subsum;

    return S;
}

int normalize(Double &a, Double &b) {
    // Normalize the input a and b so that 0 <= b <= 1/4 and 0 <= a <= 1
    // If we have to use the transformation b -> -b, a -> -a, return 1,
    // otherwise return 0

    // A return of -1 would indicate some strange (impossible) error.

    MPFR_DECL_INIT(mp_a, 100);
    MPFR_DECL_INIT(mp_b, 100);

    mpfr_set_d(mp_a, to_double(a), GMP_RNDN);
    mpfr_set_d(mp_b, to_double(b), GMP_RNDN);

    int return_value = normalize(mp_a, mp_b);

    a = mpfr_get_d(mp_a, GMP_RNDN);
    b = mpfr_get_d(mp_b, GMP_RNDN);
    return return_value;
}

int normalize(mpfr_t a, mpfr_t b) {
    // Normalize the input a and b so that 0 <= b <= 1/4 and 0 <= a <= 1
    // If we have to use the transformation b -> -b, a -> -a, return 1,
    // otherwise return 0

    // A return of -1 would indicate some strange (impossible) error.
    MPFR_DECL_INIT(t, mpfr_get_prec(a));


    //mpfr_t t;
    //mpfr_init2(t, mpfr_get_prec(a));
    mpfr_floor(t, a);
    mpfr_sub(a, a, t, GMP_RNDN);

    mpfr_floor(t, b);
    mpfr_sub(b, b, t, GMP_RNDN);

    //mpfr_clear(t);
    //a = a - floor(a);
    //b = b - floor(b);

    //if(b <= .25) {
    //    return 0;
    //}
    if(mpfr_cmp_d(b, .25) <= 0) {
        return 0;
    }
                                                        //if(b < .75) {
    if(mpfr_cmp_d(b, .75) < 0) {
        mpfr_sub_d(a, a, .5, GMP_RNDN);                         // a -= .5;
        mpfr_sub_d(b, b, .5, GMP_RNDN);                         // b -= .5;

        if(mpfr_cmp_ui(b, 0) >= 0 ) {                           // if(b >= 0)
            if(mpfr_cmp_ui(a, 0) < 0) {                         //      if(a < 0) {
                mpfr_add_ui(a, a, 1, GMP_RNDN);                 //          a = a + 1.0;
            }                                                   //      }
            return 0;                                           //
        }                                                       // }
        else {
            mpfr_neg(a, a, GMP_RNDN);
            mpfr_neg(b, b, GMP_RNDN);
            //mpfr_mul_si(a, a, -1, GMP_RNDN);                    // a = -a
            //mpfr_mul_si(b, b, -1, GMP_RNDN);                    // b = -b
            if(mpfr_cmp_ui(a, 0) < 0) {                         // if(a < 0) {
                mpfr_add_ui(a, a, 1, GMP_RNDN);                 //      a = a + 1
            }                                                   // }
            return 1;
        }
        return 0;
    }
    else {
        mpfr_ui_sub(b, 1, b, GMP_RNDN);                         // b -= 1
        mpfr_neg(a, a, GMP_RNDN);
        //mpfr_mul_si(a, a, -1, GMP_RNDN);                        // a = -a
        //mpfr_mul_si(b, b, -1, GMP_RNDN);                        // b = -b
        if(mpfr_cmp_ui(a, 0) < 0) {                             // if(a < 0) {
            mpfr_add_ui(a, a, 1, GMP_RNDN);                     //      a++
        }                                                       // }
        return 1;
    }

    return -1; // this point should not be reachable.
}


