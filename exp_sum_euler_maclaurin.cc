#include "theta_sums.h"
#include "precomputed_tables.h"

#include <cmath>
#include <iostream>

using namespace std;


Complex compute_exponential_sum_via_Euler_Maclaurin(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Double epsilon) {
    Complex S = (Complex)0;
 
    theta_cache * cache = build_theta_cache(mp_a, mp_b, j, K);

    //cout << "Warning: Euler Maclaurin case not implemented yet. Evaluating directly." << endl;
    //return direct_exponential_sum_evaluation2(mp_a, mp_b, j, 0, K);

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    check_condition(a >= -.5 && a <= .5, "In compute_exponential_sum_via_Euler_Maclaurin(), a should be between +-.5, but it isn't");
    check_condition(2 * b * K < .25, "In compute_exponential_sum_via_Euler_Maclaurin(), b should be less than .25, but it isn't");

    
    Complex C11 = compute_C11(mp_a, mp_b, K);
    Complex C12 = compute_C12(mp_a, mp_b, K);

    S = S + IC0(K, j, a, b, C11, C12, mp_a, mp_b, cache, epsilon/2);

    //cout << "IC0: " << S << endl;

    Complex z = (Double)1.0/(ExpA(mp_a, K) * ExpB(mp_b, K));
    if(j == 0) {
        z += 1.0;
    }
    S = S + (Complex).5 * z;
    
    //initialize_power_arrays(21, a, b);

    Double error = 2 * epsilon + 1;

    Complex * p;
    Complex * p_prev = new Complex[j + 1];

    for(int s = 0; s <= j - 1; s++) {
        p_prev[s] = 0;
    }
    p_prev[j] = 1;


    int r = 1;
    while(error > epsilon/2) {
        /*
        if(r >= 12) {
            cout << "Warning: not enough derivatives.  Error is" << error << endl;
            break;
        }
        */

    
        if(r > 1) {
            p = new Complex[2 * r - 2 + 1 + j];
            g_derivative_polynomial(2 * r - 2 + j, p, p_prev, a, b);
            delete [] p_prev;
            p_prev = p;
        }
        p = new Complex[2 * r - 1 + 1 + j];
        g_derivative_polynomial(2 * r - 1 + j, p, p_prev, a, b);
        delete [] p_prev;
        p_prev = p;

        Complex derivative_at_K = (Complex)0;
        Double K_power = 1;
        for(int k = 0; k <= 2 * r - 1 + j; k++) {
            derivative_at_K = derivative_at_K + K_power * p[k];
            K_power *= (Double)K;
        }

        derivative_at_K *= (Complex)1.0/(ExpA(mp_a, K) * ExpB(mp_b, K)); //exp(2 * PI * I * (alpha + b));
     //   cout << (2 * r - 1) << "th derivative at K: " << derivative_at_K << endl;
        Complex derivative_at_0 = p[0];
     //   cout << (2 * r - 1) << "th derivative at 0: " << derivative_at_0 << endl;

        //Complex z2 = bernoulli_table[2 * r]/factorial(2 * r) * (z * g_derivative_at_K_without_exponential_factor(2 * r - 1, K) - g_derivative_at_0(2 * r - 1));
        Complex z2 = bernoulli_table[2 * r]/factorial(2 * r) * (derivative_at_K - derivative_at_0) * pow(K, -j);
        S = S + z2;
        error = abs(z2);
        r = r + 1;
    }

    delete [] p;

    //cout << "Using Euler-Maclaurin summation, computed F(" << a << ", " << b << ", " << K << ") = " << S << endl;
    //cout << "Using direct evaluation, computed         F(" << a << ", " << b << ", " << K << ") = " << direct_exponential_sum_evaluation(mp_a, mp_b, 0, K) << endl;

    free_theta_cache(cache);

    return S;
}



Complex compute_exponential_sums_for_small_b(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon) {
    //
    //
    //
    
    // First we compute the last few terms so that we may assume in the following that K is
    // a multiple of 8

    //cout << "Warning: Euler-Maclauring case not implemented yet. Using direct evaluation." << endl;
    //return compute_exponential_sums_directly(mp_a, mp_b, j, K, v, epsilon);


    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

//    check_condition( 0 < a + 2 * b * K && a + 2 * b * K < 2 , "Warning: in compute_exponential_sum_for_small_b(), a + 2bK was not between 0 and 2");

    Complex S = (Complex)0;
    
    for(int l = 0; l <= j; l++) {
        S = S + v[l] * direct_exponential_sum_evaluation2(mp_a, mp_b, l, K - (K % 8), K);
    }
    int K2 = K - (K % 8);

    //cout << S << endl;

    Complex S2 = 0;

    //S = S + ((Double)1.0)/(ExpA(mp_a, K) * ExpB(mp_b, K));

    //S2 = S2 + w/(ExpA(mp_a, K2) * ExpB(mp_b, K2));

    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(mp_a));

    for(int m = 0; m < 8; m++) {
        Complex dm = (Complex)1.0/(ExpA(mp_a, m * K2/8) * ExpB(mp_b, m * K2/8));
        mpfr_mul_si(tmp, mp_b, K2, GMP_RNDN); // tmp = b * K
        mpfr_mul_si(tmp, tmp, m, GMP_RNDN);  // now tmp = b * K * m
        mpfr_div_si(tmp, tmp, 4, GMP_RNDN);  // now tmp = bKm/4
        mpfr_add(tmp, tmp, mp_a, GMP_RNDN);  // now tmp = a + bKm/4
 
        mpfr_frac(tmp, tmp, GMP_RNDN);       // now tmp = {a + bmK/4}
        if(mpfr_cmp_d(tmp, .5) > 0) {
            mpfr_sub_ui(tmp, tmp, 1., GMP_RNDN);
        }


        Complex Z[j + 1];
        for(int l = 0; l <= j; l++) {
            Complex z = 0;
            for(int s = l; s <= j; s++) {
                z = z + v[s] * binomial_coefficient(s, l) * pow(m/8.0, s - l) * pow(K2, s - l) * pow(K, -s);
            }
            z *= pow((K2/8 - 1), l);
            Z[l] = z;
        }


        Complex z = 0;
        for(int l = 0; l <= j; l++) {
            z = z + Z[l] * compute_exponential_sum_via_Euler_Maclaurin(tmp, mp_b, l, K2/8 - 1, epsilon);
        }

        //Complex z = dm * compute_exponential_sum_via_Euler_Maclaurin(tmp, mp_b, K/8 - 1, epsilon);

        z = z * dm;

        //cout << "-----"  << dm << "      "  << z << endl;
        //cout << "-----       " << direct_exponential_sum_evaluation(a, b, m * K/8, (m + 1) * K/8 - 1, epsilon);

        //S = S + dm * compute_exponential_sum_via_Euler_Maclaurin(tmp, mp_b, K/8 - 1, epsilon);
        S2 = S2 + z;
    }

    S = S + S2;

    mpfr_clear(tmp);

    return S;

}


