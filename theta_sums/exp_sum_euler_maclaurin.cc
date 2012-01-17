#include "theta_sums.h"
#include "precomputed_tables.h"

#include <cmath>
#include <iostream>

using namespace std;


Complex compute_exponential_sum_via_Euler_Maclaurin(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Double epsilon) {
    Complex S = (Complex)0;
 
    theta_cache * cache = build_theta_cache(mp_a, mp_b, j + 2, K);

    //cout << "Warning: Euler Maclaurin case not implemented yet. Evaluating directly." << endl;
    //return direct_exponential_sum_evaluation2(mp_a, mp_b, j, 0, K);

    Double a = cache->a;
    Double b = cache->b;

    //cout << "|a| + 2 * bK = " << abs(a) << " + 2 * " << b << " * " << K << " = " << abs(a) + 2 * b * K << endl;

    //check_condition(a >= -.5 && a <= .5, "In compute_exponential_sum_via_Euler_Maclaurin(), a should be between +-.5, but it isn't");
    //check_condition(2 * b * K < .25, "In compute_exponential_sum_via_Euler_Maclaurin(), b should be less than .25, but it isn't");

    
    //Complex C11 = compute_C11(mp_a, mp_b, K);
    Complex C12 = compute_C12(mp_a, mp_b, K);
    Complex C11 = I * cache->ExpABK;

    S = S + IC0(j, mp_a, mp_b, cache, epsilon/2);

    //cout << "IC0: " << S << endl;

    //Complex z = (Double)1.0/(ExpA(mp_a, K) * ExpB(mp_b, K));
    Complex z = cache->ExpABK;
    if(j == 0) {
        z += 1.0;
    }
    S = S + (Complex).5 * z;
    
    //initialize_power_arrays(21, a, b);


    Complex * p;
    Complex * p_prev = new Complex[j + 1];

    for(int s = 0; s <= j - 1; s++) {
        p_prev[s] = 0;
    }
    p_prev[j] = 1;

    int r = 1;

    //epsilon *= K_power(j, cache);
    Double error = 2 * epsilon + 1;

    Complex S2 = 0.0;


    while(error > epsilon/2) {
    //while(r < 75) {
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
        Double _K_power = K_power(-j, cache);
        Double derivative_estimate = 0.0;
        
        derivative_at_K = p[2 * r - 1 + j];
        derivative_estimate = abs(p[2 * r - 1 + j]);

        for(int k = 2 * r - 1 + j - 1; k >= 0; k--) {
            derivative_at_K = derivative_at_K * (Double)K + p[k];
            derivative_estimate += derivative_estimate * K + abs(p[k]);
        }

        derivative_at_K *= _K_power;
        derivative_estimate *= _K_power;

        //for(int k = 0; k <= 2 * r - 1 + j; k++) {
        //    derivative_at_K = derivative_at_K + _K_power * p[k];
        //    derivative_estimate += _K_power * abs(p[k]);
        //    //cout << k << ": " << _K_power << " * " << p[k] << " = " << _K_power * p[k] << endl;
        //    _K_power *= (Double)K;
        //}


        //derivative_at_K *= (Complex)1.0/(ExpA(mp_a, K) * ExpB(mp_b, K)); //exp(2 * PI * I * (alpha + b));
        //derivative_at_K *= (Complex)1.0/(ExpA(mp_a, K) * ExpB(mp_b, K));
        derivative_at_K *= cache->ExpABK;

     //   cout << (2 * r - 1) << "th derivative at K: " << derivative_at_K << endl;
        Complex derivative_at_0 = p[0] * K_power(-j, cache);
     //   cout << (2 * r - 1) << "th derivative at 0: " << derivative_at_0 << endl;

        //Complex z2 = bernoulli_table[2 * r]/factorial(2 * r) * (z * g_derivative_at_K_without_exponential_factor(2 * r - 1, K) - g_derivative_at_0(2 * r - 1));
        //Complex z2 = bernoulli_table[2 * r]/factorial(2 * r) * (derivative_at_K - derivative_at_0) * pow(K, -j);
        //Complex z2 = bernoulli_over_factorial(2*r) * (derivative_at_K - derivative_at_0) * pow(K, -j);
        Complex z2 = bernoulli_over_factorial(2*r) * (derivative_at_K - derivative_at_0);
        //cout << bernoulli_over_factorial(2 * r) << " * (" << derivative_at_K << " - " << derivative_at_0 << ") = " << z2 << endl;
        //cout << z2 << " = " << bernoulli_over_factorial(2 * r) << " * " << (derivative_at_K - derivative_at_0) << ". deriv estimate: " << pow((j + (Double)r)/K + 2 * PI * (abs(a) + abs(2 * b * K)), r)<< endl;
        //cout << r << ": estimate = " << bernoulli_over_factorial(2 * r) * pow((j + (Double)(2 * r))/K + 2 * PI * (abs(a) + abs(2 * b * K)), (2*r))<< endl;
        //cout << r << ": estimate2 = " << bernoulli_over_factorial(2 * r) * derivative_estimate << endl;
        S2 = S2 + z2;
        error = abs(bernoulli_over_factorial(2 * r)) * derivative_estimate;
        //cout << error/K_power(j, cache) << endl;
        r = r + 1;
    }

    //cout << "Number of corrections terms = " << r << endl;

    S = S + S2;

    if(r > 1)
        delete [] p;
    else
        delete [] p_prev;

    //cout << "Using Euler-Maclaurin summation, computed F(" << a << ", " << b << ", " << K << ") = " << S << endl;
    //cout << "Using direct evaluation, computed         F(" << a << ", " << b << ", " << K << ") = " << direct_exponential_sum_evaluation(mp_a, mp_b, 0, K) << endl;

    free_theta_cache(cache);

    return S;
}



Complex compute_exponential_sums_for_small_b(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon) {
    //
    //
    //
    
    //cout << "--Using Euler-Maclaurin summation." << endl;

    // First we compute the last few terms so that we may assume in the following that K is
    // a multiple of 8

    //cout << "Warning: Euler-Maclauring case not implemented yet. Using direct evaluation." << endl;
    //return compute_exponential_sums_directly(mp_a, mp_b, j, K, v, epsilon);

    if(FAKE_EULER_MACLAURIN)
        return 0.0;

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

//    check_condition( 0 < a + 2 * b * K && a + 2 * b * K < 2 , "Warning: in compute_exponential_sum_for_small_b(), a + 2bK was not between 0 and 2");

    
    {
        Double new_a = a;
        if(new_a > .5) {
            new_a = new_a - 1.0;
        }
        Double z = abs(new_a) + 2 * b * K;
        if(z <= .5 + 1.0/6) {
            mpfr_t mp_new_a;
            mpfr_init2(mp_new_a, mpfr_get_prec(mp_a));
            mpfr_set(mp_new_a, mp_a, GMP_RNDN);
            if(mpfr_cmp_d(mp_new_a, .5) > 0) {
                mpfr_sub_d(mp_new_a, mp_new_a, 1.0, GMP_RNDN);
            }
            Complex S = 0;
            for(int l = 0; l <= j; l++) {
                S = S + v[l] * compute_exponential_sum_via_Euler_Maclaurin(mp_new_a, mp_b, l, K, epsilon/abs(v[l]));
            }
            mpfr_clear(mp_new_a);
            return S;
        }
        //    cout << "a = " << a << "   2bK = " << 2 * b * K << "; |a| + 2bK = " << abs(new_a) + 2 * b * K << "***" << endl;
        //else
        //    cout << "a = " << a << "   2bK = " << 2 * b * K << "; |a| + 2bK = " << abs(new_a) + 2 * b * K << endl;

    }
    
    Complex S = (Complex)0;
    
    int number_of_divisions = ceil(6 * 2 * b * K);
    //number_of_divisions = 8;
    //cout.flush();

    for(int l = 0; l <= j; l++) {
        S = S + v[l] * direct_exponential_sum_evaluation2(mp_a, mp_b, l, K - (K % number_of_divisions), K);
    }
    int K2 = K - (K % number_of_divisions);

    //cout << S << endl;

    Complex S2 = 0;

    //S = S + ((Double)1.0)/(ExpA(mp_a, K) * ExpB(mp_b, K));

    //S2 = S2 + w/(ExpA(mp_a, K2) * ExpB(mp_b, K2));

    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(mp_a));

    for(int m = 0; m < number_of_divisions; m++) {
        Complex dm = (Complex)1.0/(ExpA(mp_a, m * K2/number_of_divisions) * ExpB(mp_b, m * K2/number_of_divisions));
        mpfr_mul_si(tmp, mp_b, 2 * K2, GMP_RNDN); // tmp = 2 * b * K
        mpfr_mul_si(tmp, tmp, m, GMP_RNDN);  // now tmp =  2 * b * K * m
        mpfr_div_si(tmp, tmp, number_of_divisions, GMP_RNDN);  // now tmp = 2 * bKm/number_of_divisions
        mpfr_add(tmp, tmp, mp_a, GMP_RNDN);  // now tmp = a + bKm/number_of_divisions
 
        mpfr_frac(tmp, tmp, GMP_RNDN);       // now tmp = {a + bmK/number_of_divisions}
        if(mpfr_cmp_d(tmp, .5) > 0) {
            mpfr_sub_ui(tmp, tmp, 1., GMP_RNDN);
        }


        Complex Z[max_j + 1];
        for(int l = 0; l <= j; l++) {
            Complex z = 0;
            for(int s = l; s <= j; s++) {
                z = z + v[s] * binomial_coefficient(s, l) * pow(m/(Double)number_of_divisions, s - l) * pow(K2, s - l) * pow(K, -s);
            }
            z *= pow((K2/number_of_divisions - 1), l);
            Z[l] = z;
            //cout << l << ": " << z << endl;
        }


        Complex z = 0;
        for(int l = 0; l <= j; l++) {
            z = z + Z[l] * compute_exponential_sum_via_Euler_Maclaurin(tmp, mp_b, l, K2/number_of_divisions - 1, epsilon/(number_of_divisions * (j+1) * abs(Z[l])));
            //z = z + Z[l] * direct_exponential_sum_evaluation2(tmp, mp_b, l, 0, K2/8 - 1);
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


