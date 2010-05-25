#include "theta_sums.h"
#include "precomputed_tables.h"

#include <iostream>
#include <cmath>



using namespace std;

Complex compute_exponential_sum_using_theta_algorithm(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon);

Complex direct_exponential_sum_evaluation(Double alpha, Double beta, int m, int M, int working_precision) {
    // return the sum 
    //
    // Sum_{n=0}^N exp(2 pi i alpha t + 2 pi i beta t^2)
    //
    // computed by direct evaluation
    
    Complex S = (Complex)0.0;

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

    for(int n = m; n <= M; n++) {
        S = S + EXP( (Complex)2.0 * PI * I * (Double)n * (alpha + beta * (Double)n) );
    }

    return S;

}

Complex direct_exponential_sum_evaluation(mpfr_t a, mpfr_t b, int m, int M) {
    mpfr_t real_part, imaginary_part;
    mpfr_t t;
    mpfr_t t2;
    mpfr_init2(real_part, mpfr_get_prec(a));
    mpfr_init2(imaginary_part, mpfr_get_prec(a));
    mpfr_set_ui(real_part, 0, GMP_RNDN);
    mpfr_set_ui(imaginary_part, 0, GMP_RNDN);

    mpfr_init2(t, mpfr_get_prec(a));
    mpfr_init2(t2, mpfr_get_prec(a));

    for(int k = m; k <= M; k++) {
        mpfr_mul_si(t, a, k, GMP_RNDN);         // t = ak
        mpfr_mul_si(t2, b, k, GMP_RNDN);        // t2 = bk
        mpfr_mul_si(t2, t2, k, GMP_RNDN);       // now t2 = bk^2
        mpfr_add(t2, t, t2, GMP_RNDN);          // now t2 = ak + bk^2
        mpfr_const_pi(t, GMP_RNDN);             // t = pi
        mpfr_mul_2ui(t, t, 1, GMP_RNDN);        // now t = 2pi
        mpfr_mul(t, t, t2, GMP_RNDN);           // now t = 2pi(ak + bk^2)
        mpfr_sin_cos(t, t2, t, GMP_RNDN);       // t = sin(2 pi(ak + bk^2)), t2 = cos(2 pi (ak + bk^2))
        mpfr_add(real_part, real_part, t2, GMP_RNDN);
        mpfr_add(imaginary_part, imaginary_part, t, GMP_RNDN);
    }
    Complex S(mpfr_get_d(real_part, GMP_RNDN), mpfr_get_d(imaginary_part, GMP_RNDN));
    mpfr_clear(real_part);
    mpfr_clear(imaginary_part);
    mpfr_clear(t);
    mpfr_clear(t2);
    return S;
}

Complex compute_exponential_sum_via_Euler_Maclaurin(Double a, Double b, int K, Double epsilon) {
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);
    mpfr_set_d(mp_a, to_double(a), GMP_RNDN);
    mpfr_set_d(mp_b, to_double(b), GMP_RNDN);
    
    Complex S = compute_exponential_sum_via_Euler_Maclaurin(mp_a, mp_b, K, epsilon);

    mpfr_clear(mp_a);
    mpfr_clear(mp_b);

    return S;
}

Complex compute_exponential_sum_via_Euler_Maclaurin(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon) {
    Complex S = (Complex)0;
    

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    check_condition(a >= -.5 && a <= .5, "In compute_exponential_sum_via_Euler_Maclaurin(), a should be between +-.5, but it isn't");
    check_condition(2 * b * K < .25, "In compute_exponential_sum_via_Euler_Maclaurin(), b should be less than .25, but it isn't");

    
    Complex C11 = compute_C11(mp_a, mp_b, K);
    Complex C12 = compute_C12(mp_a, mp_b, K);

    S = S + IC0(K, a, b, C11, C12, mp_a, mp_b, epsilon/2);

    //cout << "IC0: " << S << endl;

    Complex z = (Double)1.0/(ExpA(mp_a, K) * ExpB(mp_b, K));
    
    S = S + (Complex).5 * (z + (Double)1);
    
    //initialize_power_arrays(21, a, b);

    Double error = 2 * epsilon + 1;

    Complex * p;
    Complex * p_prev = new Complex[0];

    int r = 1;
    while(error > epsilon/2) {
        /*
        if(r >= 12) {
            cout << "Warning: not enough derivatives.  Error is" << error << endl;
            break;
        }
        */

        p = new Complex[2 * r - 2 + 1];
        g_derivative_polynomial(2 * r - 2, p, p_prev, a, b);
        delete [] p_prev;
        p_prev = p;
        p = new Complex[2 * r - 1 + 1];
        g_derivative_polynomial(2 * r - 1, p, p_prev, a, b);
        delete [] p_prev;
        p_prev = p;

        Complex derivative_at_K = (Complex)0;
        Double K_power = 1;
        for(int k = 0; k <= 2 * r - 1; k++) {
            derivative_at_K = derivative_at_K + K_power * p[k];
            K_power *= (Double)K;
        }

        derivative_at_K *= (Complex)1.0/(ExpA(mp_a, K) * ExpB(mp_b, K)); //exp(2 * PI * I * (alpha + b));
     //   cout << (2 * r - 1) << "th derivative at K: " << derivative_at_K << endl;
        Complex derivative_at_0 = p[0];
     //   cout << (2 * r - 1) << "th derivative at 0: " << derivative_at_0 << endl;

        //Complex z2 = bernoulli_table[2 * r]/factorial(2 * r) * (z * g_derivative_at_K_without_exponential_factor(2 * r - 1, K) - g_derivative_at_0(2 * r - 1));
        Complex z2 = bernoulli_table[2 * r]/factorial(2 * r) * (derivative_at_K - derivative_at_0);

        S = S + z2;
        error = abs(z2);
        r = r + 1;
    }

    delete [] p;

    //cout << "Using Euler-Maclaurin summation, computed F(" << a << ", " << b << ", " << K << ") = " << S << endl;
    //cout << "Using direct evaluation, computed         F(" << a << ", " << b << ", " << K << ") = " << direct_exponential_sum_evaluation(mp_a, mp_b, 0, K) << endl;


    return S;
}

Complex compute_exponential_sum_via_Euler_Maclaurin(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Double epsilon) {
    Complex S = (Complex)0;
    

    //cout << "Warning: Euler Maclaurin case not implemented yet. Evaluating directly." << endl;
    //return direct_exponential_sum_evaluation2(mp_a, mp_b, j, 0, K);

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    check_condition(a >= -.5 && a <= .5, "In compute_exponential_sum_via_Euler_Maclaurin(), a should be between +-.5, but it isn't");
    check_condition(2 * b * K < .25, "In compute_exponential_sum_via_Euler_Maclaurin(), b should be less than .25, but it isn't");

    
    Complex C11 = compute_C11(mp_a, mp_b, K);
    Complex C12 = compute_C12(mp_a, mp_b, K);

    S = S + IC0(K, j, a, b, C11, C12, mp_a, mp_b, epsilon/2);

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


    return S;
}




Complex compute_exponential_sum_for_small_b(Double a, Double b, int K, Double epsilon) {
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);
    mpfr_set_d(mp_a, to_double(a), GMP_RNDN);
    mpfr_set_d(mp_b, to_double(b), GMP_RNDN);
    
    Complex S = compute_exponential_sum_for_small_b(mp_a, mp_b, K, epsilon);

    mpfr_clear(mp_a);
    mpfr_clear(mp_b);

    return S;
}

Complex compute_exponential_sum_for_small_b(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon) {
    //
    //
    //
    
    // First we compute the last few terms so that we may assume in the following that K is
    // a multiple of 8


    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

//    check_condition( 0 <= a + 2 * b * K && a + 2 * b * K < 2 , "Warning: in compute_exponential_sum_for_small_b(), a + 2bK was not between 0 and 2");

    Complex S = (Complex)0;
    
    S = S + direct_exponential_sum_evaluation(mp_a, mp_b, K - (K % 8) + 1, K);
    K = K - (K % 8);

    //cout << S << endl;

    S = S + ((Double)1.0)/(ExpA(mp_a, K) * ExpB(mp_b, K));

    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(mp_a));

    for(int m = 0; m < 8; m++) {
        Complex dm = (Complex)1.0/(ExpA(mp_a, m * K/8) /*ExpA(mp_b, m * K/4)*/ * ExpB(mp_b, m * K/8));
        mpfr_mul_si(tmp, mp_b, K, GMP_RNDN); // tmp = b * K
        mpfr_mul_si(tmp, tmp, m, GMP_RNDN);  // now tmp = b * K * m
        mpfr_div_si(tmp, tmp, 4, GMP_RNDN);  // now tmp = bKm/4
        mpfr_add(tmp, tmp, mp_a, GMP_RNDN);  // now tmp = a + bKm/4
 
        mpfr_frac(tmp, tmp, GMP_RNDN);       // now tmp = {a + bmK/4}
        if(mpfr_cmp_d(tmp, .5) > 0) {
            mpfr_sub_ui(tmp, tmp, 1., GMP_RNDN);
        }

        Complex z = dm * compute_exponential_sum_via_Euler_Maclaurin(tmp, mp_b, K/8 - 1, epsilon);

        //cout << "-----"  << dm << "      "  << z << endl;
        //cout << "-----       " << direct_exponential_sum_evaluation(a, b, m * K/8, (m + 1) * K/8 - 1, epsilon);

        //S = S + dm * compute_exponential_sum_via_Euler_Maclaurin(tmp, mp_b, K/8 - 1, epsilon);
        S = S + z;
    }

    mpfr_clear(tmp);

    return S;

}

Complex compute_exponential_sum_using_theta_algorithm(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon) {
    //
    // Compute the exponential sum
    //
    // sum_{k=0}^K exp(2 pi i a k + 2 pi i b k^2)
    //
    // using Hiary's "theta sum algorithm".
    //
    

    // The sum is split into S1 + S2 + "boundary terms", which we compute shortly.
    //
    // First we compute some constants for later use, trying to keep notation consistant
    // which Hiary's paper.

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    Complex C0(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);      // (This one is not in the paper.)
    Complex C1 = I/(ExpA(mp_a, K) * ExpB(mp_b, K));
    Complex C2 = I;
    Complex C3 = (Complex)1.0/(C0 * ExpA(mp_a, K));
    Complex C4 = -I;
    Complex C5 = -C1;
    Complex C8 = -I * ExpB(mp_b, K);
    Complex CF = ExpAB(mp_a, mp_b);

    int q = to_int(a + 2 * b * K); // note that a and b are both positive, so this will do the right thing.
    Double w = a + 2 * b * K - (Double)q;
    
    int p = to_int(ceil(a));
    Double w1 = ceil(a) - a;

    int p1 = q - p;

    //----------------------------------------------------------------------------------------
    // We compute S1 first. This is the most complicated part, and involves a recursive call
    // which computes an exponential sum of a shorter length
    //----------------------------------------------------------------------------------------

    if(verbose::S1) {
        cout << "Inside S1(): q = " << q << endl;
        cout << "             w = " << w << endl;
        cout << "       a + 2bK = " << a + (Double) 2 * b * K << endl;

    }

    Complex s1 = (Complex)0;

    s1 = s1 - C1*( JBulk(w, b, p1, K, epsilon/12) + IC7(K, w, b, epsilon/12) - IC1c(K, w, b, C8, epsilon/12));
    s1 = s1 - C2*( JBulk(w1, b, p1, K, epsilon/12) + IC7(K, w1, b, epsilon/12) );
    s1 = s1 - C3 * (Complex)exp(-(Complex)2 * (Complex)PI * w * (Complex)K) * (Complex)sqrt(2.0) * IC9E(K, w, b, (epsilon/12) * exp((Double)2 * PI * w * K));

    if(verbose::S1) {
        cout << "JBulk_term1 = " << -C1*( JBulk(w, b, p1, K, epsilon/12) ) << endl;
        cout << "JBulk_term2 = " << - C2*( JBulk(w1, b, p1, K, epsilon/12)) << endl;
        cout << "IC7_term1 = " <<  -C1 * IC7(K, w, b, epsilon/12) << endl;
        cout << "IC7_term2 = " <<  -C2 * IC7(K, w1, b, epsilon/12) << endl;
        cout << "IC9E_term = " << -C3 * (Complex)exp(-(Complex)2 * (Complex)PI * w * (Complex)K) * (Complex)sqrt(2.0) * IC9E(K, w, b, (epsilon/12) * exp((Double)2 * PI * w * K)) << endl;
        cout << "IC1c_term = " << C1 * IC1c(K, w, b, C8, epsilon/12) << endl;
    }


    //Complex z = exp(-I * PI * (a * a/(2.0 * b) - .25))/sqrt((Double)2 * b);
    Complex z = exp(I * PI * .25) * CF/sqrt( (Double)2 * b );


    mpfr_t a1, b1;
    mpfr_init2(a1, mpfr_get_prec(mp_a));
    mpfr_init2(b1, mpfr_get_prec(mp_b));
    mpfr_div(a1, mp_a, mp_b, GMP_RNDN);
    mpfr_div_ui(a1, a1, 2, GMP_RNDN);  // a = a/(2b);

    mpfr_set_d(b1, -.25, GMP_RNDN);
    mpfr_div(b1, b1, mp_b, GMP_RNDN);

    s1 = s1 + z * compute_exponential_sum(a1, b1, q, (epsilon/12) * sqrt((Double)2 * b));

    if(p == 1) {
        s1 = s1 - z;
    }

    if(verbose::S1) {
        cout << "Computed S1 = " <<  s1 << endl;
    }

    mpfr_clear(a1);
    mpfr_clear(b1);

    //-----------------------------------------------------------------------
    // Now we compute S2
    //-----------------------------------------------------------------------

    Complex s2 = -C5 * IC9H(1 - w, b, epsilon/12) - C4 * IC9H(1 - w1, b, epsilon/12);
    s2 = s2 + C5 * JBoundary(2 * b * K - w1, 1 - w, b, K, epsilon/12);
    s2 = s2 + C4 * JBoundary(2 * b * K - w, 1 - w1, b, K, epsilon/12);

    if(verbose::S2) {
        cout << "Computed S2 = " << s2 << endl;
    }

    //-----------------------------------------------------------------------
    // Then the boundary terms, which are simple.
    //-----------------------------------------------------------------------


    Complex boundary_terms =  (Complex).5 + (Complex).5 / (ExpA(mp_a, K) * ExpB(mp_b, K));

    Complex S = s1 + s2 + boundary_terms;

    return S;
}

Complex compute_exponential_sum(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon, int method) {
    //
    //
    //


    stats::exponential_sum_called++;

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
    if(verbose::compute_exponential_sum) {
        cout << "compute_exponential_sum() called with a = " << a << ", b = " << b << ", K = " << K << endl;
    }
    int conjugate = normalize(mp_a, mp_b);
    a = mpfr_get_d(mp_a, GMP_RNDN);
    b = mpfr_get_d(mp_b, GMP_RNDN);
    if(verbose::compute_exponential_sum) {
        cout << "After normalization, a = " << a << ", b = " << b << ", K = " << K << endl;
    }

   
    int q = to_int(a + 2 * b * K); // note that a and b will always be positive, so this will do the right thing.
    int p = to_int(ceil(a));

    if(method == 0) {
        if(K <= 2 * pow((-LOG(epsilon)/(2 * PI)), 2) || K <= Kmin) {
            method = 1;
        }
        //else if(2.0 * b * K < 1 && b > pow((-log(epsilon))/((Double)K/(Double)8), 2)) {
        //else if(2.0 * b * K < 1) {
        else if(q <= p) {
            stats::exponential_sum_euler_maclaurin++;
            method = 3;
        }
        else {
            method = 2;
        }
    }

    if(verbose::compute_exponential_sum) {
        cout << "In compute_exponential_sum(), using method " << method << endl;
    }
   
    Complex S = (Complex)0;

    if(method == 1) {
        S = direct_exponential_sum_evaluation(mp_a, mp_b, 0, K);
    }
    else if(method == 3) {
        S = compute_exponential_sum_for_small_b(mp_a, mp_b, K, epsilon);
    }
    else {
        //S = S1(K, mp_a, mp_b, epsilon) + S2(K, mp_a, mp_b, epsilon)  + .5 + .5 / (ExpA(mp_a, K) * ExpB(mp_b, K));
        S = compute_exponential_sum_using_theta_algorithm(mp_a, mp_b, K, epsilon);
    }
    if(conjugate) {
        S = conj(S);
    }

    if(verbose::compute_exponential_sum) {
        cout << "Computed exponential_sum(";
        cout << a << ", ";
        cout << b << ", ";
        cout << K << ") = ";
        cout << S << " using method " << method << endl;

        cout << direct_exponential_sum_evaluation(mp_a, mp_b, 0, K) << endl;
    }

    return S;
}

Complex compute_exponential_sum(Double a, Double b, int K, Double epsilon, int method) {
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);
    mpfr_set_d(mp_a, to_double(a), GMP_RNDN);
    mpfr_set_d(mp_b, to_double(b), GMP_RNDN);

    Complex s = compute_exponential_sum(mp_a, mp_b, K, epsilon, method);

    mpfr_clear(mp_a);
    mpfr_clear(mp_b);
    return s;
}

int normalize(Double &a, Double &b) {
    // Normalize the input a and b so that 0 <= b <= 1/4 and 0 <= a <= 1
    // If we have to use the transformation b -> -b, a -> -a, return 1,
    // otherwise return 0

    // A return of -1 would indicate some strange (impossible) error.
    
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);

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

    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(a));
    mpfr_floor(t, a);
    mpfr_sub(a, a, t, GMP_RNDN);

    mpfr_floor(t, b);
    mpfr_sub(b, b, t, GMP_RNDN);

    mpfr_clear(t);
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
            mpfr_mul_si(a, a, -1, GMP_RNDN);                    // a = -a
            mpfr_mul_si(b, b, -1, GMP_RNDN);                    // b = -b
            if(mpfr_cmp_ui(a, 0) < 0) {                         // if(a < 0) {
                mpfr_add_ui(a, a, 1, GMP_RNDN);                 //      a = a + 1
            }                                                   // }
            return 1;
        }
        return 0;
    }
    else {
        mpfr_sub_ui(b, b, 1, GMP_RNDN);                         // b -= 1
        mpfr_mul_si(a, a, -1, GMP_RNDN);                        // a = -a
        mpfr_mul_si(b, b, -1, GMP_RNDN);                        // b = -b
        if(mpfr_cmp_ui(a, 0) < 0) {                             // if(a < 0) {
            mpfr_add_ui(a, a, 1, GMP_RNDN);                     //      a++
        }                                                       // }
        return 1;
    }
    
    return -1; // this point should not be reachable.
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






Complex compute_exponential_sums_using_theta_algorithm(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon) {
    //
    // Compute the exponential sum
    //
    // sum_{k=0}^K exp(2 pi i a k + 2 pi i b k^2)
    //
    // using Hiary's "theta sum algorithm".
    //
    

    // The sum is split into S1 + S2 + "boundary terms", which we compute shortly.
    //
    // First we compute some constants for later use, trying to keep notation consistant
    // which Hiary's paper.

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    Complex C_AK_inverse = ExpA(mp_a, K);
    Complex C_BK_inverse = ExpB(mp_b, K);
    Complex C_AK = 1.0/C_AK_inverse;
    Complex C_BK = 1.0/C_BK_inverse;
    
    Complex C_ABK = C_AK * C_BK;

    //Complex C1 = I/(ExpA(mp_a, K) * ExpB(mp_b, K));
    Complex C1 = I * C_ABK;
    //Complex C2 = I * pow(-I, j);
    //Complex C3 = exp((j + 1.0) * I * PI/4.0)/ExpA(mp_a, K);
    //Complex C4 = pow(-I, j + 1);
    Complex C5 = -C1;
    //Complex C6 = pow(I, j + 1);
    Complex C7 = -C5;
    //Complex C8 = -I * ExpB(mp_b, K);
    Complex C8 = -I * C_BK_inverse;
    Complex CF = ExpAB(mp_a, mp_b);

    int q = to_int(a + 2 * b * K); // note that a and b are both positive, so this will do the right thing.
    Double w = a + 2 * b * K - (Double)q;
    
    int p = to_int(ceil(a));
    Double w1 = ceil(a) - a;

    int p1 = q - p;

    //----------------------------------------------------------------------------------------
    // We compute S1 first. This is the most complicated part, and involves a recursive call
    // which computes an exponential sum of a shorter length
    //----------------------------------------------------------------------------------------

    if(verbose::S1) {
        cout << "Inside S1(): q = " << q << endl;
        cout << "             w = " << w << endl;
        cout << "       a + 2bK = " << a + (Double) 2 * b * K << endl;
        cout << "            w1 = " << w1 << endl;
        cout << "             p = " << p << endl;
        cout << "            p1 = " << p1 << endl;
        cout << "            C1 = " << C1 << endl;
        cout << "            C5 = " << C5 << endl;
        cout << "            C7 = " << C7 << endl;
        cout << "            C8 = " << C8 << endl;
    }

    Complex S1 = 0;

    Complex Z[j + 1];
    for(int l = 0; l <= j; l++) {
        Complex z = 0;
        for(int s = l; s <= j; s++) {
            z = z + v[s] * binomial_coefficient(s, l);
        }
        z = z * pow(I, l);
        Z[l] = z;
    }

    if(verbose::S1) {
        cout << "Z == ";
        for(int l = 0; l <= j; l++ ) {
            cout << Z[l] << " ";
        }
        cout << endl;

        cout << "v == ";
        for(int l = 0; l <= j; l++ ) {
            cout << v[l] << " ";
        }
        cout << endl;
    }

    Double Z_epsilon[j + 1];
    Double V_epsilon[j + 1];
    for(int l = 0; l <= j; l++) {
        Z_epsilon[l] = epsilon/(12.0 * abs(Z[l]) * (j + 1.0));
        V_epsilon[l] = epsilon/(12.0 * abs(v[l]) * (j + 1.0));
    }

    Complex JBulk_term1 = 0;
    for(int l = 0; l <= j; l++) {
        JBulk_term1 += Z[l] * JBulk(w, b, l, p1, K, Z_epsilon[l] );        //----------
    }
    JBulk_term1 *= -C1;

    S1 = S1 + JBulk_term1;

    Complex IC7_term1 = 0;
    for(int l = 0; l <= j; l++) {
        IC7_term1 += Z[l] * IC7(K, l, w, b, Z_epsilon[l]);                    //---------
    }

    IC7_term1 *= -C1;
    S1 = S1 + IC7_term1;

    Complex IC1c_term = 0;
    Complex IC9E_term = 0;
    
    if(2.0 * PI * w * K <= -log(epsilon) + (j + 1) * log(2.0)) {
        for(int l = 0; l <= j; l++) {
            IC1c_term += Z[l] * IC1c(K, l, w, b, C8, Z_epsilon[l]);  //---------
        }

        IC1c_term *= C1;

        Complex Z2[j + 1];
        for(int l = 0; l <= j; l++) {
            Complex z = 0;
            for(int s = l; s <= j; s++) {
                z = z + v[s] * binomial_coefficient(s, l) * pow(2, (s + 1.0)/2.0)  * exp(I * PI * (s + 1.0)/4.0);
            }
            Z2[l] = z;
        }

        for(int l = 0; l <= j; l++) {
            IC9E_term += Z2[l] * IC9E(K, l, w, b, epsilon * exp(2.0 * PI * w * K)/(12 * abs(Z2[l]) * (j + 1)) ); //----------
        }
        //IC9E_term = -IC9E_term * exp(-2.0 * PI * w * K) / ExpA(mp_a, K);
        IC9E_term = -IC9E_term * exp(-2.0 * PI * w * K) * C_AK_inverse;
    }
    else if(verbose::S1) {
        cout << "Not computing IC9E or IC1c terms because they are going to be 0." << endl;
    }

    S1 = S1 + IC9E_term;
    S1 = S1 + IC1c_term;

    Complex JBulk_term2 = 0;
    for(int l = 0; l <= j; l++) {
        JBulk_term2 += pow(-1, l) * pow(I, l + 1) * v[l] * JBulk(w1, b, l, p1, K, V_epsilon[l]);  //-----------
    }
    JBulk_term2 *= -1;

    S1 = S1 + JBulk_term2;

    Complex IC7_term2 = 0;
    for(int l = 0; l <= j; l++) {
        IC7_term2 += pow(-1, l) * pow(I, l + 1) * v[l] * IC7(K, l, w1, b, V_epsilon[l] );          //-----------
    }
    IC7_term2 *= -1;

    S1 = S1 + IC7_term2;

    if(verbose::S1) {
        cout << "JBulk_term1 = " << JBulk_term1 << endl;
        cout << "JBulk_term2 = " << JBulk_term2 << endl;
        cout << "IC7_term1 = " << IC7_term1 << endl;
        cout << "IC7_term2 = " << IC7_term2 << endl;
        cout << "IC9E_term = " << IC9E_term << endl;
        cout << "IC1c_term = " << IC1c_term << endl;

        cout << endl;
        cout << -IC7_term1 - IC1c_term << endl;
        cout << endl;
    }

    //s1 = s1 - C1*( JBulk(w, b, p1, K, epsilon/12) + IC7(K, w, b, epsilon/12) - IC1c(K, w, b, C8, epsilon/12));
    //s1 = s1 - C2*( JBulk(w1, b, p1, K, epsilon/12) + IC7(K, w1, b, epsilon/12) );
    //s1 = s1 - C3 * (Complex)exp(-(Complex)2 * (Complex)PI * w * (Complex)K) * (Complex)sqrt(2.0) * IC9E(K, w, b, (epsilon/12) * exp((Double)2 * PI * w * K));

    //Complex z = exp(-I * PI * (a * a/(2.0 * b) - .25))/sqrt((Double)2 * b);
    
    mpfr_t a1, b1;
    mpfr_init2(a1, mpfr_get_prec(mp_a));
    mpfr_init2(b1, mpfr_get_prec(mp_b));
    mpfr_div(a1, mp_a, mp_b, GMP_RNDN);
    mpfr_div_ui(a1, a1, 2, GMP_RNDN);  // a = a/(2b);

    mpfr_set_d(b1, -.25, GMP_RNDN);
    mpfr_div(b1, b1, mp_b, GMP_RNDN);

    Complex v2[j+1];

    for(int l = 0; l <= j; l++) {
        v2[l] = 0;
        for(int s = l; s <= j; s++) {
            v2[l] += v[s] * w_coefficient(mp_a, mp_b, K, l, s, CF);
        }
    }

    //s1 = s1 + z * compute_exponential_sum(a1, b1, q, (epsilon/12) * sqrt((Double)2 * b));
    S1 = S1 + compute_exponential_sums(a1, b1, j, q, v2, epsilon);

    //cout << "--------------------" << abs(v2[0]) * .5 << endl;

    if(p == 1) {
        S1 = S1 - v2[0];
    }

    if(verbose::S1) {
        cout << "Computed S1 = " <<  S1 << endl;
    }

    mpfr_clear(a1);
    mpfr_clear(b1);

    //-----------------------------------------------------------------------
    // Now we compute S2
    //-----------------------------------------------------------------------

    Complex S2 = 0;

    Complex IC9H_term1 = 0;
    for(int l = 0; l <= j; l++) {
        IC9H_term1 += pow(-1, l) * Z[l] * IC9H(K, l, 1 - w, b, Z_epsilon[l] );                           //-----------------
    }
    IC9H_term1 *= -C5;
    S2 = S2 + IC9H_term1;

    Complex IC9H_term2 = 0;
    for(int l = 0; l <= j; l++) {
        IC9H_term2 += pow(I, l + 1) * v[l] * IC9H(K, l, 1 - w1, b, V_epsilon[l] );                       //------------------
    }
    //IC9H_term2 *= C6;
    S2 = S2 + IC9H_term2;

    Complex JBoundary_term1 = 0;
    for(int l = 0; l <= j; l++) {
        JBoundary_term1 += pow(1, l) * Z[l] * JBoundary(2 * b * K - w1, 1 - w, b, l, K, Z_epsilon[l] );     //---------------
    }
    JBoundary_term1 *= C5;
    S2 = S2 + JBoundary_term1;

    Complex JBoundary_term2 = 0;
    for(int l = 0; l <= j; l++) {
        JBoundary_term2 += pow(-I, l + 1) * v[l] * JBoundary(2 * b * K - w, 1 - w1, b, l, K, V_epsilon[l] );     //-----------
    }
    //JBoundary_term2 *= C4;
    S2 = S2 + JBoundary_term2;


    //Complex s2 = -C5 * IC9H(1 - w, b, epsilon/12) - C4 * IC9H(1 - w1, b, epsilon/12);
    //s2 = s2 + C5 * JBoundary(2 * b * K - w1, 1 - w, b, K, epsilon/12);
    //s2 = s2 + C4 * JBoundary(2 * b * K - w, 1 - w1, b, K, epsilon/12);

    if(verbose::S2) {
        cout << "Computed S2 = " << S2 << endl;
        cout << "JBoundary_term1 = " << JBoundary_term1 << endl;
        cout << "JBoundary_term2 = " << JBoundary_term2 << endl;
        cout << "IC9H_term1 = " << IC9H_term1 << endl;
        cout << "IC9H_term2 = " << IC9H_term2 << endl;
    }

    //-----------------------------------------------------------------------
    // Then the boundary terms, which are simple.
    //-----------------------------------------------------------------------

    Complex boundary_terms = 0;
    for(int l = 0; l <= j; l++) {
        boundary_terms += v[l];
    }
    boundary_terms = boundary_terms * .5 * C_ABK;
    boundary_terms += .5 * v[0];

    //cout << "Boundary terms = " << boundary_terms << endl;

    Complex S = S1 + S2 + boundary_terms;

    return S;


}

Complex compute_exponential_sums_directly(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon) {
    Complex S = 0;

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    for(int l = 0; l <= j; l++) {
        if(v[l] != 0.0) {
            if(K > 200) {
                S = S + v[l] * direct_exponential_sum_evaluation2(mp_a, mp_b, l, 0, K);
            }
            else
                S = S + v[l] * direct_exponential_sum_evaluation2(a, b, l, 0, K);
        }
    }

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

Complex compute_exponential_sums(Double a, Double b, int j, int K, Complex * v, Double epsilon, int method) {
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);
    mpfr_set_d(mp_a, a, GMP_RNDN);
    mpfr_set_d(mp_b, b, GMP_RNDN);

    Complex S = compute_exponential_sums(mp_a, mp_b, j, K, v, epsilon, method);

    mpfr_clear(mp_a);
    mpfr_clear(mp_b);
    return S;
}

Complex compute_exponential_sums(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int method) {
    //
    // compute the linear combination of exponential sums
    //
    // sum_{i=0}^j v[i] * 1/K^j sum_{k=0}^K k^j exp(2 pi i a k + 2 pi i b k^2)
    //

    stats::exponential_sum_called++;

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
    if(verbose::compute_exponential_sum) {
        cout << "compute_exponential_sums() called with a = " << a << ", b = " << b << ", K = " << K << ", j = " << j << ", epsilon = " << epsilon << endl;
        cout << "   v = ";
        for(int l = 0; l <= j; l++) {
            cout << v[l] << " " << endl;
        }
    }
    int conjugate = normalize(mp_a, mp_b);
    a = mpfr_get_d(mp_a, GMP_RNDN);
    b = mpfr_get_d(mp_b, GMP_RNDN);
    if(verbose::compute_exponential_sum) {
        cout << "After normalization, a = " << a << ", b = " << b << ", K = " << K << endl;
    }
   
    int q = to_int(a + 2 * b * K); // note that a and b will always be positive, so this will do the right thing.
    int p = to_int(ceil(a));

    if(method == 0) {
        if(K <= 2 * pow((-LOG(epsilon)/(2 * PI)), 2) || K <= Kmin || K <= 10 * (j + 1)) {
            method = 1;
        }
        //else if(2.0 * b * K < 1 && b > pow((-log(epsilon))/((Double)K/(Double)8), 2)) {
        //else if(2.0 * b * K < 1) {
        else if(q <= p) {
            stats::exponential_sum_euler_maclaurin++;
            method = 3;
        }
        else {
            method = 2;
        }
    }

    if(verbose::compute_exponential_sum) {
        cout << "In compute_exponential_sum(), using method " << method << endl;
    }
   
    Complex S = (Complex)0;

    Complex v2[j + 1];

    if(conjugate) {
        for(int l = 0; l <= j; l++) {
            v2[l] = conj(v[l]);
        }
    }
    else {
        for(int l = 0; l <= j; l++) {
            v2[l] = v[l];
        }
    }

    if(method == 1) {
        // direct evaluation

        S = compute_exponential_sums_directly(mp_a, mp_b, j, K, v2, epsilon);
    }
    else if(method == 3) {
        //S = compute_exponential_sum_for_small_b(mp_a, mp_b, K, epsilon);
        S =  compute_exponential_sums_for_small_b(mp_a, mp_b, j, K, v2, epsilon);
    }
    else {
        //S = S1(K, mp_a, mp_b, epsilon) + S2(K, mp_a, mp_b, epsilon)  + .5 + .5 / (ExpA(mp_a, K) * ExpB(mp_b, K));
        S = compute_exponential_sums_using_theta_algorithm(mp_a, mp_b, j, K, v2, epsilon);
    }
    if(conjugate) {
        S = conj(S);
    }

    if(verbose::compute_exponential_sum) {
        cout << "Computed exponential_sum(";
        cout << a << ", ";
        cout << b << ", ";
        cout << j << ", ";
        cout << K << ") = ";
        cout << S << " using method " << method << endl;

        cout << v[0] * direct_exponential_sum_evaluation2(a, b, j, 0, K) << endl;
    }

    return S;

    

}

