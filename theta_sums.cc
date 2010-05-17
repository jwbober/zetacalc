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

Double sum_of_offset_inverse_powers(Double a, int m, int M, int j, Double epsilon, int method) {
    //
    // Compute the sum_{k=m}^M (k + a)^(-j)
    //
    // method = 0 corresponds to a general method which may evaluate the first few terms
    //    directly and then use Euler-Maclaurin summation for the tail.
    // method = anything else uses direct evaluation, and is just there for debugging purposes.
    
    Double S = 0;

    //cout << a << "  " << m << "  " << M << "  " << j << "  " << epsilon << endl;

    if(method != 0) {

        for(int k = m; k <= M; k++) {
            S += pow( k + a, -j);
        }

        return S;
    }

    int new_m = min(M + 1, max(m + 1, to_int(ceil(-LOG(epsilon)) + 1)));
    for(int k = m; k < new_m; k++) {
        S += pow( k + a, -j );
    }

    m = new_m;

    if(j == 1) {
        S = S + LOG(M + a) - LOG(m + a) + .5 * ((Double)1/(Double)(m + a) + (Double)1/(Double)(M + a));
    }
    else {
        S = S + (Double)1.0/(Double)(j - 1) * ( pow(m + a, 1 - j) - pow(M + a, 1 - j) ) + .5 * ( pow(m + a, -j) + pow(M + a, -j) );
    }

    Double error = epsilon + 1;

    int r = 1;
    while(error > epsilon) {
        Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( pow(m + a, -(j + 2*r - 1)) - pow(M + a, -(j + 2 * r - 1)));

        error = abs(z);
        //cout << error << endl;
        S = S + z;
        r++;
    }

    return S;
}

Double infinite_sum_of_differenced_inverse_powers(Double a1, Double a2, int m, int j, Double epsilon) {
    //
    // Return sum_{k=m}^infty 1/(m + a1)^j - 1/(m + a2)^j to withing precision epsilon.
    //
    // Computed using Euler-Maclauring summation.
    //
    Double S = 0;

    int new_m = max(m + 1, to_int(ceil(-LOG(epsilon)) + 1));
    for(int k = m; k < new_m; k++) {
        S += pow( k + a1, -j ) - pow(k + a2, -j);
    }

    m = new_m;

    if(j == 1) {
        //S = S + log(M + a) - log(m + a) + .5 * ((Double)1/(Double)(m + a) + (Double)1/(Double)(M + a));
        S = S + LOG(m + a2) - LOG(m + a1) + .5*( (Double)1/(m + a1) - (Double)1/(m + a2));
    }
    else {
        S = S + (Double)1.0/(Double)(j - 1) * ( pow(m + a1, 1 - j) - pow(m + a2, 1 - j) ) + .5 * ( pow(m + a1, -j) - pow(m + a2, -j) );
    }

    Double error = epsilon + 1;

    int r = 1;
    Double m_plus_a1_power = pow(m + a1, -(j + 1));
    Double m_plus_a2_power = pow(m + a2, -(j + 1));

    while(error > epsilon) {
        //Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( pow(m + a1, -(j + 2*r - 1)) - pow(m + a2, -(j + 2 * r - 1)));
        Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( m_plus_a1_power - m_plus_a2_power);

        m_plus_a1_power = m_plus_a1_power / (m + a1);
        m_plus_a1_power = m_plus_a1_power / (m + a1);
        m_plus_a2_power = m_plus_a2_power / (m + a2);
        m_plus_a2_power = m_plus_a2_power / (m + a2);

        error = abs(z);
        //cout << error << endl;
        //cout << a1 << "   " << a2 << "   " << m << "   " << j << endl;
        S = S + z;
        r++;
    }

    return S;

}


Complex S1(int K, mpfr_t mp_a, mpfr_t mp_b, Double epsilon) {
    
    
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    Complex C0(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
    Complex C4 = -I;
    Complex C5 = -I/(ExpA(mp_a, K) * ExpB(mp_b, K));
    Complex C1 = -C5;
    Complex C2 = I;
    Complex C3 = (Complex)((Double)1.0/(C0 * ExpA(mp_a, K)));
    int p = to_int(ceil(a));

    Complex C8 = -I * ExpB(mp_b, K);

    //mpfr_t tmp1;
    //mpfr_init2(tmp1, mpfr_get_prec(mp_a));
    //mpfr_mul_si(tmp1, mp_b, 2, GMP_RNDN);
    //mpfr_mul_si(tmp1, tmp1, K, GMP_RNDN); // now tmp1 = 2bK
    //mpfr_add(tmp1, tmp1, mp_a, GMP_RNDN);
    //mpfr_frac(tmp1, tmp1, GMP_RNDN);

    //Double w = mpfr_get_d(tmp1, GMP_RNDN); // now w = frac(a + 2bK)
    //int q = round(a + 2 * b * K - w);

    int q = to_int(a + 2 * b * K);
    Double w = a + 2 * b * K - (Double)q;

    Double w1 = ceil(a) - a;

    int p1 = q - p;

    if(verbose::S1) {
        cout << "Inside S1(): q = " << q << endl;
        cout << "             w = " << w << endl;
        cout << "       a + 2bK = " << a + (Double) 2 * b * K << endl;

    }

    Complex S = (Complex)0;

    //cout << w << "   " << b << "   " << "   " << p1 << "   " << K << endl;
    S = S - C1*( JBulk(w, b, p1, K, epsilon) + IC7(K, w, b, epsilon) - IC1c(K, w, b, C8, epsilon));
    S = S - C2*( JBulk(w1, b, p1, K, epsilon) + IC7(K, w1, b, epsilon) );
    
    Complex x = C3 * exp(-(Double)2 * PI * w * (Double)K) * (Complex)sqrt(2.0) * IC9E(K, w, b, epsilon * exp((Double)2 * PI * w * K));

//    S = S - C3 * exp(-(Double)2 * PI * w * K) * sqrt(2.0) * IC9E(K, w, b, epsilon * exp((Double)2 * PI * w * K));
    //cout << x << endl;
    S = S - x;

    Complex z = exp(-I * PI * (a * a/(2.0 * b) - .25))/sqrt((Double)2 * b);

    //Double a1 = a/((Double)2 * b);
    //Double b1 = -1/((Double)4 * b);

    mpfr_t a1, b1;
    mpfr_init2(a1, mpfr_get_prec(mp_a));
    mpfr_init2(b1, mpfr_get_prec(mp_b));
    mpfr_div(a1, mp_a, mp_b, GMP_RNDN);
    mpfr_div_ui(a1, a1, 2, GMP_RNDN);  // a = a/(2b);

    mpfr_set_d(b1, -.25, GMP_RNDN);
    mpfr_div(b1, b1, mp_b, GMP_RNDN);

    //cout << q << endl;
    //S = S + z * direct_exponential_sum_evaluation(a1, b1, q);
    S = S + z * compute_exponential_sum(a1, b1, q, epsilon * sqrt((Double)2 * b));

    if(p == 1) {
        S = S - z;
    }

    //mpfr_clear(tmp1);
    mpfr_clear(a1);
    mpfr_clear(b1);

    if(verbose::S1) {
        cout << "computed S1 = " << S << endl;
    }

    return S;
}

Complex S2(int K, mpfr_t mp_a, mpfr_t mp_b, Double epsilon) {
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    Complex C4 = -I;
    Complex C5 = -I/(ExpA(mp_a, K) * ExpB(mp_b, K));
    int q = to_int(a + 2 * b * K);
    Double w = a + 2 * b * K - (Double)q;

    Double w1 = ceil(a) - a;

    Complex z = -C5 * IC9H(1 - w, b, epsilon/4) - C4 * IC9H(1 - w1, b, epsilon/4);
    z = z + C5 * JBoundary(2 * b * K - w1, 1 - w, b, K, epsilon/4);
    z = z + C4 * JBoundary(2 * b * K - w, 1 - w1, b, K, epsilon/4);

    if(verbose::S2) {
        cout << "computed S2 = " << z << endl;
    }

    return z;
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

    check_condition( 0 < a + 2 * b * K && a + 2 * b * K < 2 , "Warning: in compute_exponential_sum_for_small_b(), a + 2bK was not between 0 and 2");

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


    /*
    a = a - floor(a);
    b = b - floor(b);

    if(b <= .25) {
        return 0;
    }
    if(b < .75) {
        a -= .5;
        b -= .5;
        if(b >= 0) {
            if(a < 0) {
                a = a + 1.0;
            }
            return 0;
        }
        else {
            a = -a;
            b = -b;
            if(a < 0) {
                a = a + 1.0;
            }
            return 1;
        }
        return 0;
    }
    else {
        b -= 1.0;
        a = -a;
        b = -b;
        if(a < 0) {
            a = a + 1.0;
        }
        return 1;
    }
    
    return -1; // this point should not be reachable.
    */
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
        mpfr_sub_d(a, a, .5, GMP_RNDN);  //a -= .5;
        mpfr_sub_d(b, b, .5, GMP_RNDN);  //b -= .5;
        
        if(mpfr_cmp_ui(b, 0) >= 0 ) {  //if(b >= 0) 
            if(mpfr_cmp_ui(a, 0) < 0) { //  if(a < 0) {
                mpfr_add_ui(a, a, 1, GMP_RNDN); // a = a + 1.0;
            }
            return 0;
        }
        else {
            mpfr_mul_si(a, a, -1, GMP_RNDN); // a = -a
            mpfr_mul_si(b, b, -1, GMP_RNDN); // b = -b
            if(mpfr_cmp_ui(a, 0) < 0) { //if(a < 0) {
                mpfr_add_ui(a, a, 1, GMP_RNDN); // a = a + 1
            }
            return 1;
        }
        return 0;
    }
    else {
        mpfr_sub_ui(b, b, 1, GMP_RNDN); // b -= 1
        mpfr_mul_si(a, a, -1, GMP_RNDN); // a = -a
        mpfr_mul_si(b, b, -1, GMP_RNDN); // b = -b
        if(mpfr_cmp_ui(a, 0) < 0) {  // if(a < 0) {
            mpfr_add_ui(a, a, 1, GMP_RNDN); // a++
        }
        return 1;
    }
    
    return -1; // this point should not be reachable.
}

Complex theta_sum2(Double a, Double b, int K, Double epsilon) {
    Double h = 1/(Double)K;
    Complex z1 = compute_exponential_sum(a + h, b, K, epsilon);
    Complex z2 = compute_exponential_sum(a - h, b, K, epsilon);

    return (Double)1/((Double)2 * PI * I * (Double)K) * (z1 - z2)/(2 * h);
}
