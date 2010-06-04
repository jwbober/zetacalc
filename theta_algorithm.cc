#include "theta_sums.h"
#include "precomputed_tables.h"

#include <cmath>
#include <iostream>


using namespace std;


Complex compute_exponential_sums_using_theta_algorithm(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon) {
    //
    // Compute the exponential sum
    //
    // sum_{k=0}^K exp(2 pi i a k + 2 pi i b k^2)
    //
    // using Hiary's "theta sum algorithm".
    //
    
    theta_cache * cache = build_theta_cache(mp_a, mp_b, j, K);

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
    Complex CF;
    if(53 + 2 + log2(epsilon * b) < 0)
        CF = ExpAB(mp_a, mp_b);
    else
        CF = exp(-I * PI * a * a * .5/b);

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
        z = z * I_power(l);
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
        IC7_term1 += Z[l] * IC7(K, l, w, b, cache, Z_epsilon[l]);                    //---------
    }

    IC7_term1 *= -C1;
    S1 = S1 + IC7_term1;

    Complex IC1c_term = 0;
    Complex IC9E_term = 0;
    
    if(2.0 * PI * w * K <= -log(epsilon) + (j + 1) * log(2.0)) {
        for(int l = 0; l <= j; l++) {
            IC1c_term += Z[l] * IC1c(K, l, w, b, C8, cache, Z_epsilon[l]);  //---------
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
            IC9E_term += Z2[l] * IC9E(K, l, w, b, cache, epsilon * exp(2.0 * PI * w * K)/(12 * abs(Z2[l]) * (j + 1)) ); //----------
        }
        //IC9E_term = -IC9E_term * exp(-2.0 * PI * w * K) / ExpA(mp_a, K);
        IC9E_term = -IC9E_term * exp(-2.0 * PI * w * K) * C_AK;
    }
    else if(verbose::S1) {
        cout << "Not computing IC9E or IC1c terms because they are going to be 0." << endl;
    }

    S1 = S1 + IC9E_term;
    S1 = S1 + IC1c_term;

    Complex JBulk_term2 = 0;
    for(int l = 0; l <= j; l++) {
        JBulk_term2 += (Double)minus_one_power(l) * I_power(l+1) * v[l] * JBulk(w1, b, l, p1, K, V_epsilon[l]);  //-----------
    }
    JBulk_term2 *= -1;

    S1 = S1 + JBulk_term2;

    Complex IC7_term2 = 0;
    for(int l = 0; l <= j; l++) {
        IC7_term2 += (Double)minus_one_power(l) * I_power(l+1) * v[l] * IC7(K, l, w1, b, cache, V_epsilon[l] );          //-----------
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

    int N = j;
    if(N == 0) {
        N = N + 1;
    }

    Double a_powers[N+1];
    Double sqrt_b_powers[N+2];
    Double q_powers[N+1];
    Double K_powers[N + 1];

    a_powers[0] = 1;
    sqrt_b_powers[0] = 1;
    Double sqrt_b_inverse = 1.0/sqrt(b);
    q_powers[0] = 1;
    K_powers[0] = 1;
    Double K_inverse = 1.0/K;
    for(int k = 1; k <= N; k++) {
        a_powers[k] = a * a_powers[k-1];
        sqrt_b_powers[k] = sqrt_b_inverse * sqrt_b_powers[k-1];
        q_powers[k] = q * q_powers[k-1];
        K_powers[k] = K_inverse * K_powers[k-1];
    }
    sqrt_b_powers[N + 1] = sqrt_b_powers[N] * sqrt_b_inverse;

    for(int l = 0; l <= j; l++) {
        v2[l] = 0;
        for(int s = l; s <= j; s++) {
            v2[l] += v[s] * w_coefficient(a_powers, sqrt_b_powers, q_powers, K_powers, l, s, CF);
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
        IC9H_term1 += (Double)minus_one_power(l) * Z[l] * IC9H(K, l, 1 - w, b, cache, Z_epsilon[l] );                           //-----------------
    }
    IC9H_term1 *= -C5;
    S2 = S2 + IC9H_term1;

    Complex IC9H_term2 = 0;
    for(int l = 0; l <= j; l++) {
        IC9H_term2 += I_power(l + 1) * v[l] * IC9H(K, l, 1 - w1, b, cache, V_epsilon[l] );                       //------------------
    }
    //IC9H_term2 *= C6;
    S2 = S2 + IC9H_term2;

    Complex JBoundary_term1 = 0;
    for(int l = 0; l <= j; l++) {
        JBoundary_term1 += Z[l] * JBoundary(2 * b * K - w1, 1 - w, b, l, K, Z_epsilon[l] );     //---------------
    }
    JBoundary_term1 *= C5;
    S2 = S2 + JBoundary_term1;

    Complex JBoundary_term2 = 0;
    for(int l = 0; l <= j; l++) {
        JBoundary_term2 += minus_I_power(l + 1) * v[l] * JBoundary(2 * b * K - w, 1 - w1, b, l, K, V_epsilon[l] );     //-----------
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

    free_theta_cache(cache);

    return S;


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


