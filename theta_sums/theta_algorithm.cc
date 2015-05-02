#include "theta_sums.h"
#include "precomputed_tables.h"
#include "log.h"

#include <cmath>
#include <iostream>

using namespace std;


Complex compute_exponential_sums_using_theta_algorithm(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int _Kmin) {
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
    // which Hiary's paper. This "cache" variable contains a bunch of commonly used
    // quantities, like powers of a, b, and K, which we pass to pretty much every function.

    theta_cache * cache = build_theta_cache(mp_a, mp_b, j, K);

    Double a = cache->a;
    Double b = cache->b;

    Complex C_BK_inverse = cache->ExpBK_inverse;
    Complex C_AK = cache->ExpAK;

    Complex C_ABK = cache->ExpABK;

    Complex C1 = I * C_ABK;
    Complex C5 = -C1;
    Complex C8 = -I * C_BK_inverse;

    int q = cache->q;
    Double w = a + 2 * b * K - (Double)q;
    
    int p = to_int(ceil(a));
    Double w1 = ceil(a) - a;

    int p1 = q - p;

    //----------------------------------------------------------------------------------------
    // We compute S1 first. This is the most complicated part, and involves a recursive call
    // which computes an exponential sum of a shorter length
    //----------------------------------------------------------------------------------------

    Complex S1 = 0;

    Complex Z[max_j + 1];
    for(int l = 0; l <= j; l++) {
        Complex z = 0;
        for(int s = l; s <= j; s++) {
            z = z + v[s] * binomial_coefficient(s, l);
        }
        z = z * I_power(l);
        if( (j + 1) * abs(z) * K * 12 < epsilon) {
            Z[l] = 0;
        }
        else {
            Z[l] = z;
        }
    }

    Double Z_epsilon[max_j + 1];
    Double V_epsilon[max_j + 1];

    complex<double> * X = new complex<double>[j + 1];

    {
        Double x = 1.0/(12 * (j + 1.0));
        for(int l = 0; l <= j; l++) {
            Z_epsilon[l] = x * epsilon/abs(Z[l]);
            V_epsilon[l] = x * epsilon/abs(v[l]);
        }
    }
    Complex JBulk_term1 = 0;
    for(int l = 0; l <= j; l++) {
        Complex x = 0;
        if(Z[l] != 0.0)  {
            //x = JBulk(w, b, l, p1, K, cache, Z_epsilon[l]) + IC7(K, l, w, b, cache, Z_epsilon[l]);
            x = IC7(K, l, w, b, cache, Z_epsilon[l]);
            JBulk_term1 += Z[l] * x;
        }
    }

    JBulk(X, w, b, j, p1, K, cache, Z_epsilon);
    for(int l = 0; l <= j; l++) {
        JBulk_term1 += Z[l] * X[l];
    }

    JBulk_term1 *= -C1;

    S1 = S1 + JBulk_term1;

    Complex IC1c_term = 0;
    Complex IC9E_term = 0;
    
    if(2.0 * PI * w * K <= -fastlog(epsilon) + (j + 1) * log(2.0)) {
        for(int l = 0; l <= j; l++) {
            Complex x = 0.0;
            if(Z[l] != 0.0)
                x = Z[l] * IC1c(K, l, w, b, C8, cache, Z_epsilon[l]);  //---------
            
            IC1c_term += x;
        
        }

        IC1c_term *= C1;

        Complex Z2[max_j + 1];
        for(int l = 0; l <= j; l++) {
            Complex z = 0;
            for(int s = l; s <= j; s++) {
                //z = z + v[s] * binomial_coefficient(s, l) * pow(2, (s + 1.0)/2.0)  * exp(I * PI * (s + 1.0)/4.0);
                z = z + v[s] * binomial_coefficient(s, l) * pow(2, (s + 1.0)/2.0) * exp_i_pi4(s + 1);       // TODO get rid of this pow
            }
            if( (j + 1) * abs(z) * K * 12 < epsilon) {
                Z2[l] = 0;
            }
            else {
                Z2[l] = z;
            }
        }



        for(int l = 0; l <= j; l++) {
            Complex x = 0;
            if(Z2[l] != 0.0)
                x = Z2[l] * IC9E(K, l, w, b, cache, epsilon * EXP(2.0 * PI * w * K)/(12 * abs(Z2[l]) * (j + 1)) ); //----------

            IC9E_term += x;

        }
        //IC9E_term = -IC9E_term * exp(-2.0 * PI * w * K) / ExpA(mp_a, K);
        IC9E_term = -IC9E_term * EXP(-2.0 * PI * w * K) * C_AK;

        S1 = S1 + IC9E_term;
        S1 = S1 + IC1c_term;

    }


    Complex JBulk_term2 = 0;

    JBulk(X, w1, b, j, p1, K, cache, V_epsilon);

    for(int l = 0; l <= j; l++) {
        Complex x = 0;
        if(v[l] != 0.0)
            //x = (Double)minus_one_power(l) * I_power(l+1) * v[l] * (JBulk(w1, b, l, p1, K, cache, V_epsilon[l]) + IC7(K, l, w1, b, cache, V_epsilon[l] )); 
            x = (Double)minus_one_power(l) * I_power(l+1) * v[l] * (X[l] + IC7(K, l, w1, b, cache, V_epsilon[l] ));

        JBulk_term2 += x;


    }
    JBulk_term2 *= -1;

    S1 = S1 + JBulk_term2;

    //s1 = s1 - C1*( JBulk(w, b, p1, K, epsilon/12) + IC7(K, w, b, epsilon/12) - IC1c(K, w, b, C8, epsilon/12));
    //s1 = s1 - C2*( JBulk(w1, b, p1, K, epsilon/12) + IC7(K, w1, b, epsilon/12) );
    //s1 = s1 - C3 * (Complex)exp(-(Complex)2 * (Complex)PI * w * (Complex)K) * (Complex)sqrt(2.0) * IC9E(K, w, b, (epsilon/12) * exp((Double)2 * PI * w * K));

    //Complex z = exp(-I * PI * (a * a/(2.0 * b) - .25))/sqrt((Double)2 * b);
    
    mpfr_t a1, b1;
    mpfr_init2(a1, mpfr_get_prec(mp_a));
    mpfr_init2(b1, mpfr_get_prec(mp_b));
    mpfr_ui_div(b1, 1, mp_b, GMP_RNDN);
    mpfr_div_2ui(b1, b1, 1, GMP_RNDN);
    mpfr_mul(a1, mp_a, b1, GMP_RNDN);

    //mpfr_mul_d(b1, b1, -.5, GMP_RNDN);
    mpfr_div_2ui(b1, b1, 1, GMP_RNDN);
    mpfr_neg(b1, b1, GMP_RNDN);

    Complex v2[max_j+1];
    compute_subsum_coefficients(v2, v, cache);



    //s1 = s1 + z * compute_exponential_sum(a1, b1, q, (epsilon/12) * sqrt((Double)2 * b));
    //mpfr_t A1, B1;
    //mpfr_init2(A1, mpfr_get_prec(a1));
    //mpfr_set(A1, a1, GMP_RNDN);
    //mpfr_init2(B1, mpfr_get_prec(b1));
    //mpfr_set(B1, b1, GMP_RNDN);
    Complex subsum = compute_exponential_sums(a1, b1, j, q, v2, epsilon, _Kmin);
    //complex<double> subsum2 = compute_exponential_sums(A1, B1, j, q, v2, epsilon, _Kmin, 4);
    //cout << mpfr_get_d(A1, GMP_RNDN) << " " << mpfr_get_d(B1, GMP_RNDN) << " " << j << " " << q << " " << epsilon << " ";
    //for(int l = 0; l <= j; l++) {
    //    cout << "\\(" << v2[l].real() << "," << v[l].imag() << "\\)" << " ";
    //}
    //cout << endl;
    //cout << "*" << subsum - subsum2 << endl;
    S1 = S1 + subsum;
    //cout << "--------------------" << abs(v2[0]) * .5 << endl;
    if(p == 1) {
        S1 = S1 - v2[0];
    }

    mpfr_clear(a1);
    mpfr_clear(b1);

    //-----------------------------------------------------------------------
    // Now we compute S2
    //-----------------------------------------------------------------------

    Complex S2 = 0;

    Complex IC9H_term1 = 0;

    JBoundary(X, 2 * b * K - w1, 1 - w, b, j, K, cache, Z_epsilon);
    for(int l = 0; l <= j; l++) {
        Complex x = 0.0;
        if(Z[l] != 0.0) {
            x = Z[l] * ((Double)minus_one_power(l) * IC9H(K, l, 1 - w, b, cache, Z_epsilon[l] ) - X[l] ); 
        }
        IC9H_term1 += x;
    }
    IC9H_term1 *= -C5;
    S2 = S2 + IC9H_term1;

    Complex JBoundary_term2 = 0;

    JBoundary(X, 2 * b * K - w, 1 - w1, b, j, K, cache, V_epsilon);
    for(int l = 0; l <= j; l++) {
        Complex x = 0.0;
        if(v[l] != 0.0)
            //x = I_power(l + 1) * v[l] * ((Double)minus_one_power(l + 1) * JBoundary(2 * b * K - w, 1 - w1, b, l, K, cache, V_epsilon[l] ) + IC9H(K, l, 1 - w1, b, cache, V_epsilon[l] )) ;     //-----------
            x = I_power(l + 1) * v[l] * ((Double)minus_one_power(l + 1) * X[l] + IC9H(K, l, 1 - w1, b, cache, V_epsilon[l] )) ;     //-----------

        JBoundary_term2 += x;
    }
    S2 = S2 + JBoundary_term2;

    //-----------------------------------------------------------------------
    // Then the boundary terms, which are simple.
    //-----------------------------------------------------------------------

    Complex boundary_terms = 0;
    for(int l = 0; l <= j; l++) {
        boundary_terms += v[l];
    }
    boundary_terms = .5 * (boundary_terms * C_ABK + v[0]);

    Complex S = S1 + S2 + boundary_terms;

    free_theta_cache(cache);
    delete [] X;

    return S;


} // end of compute_exponential_sums_using_theta_algorithm()

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


