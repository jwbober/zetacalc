#include <cmath>
#include <iostream>
#include <complex>
#include "mpfr.h"


#include "misc.h"
#include "precomputed_tables.h"
#include "log.h"

using namespace std;

Double infinite_sum_of_differenced_inverse_powers(Double a1, Double a2, int m, int j, Double epsilon) {
    //
    // Return sum_{k=m}^infty 1/(m + a1)^j - 1/(m + a2)^j to withing precision epsilon.
    //
    // Computed using Euler-Maclaurin summation.
    //
    // We assume and a1 and a2 are both positive.
    Double S = 0;

    int p = (int)ceil( (-fastlog2(epsilon) + .61 + j * log(2 * PI) - j * (fastlog(j) + 1))/2.0 );
    int new_m = max(m + 1, (int)ceil((j + 2 * p - 1)/(2 * PI)) + 1);

    //int new_m = max(m + 1, to_int(ceil(-fastlog(epsilon)) + 1));

    for(int k = m; k < new_m; k++) {
        S += pow( k + a1, -j ) - pow(k + a2, -j);
    }

//    cout << "Length of direct computation: " << new_m - m << endl;

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

    Double m_plus_a1_pow_minus_2 = 1.0/((m + a1) * (m + a1));
    Double m_plus_a2_pow_minus_2 = 1.0/((m + a2) * (m + a2));

    while(error > epsilon) {
        //Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( pow(m + a1, -(j + 2*r - 1)) - pow(m + a2, -(j + 2 * r - 1)));
        Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( m_plus_a1_power - m_plus_a2_power);

        //m_plus_a1_power = m_plus_a1_power / (m + a1);
        //m_plus_a1_power = m_plus_a1_power / (m + a1);
        //m_plus_a2_power = m_plus_a2_power / (m + a2);
        //m_plus_a2_power = m_plus_a2_power / (m + a2);

        m_plus_a1_power *= m_plus_a1_pow_minus_2;
        m_plus_a2_power *= m_plus_a2_pow_minus_2;

        error = abs(z);
        //cout << error << endl;
        //cout << a1 << "   " << a2 << "   " << m << "   " << j << endl;
        S = S + z;
        r++;
    }

//    cout << "Number of correction terms: " << r - 1 << endl;
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


    /* The following code does not seem good to put it. 
     * It has just about no effect, but might slow things down
     * just a little bit. Probably in the cases where we
     * could return zero immediately we end up computing
     * no terms in the sum and the one correction term,
     * which doesn't take very long.
     *
    static int returning_zero = 0;

    if(j == 1) {
        if(epsilon > a * (fastlog(M) - fastlog(m) + 1)) {
            return 0.0;
        }
    }
    else {
        if(epsilon > a * (pow(m, -(j-1)) - pow(M, -(j-1)))) {
            return 0.0;
        }
    }
*/
    //cout << a << "  " << m << "  " << M << "  " << j << "  " << epsilon << endl;

    if(j == 1 && M == -1) {
        cout << "Warning: sum does not converge." << endl;
        return 1.0/0.0;
    }

    if(method != 0) {

        for(int k = m; k <= M; k++) {
            S += pow( k + a, -j);
        }

        return S;
    }
    
    // We calculate a few terms directly before we use Euler-Maclaurin summation,
    // so that the Euler-Maclaurin summation can calculate a good enough error term.

    int p = (int)ceil( (-fastlog2(epsilon) + .61 + j * log(2 * PI) - j * (fastlog(j) + 1))/2.0 );
    int new_m = max(m + 1, (int)ceil((j + 2 * p - 1)/(2 * PI)) + 1);

    //int new_m = max(m + 1, to_int(ceil(-fastlog(epsilon) + 1)));
    new_m = max(new_m, 3);
    if(M != -1) {
        new_m = min(M + 1, new_m);
    }

    //int new_m = min(M + 1, max(m + 1, to_int(ceil(-LOG(epsilon)) + 1)));
    for(int k = m; k < new_m; k++) {
        S += pow( k + a, -j );
    }

    m = new_m;

    Double m_plus_a_power;
    Double M_plus_a_power;

    Double one_over_m_plus_a = 1.0/(m + a);
    Double one_over_M_plus_a;
    if(M != -1)
        one_over_M_plus_a = 1.0/(M + a);
    else {
        one_over_M_plus_a = 0;
        M_plus_a_power = 0;
    }

    if(j == 1) {
        m_plus_a_power = one_over_m_plus_a;
        M_plus_a_power = one_over_M_plus_a;
        S = S + LOG(M + a) - LOG(m + a) + .5 * (m_plus_a_power + M_plus_a_power);
    }
    else {
        m_plus_a_power = pow(m + a, 1 - j);
        if(M != -1) {
            M_plus_a_power = pow(M + a, 1 - j);
            S = S + 1.0/(j - 1.0) * (m_plus_a_power - M_plus_a_power);
            M_plus_a_power *= one_over_M_plus_a;
            m_plus_a_power *= one_over_m_plus_a;
            S = S + .5 * (m_plus_a_power + M_plus_a_power);
            //S = S + (Double)1.0/(Double)(j - 1) * ( pow(m + a, 1 - j) - pow(M + a, 1 - j) ) + .5 * ( pow(m + a, -j) + pow(M + a, -j) );
        }
        else { // M == -1
            S = S + 1.0/(j - 1.0) * m_plus_a_power;
            m_plus_a_power *= one_over_m_plus_a;
            S = S + .5 * m_plus_a_power;
            //S = S + (Double)1.0/(Double)(j - 1) * pow(m + a, 1 - j) + .5 * pow(m + a, -j);
        }
    }

    Double error = epsilon + 1;

    Double one_over_m_plus_a_squared = one_over_m_plus_a * one_over_m_plus_a;
    Double one_over_M_plus_a_squared = one_over_M_plus_a * one_over_M_plus_a;
    m_plus_a_power *= one_over_m_plus_a;
    M_plus_a_power *= one_over_M_plus_a;

    int r = 1;
    if(M != -1) {
        while(error > epsilon) {
            //Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( pow(m + a, -(j + 2*r - 1)) - pow(M + a, -(j + 2 * r - 1)));
            Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( m_plus_a_power - M_plus_a_power);

            m_plus_a_power *= one_over_m_plus_a_squared;
            M_plus_a_power *= one_over_M_plus_a_squared;

            error = abs(z);
        //cout << error << endl;
            S = S + z;
            r++;
        }
    }
    else { // M == -1
        while(error > epsilon) {
            Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * m_plus_a_power;
            m_plus_a_power *= one_over_m_plus_a_squared;

            error = abs(z);
            S = S + z;
            r++;
        }
    }

    return S;
}

Complex ExpA(mpfr_t A, int K) {
    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(A));
    mpfr_mul_si(tmp, A, K, GMP_RNDN);
    mpfr_frac(tmp, tmp, GMP_RNDN);
    Complex S = exp(-2.0 * PI * I * mpfr_get_d(tmp, GMP_RNDN));
    mpfr_clear(tmp);
    return S;
}

Complex ExpAK(mpfr_t A, int K) {
    // Return exp(2 pi i A K)
    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(A));
    mpfr_mul_si(tmp, A, K, GMP_RNDN);
    mpfr_frac(tmp, tmp, GMP_RNDN);
    Complex S = exp(2.0 * PI * I * mpfr_get_d(tmp, GMP_RNDN));
    mpfr_clear(tmp);
    return S;
}

Complex ExpB(mpfr_t B, int K) {
    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(B));
    mpfr_mul_si(tmp, B, K, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, K, GMP_RNDN);
    mpfr_frac(tmp, tmp, GMP_RNDN);
    Complex S = exp(-2.0 * PI * I * mpfr_get_d(tmp, GMP_RNDN));
    mpfr_clear(tmp);
    return S;
}

Complex ExpBK(mpfr_t B, int K) {
    // Return exp(2 pi i B K^2)
    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(B));
    mpfr_mul_si(tmp, B, K, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, K, GMP_RNDN);
    mpfr_frac(tmp, tmp, GMP_RNDN);
    Complex S = exp(2.0 * PI * I * mpfr_get_d(tmp, GMP_RNDN));
    mpfr_clear(tmp);
    return S;
}


Complex ExpAB(mpfr_t A, mpfr_t B) {
    // Return exp(-2 pi i A^2/4B)
    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(A));
    mpfr_mul(tmp, A, A, GMP_RNDN);
    mpfr_mul_d(tmp, tmp, .25, GMP_RNDN);
    mpfr_div(tmp, tmp, B, GMP_RNDN);
    
    Complex S = exp(-2.0 * PI * I * mpfr_get_d(tmp, GMP_RNDN));
    mpfr_clear(tmp);
    return S;

}
Complex ExpABK(mpfr_t A, mpfr_t B, int K) {
    //
    //  return exp(2 pi i A K + 2 pi i B K^2)
    //

    mpfr_t real_part;
    mpfr_t imag_part;
    mpfr_t t1;
    mpfr_t t2;

    mpfr_init2(real_part, 53);
    mpfr_init2(imag_part, 53);
    
    mpfr_init2(t1, mpfr_get_prec(A));
    mpfr_init2(t2, mpfr_get_prec(A));

    mpfr_mul_si(t2, B, K, GMP_RNDN);    // t2 = bK
    mpfr_add(t2, t2, A, GMP_RNDN);      // t2 = a + bK
    mpfr_mul_si(t2, t2, K, GMP_RNDN);   // t2 = aK + bK^2

    mpfr_const_pi(t1, GMP_RNDN);        // t1 = pi
    mpfr_mul(t2, t2, t1, GMP_RNDN);     // t2 = 2 pi (aK + bK^2)
    mpfr_mul_si(t2, t2, 2, GMP_RNDN);   // t2 = 2 pi (aK + bK^2)

    mpfr_sin_cos(imag_part, real_part, t2, GMP_RNDN);

    Complex S(mpfr_get_d(real_part, GMP_RNDN), mpfr_get_d(imag_part, GMP_RNDN));

    mpfr_clear(real_part);
    mpfr_clear(imag_part);
    mpfr_clear(t1);
    mpfr_clear(t2);

    return S;
}
