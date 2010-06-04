#include "theta_sums.h"
#include "precomputed_tables.h"

#include "w_coefficient.h"

using namespace std;

namespace wstats {
    const bool stats = true;
    int z_was_zero = 0;
    int z_wasnt_zero = 0;
};

void print_w_coefficient_stats() {
    cout << endl;
    cout << "w_coefficient statistics: " << endl;
    cout << "       Number of times z was 0: " << wstats::z_was_zero << endl;
    cout << "       Number of times z wasn't 0: " << wstats::z_wasnt_zero << endl;
}

Complex w_coefficient(Double * a_powers, Double * b_powers, Double * q_powers, Double * K_powers, int s, int j, Complex CF) {
    //
    // a_powers should be an array of length at least j - s + 1, filled with 1, a, a^2, ..., a^{j - s}
    // b_powers should be an array of length at least j + 2, filled with 1, b^{-1/2}, 1/b, b^{-3/2}, b^{-2}, ..., b^{-j - 1)/2
    // q_powers really only needs to satisfy q_powers[s] = floor(a + 2 b K)^s
    // K_powers really only needs to satisfy K_powers[j] = K^{-j}
    // CF should equal ExpAB(mp_a, mp_b)
    //
    //

    if (BUILTIN_EXPECT( s > j || a_powers[1] + 2.0 /(b_powers[2] * K_powers[1]) <= 0, 0 )) {
        cout << "Warning: w_coefficient called with bad input." << endl;
        cout << "Got called with:" << endl;
        cout << "s = " << s << endl;
        cout << "j = " << j << endl;
        cout << "a = " << a_powers[1] << endl;
        cout << "b = " << 1.0/b_powers[2] << endl;
        cout << "K = " << K_powers[1] << endl;
        return 0.0/0.0;
    }

    Complex z = CF * q_powers[s] * b_powers[j +1] * b_powers[s] * K_powers[j] * A[s][j];

//    if(stats::stats || wstats::stats) {
//        if(z == 0.0) {
//            wstats::z_was_zero++;
//        }
//        else {
//            wstats::z_wasnt_zero++;
//        }
//    }

//    Complex z = CF * pow(floor(a + 2 * b * K), s)
//                  * pow(b, -(j + 1)/2.0 - s/2.0)
//                   * pow(K, -j) * A[s][j];
//                   * pow(PI, -(j + 1)/2.0) 
//                   * pow(2 * PI, s/2.0)
//                   * pow(2.0, -3.0 * j/2.0 - 1)
//                   * factorial(j)
//                   * sqrt(2 * PI)
//                   * exp(PI * I / 4.0 + (j - s) * 3.0 * PI * I / 4.0)
//                   / factorial(s);

    Complex S = 0;
    int l = 1;
    if( (j - s) % 2 == 0) {
        l = 0;
    }
    for(; l <= j - s; l+=2)
        S = S + a_powers[l] * b_powers[l] * B[s][j][l];

    S = S * z;

    return S;
}


Complex w_coefficient_slow(mpfr_t mp_a, mpfr_t mp_b, int K, int s, int j, Complex CF) {
    //
    // note: should be passed with CF = ExpAB(mp_a, mp_b);
    //
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    if(s > j || a + 2.0 * b * K <= 0) {
        cout << "Warning: w_coefficient called with bad input." << endl;
        return 0.0/0.0;
    }

    Complex z = CF * pow(floor(a + 2 * b * K), s)
                   * pow(b * PI, -(j + 1)/2.0) 
                   * pow(2 * PI / b, s/2.0)
                   * pow(K, -j)
                   * pow(2.0, -3.0 * j/2.0 - 1)
                   * factorial(j)
                   * sqrt(2 * PI)
                   * exp(PI * I / 4.0 + (j - s) * 3.0 * PI * I / 4.0)
                   / factorial(s);

    Complex S = 0;
    for(int l = 0; l <= j - s; l++) {
        if( (j - s - l) % 2 == 0 ) {
            Double sign = 0;
            if( ((j + l - s)/2) % 2 == 0 )
                sign = 1;
            else
                sign = -1;

            S = S + sign * (  pow(a, l) * pow(2.0 * PI / b, l/2.0) *
                              exp(-3.0 * PI * I * (Double)l/4.0) /
                              ( factorial(l) * factorial( (j - s - l)/2 ) ) );
        }
    }

    S = S * z;

    return S;
}

Double infinite_sum_of_differenced_inverse_powers(Double a1, Double a2, int m, int j, Double epsilon) {
    //
    // Return sum_{k=m}^infty 1/(m + a1)^j - 1/(m + a2)^j to withing precision epsilon.
    //
    // Computed using Euler-Maclaurin summation.
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


Double sum_of_offset_inverse_powers(Double a, int m, int M, int j, Double epsilon, int method) {
    //
    // Compute the sum_{k=m}^M (k + a)^(-j)
    //
    // method = 0 corresponds to a general method which may evaluate the first few terms
    //    directly and then use Euler-Maclaurin summation for the tail.
    // method = anything else uses direct evaluation, and is just there for debugging purposes.
    
    Double S = 0;

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

    int new_m = max(m + 1, to_int(ceil(-LOG(epsilon) + 1)));
    new_m = max(new_m, 3);
    if(M != -1) {
        new_m = min(M + 1, new_m);
    }

    //int new_m = min(M + 1, max(m + 1, to_int(ceil(-LOG(epsilon)) + 1)));
    for(int k = m; k < new_m; k++) {
        S += pow( k + a, -j );
    }

    m = new_m;

    if(j == 1) {
        S = S + LOG(M + a) - LOG(m + a) + .5 * ((Double)1/(Double)(m + a) + (Double)1/(Double)(M + a));
    }
    else {
        if(M != -1) {
            S = S + (Double)1.0/(Double)(j - 1) * ( pow(m + a, 1 - j) - pow(M + a, 1 - j) ) + .5 * ( pow(m + a, -j) + pow(M + a, -j) );
        }
        else { // M == -1
            S = S + (Double)1.0/(Double)(j - 1) * pow(m + a, 1 - j) + .5 * pow(m + a, -j);
        }
    }

    Double error = epsilon + 1;

    int r = 1;
    if(M != -1) {
        while(error > epsilon) {
            Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) * ( pow(m + a, -(j + 2*r - 1)) - pow(M + a, -(j + 2 * r - 1)));

            error = abs(z);
        //cout << error << endl;
            S = S + z;
            r++;
        }
    }
    else { // M == -1
        while(error > epsilon) {
            Double z = bernoulli_table[2 * r] / ( factorial(2 * r) * factorial(j - 1) ) * factorial(2 * r - 2 + j) *  pow(m + a, -(j + 2*r - 1));

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
