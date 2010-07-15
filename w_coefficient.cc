#include <cmath>
#include <iostream>
#include <complex>
#include "mpfr.h"

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

Complex w_coefficient(Double * b_powers, Double * q_powers, int s, int j, theta_cache * cache, Complex * inner_sums) {
    //
    // a_powers should be an array of length at least j - s + 1, filled with 1, a, a^2, ..., a^{j - s}
    // b_powers should be an array of length at least j + 2, filled with 1, b^{-1/2}, 1/b, b^{-3/2}, b^{-2}, ..., b^{-j - 1)/2
    // q_powers really only needs to satisfy q_powers[s] = floor(a + 2 b K)^s
    // K_powers really only needs to satisfy K_powers[j] = K^{-j}
    // CF should equal ExpAB(mp_a, mp_b)
    //
    //

    /*
    if (__builtin_expect( s > j || a_powers[1] + 2.0 /(b_powers[2] * cache->K) <= 0, 0 )) {
        cout << "Warning: w_coefficient called with bad input." << endl;
        cout << "Got called with:" << endl;
        cout << "s = " << s << endl;
        cout << "j = " << j << endl;
        cout << "a = " << a_powers[1] << endl;
        cout << "b = " << 1.0/b_powers[2] << endl;
        cout << "K = " << cache->K << endl;
        return 0.0/0.0;
    }
    */

    Complex z = cache->ExpAB * q_powers[s] * b_powers[j +1] * b_powers[s] * K_power(-j, cache) * A[s][j];
    z = z * inner_sums[j - s];
    return z;
    //Complex S = 0;
    //int l = 1;
    //if( (j - s) % 2 == 0) {
    //    l = 0;
    //}
    //for(; l <= j - s; l+=2)
    //    S = S + a_powers[l] * b_powers[l] * B[j-s][l];
    //
    //S = S * z;

    //return S;
}

void compute_subsum_coefficients(Complex * v2, Complex * v, theta_cache * cache) {
 
    int j = cache->j;
    Double a = cache->a;
    Double b = cache->b;
    Double q = cache->q;
    Complex CF = cache->ExpAB;


    int N = j;
    if(N == 0) {
        N = N + 1;
    }

    Double a_powers[N+1];
    Double sqrt_b_powers[N+2];
    Double q_powers[N+1];

    a_powers[0] = 1;
    sqrt_b_powers[0] = 1;
    Double sqrt_b_inverse = 1.0/sqrt(b);
    q_powers[0] = 1;
    for(int k = 1; k <= N; k++) {
        a_powers[k] = a * a_powers[k-1];
        sqrt_b_powers[k] = sqrt_b_inverse * sqrt_b_powers[k-1];
        q_powers[k] = q * q_powers[k-1];
    }
    sqrt_b_powers[N + 1] = sqrt_b_powers[N] * sqrt_b_inverse;

    Complex inner_sums[j + 1];
    for(int k = 0; k <= j; k++) {
        Complex S = 0.0;
        int l;
        if(k % 2 == 0)
            l = 0;
        else
            l = 1;
        for(; l <= k; l+=2) {
            S = S + a_powers[l] * sqrt_b_powers[l] * B[k][l];
        }
        inner_sums[k] = S;
    }


    for(int l = 0; l <= j; l++) {
        v2[l] = 0;
        for(int s = l; s <= j; s++) {
            if(v[s] != 0.0) {
                //Complex z = w_coefficient(a_powers, sqrt_b_powers, q_powers, K_powers, l, s, CF);
                //cout << "w_coefficient(" << a << ", " << b << ", " << l << ", " << s << ",) = " << z << endl;
                //v2[l] += v[s] * w_coefficient(sqrt_b_powers, q_powers, l, s, cache, inner_sums);

                Complex z = sqrt_b_powers[s +1] * K_power(-s, cache) * A[l][s] * inner_sums[s - l];
                v2[l] += v[s] * z;

                //v2[l] += v[s] * w_coefficient_slow(mp_a, mp_b, K, l, s, CF);
            }
        }
        v2[l] *= cache->ExpAB * q_powers[l] * sqrt_b_powers[l];
    }


}

Complex IC8(int K, int j, mpfr_t mp_a, mpfr_t mp_b, theta_cache * cache) {
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    Complex z = ExpAB(mp_a, mp_b);

    //z = z * pow(2.0, -3.0 * j/2.0 - 1) * pow(b * PI, -(j + 1)/2.0) *
    //            K_power(-j, cache) * factorial(j) * sqrt(2 * PI) * exp(PI * I / 4.0 + j * 3.0 * PI * I / 4.0);

    z = z * pow(b, -(j + 1)/2.0) * K_power(-j, cache) * A[0][j];

    Complex S = 0;
    for(int l = 0; l <= j; l++) {
        if( (j - l) % 2 == 0 ) {
            //S = S + sign * (  pow(a, l) * exp(-3.0 * PI * I * (Double)l/4.0) * pow(2.0 * PI / b, l/2.0)/(factorial(l) * factorial( (j - l)/2 ) ) );
            S = S + pow(a, l) * pow(b, -l/2.0) * B[j][l];
        }
    }

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


