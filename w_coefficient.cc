#include <cmath>
#include <iostream>
#include <complex>
#include "mpfr.h"

#include "misc.h"
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

    if (__builtin_expect( s > j || a_powers[1] + 2.0 /(b_powers[2] * K_powers[1]) <= 0, 0 )) {
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


