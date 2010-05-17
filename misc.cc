#include "theta_sums.h"
#include "precomputed_tables.h"

using namespace std;

Complex w_coefficient(mpfr_t mp_a, mpfr_t mp_b, int K, int s, int j) {
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    if(s > j || a + 2.0 * b * K <= 0) {
        cout << "Warning: w_coefficient called with bad input." << endl;
        return 0.0/0.0;
    }

    Complex z = ExpAB(mp_a, mp_b);

    z = z * pow(floor(a + 2 * b * K), s) * pow(2.0, -3.0 * j/2.0 - 1) * pow(b * PI, -(j + 1)/2.0) *
                pow(K, -j) * factorial(j) * sqrt(2 * PI) * exp(PI * I / 4.0 + (j - s) * 3.0 * PI * I / 4.0) * pow(2 * PI / b, s/2.0) / factorial(s);

    Complex S = 0;
    for(int l = 0; l <= j - s; l++) {
        if( (j - s - l) % 2 == 0 ) {
            Double sign = 0;
            if( ((j + l - s)/2) % 2 == 0 )
                sign = 1;
            else
                sign = -1;

            S = S + sign * (  pow(a, l) * exp(-3.0 * PI * I * (Double)l/4.0) * pow(2.0 * PI / b, l/2.0)/(factorial(l) * factorial( (j - s - l)/2 ) ) );
        }
    }

    S = S * z;

    return S;
}
