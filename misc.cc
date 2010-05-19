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


