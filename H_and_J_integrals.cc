#include "theta_sums.h"
#include "precomputed_tables.h"

#include <iostream>
#include <cmath>

using namespace std;

Complex H_Integral_0(int j, Double a, int M, Double epsilon) {
    //
    // Compute the integral int_0^1 t^j (1 - exp(-2 pi M t))/(exp(2 pi t) - t)dt,
    // which can be written (by expanding the geometric sum) as
    //
    // sum_{m=1}^M int_0^1 t^j exp(-2 pi (a + m) t) dt
    //
    //
    // In terms of our function H(j, a, epsilon), this is
    //
    // sum_{m=1}^M H(j, a + m, epsilon/(2 j + 1)
    //
    // We don't want to compute all of these terms directly, however, because
    // M can be large. We separate this into the case of small m, where
    // we compute directly, and then large m, where the antiderivative
    // can be approximated simply. We compute "simple" approximation
    // with Euler-Maclaurin summation.
    //
    // small m is just m < C = max(j, ceil(-log(epsilon)/2pi)
    //

    stats::H_Integral_0++;

    Complex S = (Complex)0;

    //int C = min(M, max(j, (int)ceil(-log(epsilon)/(2 * PI) ) ));
    int C = min(M, max(to_int(ceil(j/(2 * PI * E))), to_int(ceil(-LOG(epsilon)/(2 * PI))) ));
//    int C = M;

    for(int m = 1; m <= C; m++) {
        Complex z = H(j, a + m, epsilon/(C + 1));
//        cout << z << endl;
        S = S + z;
    }

    if(C == M) {
        return S;
    }

    S = S + factorial(j)/two_pi_power(j + 1) * sum_of_offset_inverse_powers(a, C + 1, M, j+1, epsilon * two_pi_power(j + 1)/factorial(j));

    return S;

}


Complex J_Integral_0(Double a, Double b, int M, Double epsilon) {
    //
    // Compute the integral int_0^1 exp(-2 pi a t - 2 pi i b t^2) (1 - exp(-2 pi M t))/(exp(2 pi t) - 1) dt
    //
    // We can do this by expanding the denominator in a geometric series
    // and doing a taylor expansion in b, which puts this integral in the form
    //
    // sum_{r=0}^\infty ((2 pi i b)^r)/r! \sum_{m=1}^M int_0^1 t^j exp(-2 pi(m + a)t)dt
    stats::J_Integral_0++;

    Complex S = (Complex)0;

    int r = 0;
    Double error = 2 * epsilon;
    int sign = -1;
    
    if(b == 0) {
        return H_Integral_0(0, a, M, epsilon);
    }
    
    int N = max(1, to_int(ceil(-LOG(epsilon))));   // an approxmation of the number of terms
                                                // we compute in the taylor expansion

    Complex Ib_power = (Complex)1.0/(I * b);
    while(error > epsilon) {
        //cout << r << endl;
        sign = -sign;
        Ib_power *= (I * b);
        Complex z = sign * two_pi_over_factorial_power(r) * Ib_power;
        error = abs(z);
        z *= H_Integral_0(2 * r, a, M, epsilon/(error *N));  // TODO: We really don't need to compute this term to
                                                                                // this much precision usually. We really should figure
                                                                                // our how many terms we are going to compute
                                                                                // and then divide appropriately.
//        cout << r << endl;
        S = S + z;
        r++;
    }

    return S;
}


Complex J_Integral_1(Double a, Double b, int M, int K, Double epsilon) {
    //
    // Compute the integral int_1^K exp(-2 pi a t - 2 pi i b t^2)(1 - exp(-2 pi M t))/(exp(2 pi t) - 1) dt
    //
    // If M == -1, then we treat it as positive infinity.
    //
    
    // We truncate the integral at L, where L is given by
    stats::J_Integral_1++;

    int L = min(K, max(1, to_int(ceil(-LOG(epsilon)/(2 * PI * (1 + a))))));

    // Now we compute the integral as a sum over unit length integrals.
    // On each unit length interval, we do a change of variables to get
    // the range of integration between 0 and 1, and reduce the computation
    // to a sum of exponentials times G()

    Complex S = (Complex)0;

    //cout << L << endl;

    for(Double n = (Double)1; n <= L - 1; n = n + 1) {
        int end_point;
        if(M == -1) {
            end_point = to_int(ceil(-LOG(epsilon)/(2 * PI * n)));
        }
        else {
            end_point = min(M, to_int(ceil(-LOG(epsilon)/(2 * PI * n) ) ));
        }
        for(Double m = (Double)1; m <= end_point; m = m + 1) {
   //         cout << n << " " << m << endl;
            Complex z =  G(I*(m + a + (Complex)2.0 * I * b * n), -b, epsilon * exp( 2 * PI * (m + a) * n)/(end_point * Double(L - 1)));
   //         cout << z << endl;
            z *= EXP(-2.0 * PI * n * (m + a + I * b * n));
   //         cout << z << endl;
            S = S + z;
        }
    }

    return S;

}

Complex H_Integral_2(int j, Double a1, Double a2, Double epsilon) {
    stats::H_Integral_2++;

    Complex S = (Complex)0;

    int C = max(to_int(ceil(j/(2*PI*E))), to_int(ceil(-LOG(epsilon)/(2 * PI) ) ));
//    int C = M;

    for(int m = 1; m <= C; m++) {
        Complex z = H(j, m + a1, epsilon/(C + 1)) - H(j, m + a2, epsilon/(C + 1));
//        cout << z << endl;
        S = S + z;
    }

    S = S + factorial(j)/two_pi_power(j + 1) * infinite_sum_of_differenced_inverse_powers(a1, a2, C + 1, j+1, epsilon);

    return S;

}

Complex J_Integral_2(Double a1, Double a2, Double b, Double epsilon) {
    //
    // Compute the integral int_0^1 exp(-2 pi i b t^2) (exp(-2 pi a1 t) - exp(-2 pi a2 t))(exp(2 pi t) - 1) dt,
    //
    // which is equal to 
    //
    // lim_{M --> oo} J_Integral_0(a1, b, M) - J_Integral_0(a2, b, M)
    //
    stats::J_Integral_2++;

    Complex S = (Complex)0;

    int r = 0;
    Double error = 2 * epsilon;
    int sign = -1;
    
    if(b == 0) {
        return H_Integral_2(0, a1, a2, epsilon);
    }
   
    int N = max(1, to_int(ceil(-LOG(epsilon))));

    Complex Ib_power = (Complex)1.0/(I * b);
    while(error > epsilon) {
        sign = -sign;
        Ib_power *= (I * b);
        Complex z = sign * two_pi_over_factorial_power(r) * Ib_power;
        error = abs(z);
        //cout << error << endl;
        z *= H_Integral_2(2 * r, a1, a2, epsilon/(error * N));  // TODO: We really don't need to compute this term to
                                                                                // this much precision usually. We really should figure
                                                                                // our how many terms we are going to compute
                                                                                // and then divide appropriately.
//        cout << r << endl;
        S = S + z;
        r++;
    }

    return S;

}


