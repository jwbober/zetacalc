#include "theta_sums.h"

#include <iostream>
#include <cmath>



using namespace std;
Complex compute_exponential_sums(Double a, Double b, int j, int K, Complex * v, Double epsilon, int method) {
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);
    mpfr_set_d(mp_a, a, GMP_RNDN);
    mpfr_set_d(mp_b, b, GMP_RNDN);

    Complex S = compute_exponential_sums(mp_a, mp_b, j, K, v, epsilon, method);

    mpfr_clear(mp_a);
    mpfr_clear(mp_b);
    return S;
}

Complex compute_exponential_sums(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int method) {
    //
    // compute the linear combination of exponential sums
    //
    // sum_{i=0}^j v[i] * 1/K^j sum_{k=0}^K k^j exp(2 pi i a k + 2 pi i b k^2)
    //

    stats::exponential_sum_called++;

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
    if(verbose::compute_exponential_sum) {
        cout << "compute_exponential_sums() called with a = " << a << ", b = " << b << ", K = " << K << ", j = " << j << ", epsilon = " << epsilon << endl;
        cout << "   v = ";
        for(int l = 0; l <= j; l++) {
            cout << v[l] << " " << endl;
        }
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
        if(K <= 2 * pow((-LOG(epsilon)/(2 * PI)), 2) || K <= Kmin || K <= 10 * (j + 1)) {
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

    Complex v2[j + 1];

    if(conjugate) {
        for(int l = 0; l <= j; l++) {
            v2[l] = conj(v[l]);
        }
    }
    else {
        for(int l = 0; l <= j; l++) {
            v2[l] = v[l];
        }
    }

    if(method == 1) {
        // direct evaluation

        S = compute_exponential_sums_directly(mp_a, mp_b, j, K, v2, epsilon);
    }
    else if(method == 3) {
        //S = compute_exponential_sum_for_small_b(mp_a, mp_b, K, epsilon);
        S =  compute_exponential_sums_for_small_b(mp_a, mp_b, j, K, v2, epsilon);
    }
    else {
        //S = S1(K, mp_a, mp_b, epsilon) + S2(K, mp_a, mp_b, epsilon)  + .5 + .5 / (ExpA(mp_a, K) * ExpB(mp_b, K));
        S = compute_exponential_sums_using_theta_algorithm(mp_a, mp_b, j, K, v2, epsilon);
    }
    if(conjugate) {
        S = conj(S);
    }

    if(verbose::compute_exponential_sum) {
        cout << "Computed exponential_sum(";
        cout << a << ", ";
        cout << b << ", ";
        cout << j << ", ";
        cout << K << ") = ";
        cout << S << " using method " << method << endl;

        cout << v[0] * direct_exponential_sum_evaluation2(a, b, j, 0, K) << endl;
    }

    return S;

}

