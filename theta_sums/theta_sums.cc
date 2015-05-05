#include "theta_sums.h"

#include <iostream>
#include <cmath>



using namespace std;
Complex compute_exponential_sums(Double a, Double b, int j, int K, Complex * v, Double epsilon, int _Kmin, int method) {
    int working_precision = int(2 * log2(K) + 53 + 5);
    
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, working_precision);
    mpfr_init2(mp_b, working_precision);
    mpfr_set_d(mp_a, a, GMP_RNDN);
    mpfr_set_d(mp_b, b, GMP_RNDN);

    Complex S = compute_exponential_sums(mp_a, mp_b, j, K, v, epsilon, _Kmin, method);

    mpfr_clear(mp_a);
    mpfr_clear(mp_b);
    return S;
}

Complex compute_exponential_sums(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int _Kmin, int method) {
    //
    // compute the linear combination of exponential sums
    //
    // sum_{i=0}^j v[i] * 1/K^j sum_{k=0}^K k^j exp(2 pi i a k + 2 pi i b k^2)
    //
    

    if(_Kmin == 0) {
        _Kmin = Kmin;
    }

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
    int conjugate = normalize(mp_a, mp_b);
    a = mpfr_get_d(mp_a, GMP_RNDN);
    b = mpfr_get_d(mp_b, GMP_RNDN);
   
    int q = to_int(a + 2 * b * K); // note that a and b will always be positive, so this will do the right thing.
    int p = to_int(ceil(a));

    //cout << "compute exponential sums called with" << endl
    //    << "a = " << a << endl
    //    << "b = " << b << endl
    //    << "j = " << j << endl
    //    << "K = " << K << endl;

    //cout << "v =";
    //for(int k = 0; k <= j; k++) {cout << " " << v[k];}
    //cout << endl << "method = " << method << endl;

    if(method == 0) {
        //if(K <= 2 * pow((-LOG(epsilon)/(2 * PI)), 2) || K <= _Kmin || K <= 5 * (j + 1)) {
        if(K <= _Kmin || K <= 5 * (j + 1)) {
            method = 1;
        }
        //else if(2.0 * b * K < 1 && b > pow((-log(epsilon))/((Double)K/(Double)8), 2)) {
        //else if(2.0 * b * K < 1) {
        else if(q <= 3) {
            if(K <= 1000) {method = 1;} // we need some bound here, but I'm
                                        // sure exactly what it should be.
            else {method = 3;}
        }
        else {
            method = 2;
        }
        //cout << "selecting method = " << method << endl;
    }


    Complex S = (Complex)0;

    Complex v2[max_j + 1];

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

    for(int l = 0; l <= j; l++) {
        if(abs(v2[l]) * (j + 1) * K < epsilon) {
            v2[l] = 0;
        }
    }

    if(method == 1 || method > 3) {
        // direct evaluation

        S = compute_exponential_sums_directly(mp_a, mp_b, j, K, v2, epsilon, method);
    }
    else if(method == 3) {
        //S = compute_exponential_sum_for_small_b(mp_a, mp_b, K, epsilon);
        S = compute_exponential_sums_for_small_b(mp_a, mp_b, j, K, v2, epsilon);
    }
    else {
        //S = S1(K, mp_a, mp_b, epsilon) + S2(K, mp_a, mp_b, epsilon)  + .5 + .5 / (ExpA(mp_a, K) * ExpB(mp_b, K));
        S = compute_exponential_sums_using_theta_algorithm(mp_a, mp_b, j, K, v2, epsilon, _Kmin);
    }
    if(conjugate) {
        S = conj(S);
    }

    return S;

}

