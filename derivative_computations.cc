#include <iostream>

#include "theta_sums.h"

void g_derivative_polynomial(int n, Complex * p, Complex * q, const Complex a, const Complex b) {
    // q should be an array of length n with coefficients for the n-1 derivative
    // polynomial.
    //
    // p should be an array of length n + 1 which will be filled with the coefficients
    // of the new derivative polynomial
    if(n == 0) {
        p[0] = Complex(1, 0);
        return;
    }
    if(n == 1) {
        p[0] = (Double)2 * PI * I * a;
        p[1] = (Double)4 * PI * I * b;
        return;
    }
    p[0] = (Double)2 * PI * I * a * q[0] + q[1];
    p[n-1] = (Double)4 * PI * I * b * q[n-2] + (Double)2 * PI * I * a * q[n-1];
    p[n] = (Double)4 * PI * I * b * q[n-1];
    for(int k = 1; k < n - 1; k++) {
        p[k] = (Double)4 * PI * I * b * q[k-1] + (Double)2 * PI * I * a * q[k] + (Double)(k + 1) * q[k + 1];
    }
}
/*
Complex g_derivative_at_1(int n, Complex a, Complex b) {
    Complex * pk;
    Complex * pkminus1 = new Complex[0];
    for(int k = 0; k <= n; k++) {
        pk = new Complex[k+1];
        g_derivative_polynomial(k, pk, pkminus1, a, b);
        delete [] pkminus1;
        pkminus1 = pk;
    }
    Complex S = 0;
    for(int k = 0; k <= n; k++) {
        S = S + pk[k];
    }
    S = S * exp((Double)2 * PI * I * (a + b));
    delete [] pk;
    return S;
}
*/
