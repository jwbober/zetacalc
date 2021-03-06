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

    Complex twopi_i_a = 2.0 * PI * I * a;
    Complex fourpi_i_b = 4.0 * PI * I * b;

    if(n == 1) {
        p[0] = twopi_i_a;
        p[1] = fourpi_i_b;
        return;
    }
    p[0] = twopi_i_a * q[0] + q[1];
    p[n-1] = fourpi_i_b * q[n-2] + twopi_i_a * q[n-1];
    p[n] = fourpi_i_b * q[n-1];
    for(int k = 1; k < n - 1; k++) {
        p[k] = fourpi_i_b * q[k-1] + twopi_i_a * q[k] + (Double)(k + 1) * q[k + 1];
    }
}


void g_derivative_polynomial_R(int n, Complex * p, Complex * q, const Complex a, const Double b) {
    // q should be an array of length n with coefficients for the n-1 derivative
    // polynomial.
    //
    // p should be an array of length n + 1 which will be filled with the coefficients
    // of the new derivative polynomial
    //
    // here b is real.
    if(n == 0) {
        p[0] = Complex(1, 0);
        return;
    }

    Complex twopi_i_a = 2.0 * PI * I * a;
    //Complex fourpi_i_b = 4.0 * PI * I * b;
    Double fourpi_b = 4.0 * PI * b;

    if(n == 1) {
        p[0] = twopi_i_a;
        p[1] = Complex(0.0, fourpi_b);
        return;
    }
    p[0] = twopi_i_a * q[0] + q[1];
    //p[n-1] = fourpi_i_b * q[n-2] + twopi_i_a * q[n-1];
    p[n-1] = Complex( -fourpi_b * imag(q[n-2]), fourpi_b * real(q[n-2])) + twopi_i_a * q[n-1];
    //p[n] = fourpi_i_b * q[n-1];
    p[n] = Complex( -fourpi_b * imag(q[n-1]), fourpi_b * imag(q[n-1]));
    for(int k = 1; k < n - 1; k++) {
        //p[k] = fourpi_i_b * q[k-1] + twopi_i_a * q[k] + (Double)(k + 1) * q[k + 1];
        p[k] = Complex(-fourpi_b * imag(q[k-1]), fourpi_b * real(q[k-1])) + twopi_i_a * q[k] + (Double)(k + 1) * q[k + 1];
    }
}

void g_derivative_polynomial_I(int n, Complex * p, Complex * q, const Complex a, const Double b) {
    // q should be an array of length n with coefficients for the n-1 derivative
    // polynomial.
    //
    // p should be an array of length n + 1 which will be filled with the coefficients
    // of the new derivative polynomial
    //
    // here we treat b as purely imaginary, with an implicit I in front of it.
    if(n == 0) {
        p[0] = Complex(1, 0);
        return;
    }

    Complex twopi_i_a = 2.0 * PI * I * a;
    Double fourpi_i_b = -4.0 * PI * b;

    if(n == 1) {
        p[0] = twopi_i_a;
        p[1] = fourpi_i_b;
        return;
    }
    p[0] = twopi_i_a * q[0] + q[1];
    p[n-1] = fourpi_i_b * q[n-2] + twopi_i_a * q[n-1];
    p[n] = fourpi_i_b * q[n-1];
    for(int k = 1; k < n - 1; k++) {
        p[k] = fourpi_i_b * q[k-1] + twopi_i_a * q[k] + (Double)(k + 1) * q[k + 1];
    }
}

void g_derivative_polynomial_I_over_twopi(int n, Complex * p, Complex * q, const Complex a) {
    // q should be an array of length n with coefficients for the n-1 derivative
    // polynomial.
    //
    // p should be an array of length n + 1 which will be filled with the coefficients
    // of the new derivative polynomial
    //
    // here b is specialized to be I/2pi
    if(n == 0) {
        p[0] = Complex(1, 0);
        return;
    }

    Complex twopi_i_a = 2.0 * PI * I * a;
    Double fourpi_i_b = -2.0;

    if(n == 1) {
        p[0] = twopi_i_a;
        p[1] = fourpi_i_b;
        return;
    }
    p[0] = twopi_i_a * q[0] + q[1];
    p[n-1] = fourpi_i_b * q[n-2] + twopi_i_a * q[n-1];
    p[n] = fourpi_i_b * q[n-1];
    for(int k = 1; k < n - 1; k++) {
        p[k] = fourpi_i_b * q[k-1] + twopi_i_a * q[k] + (Double)(k + 1) * q[k + 1];
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
