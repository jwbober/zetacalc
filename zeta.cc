#include "theta_sums.h"
#include "zeta.h"

using namespace std;

Complex zeta_block_mpfr(mpfr_t v, int K, mpfr_t t);

void compute_taylor_coefficients(mpfr_t t, Complex Z[13]) {
    Double tt1 = mpfr_get_d(t, GMP_RNDN);
    Double tt2 = tt1 * tt1;
    Double tt3 = tt2 * tt1;
    Double tt4 = tt3 * tt1;

    Z[0] = 1.0;
    Z[1] = -.5;
    Z[2] =  .375;
    Z[3] = I * tt1 / 3.0;
    Z[4] = -I * tt1 * 5.0/12.0;
    Z[5] = I * tt1 * 9.0/20.0;
    Z[6] = - tt2 / 18.0;
    Z[7] =  tt2 / 9.0;
    Z[8] = - tt2 * 77.0/480.0;
    Z[9] = -I * tt3 / 162.0;
    Z[10] = I  * tt3 * 11.0/648.0;
    Z[11] = -I * tt3 * 133.0/4320.0;
    Z[12] = tt4 / 1944.0;

}

Complex zeta_block(mpfr_t v, int K, mpfr_t t, Complex ZZ[13], int method) {
    //
    // Compute sum_{n=v}^{v+K} exp(it log n)/sqrt(n)
    //
    // We assume that K is significantly smaller than v -- something
    // like K/V = t^{1/3}. We do a taylor expansion on exp(i t log n)
    // using just two terms, which is where this assumption is used.
    //

    if(method == 2) {
        return zeta_block_mpfr(v, K, t);
    }

    Double vv = mpfr_get_d(v, GMP_RNDN);

    Double w = K/vv;
    Double w_power = 1;

    Complex Z[13];

    for(int l = 0; l < 13; l++) {
        Z[l] = ZZ[l] * w_power;
        w_power *= w;
    }

    int j = 9;

    // Compute Z[l]
 
    mpfr_t a, b, x;

    int precision = mpfr_get_prec(t);

    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(x, precision);

    mpfr_const_pi(x, GMP_RNDN);             // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);         // x = 2 pi
    mpfr_mul(x, x, v, GMP_RNDN);            // x = 2 pi v
    mpfr_div(a, t, x, GMP_RNDN);            // a = t / (2 pi v)

//    mpfr_mul_si(x, x, -2, GMP_RNDN);        // x = -4 pi v
//    mpfr_mul(x, x, v, GMP_RNDN);            // x = -4 pi v^2
//    mpfr_div(b, t, x, GMP_RNDN);            // b = -t/ (4 pi v^2)

    mpfr_mul_si(b, v, -2, GMP_RNDN);
    mpfr_div(b, a, b, GMP_RNDN);

    Complex S = compute_exponential_sums(a, b, j, K, Z, exp(-20));

    // we don't need the previous values of a and b anymore, so
    // we can erase them.

    mpfr_log(x, v, GMP_RNDN);               // x = log v
    mpfr_mul(x, x, t, GMP_RNDN);            // x = t log v

    mpfr_const_pi(a, GMP_RNDN);
    mpfr_mul_si(a, a, 2, GMP_RNDN);
    mpfr_fmod(x, x, a, GMP_RNDN);
    
    Complex z = exp(I * mpfr_get_d(x, GMP_RNDN));
    z = z / sqrt(mpfr_get_d(v, GMP_RNDN));
    S = S * z;

//    mpfr_sin_cos(b, a, x, GMP_RNDN);        // a + ib = exp(i t log v)
//    mpfr_sqrt(x, v, GMP_RNDN);              // x = sqrt(v)
    
//    mpfr_div(a, a, x, GMP_RNDN);
//    mpfr_div(b, b, x, GMP_RNDN);            // a + ib is now (exp(i t log v)/sqrt(v))

//    Complex z(mpfr_get_d(a, GMP_RNDN), mpfr_get_d(b, GMP_RNDN));
//    S = S * z;

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(x);

    return S;
}

Complex zeta_block_mpfr(mpfr_t v, int K, mpfr_t t) {
    mpfr_t real_part, imaginary_part, a, b, c;
    mpfr_t x, y;

    int precision = mpfr_get_prec(t);

    Complex S = 0.0;
    
    mpfr_init2(x, 53);
    mpfr_init2(y, 53);

    mpfr_set_si(real_part, 0, GMP_RNDN);
    mpfr_set_si(imaginary_part, 0, GMP_RNDN);

    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(c, precision);

    for(int k = 0; k <= K; k++) {
        mpfr_add_si(c, v, k, GMP_RNDN);     // c = n (= v + k)
        mpfr_log(a, c, GMP_RNDN);           // a = log(n)
        mpfr_mul(a, a, t, GMP_RNDN);        // a = t log n
        mpfr_const_pi(b, GMP_RNDN);
        mpfr_mul_si(b, b, 2, GMP_RNDN);
        mpfr_fmod(a, a, b, GMP_RNDN);
        Complex z = exp(I * mpfr_get_d(a, GMP_RNDN));
        mpfr_sqrt(c, c, GMP_RNDN);
        z = z/mpfr_get_d(c, GMP_RNDN);
        S = S + z;
//        mpfr_sin_cos(y, x, a, GMP_RNDN);    // x + iy = exp(i t log n)
//        mpfr_sqrt(c, c, GMP_RNDN);          // c = sqrt(n)
//        mpfr_div(x, x, c, GMP_RNDN);
//        mpfr_div(y, y, c, GMP_RNDN);
//        S = S + Complex(mpfr_get_d(x, GMP_RNDN), mpfr_get_d(y, GMP_RNDN));
    }
    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(c);
    mpfr_clear(x);
    mpfr_clear(y);

    return S;
}

Complex initial_zeta_sum(mpfr_t M, mpfr_t t) {
    int precision = mpfr_get_prec(t);
    
    mpfr_t a, b, c, n, real_part, imaginary_part;
    mpfr_init2(n, precision);
    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(c, precision);
    mpfr_init2(real_part, precision);
    mpfr_init2(imaginary_part, precision);

    mpfr_set_si(real_part, 0, GMP_RNDN);
    mpfr_set_si(imaginary_part, 0, GMP_RNDN);

    Complex S1 = 0;

    for(mpfr_set_si(n, 1, GMP_RNDN); mpfr_cmp(n, M) <= 0; mpfr_add_ui(n, n, 1, GMP_RNDN)) {
        mpfr_log(a, n, GMP_RNDN);           // a = log(n)
        mpfr_mul(a, a, t, GMP_RNDN);        // a = t log n

        mpfr_const_pi(b, GMP_RNDN);
        mpfr_mul_si(b, b, 2, GMP_RNDN);
        mpfr_fmod(a, a, b, GMP_RNDN);

        Double x = mpfr_get_d(a, GMP_RNDN);
        S1 = S1 + exp(I * x)/sqrt(mpfr_get_d(n, GMP_RNDN));

        unsigned int nn = mpfr_get_ui(n, GMP_RNDN);
        if(nn % 10000 == 0) {
            cout << nn << endl;
        }

//        mpfr_sin_cos(b, a, a, GMP_RNDN);    // a + ib = exp(i t log n)
//        mpfr_sqrt(c, n, GMP_RNDN);          // c = sqrt(n)
//        mpfr_div(a, a, c, GMP_RNDN);
//        mpfr_div(b, b, c, GMP_RNDN);
//        mpfr_add(real_part, real_part, a, GMP_RNDN);
//        mpfr_add(imaginary_part, imaginary_part, b, GMP_RNDN);
    }
    
    Complex S(mpfr_get_d(real_part, GMP_RNDN), mpfr_get_d(imaginary_part, GMP_RNDN));

    mpfr_clear(n);
    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(c);
    mpfr_clear(real_part);
    mpfr_clear(imaginary_part);

    return S1;
}
Complex zeta_sum(mpfr_t t) {
    //
    // 
    //

    int precision = mpfr_get_prec(t);

    mpfr_t M, n1, M1;
    int m;

    mpfr_init2(M, precision);
    mpfr_init2(n1, precision);
    mpfr_init2(M1, precision);

    mpfr_cbrt(M, t, GMP_RNDN);
    mpfr_ceil(M, M);
    mpfr_mul_si(M, M, 3, GMP_RNDN);     // now M = 3 ceil(t^{1/3})

    mpfr_const_pi(n1, GMP_RNDN);
    mpfr_mul_si(n1, n1, 2, GMP_RNDN);
    mpfr_div(n1, t, n1, GMP_RNDN);
    mpfr_sqrt(n1, n1, GMP_RNDN);        // n1 = sqrt(t/2pi)

    mpfr_div(M1, n1, M, GMP_RNDN);
    mpfr_floor(M1, M1);
    mpfr_log2(M1, M1, GMP_RNDN);
    mpfr_floor(M1, M1);                 // m = floor(log_2(floor(n1/M)))

    m = mpfr_get_si(M1, GMP_RNDN);
    
    mpfr_set_si(M1, 2, GMP_RNDN);
    mpfr_pow_si(M1, M1, m, GMP_RNDN);
    mpfr_mul(M1, M1, M, GMP_RNDN);
    mpfr_sub(M1, n1, M1, GMP_RNDN);     // M1 = n1 - 2^m M

    cout << "m = " << m << endl;

    Double a;
    a = mpfr_get_d(M, GMP_RNDN);
    cout << "M = " << a << endl;
    a = mpfr_get_d(n1, GMP_RNDN);
    cout << "n1 = " << a << endl;
    a = mpfr_get_d(M1, GMP_RNDN);
    cout << "M1 = " << a << endl;

    cout << initial_zeta_sum(M, t) << endl;

    mpfr_clear(M);
    mpfr_clear(n1);
    mpfr_clear(M1);

    return 0.0;

}
