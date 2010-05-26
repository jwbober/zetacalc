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

    for(int k = 0; k <= K-1; k++) {
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

Complex zeta_block_d(mpfr_t v, int K, mpfr_t t, Double epsilon) {
    //
    // This routine calculates the sum
    //
    // sum_{n=v}^{v + K - 1} n^{.5 + it) = sum_{n=v}^{v + K} exp(it log n)/sqrt(n)
    //
    // to a nominal precision of epsilon*K/sqrt(v)
    // 
    // To deal with precision issues in the calculation of the exponential
    // and to get significant speedup, we do a change of variables to
    // write this sum as 
    //
    // exp(it log v)/sqrt(v) sum_{k=0}^{K-1} exp(it log(1 + k/v))/sqrt(1 + k/v).
    //
    // Then we write 1/sqrt(v + k) as exp(log -.5(1 + k/v)),
    // so the sum is
    //
    // exp(i t log(1 + k/v) + .5 log(1 + k/v))
    //
    // We Taylor expand the logarithms with a few terms, and instead
    // of computing each term using mpfr, we just calculate the first
    // terms mod 2pi using mpfr, and then multiply by powers of k
    // in the inner sum using double arithmetic.
    //
    // Let x = K/v.
    // The number of terms we compute in the taylor expansion of the log is
    //
    // -(log(t) - log(epsilon))/log(x)
    //
    // We compute the initial terms using mpfr up to -log(t)/log(x).
    //
    // The number of terms of the in the taylor expansion for the square root
    // term is log(epsilon)/log(x).

    Complex S = 0;

    if(K == 0) {
        return S;
    }
    if(K <= 3) {
        return zeta_block_mpfr(v, K, t);
    }

    Double vv = mpfr_get_d(v, GMP_RNDN);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    Double x = K/vv;

    if(x > pow(tt, -.2)) {
        return zeta_block_mpfr(v, K, t);
    }

    int number_of_log_terms = (int)( (log(epsilon) - log(tt))/log(x));
    int number_of_log_terms_mpfr = (int)ceil(-log(tt)/log(x));
    int number_of_sqrt_terms = (int)( log(epsilon)/log(x) );

    //number_of_log_terms = 20;
    //number_of_log_terms_mpfr = number_of_log_terms;
    //number_of_sqrt_terms = 20;

    if(verbose::zeta_block_d) {
        cout << "zeta_block_d() called with " << endl;
        cout << "                          v = " << vv << endl;
        cout << "                          K = " << K << endl;
        cout << "                          t = " << tt << endl;
        cout << "                        K/v = " << x << endl;

        cout << "   Number of terms in log taylor expansion is " << number_of_log_terms << endl;
        cout << "                  Number of terms using mpfr: " << number_of_log_terms_mpfr << endl;
        cout << "   Number of terms in sqrt taylor expansion is " << number_of_sqrt_terms << endl;
    }

    Double aa = 53 + ceil((-log(tt)/log(vv))) * log2(K) + 100;
    if(aa < -log2(epsilon * sqrt(vv)/K)) {
        if(verbose::zeta_block_d) {
            cout << " --Evaluating directly using mpfr" << endl;
            cout << "    aa was " << aa << endl;
        }
        return zeta_block_mpfr(v, K, t);
    }


//    cout << "here" << endl;
    Double a[number_of_log_terms + 1];
    Double b[number_of_sqrt_terms + 1];
    mpfr_t mp_v_power, z, twopi;
    mpfr_init2(mp_v_power, mpfr_get_prec(v));
    mpfr_init2(z, mpfr_get_prec(v));
    mpfr_init2(twopi, mpfr_get_prec(v));
    
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_set_si(mp_v_power, 1, GMP_RNDN);
    Double v_power = 1;
    int sign = 1;
    for(int l = 1; l <= number_of_log_terms; l++) {
        if(l <= number_of_log_terms_mpfr) {
            mpfr_mul(mp_v_power, mp_v_power, v, GMP_RNDN);
            v_power = v_power * vv;
            mpfr_div(z, t, mp_v_power, GMP_RNDN);
            mpfr_div_si(z, z, l, GMP_RNDN);
            mpfr_fmod(z, z, twopi, GMP_RNDN);
            a[l] = sign * mpfr_get_d(z, GMP_RNDN);
        }
        else {
            v_power = v_power * vv;
            a[l] = sign/(l * v_power);
        }
        sign = -sign;
    }


    Double s = 1;
    for(int l = 1; l <= number_of_sqrt_terms; l++) {
        s = s * (-.5/vv);
        b[l] = s/l;
    }

    for(int k = 0; k <= K-1; k++) {
        Double k_power = 1;
        Double x = 0;
        Double y = 0;
        for(int l = 1; l <= number_of_log_terms; l++) {
            k_power = k_power * k;
            y = y + a[l] * k_power;
            if(l <= number_of_sqrt_terms)
                x = x + b[l] * k_power;
        }
        S = S + exp(x + I * y);
    }
    S = S / sqrt(vv);

    mpfr_log(z, v, GMP_RNDN);
    mpfr_mul(z, z, t, GMP_RNDN);
    mpfr_fmod(z, z, twopi, GMP_RNDN);
    S = S * exp(I * mpfr_get_d(z, GMP_RNDN));
    mpfr_clear(mp_v_power);
    mpfr_clear(z);
    mpfr_clear(twopi);

    if(verbose::zeta_block_d >= 2) {
        cout << "Computed zeta_block_d = " << S << endl;
        cout << "Answer should be: " << zeta_block_mpfr(v, K, t) << endl;
    }

    return S;
}

Complex zeta_block_d_stupid(mpfr_t v, int K, mpfr_t t) {
    Double vv = mpfr_get_d(v, GMP_RNDN);
    Double tt = mpfr_get_d(t, GMP_RNDN);

    Complex S = 0;
    for(int l = 0; l <= K-1; l++) {
        Double n = vv + l;
        S = S + pow(n, -.5 + I * tt);
    }
    return S;
}

Complex initial_zeta_sum_mpfr(mpfr_t M, mpfr_t t) {
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
//        if(nn % 10000 == 0) {
//            cout << nn << endl;
//        }

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

inline void printmp(mpfr_t x) {
    cout << mpfr_get_d(x, GMP_RNDN) << endl;
}

Complex initial_zeta_sum(mpfr_t M, mpfr_t t, Double epsilon) {
    Complex S = 0;
    
    int m = 9;
    int mm = pow(2, m);
    mpfr_t M2, M3;
    mpfr_t x;
    mpfr_t r;

    mpfr_init2(M2, mpfr_get_prec(M));
    mpfr_div_si(M2, M, mm, GMP_RNDN);
    mpfr_floor(M2, M2);

    mpfr_init2(M3, mpfr_get_prec(M));
    mpfr_mul_si(M3, M2, mm, GMP_RNDN);
    mpfr_sub(M3, M, M3, GMP_RNDN);

    printmp(M);
    printmp(M2);
    printmp(M3);

    mpfr_init2(x, mpfr_get_prec(M));
    mpfr_init2(r, mpfr_get_prec(M));

    mpfr_sub_si(x, M2, 1, GMP_RNDN);
    S = S + initial_zeta_sum_mpfr(x, t);


    for(int l = 0; l <= m-1; l++) {
        int K = pow(2, l);
        for(mpfr_set_si(r, 0, GMP_RNDN); mpfr_cmp(r, M2) < 0; mpfr_add_si(r, r, 1, GMP_RNDN)) {
            mpfr_add(x, M2, r, GMP_RNDN);
            mpfr_mul_si(x, x, K, GMP_RNDN);
            S = S + zeta_block_d(x, K, t, epsilon);
            //S = S + zeta_block_mpfr(x, K, t);
        }
    }

    mpfr_t R;
    mpfr_init2(R, mpfr_get_prec(M));


    mpfr_add_si(R, M3, 1, GMP_RNDN);
    mpfr_div_ui(R, R, mm, GMP_RNDN);
    mpfr_floor(R, R);

    printmp(R);

    int K = mm;
    for(mpfr_set_si(r, 0, GMP_RNDN); mpfr_cmp(r, R) < 0; mpfr_add_si(r, r, 1, GMP_RNDN)) {
        mpfr_add(x, M2, r, GMP_RNDN);
        mpfr_mul_si(x, x, K, GMP_RNDN);
        S = S + zeta_block_d(x, K, t, epsilon);
        //S = S + zeta_block_mpfr(x, K, t);
    }

    mpfr_add(x, M2, R, GMP_RNDN);
    mpfr_mul_ui(x, x, mm, GMP_RNDN);

    mpfr_mul_si(R, R, mm, GMP_RNDN);
    mpfr_sub(M3, M3, R, GMP_RNDN);
    mpfr_add_si(M3, M3, 1, GMP_RNDN);

    S = S + zeta_block_d(x, mpfr_get_si(M3, GMP_RNDN), t, epsilon);
    //S = S + zeta_block_mpfr(x, mpfr_get_si(M3, GMP_RNDN), t);

    mpfr_clear(M2);
    mpfr_clear(M3);
    mpfr_clear(x);
    mpfr_clear(r);
    return S;
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
    mpfr_floor(n1, n1);

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

    cout << initial_zeta_sum(M, t, exp(-20)) << endl;

    mpfr_clear(M);
    mpfr_clear(n1);
    mpfr_clear(M1);

    return 0.0;

}
