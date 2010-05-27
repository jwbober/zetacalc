#include "theta_sums.h"
#include "zeta.h"

using namespace std;

namespace zeta_stats {
    int zeta_block_d = 0;
    int zeta_block_d_using_mpfr = 0;
    int zeta_block_d_using_mpfr_x_large = 0;
};

void print_zeta_stats() {
    cout << "zeta_block_d() called " << zeta_stats::zeta_block_d << " times." << endl;
    cout << "zeta_block_d() used zeta_block_mpfr() " << zeta_stats::zeta_block_d_using_mpfr << " times." << endl;
    cout << "       " << zeta_stats::zeta_block_d_using_mpfr_x_large << " times because K/v was too big." << endl;
}


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
        mpz_t vv;
        mpz_init(vv);
        mpfr_get_z(vv, v, GMP_RNDN);
        Complex S = zeta_block_mpfr(vv, K, t);
        mpz_clear(vv);
        return S;
    }

    Double vv = mpfr_get_d(v, GMP_RNDN);

    Double w = (K-1)/vv;
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

    Complex S = compute_exponential_sums(a, b, j, K-1, Z, exp(-20));

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

    if(verbose::zeta_block) {
        cout << "zeta block returning  " << S << endl;
        //cout << "using mpfr, answer is " << zeta_block_mpfr(v, K, t) << endl;
        cout << "   v = " << mpfr_get_d(v, GMP_RNDN) << endl;
        cout << "   K = " << K << endl;
    }

    return S;
}

Complex zeta_block_mpfr(mpfr_t v, unsigned int K, mpfr_t t) {
    mpz_t vv;
    mpz_init(vv);

    mpfr_get_z(vv, v, GMP_RNDN);

    Complex S = zeta_block_mpfr(vv, K, t);

    mpz_clear(vv);
    return S;
}

Complex zeta_block_mpfr(mpz_t v, unsigned int K, mpfr_t t) {
    //
    // Compute the sum_{k=v}^{v + K - 1) n^{-.5 + it}
    //
    // We use machine doubles for the terms in the sum, but we use
    // mpfr to accurately calculate the quantity t log n mod 2 pi
    // for each term in the sum.
 
    //
    // First we figure out how much precision we will need from mpfr.
    //
    // We want to accurately calculate t log n mod 2pi to 53 bits, which
    // means that we need to compute t log n to 53 + log_2(t log n) bits.
    //
    // For safety we add an extra 2 bits.
    //

    mpfr_t x, y, n, twopi;
    mpfr_init2(x, 53);
    mpfr_init2(y, 53);

    mpfr_log2(x, t, GMP_RNDN);                          // Now x = log2(t)
    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log(y, y, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // And y = log2(log n) (We calculate these quantities
                                                        // to low precision because we are only interested
                                                        // in their size, really. There is probably 
                                                        // a clever way to do faster using less precision.

    mpfr_add(x, x, y, GMP_RNDN);
    int mod_precision = mpfr_get_ui(x, GMP_RNDN) + 55;  // This is the precision that we need when
                                                        // need to calculate the quantity t log n
                                                        // when we mod by 2 pi

    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // y = log2(v + K) now.
    int n_precision = mpfr_get_ui(y, GMP_RNDN) + 2;     // This is the precision that we need to exactly
                                                        // represent the largest integer that will occur
                                                        // in the summation.
 
    mpfr_init2(n, n_precision);

    mpfr_init2(twopi, mod_precision);
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_clear(x);
    mpfr_init2(x, mod_precision);

    Complex S = 0.0;
    
    mpfr_set_z(n, v, GMP_RNDN);             // The summation starts with n = v
    for(unsigned int k = 0; k <= K-1; k++) {
        mpfr_log(x, n, GMP_RNDN);           // x = log(n)
        mpfr_mul(x, x, t, GMP_RNDN);        // a = t log n
        mpfr_fmod(x, x, twopi, GMP_RNDN);
        Complex z = exp(I * mpfr_get_d(x, GMP_RNDN));
        z = z/sqrt(mpfr_get_d(n, GMP_RNDN));
        S = S + z;
    
        mpfr_add_ui(n, n, 1, GMP_RNDN);     // now on the next iteration, n will be v + k + 1
    }

    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(n);
    mpfr_clear(twopi);

    return S;
}

Complex zeta_block_d(mpz_t v, int K, mpfr_t t, Double epsilon) {
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

    zeta_stats::zeta_block_d++;

    Complex S = 0;

    if(K == 0) {
        return S;
    }
    if(K <= 1) {
        zeta_stats::zeta_block_d_using_mpfr++;
        return zeta_block_mpfr(v, K, t);
    }

    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    Double x = K/vv;

    if(x > pow(tt, -1.0/6.0)) {
        zeta_stats::zeta_block_d_using_mpfr_x_large++;
        zeta_stats::zeta_block_d_using_mpfr++;
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
        zeta_stats::zeta_block_d_using_mpfr++;
        return zeta_block_mpfr(v, K, t);
    }


//    cout << "here" << endl;
    Double a[number_of_log_terms + 1];
    Double b[number_of_sqrt_terms + 1];
    mpfr_t mp_v_power, z, twopi;
    mpfr_init2(mp_v_power, mpfr_get_prec(t));       // TODO: do a better job selecting the precision here.
    mpfr_init2(z, mpfr_get_prec(t));
    mpfr_init2(twopi, mpfr_get_prec(t));
    
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_set_si(mp_v_power, 1, GMP_RNDN);
    Double v_power = 1;
    int sign = 1;
    for(int l = 1; l <= number_of_log_terms; l++) {
        if(l <= number_of_log_terms_mpfr) {
            mpfr_mul_z(mp_v_power, mp_v_power, v, GMP_RNDN);
            v_power = v_power * vv;
            mpfr_div(z, t, mp_v_power, GMP_RNDN);
            mpfr_div_si(z, z, l, GMP_RNDN);
            mpfr_fmod(z, z, twopi, GMP_RNDN);
            a[l] = sign * mpfr_get_d(z, GMP_RNDN);
        }
        else {
            v_power = v_power * vv;
            a[l] = tt * sign/(l * v_power);
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

    mpfr_set_z(z, v, GMP_RNDN);
    mpfr_log(z, z, GMP_RNDN);
    mpfr_mul(z, z, t, GMP_RNDN);
    mpfr_fmod(z, z, twopi, GMP_RNDN);
    S = S * exp(I * mpfr_get_d(z, GMP_RNDN));

    if(verbose::zeta_block_d >= 2) {
        cout << "Computed zeta_block_d = " << S << endl;
        Complex z = zeta_block_mpfr(v, K, t);
    }

    if(0) {
        Complex z1 = zeta_block_mpfr(v, K, t);
        Double logerror = log(abs(z1 - S));
        if(logerror > -19) {
            cout << "zeta_block_d() called with " << endl;
            cout << "                          v = " << vv << endl;
            cout << "                          K = " << K << endl;
            cout << "                          t = " << tt << endl;
            cout << "                        K/v = " << x << endl;

            cout << "   Number of terms in log taylor expansion is " << number_of_log_terms << endl;
            cout << "                  Number of terms using mpfr: " << number_of_log_terms_mpfr << endl;
            cout << "   Number of terms in sqrt taylor expansion is " << number_of_sqrt_terms << endl;
            cout << "Computed zeta_block_d = " << S << endl;
            cout << "      Answer should be: " << z1 << endl;
            cout << "   log of difference is " << logerror << endl;
            cout << "  The logarithm taylor coefficients we computed were: " << endl;
            for(int n = 1; n <= number_of_log_terms; n++) {
                cout << a[n] << "  ";
            }
            cout << endl;
            cout << " Now computing them again, all with mpfr... I get:" << endl;
            mpfr_set_si(mp_v_power, 1, GMP_RNDN);
            int sign = 1;
            for(int l = 1; l <= number_of_log_terms; l++) {
                mpfr_mul_z(mp_v_power, mp_v_power, v, GMP_RNDN);
                mpfr_div(z, t, mp_v_power, GMP_RNDN);
                mpfr_div_si(z, z, l, GMP_RNDN);
                mpfr_fmod(z, z, twopi, GMP_RNDN);
                a[l] = sign * mpfr_get_d(z, GMP_RNDN);
                sign = -sign;
            }
            for(int n = 1; n <= number_of_log_terms; n++) {
                cout << a[n] << "  ";
            }

            cout << endl;
        }
    }

    mpfr_clear(mp_v_power);
    mpfr_clear(z);
    mpfr_clear(twopi);

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

Complex initial_zeta_sum_mpfr(mpz_t M, mpfr_t t) {
    mpfr_t x, y;
    mpfr_init2(x, 53);
    mpfr_init2(y, 53);
    mpfr_log2(x, t, GMP_RNDN);
    mpfr_set_z(y, M, GMP_RNDN);
    mpfr_log(y, y, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);
    mpfr_add(x, x, y, GMP_RNDN);
    int mod_precision = mpfr_get_ui(x, GMP_RNDN) + 55;  // This is the precision that we need when
                                                        // need to calculate the quantity t log n
                                                        // when we mod by 2 pi

    mpfr_set_z(y, M, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);
    int n_precision = mpfr_get_ui(y, GMP_RNDN) + 2;     // This is the precision that we need to exactly
                                                        // represent the largest integer that will occur
                                                        // in the summation.
    
    if(verbose::initial_zeta_sum_mpfr) {
        cout << "In initial_zeta_sum_mpfr using " << mod_precision << " bits of precision for computation." << endl;
    }

    //int precision = mpfr_get_prec(t);
    
    mpfr_t nn, twopi;
    mpfr_init2(nn, n_precision);
    mpfr_init2(twopi, mod_precision);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_clear(x);
    mpfr_init2(x, mod_precision);

    mpz_t n;
    mpz_init(n);

    Complex S1 = 0;

    for(mpz_set_si(n, 1); mpz_cmp(n, M) <= 0; mpz_add_ui(n, n, 1)) {
        mpfr_set_z(nn, n, GMP_RNDN);
        mpfr_log(x, nn, GMP_RNDN);           // x = log(n)
        mpfr_mul(x, x, t, GMP_RNDN);         //  x = t log n

        mpfr_fmod(y, x, twopi, GMP_RNDN);

        Double z = mpfr_get_d(y, GMP_RNDN);
        S1 = S1 + exp(I * z)/sqrt(mpz_get_d(n));

//        unsigned int nn = mpfr_get_ui(n, GMP_RNDN);
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
    
    mpz_clear(n);
    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(twopi);

    return S1;
}

inline void printmp(mpfr_t x) {
    cout << mpfr_get_d(x, GMP_RNDN) << endl;
}

Complex initial_zeta_sum(mpz_t M, mpfr_t t, Double epsilon) {
    Complex S, S1, S2, S3, S4;
    S = 0.0;
    S1 = 0.0;
    S2 = 0.0;
    S3 = 0.0;
    S4 = 0.0;
    
    unsigned int m = 2;
    unsigned int mm = pow(2, m);
    unsigned int K = 0;
    mpz_t M2, M3, r, v, R;
    mpfr_t x;

    mpz_init(M2);
    mpz_init(M3);
    mpz_init(r);
    mpz_init(v);
    mpz_init(R);

    mpfr_init2(x, mpfr_get_prec(t));            //
    mpfr_cbrt(x, t, GMP_RNDN);                  //
    mpfr_div_si(x, x, mm, GMP_RNDN);            //
    mpfr_get_z(M2, x, GMP_RNDD);                // M2 = floor( t^(1/3)/2^m )


    mpz_mul_ui(M3, M2, mm);
    mpz_sub(M3, M, M3);                         // M3 = M - 2^m M2

    mpz_add_ui(R, M3, 1);                       // R = M3 + 1
    mpz_div_ui(R, R, mm);                       // R = floor((M3 + 1)/2^m)

    if(verbose::initial_zeta_sum) {
        cout << "In initial_zeta_sum(): " << endl;
        cout << "        t = " << mpfr_get_d(t, GMP_RNDN) << endl;
        cout << "        M = " << M << endl;
        cout << "       M2 = " << M2 << endl;
        cout << "       M3 = " << M3 << endl;
        cout << "        R = " << R << endl;
    }


    mpz_sub_ui(r, M2, 1);                       //
    S1 = initial_zeta_sum_mpfr(r, t);           // We first compute the sum up to M2 - 1 using mpfr




    // Next we compute the sum in block sizes of length
    // K, where K starts at 1 and goes up by powers of
    // 2 until K = 2^(m-1)
    //
    // (Really we should start at K = 2 instead of K = 1, for efficiency)
    
    S2 = 0;
    K = 1;
    mpz_set(v, M2);
    if(verbose::initial_zeta_sum) {
        cout << "   In initial_zeta_sum(), to start with, v = " << v << endl;
    }
    for(unsigned int l = 0; l <= m-1; l++) {
//        mpz_mul_ui(v, M2, K);        // v = M2 * 2^l
        for(mpz_set_ui(r, 0); mpz_cmp(r, M2) < 0; mpz_add_ui(r, r, 1)) {
            S2 = S2 + zeta_block_d(v, K, t, epsilon);
            //S2 = S2 + zeta_block_mpfr(v, K, t);
            mpz_add_ui(v, v, K);    // v = (M2 + r + 1) * 2^l
        }
        K = K * 2;
    }

    if(verbose::initial_zeta_sum) {
        cout << "                     After first loop,   v = " << v << endl;
        cout << "                                         K = " << K << endl;
    }

    // Now we compute some remainder terms with block size mm = 2^m

    S3 = 0;
    //K = mm;
    //mpz_mul_ui(v, M2, K);            // v = M2 * 2^m to begin with
    for(mpz_set_ui(r, 0); mpz_cmp(r, R) < 0; mpz_add_ui(r, r, 1)) {
        S3 = S3 + zeta_block_d(v, K, t, epsilon);
        //S3 = S3 + zeta_block_mpfr(v, K, t);
        mpz_add_ui(v, v, K);        // now v = (M2 + r + 1) * 2^m, so it will have the correct value on the next iteration of the loop
    }

    if(verbose::initial_zeta_sum) {
        cout << "                    After second loop,   v = " << v << endl;
    }

    // Finally we compute one more block, from v = (M2 + R)2^m
    // (which already has the correct value), to the end, which
    // is M, unless we have already computed the whole sum correctly.

    mpz_sub(M3, M, v);
    mpz_add_ui(M3, M3, 1);      // to include M, the number of terms we have to compute is M - v + 1
    K = mpz_get_ui(M3);

    if(verbose::initial_zeta_sum) {
        cout << "                        for final block, K = " << K << endl;
    }

    S4 = zeta_block_d(v, K, t, epsilon);
    //S4 = zeta_block_mpfr(v, K, t);

    mpz_clear(M2);
    mpz_clear(M3);
    mpz_clear(r);
    mpz_clear(v);
    mpz_clear(R);
    mpfr_clear(x);

    S = S1 + S2 + S3 + S4;

    return S;
}

Complex zeta_sum_basic(mpfr_t t) {
    mpfr_t x;
    mpz_t z;

    mpfr_init2(x, mpfr_get_prec(t));
    mpz_init(z);

    mpfr_const_pi(x, GMP_RNDN);                 // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);             // x = 2 pi
    mpfr_div(x, t, x, GMP_RNDN);                // x = t/2pi
    mpfr_sqrt(x, x, GMP_RNDN);                  // x = sqrt(t/2pi)
    mpfr_floor(x, x);                           // x = floor(sqrt(t/2pi))

    mpfr_get_z(z, x, GMP_RNDN);

    Complex S = initial_zeta_sum(z, t, exp(-20));

    mpfr_clear(x);
    mpz_clear(z);
    return S;

}

Complex zeta_sum_mpfr(mpfr_t t) {
    mpfr_t x;
    mpz_t z;

    mpfr_init2(x, mpfr_get_prec(t));
    mpz_init(z);

    mpfr_const_pi(x, GMP_RNDN);                 // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);             // x = 2 pi
    mpfr_div(x, t, x, GMP_RNDN);                // x = t/2pi
    mpfr_sqrt(x, x, GMP_RNDN);                  // x = sqrt(t/2pi)
    mpfr_floor(x, x);                           // x = floor(sqrt(t/2pi))

    mpfr_get_z(z, x, GMP_RNDN);

    Complex S = initial_zeta_sum_mpfr(z, t);

    mpfr_clear(x);
    mpz_clear(z);
    return S;
}





Complex zeta_sum(mpfr_t t) {
    //
    // 
    //

    int precision = mpfr_get_prec(t);

    mpz_t M0, M, n1, M1, R;


    int m;
    int m0 = 2;

    mpz_init(M);
    mpz_init(M0);
    mpz_init(n1);
    mpz_init(M1);
    mpz_init(R);

    mpfr_t x, y, z;
    mpfr_init2(x, precision);
    mpfr_init2(y, precision);
    mpfr_init2(z, precision);

    mpfr_cbrt(x, t, GMP_RNDN);                  // x = t^(1/3)
    mpfr_ceil(x, x);                            // x = ceil(t^(1/3))
    mpfr_get_z(M0, x, GMP_RNDN);                // M0 = ceil(t^(1/3))
    mpz_mul_ui(M0, M0, 6);                      // M0 = 3 ceil(t^(1/3))
    mpz_mul_si(M, M0, pow(2, m0));              // now M = 3 2^m0 ceil(t^{1/3})

    mpfr_const_pi(x, GMP_RNDN);                 // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);             // x = 2 pi
    mpfr_div(x, t, x, GMP_RNDN);                // x = t/2pi
    mpfr_sqrt(x, x, GMP_RNDN);                  // x = sqrt(t/2pi)
    mpfr_floor(x, x);                           // x = floor(sqrt(t/2pi))
    mpfr_get_z(n1, x, GMP_RNDN);                // n1 = floor(sqrt(t/2pi))

    cout << "n1 = ";
    mpz_out_str(0, 10, n1);
    cout << endl;

    cout << "M = ";
    mpz_out_str(0, 10, M);
    cout << endl;

    if(mpz_cmp(n1, M) <= 0) {

        cout << "Warning. t was too small so we didn't apply the theta algorithm at all." << endl;

        Complex S = initial_zeta_sum(n1, t, exp(-20));
        mpz_clear(M);
        mpz_clear(n1);
        mpz_clear(M1);
        mpz_clear(M0);
        mpfr_clear(x);
        mpfr_clear(y);
        mpfr_clear(z);
        return S;
    }

    mpfr_set_z(x, n1, GMP_RNDN);        // x = n1
    mpfr_div_z(x, x, M0, GMP_RNDN);     // x = n1/M0
    mpfr_floor(x, x);                   // x = floor(n1/M0)
    mpfr_log2(x, x, GMP_RNDN);          // x = log2( floor(n1/M0) )
    mpfr_floor(x, x);                   // x = floor(log2(floor(n1/M0)))
    m = mpfr_get_si(x, GMP_RNDN);       // m = floor(log2(floor(n1/M0)))
    
    if(m < m0) {
        cout << "m < m0... Is it going to work?" << endl;
    }
    if(m == m0) {
        cout << "m == m0... Is it going to work?" << endl;
    }
    cout << "m = " << m << endl;
    cout << "m0 = " << m0 << endl;

    mpz_set_si(M1, 2);                  // M1 = 2;
    mpz_pow_ui(M1, M1, max(m, m0));              // M1 = 2^m
    mpz_mul(M1, M1, M0);                 // M1 = 2^m M
    mpz_sub(M1, n1, M1);                // M1 = n1 - 2^m M

    cout << "M1 = ";
    mpz_out_str(0, 10, M1);
    cout << endl;

    mpz_sub_ui(M, M, 1);
    Complex S1 = initial_zeta_sum(M, t, exp(-20));
//    Complex S1_ = initial_zeta_sum_mpfr(M, t);
    mpz_add_ui(M, M, 1);

    cout << "Computed S1 =   " << S1 << endl;
//    cout << "using mpfr, get " << S1_ << endl;
    

    mpz_t r;
    mpz_init(r);

    Complex Z[13];
    compute_taylor_coefficients(t, Z);

    Complex S2 = 0;

    for(int l = m0; l <= m - 1; l++) {
        int K = pow(2, l);
        for(mpz_set_si(r, 0); mpz_cmp(r, M0) < 0; mpz_add_ui(r, r, 1)) {
//            cout << "here" << endl;
            mpfr_set_z(x, r, GMP_RNDN);                 // x = r
            mpfr_add_z(x, x, M0, GMP_RNDN);              // x = r + M0
            mpfr_mul_si(x, x, K, GMP_RNDN);             // x = (2^l)(r + M0)
            S2 = S2 + zeta_block(x, K, t, Z, exp(-20));
//            S2 = S2 + zeta_block_mpfr(x, K, t);
        }
    }

    cout << "Computed S2 =   " << S2 << endl;
    
//    mpz_t w;
//    mpz_init(w);

//    mpz_mul_ui(w, M, pow(2, m - m0));
//    mpz_sub_ui(w, w, 1);

//    Complex S2_ = initial_zeta_sum_mpfr(w, t) - S1_;

//    cout << "using mpfr, get " << S2_ << endl;

    int K = pow(2, max(m, m0));

    mpz_add_ui(R, M1, 1);
    mpz_div_ui(R, R, pow(2, max(m, m0)));
    cout << "R = ";
    mpz_out_str(0, 10, R);
    cout << endl;

    Complex S3 = 0;

    for(mpz_set_si(r, 0); mpz_cmp(r, R) < 0; mpz_add_ui(r, r, 1)) {
        mpfr_set_z(x, r, GMP_RNDN);
        mpfr_add_z(x, x, M0, GMP_RNDN);
        mpfr_mul_si(x, x, K, GMP_RNDN);
        S3 = S3 + zeta_block(x, K, t, Z, exp(-20));
//        S3 = S3 + zeta_block_mpfr(x, K, t);
    }
    
    cout << "S3 = " << S3 << endl;

    mpfr_set_z(x, R, GMP_RNDN);                     // x = R
    mpfr_add_z(x, x, M0, GMP_RNDN);                 // x = R + M0
    mpfr_mul_si(x, x, K, GMP_RNDN);                 // x = 2^(max(m, m0))(R + M0)

    mpfr_sub_z(y, x, n1, GMP_RNDN);                 // y = x - n1
    mpfr_mul_si(y, y, -1, GMP_RNDN);                // y = n1 - x
    mpfr_add_si(y, y, 1, GMP_RNDN);                 // y = n1 - x + 1

    K = mpfr_get_si(y, GMP_RNDN);                   // K = n1 - x + 1

    cout << K << endl;

    Complex S4 = zeta_block(x, K, t, Z, exp(-20));
//    Complex S4 = zeta_block_mpfr(x, K, t);

    cout << "S4 = " << S4 << endl;

    mpfr_clear(x);
    mpfr_clear(y);

    mpz_clear(M1);

    return S1 + S2 + S3 + S4;

}
