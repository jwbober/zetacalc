#include "theta_sums.h"
#include "precomputed_tables.h"

#include "log.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>

#include <sys/mman.h>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;
Complex IC0_method1(int j, mpfr_t mp_a, mpfr_t mp_b, const theta_cache * cache, Double epsilon) {
    MPFR_DECL_INIT(mp_a2, mpfr_get_prec(mp_a));
    MPFR_DECL_INIT(mp_b2, mpfr_get_prec(mp_a));
    MPFR_DECL_INIT(tmp, mpfr_get_prec(mp_a));
    MPFR_DECL_INIT(tmp2, mpfr_get_prec(mp_a));



    int K = cache->K;

    //int N = to_int(ceil(std::max ((Double)2.0, sqrt(2.0) * abs(fastlog(epsilon)) ) ));
    int N = to_int(ceil(std::max ((Double)2.0, sqrt(2.0 * cache->b * (Double)K * (Double)K) ) ));

    //cout << "K: " << K << endl;
    //cout << "N: " << N << endl;
    //cout << "a: " << a << endl;
    //cout << "b: " << b << endl;
    //cout << "j: " << j << endl;
    //cout << "epsilon: " << epsilon << endl;


    mpfr_mul_si(mp_a2, mp_a, K, GMP_RNDN);      //mp_a2 = aK
    mpfr_div_si(mp_a2, mp_a2, N, GMP_RNDN);     //now mp_a2 = aK/N

    Double a2 = mpfr_get_d(mp_a2, GMP_RNDN);    // a2 = aK/N
    
    //cout << "a2: " << a2 << endl;

    mpfr_floor(tmp, mp_a2);                     // tmp = floor(aK/N) 
    mpfr_sub(mp_a2, mp_a2, tmp, GMP_RNDN);      // mp_a2 = {aK/N}

    mpfr_mul_si(mp_b2, mp_b, K, GMP_RNDN);      // mp_b2 = bK
    mpfr_mul_si(mp_b2, mp_b2, K, GMP_RNDN);     // mp_b2 = bK^2
    mpfr_div_si(mp_b2, mp_b2, N * N, GMP_RNDN); // mp_b2 = bK^2/N^2

    Double b2 = mpfr_get_d(mp_b2, GMP_RNDN);    // b2 = bK^2/N^2

    //cout << "b2: " << b2 << endl;

    mpfr_floor(tmp, mp_b2);                     // tmp = floor(bK^2/N^2)
    mpfr_sub(mp_b2, mp_b2, tmp, GMP_RNDN);      // mp_b2 = {bK^2/N^2}

    Double a2_mod1 = mpfr_get_d(mp_a2, GMP_RNDN);


    Complex S = (Complex)0;

    Double N_to_the_j = pow(N, j);

    Double exp_linear_term = 0;
    Double new_epsilon = epsilon * N_to_the_j * K_power(-1, cache);

    for(int n = 0; n < N; n++) {
        mpfr_mul_si(tmp, mp_b2, n * n, GMP_RNDN);
        mpfr_mul_si(tmp2, mp_a2, n, GMP_RNDN);
        mpfr_add(tmp, tmp, tmp2, GMP_RNDN);
        mpfr_frac(tmp, tmp, GMP_RNDN);
        Double x = mpfr_get_d(tmp, GMP_RNDN);
        Complex C = Complex( cos(2 * PI * x), sin(2 * PI * x) );
        Complex z;
        z = G(a2 + (Double) 2 * (Double)n * b2, b2, n, j, new_epsilon);
        S = S + C * z;
        exp_linear_term += 2 * PI * a2_mod1;
    }
    S *= K * pow(N, -(j + 1));


    return S;
}

inline Complex IC0_method2(int j, mpfr_t mp_a, mpfr_t mp_b, const theta_cache * cache, Double epsilon) {
    Double a = cache->a;
    Double b = cache->b;
    int K = cache->K;

    Complex A = IC8(K, j, mp_a, mp_b, cache);
    Complex B = IC6(K, j, a, b, mp_a, cache, epsilon/4);
    Complex C = IC5(K, j, a, b, cache, epsilon/4);
    Complex D = IC1(K, j, a, b, cache, epsilon/4);
    // Really here we are computing IC2 - IC1...
    return A - B - C - D;
}

inline Complex IC0_method3(int j, const theta_cache * cache, Double epsilon) {
    int K = cache->K;
    Double a = cache->a;
    Double b = cache->b;
    
    Complex A = IC3(K, j, a, b, cache, epsilon/2);
    Complex B = IC4(K, j, a, b, I * cache->ExpABK, cache, epsilon/2);
    return A - B;
}

Complex IC0_method4(int j, const theta_cache * cache, Double epsilon) {

    int K = cache->K;
    Double a = cache->a;
    Double b = cache->b;

    Complex A = IC3c(K, j, a, b, cache, epsilon/2);
    Complex B = IC4c(K, j, a, b, I*cache->ExpABK, cache, epsilon/2);
    return A - B;
}

Complex IC0(int j, mpfr_t mp_a, mpfr_t mp_b, const theta_cache * cache, Double epsilon) {
    Double a = cache->a;
    Double b = cache->b;
    int K = cache->K;

    Double logepsilon = max(-fastlog(epsilon), 1);
    //if(b <= logepsilon * logepsilon * K_power(-2, cache)) {
    //    return IC0_method1(j, mp_a, mp_b, cache, epsilon);
    //}
    if(-a/(2 * b) >= 0 && -a/(2 * b) <= K) {
        if(b <= logepsilon * logepsilon * K_power(-2, cache)) {
            return IC0_method1(j, mp_a, mp_b, cache, epsilon);
        }
        return IC0_method2(j, mp_a, mp_b, cache, epsilon);
    }
    else {
        //if(b <= logepsilon * logepsilon * K_power(-2, cache)) {
        if(b <= K_power(-2, cache)) {
            return IC0_method1(j, mp_a, mp_b, cache, epsilon);
        }
        if(-a/(2 * b) > K) {
            return IC0_method3(j, cache, epsilon);
        }
        else {// -a/(2b) < 0
            return IC0_method4(j, cache, epsilon);
        }
    }
}

Complex IC1(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon) {
    //
    // C11 should be passed as I * exp(2 pi I a K + 2 PI i b K^2)
    // C12 should be passed as I * exp(-2 pi(a + 2bK)K - 2 pi i b K^2)/sqrt(2 * PI * b)
    // 
    // (NEED (-log(epsilon)^2/(K^2 * 2 pi)) < b <= 1/8K )
    //

    //if(b < (LOG(epsilon) * LOG(epsilon)/(K * K))) {
    //    cout << "************Warning: b is too small in IC1" << endl;
    //}

    Complex C11 = I * cache->ExpABK;
    //cout  << "C11 = " << C11 << endl;

    Complex S = 0;
    for(int r = 0; r <= j; r++) {
        Double z = binomial_coefficient(j, r);
        Complex y = IC7(K, r, a + 2 * b * K, b, cache, epsilon /(2.0 * z * (j + 1)));
        S = S + z * I_power(r) * y;
    }

    if( a + 2 * b * K <= - LOG(epsilon * sqrt(b))/((Double)2 * PI * K) + 1) {
        Complex C12 = I * cache->ExpBK_inverse * EXP(-2 * PI * (a + 2 * b * K) * K);

        Complex S2 = 0;
        for(int l = 0; l <= j; l++) {
            Complex coeff = 0;
            for(int r = l; r <= j; r++) {
                coeff += I_power(r) * binomial_coefficient(j, r) * binomial_coefficient(r, l);
            }
            coeff *= minus_I_power(l) * pow(2 * PI * b, -(l + 1)/2.0) * pow(K, -l);
            Complex y = G((a + 2 * b * K)/sqrt(2 * PI * b) + I * sqrt(2 * b) * (Double)K/sqrt(PI), (Double)1/(2 * PI), 0, l, epsilon/(abs(C12) * (Double)2 * abs(coeff) * (j + 1.0)));
            S2 = S2 + coeff * y;
            //cout << l << ": " << y * pow( 2* PI * b, -(l + 1.0)/2.0) << endl;
        }
        S2 *= C12;
        S = S + S2;
    }

    S *= C11;

    return S;
}

Complex IC1c(int K, int j, Double a, Double b, Complex C8, const theta_cache * cache, Double epsilon) {
    //
    // Compute C8 * exp(-2 pi a K) int_0^K t^j exp(2 pi i a t - 4 pi b K t + 2 pi i b t^2),
    //
    // ( where we expect that C8 = -i exp(-2 pi i b K^2) )
    //
    // requires that 2bK >= 1

    Complex S = (Complex)0;

    int L = min(K, max(0, to_int(ceil(-LOG(epsilon) - 2 * PI * a * (Double)K/(4 * PI * b * K)) ) ));

    for(int l = 0; l <= j; l++) {
        Complex S1 = 0;
        Double z = K_power(l, cache) * inverse_binomial_coefficient(j, l);
        for(int n = 0; n <= L - 1; n++) {
            S1 = S1 + EXP(2.0 * PI * n * (I * a - 2.0 * b * (Double)K + I * b * (Double)n) ) 
                    * G(a + (Double)2.0 * I * b * (Double)K + (Double)2.0 * b * (Double)n, b, n, l, epsilon * EXP(4 * PI * b * K * (Double)n + 2 * PI * a * K) * z, 0);
        }
        S1 = S1 * minus_I_power(l)/z;
        S = S + S1;
    }

    //for(int n = 0; n <= L - 1; n++) {
    //    S = S + EXP(2.0 * PI * n * (I * a - 2.0 * b * (Double)K + I * b * (Double)n) ) * G(a + (Double)2.0 * I * b * (Double)K + (Double)2.0 * b * (Double)n, b, epsilon * exp(4 * PI * b * K * (Double)n + 2 *
    //     PI * a * K) );
    //}

    S *= EXP(-2 * PI * a * K);
    S *= C8;
    return S;
}

// IC2 not defined anywhere

// IC3  defined inline in theta_sums.h
// IC3c defined inline in theta_sums.h


Complex IC4(int K, int j, Double a, Double b, Complex C11, const theta_cache * cache, Double epsilon) {
    // needs a + 2bK <= 0
    // b is always positive, so this will be satisfied if 
    
    Complex S = 0;
    for(int r = 0; r <= j; r++) {
        Double z = binomial_coefficient(j, r);
        S = S + z * minus_I_power(r) * IC9H(K, r, -(a + 2 * b * (Double)K), b, cache, epsilon/(z * (j + 1)));
    }
    return -S * C11;
}




Complex IC4c(int K, int j, Double a, Double b, Complex C11, const theta_cache * cache, Double epsilon) {
    // called when a > 0 and b is small but not too small.
    
    Complex S = 0;
    for(int r = 0; r <= j; r++) {
        Double z = binomial_coefficient(j, r);
        S = S + z * I_power(r) * IC9H(K, r, (a + 2 * b * (Double)K), b, cache, epsilon/(z * (j + 1)));
    }
    return S * C11;
}

Complex IC5(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi i exp(i pi /4) a - 2 pi b t^2) dt 
    // assuming that a is negative (??)
    //

    return (Double)minus_one_power(j) * I_power(j + 1) * IC7(-1, j, -a, b, cache, epsilon * K_power(j, cache)) * K_power(-j, cache);
}

Complex IC6(int K, int j, Double a, Double b, mpfr_t mp_a, const theta_cache * cache, Double epsilon) {
    //
    // Compute the integral of exp(2 pi i a t + 2 pi i b t^2) over the
    // contour C_6 = {t exp(i pi/4) | sqrt{2} K < t < infinity}.
    //
    //
    // After some changes of variables this integral is equal (within epsilon) to
    //
    // (1/sqrt(2 pi b)) exp( i pi/4 + 2 pi (i - 1)a K - 4 PI b K^2) int_0^L exp(pi (i - 1) a t/sqrt(pi b) - 4 sqrt(pi b) K t - t^2
    //
    // where L is defined below.

    

    Double x = sqrt(PI/b) * a + 4 * sqrt(PI * b) * (Double)K;
    Double y = max(-LOG(epsilon * sqrt(2 * PI) * b) - (Double)2 * PI * (Double)K*(a + (Double)2 * b * (Double)K), (Double)0);

    int L = to_int(ceil(min(sqrt(y), y/x)));
    Complex z1 = 1.0/ExpA(mp_a, K);

    //Complex z = (Double)1/sqrt(2 * PI * b) * exp(I * PI/(Double)4 + 2 * PI*(I - (Double)1) * a * (Double)K - 4 * PI * b * (Double)K * Double(K));
    //Complex z = (Double)1/sqrt(2 * PI * b) * exp(I * PI/(Double)4 - 2.0 * PI* a * (Double)K - 4 * PI * b * (Double)K * Double(K)) * z1;

    //Complex z = pow(2.0 * PI * b, -((Double)j + 1.0)/2) * exp(I * PI * ((Double)j + 1)/(Double)4 - 2.0 * PI* a * (Double)K - 4 * PI * b * (Double)K * Double(K)) * z1;
    Complex z = EXP(I * PI * ((Double)j + 1)/(Double)4 - 2.0 * PI* a * (Double)K - 4 * PI * b * (Double)K * Double(K)) * z1;

    Complex alpha = (I - (Double)1)*a/(2 * sqrt(PI * b)) - 2 * sqrt(b) * (Double)K / sqrt(PI);
    Double beta = -1/(2 * PI);

    Complex S = (Complex)0;
    
    for(int r = 0; r <= j; r++) {
     
        Complex S1 = 0;

        Double z2 = K_power(r, cache) * pow(2.0 * PI * b, ((Double)r + 1)/2) / ( pow(2, ((Double)(j-r))/2) * binomial_coefficient(j, r));

        for(int n = 0; n < L; n++) {
            Complex w1 = EXP(2 * PI * alpha * (Double)n + 2 * PI * beta * (Double)n * (Double)n);
            Complex w2 = G_I(-I * alpha - (Double)2 * I * beta * (Double)n, -beta, n, r, epsilon * z2/(abs(z) * abs(w1)));
            S1 = S1 + z * w1 * w2;
        }
        S1 = S1/z2;
        S = S + S1;
    }
    S = S;

    return S;
}

Complex IC7(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon) {
    //
    // We compute C9 * K^(-j) int_0^{sqrt(2) K} t^j exp(-sqrt(2) PI a t + sqrt(2) PI i a t - 2 PI b t^2) dt
    // where C9 = exp(-I pi (j + 1)/4)
    //
    // K = -1 corresponds to infinity, in which case we don't divide by K^j.
    //
    // a, K, b, and epsilon should satisfy:
    //   
    //  (1) K > 2 * (-log(epsilon)/2pi)^2
    //  (2) 2bK >= 1
    //  (3) a >= 0

    //assert(a <= 1);

    Double z;
    int L;
    if(K == -1) {
        //z = sqrt( 2 * max(0.0, -log(epsilon * pow(b, (j + 1.0)/2.0))) );
        z = sqrt( 2 * max(0.0, - (fastlog(epsilon) + (j + 1.0)/2.0 * fastlog(b))));
    }
    else {
        //z = sqrt( 2 * max(0.0, -log(epsilon * pow(b, (j + 1.0)/2.0) * K_power(j, cache) )) );
        z = sqrt(2 * max(0.0, - (fastlog(epsilon) + (j + 1.0)/2.0 * fastlog(b) + fastlog(K_power(j, cache)))));
    }
    if(j > z) {
        L = j + 1;
    }
    else {
        L = ceil(z);
    }
    
    if(L <= 0) {
        return 0.0;
    }

    if(K != -1 && (4 * PI * b * K * K < L * L)) {
        return IC7_method1(K, j, a, b, cache, epsilon, L);
    }

    //return IC7_method1(K, j, a, b, cache, epsilon, L);
    
    if(K == -1)
        return exp_minus_i_pi4(j + 1) * root_2pi_b_power(-(j + 1), cache) * IC7star(a * root_2pi_b_power(-1, cache), j, epsilon * root_2pi_b_power(j + 1, cache));
    else
        return exp_minus_i_pi4(j + 1) * root_2pi_b_power(-(j + 1), cache) * K_power(-j, cache) * IC7star(a * root_2pi_b_power(-1, cache), j, epsilon * root_2pi_b_power(j + 1, cache) * K_power(j, cache));
        

//    return IC7_method1(K, j, a, b, cache, epsilon, L);

}

Complex IC7_method1(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon, int L) {
    Complex C9;
    C9 = exp_minus_i_pi4(j + 1);

    Complex C10(sqrt(2.0)/2.0, sqrt(2.0)/2.0);

    //
    // Note that conditions (1) and (2) above should ensure that we
    // can make a change of variables and then truncate this
    // integral at an integer.
    //
    
    //
    // TODO: When K == -1 and j > 1, we are not truncating this integral
    // as early as we should. We should think about this a little bit more.
    //

    Complex S = (Complex)0;


    Double K_to_the_j;
    if(K == -1) {
        //Double two_pi_b_power = pow(2 * PI * b, ((Double)j + 1.0)/2);
        Double two_pi_b_power = root_2pi_b_power(j + 1, cache);

        {
            Double one_over_L = 1.0/L;
            Complex alpha = C10 * a * root_2pi_b_power(-1, cache);
            Complex s1 = alpha;
            Complex y1 = EXP(2 * PI * I * s1);
            Complex y = 1.0;
            Double newepsilon = epsilon * two_pi_b_power * one_over_L;
            for(int n = 0; n <= L-1; n = n + 1) {
                Double y3 = EXP(-n * n);
                Complex z3 = y * y3;
                Complex z = 0.0;
                if( j * fastlog(n + 1) + fastlog(abs(z3)) < fastlog(newepsilon)) {          // TODO: This is a possible source of error,
                    break;                                                                  // because this isn't necessarily decreasing.
                }
                //z = G( alpha, I/( (Double)2 * PI), n, j, newepsilon/abs(z3));
                z = G_I_over_twopi(alpha, n, j, newepsilon/abs(z3));
                z = z * z3;
                S = S + z;
                alpha = alpha + I/PI;
                y = y * y1;
                
            }
        }

        S = S * C9 * root_2pi_b_power(-(j + 1), cache);

    }
    else {
        K_to_the_j = K_power(j, cache);
        //Double two_pi_b_power = pow(2 * PI * b, ((Double)j + 1.0)/2);

        Double two_pi_b_power = root_2pi_b_power(j + 1, cache);
        {
            Double one_over_L = 1.0/L;
            Complex alpha = C10 * a * root_2pi_b_power(-1, cache);
            Complex s1 = alpha;
            Complex y1 = EXP(2 * PI * I * s1);
            Complex y = 1.0;
            Double newepsilon = epsilon * two_pi_b_power * K_to_the_j * one_over_L;
 
            for(int n = 0; n <= L-1; n = n + 1) {
                //Complex z3 = exp(2 * PI * I * n * s1- n * n);
                Double y3 = EXP(-n * n);


                Complex z3 = y * y3;

                Complex z = 0.0;
                if( j * fastlog(n + 1) + fastlog(abs(z3)) < fastlog(newepsilon)) {      // TODO: This is a possible source of error,
                    break;                                                              // because this isn't necessarily decreasing.
                }
                //z = G( alpha, I/( (Double)2 * PI), n, j, newepsilon/abs(z3));
                z = G_I_over_twopi( alpha, n, j, newepsilon/abs(z3));
                z = z * z3;
                S = S + z;
                alpha = alpha + I/PI;
                y = y * y1;

            }
        }


        S = S * C9 * K_power(-j, cache) * root_2pi_b_power(-(j + 1), cache);
    }
    return S;




}

Complex IC7star(Double a, int j, Double epsilon) {
    Complex C10(sqrt(2.0)/2.0, sqrt(2.0)/2.0);

    Double z;
    int L;
    
    //z = sqrt( 2 * max(0.0, - (fastlog(epsilon) + (j + 1.0)/2.0 * fastlog(b))));
    z = sqrt( 2 * max(0, -fastlog(epsilon)) );
    if(j > z) {
        L = j + 1;
    }
    else {
        L = ceil(z);
    }
    
    L = max(0, L);

    if(L == 0) {
        return 0.0;
    }

    Complex S = (Complex)0;

    {
        Double one_over_L = 1.0/L;
        Complex alpha = C10 * a;
        Complex s1 = alpha;
        Complex y1 = EXP(2 * PI * I * s1);
        Complex y = 1.0;
        Double newepsilon = epsilon * one_over_L;
        for(int n = 0; n <= L-1; n = n + 1) {
            //Complex z3 = exp(2 * PI * I * n * s1- n * n);
            Double y3 = EXP(-n * n);
            Complex z3 = y * y3;
            Complex z = 0.0;
            if( j * fastlog(n + 1) + fastlog(abs(z3)) < fastlog(newepsilon)) {          // TODO: This is a possible source of error,
                break;                                                                  // because this isn't necessarily decreasing.
            }
            //z = G( alpha, I/( (Double)2 * PI), n, j, newepsilon/abs(z3));
            z = G_I_over_twopi(alpha, n, j, newepsilon/abs(z3));
            z = z * z3;

            S = S + z;
            alpha = alpha + I/PI;
            y = y * y1;
        }
    }

    return S;

}

/*
Complex IC8(int K, int j, mpfr_t mp_a, mpfr_t mp_b, theta_cache * cache) {
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    Complex z = ExpAB(mp_a, mp_b);

    z = z * pow(2.0, -3.0 * j/2.0 - 1) * pow(b * PI, -(j + 1)/2.0) *
                K_power(-j, cache) * factorial(j) * sqrt(2 * PI) * exp(PI * I / 4.0 + j * 3.0 * PI * I / 4.0);

    Complex S = 0;
    for(int l = 0; l <= j; l++) {
        if( (j - l) % 2 == 0 ) {
            Double sign = 0;
            if( ((j + l)/2) % 2 == 0 )
                sign = 1;
            else
                sign = -1;

            S = S + sign * (  pow(a, l) * exp(-3.0 * PI * I * (Double)l/4.0) * pow(2.0 * PI / b, l/2.0)/(factorial(l) * factorial( (j - l)/2 ) ) );
        }
    }

    S = S * z;

    return S;
}
*/

Complex IC9E(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi(a - ia + 2bK + 2ibK)t - 4 pi b t^2) dt
    // for a and b positive, and 2bK >= 1
    //
    
    
    int L = to_int(ceil(-LOG(epsilon)/(2 * PI *(a + 2 * b * K)) ));
    L = max(0, L);
    
    Complex S = (Complex)0;
    Complex c = a - I * a + (Double)2 * b * (Double)K + (Double)2 * I * b * (Double)K;
    Double K_to_the_j = K_power(j, cache);
    for(Double n = (Double)0; n <= L-1; n = n + 1) {
        
        Complex z2 = EXP(-2 * PI * (n * c + (Double)2 * b * n * n));
        Complex z = G_I( I*(c + (Double)4 * b * n), 2 * b, n, j, epsilon * K_to_the_j/((Double)L * abs(z2)));
        z = z * z2;
        S = S + z;
    }
    S = S * K_power(-j, cache);
    return S;
}


// IC9H is an inline function in the header file

