#include "theta_sums.h"
#include "precomputed_tables.h"

#include "log.h"

#include <iostream>
#include <cmath>
#include <assert.h>

using namespace std;

Complex IC0(int K, int j, Double a, Double b, Complex C11, Complex C12, mpfr_t mp_a, mpfr_t mp_b, theta_cache * cache, Double epsilon) {
    if(verbose::IC0 >= 1) {
        cout << "IC0 called with:   " << endl;
        cout << "               a = " << a << endl;
        cout << "               b = " << b << endl;
        cout << "               j = " << j << endl;
        cout << "               K = " << K << endl;
        cout << "         epsilon = " << epsilon << endl;
    }
    if(b <= (-LOG(epsilon)) * (-LOG(epsilon))/((Double)K * (Double)K)) {
        if(verbose::IC0) {
            cout << "IC0() using method 1. " << endl;
        }

        if(verbose::IC0 >= 2) {
            cout << "In ICO(), 'really small b' case:" << endl;
        }

        int N = to_int(ceil(std::max ((Double)2.0, sqrt(2.0) * abs(-LOG(epsilon))) ));
        if(verbose::IC0 >= 2) {
            cout << "N = " << N << endl;
        }

        //cout << "K: " << K << endl;
        //cout << "N: " << N << endl;
        //cout << "a: " << a << endl;
        //cout << "b: " << b << endl;
        //cout << "j: " << j << endl;
        //cout << "epsilon: " << epsilon << endl;

        mpfr_t mp_a2, mp_b2, tmp;
        mpfr_init2(mp_a2, mpfr_get_prec(mp_a));
        mpfr_init2(mp_b2, mpfr_get_prec(mp_b));
        mpfr_init2(tmp, mpfr_get_prec(mp_a));

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
        Double b2_mod1 = mpfr_get_d(mp_b2, GMP_RNDN);

        //cout << "a2 mod 1: " << a2_mod1 << endl;
        //cout << "b2 mod 1: " << b2_mod1 << endl;

        //mpfr_clear(mp_a2);
        //mpfr_clear(mp_b2);
        //mpfr_clear(tmp);

        Complex S = (Complex)0;

        Double N_to_the_j = pow(N, j);

        Double exp_linear_term = 0;
        Double new_epsilon = epsilon * N_to_the_j * K_power(-1, cache);

        mpfr_t tmp2;
        mpfr_init2(tmp2, mpfr_get_prec(mp_b2));
        for(int n = 0; n < N; n++) {
            //Complex C = exp((Double)2 * PI * I * (Double)n * (a2_mod1 + b2_mod1 * (Double)n));
            mpfr_mul_si(tmp, mp_b2, n * n, GMP_RNDN);
            mpfr_mul_si(tmp2, mp_a2, n, GMP_RNDN);
            mpfr_add(tmp, tmp, tmp2, GMP_RNDN);
            mpfr_frac(tmp, tmp, GMP_RNDN);
            Double x = mpfr_get_d(tmp, GMP_RNDN);
            //Double x = b2_mod1 * n * n;
            //Double exp_quadratic_term = 2 * PI * x;
            //Complex C = Complex( cos(exp_linear_term + exp_quadratic_term), sin(exp_linear_term + exp_quadratic_term) );
            Complex C = Complex( cos(2 * PI * x), sin(2 * PI * x) );
            //Complex z = G_method1_R(a2 + (Double) 2 * (Double)n * b2, b2, n, j, new_epsilon);
            Complex z = G(a2 + (Double) 2 * (Double)n * b2, b2, n, j, new_epsilon);
            //Complex z = G_method1(a2 + (Double) 2 * (Double)n * b2, b2, n, j, new_epsilon);
            //Complex z2 = G_via_Euler_MacLaurin(a2 + (Double) 2 * (Double)n * b2, b2, n, j, new_epsilon);
            //Complex z = G_via_Euler_MacLaurin(a2 + (Double) 2 * (Double)n * b2, b2, n, j, new_epsilon);
            //Complex z = G_via_Euler_MacLaurin(a2 + (Double) 2 * (Double)n * b2, b2, n, j, epsilon * N_to_the_j/K);

            //S = S + C * z;
            if(verbose::IC0 >= 2) {
//                cout << n << ": " << z << endl;
                cout << n << ": " << K * pow(N, -(j+1)) * C * z << endl;
                cout << "   G returned " << z << endl;
                //cout << "   G_via_EulerMaclaurin returned " << z2 << endl;
                //cout << "   difference was << " << z - z2 << endl;
            }
            //S = S + C * z * (Double)K * pow(N, -(j+1));
            S = S + C * z;
            exp_linear_term += 2 * PI * a2_mod1;
        }
        S *= K * pow(N, -(j + 1));
        if(verbose::IC0 >= 2) {
            cout << "IC0 returning S = " << S << endl;
        }
//        S = S * (Double)K;
//        S = S / (Double)(N);
//        S = S / N_to_the_j;

        return S;
    }
    if(-a/(2 * b) >= 0 && -a/(2 * b) <= K) {

        if(verbose::IC0) {
            cout << "IC0() using method 2. " << endl;
        }

        Complex A = IC8(K, j, mp_a, mp_b, cache);
        Complex B = IC6(K, j, a, b, mp_a, cache, epsilon/4);
        Complex C = IC5(K, j, a, b, cache, epsilon/4);
        //Complex D = IC1(K, j, a, b, C11, C12, cache, epsilon/4);
        Complex D = IC1(K, j, a, b, C11, C12, cache, epsilon);
        if(verbose::IC0) {
            std::cout << "   IC1:" << D << std::endl;
            std::cout << "   IC5:" << C << std::endl;
            std::cout << "   IC6:" << B << std::endl;
            std::cout << "   IC8:" << A << std::endl;
        }
        if(verbose::IC0 >= 2) {
            cout << "IC0 returning S = " << A - B - C - D << endl;
        }
        // Really here we are computing IC2 - IC1...
        return A - B - C - D;
    }
    else {
        if(-a/(2 * b) > K) {
            if(verbose::IC0) {
                cout << "IC0() using method 3. " << endl;
            }
            
            Complex A = IC3(K, j, a, b, cache, epsilon/2);
            Complex B = IC4(K, j, a, b, C11, cache, epsilon/2);
            if(verbose::IC0) {
                std::cout << "   IC3:" << A << std::endl;
                std::cout << "   IC4:" << B << std::endl;
            }
            if(verbose::IC0 >= 2) {
                cout << "IC0 returning S = " << A - B << endl;
            }
            return A - B;
        }
        else {// -a/(2b) < 0
            if(verbose::IC0) {
                cout << "IC0() using method 4. " << endl;
            }

            Complex A = IC3c(K, j, a, b, cache, epsilon/2);
            Complex B = IC4c(K, j, a, b, C11, cache, epsilon/2);
            if(verbose::IC0) {
                std::cout << "   IC3c:" << A << std::endl;
                std::cout << "   IC4c:" << B << std::endl;
            }
            if(verbose::IC0 >= 2) {
                cout << "IC0 returning S = " << A - B << endl;
            }
            return A - B;
        }
    }
}

Complex IC1(int K, int j, Double a, Double b, Complex C11, Complex C12, theta_cache * cache, Double epsilon) {
    //
    // C11 should be passed as I * exp(2 pi I a K + 2 PI i b K^2)
    // C12 should be passed as I * exp(-2 pi(a + 2bK)K - 2 pi i b K^2)/sqrt(2 * PI * b)
    // 
    // (NEED (-log(epsilon)^2/(K^2 * 2 pi)) < b <= 1/8K )
    //

    if(b < (LOG(epsilon) * LOG(epsilon)/(K * K))) {
        cout << "************Warning: b is too small in IC1" << endl;
    }

    C11 = I * cache->ExpABK;
    //cout  << "C11 = " << C11 << endl;

    Complex S = 0;
    for(int r = 0; r <= j; r++) {
        Double z = binomial_coefficient(j, r);
        Complex y = IC7(K, r, a + 2 * b * K, b, cache, epsilon /(2.0 * z * (j + 1)));
        S = S + z * I_power(r) * y;
        if(verbose::IC1) {
            cout << r << ": IC7(" << a + 2 * b * K << ", " << b << ", " << r << ", " << K << ") returned " << y << endl;
        }
    }

    if( a + 2 * b * K <= - LOG(epsilon * sqrt(b))/((Double)2 * PI * K) + 1) {
        Complex C12 = I * cache->ExpBK_inverse * exp(-2 * PI * (a + 2 * b * K) * K);

        /*
        Complex S2 = 0;
        for(int r = 0; r <= j; r++) {
            Double z = binomial_coefficient(j, r) * pow(2 * PI * b, (Double)(-(r+1)) / 2.0) * pow(K, -r);
            //S2 = S2 + z * I_power(r) * G((a + 2 * b * K)/sqrt(2 * PI * b) + I * (Double)2 * sqrt(b) * (Double)K/sqrt(2 * PI), (Double)1/(2 * PI), 0, r, sqrt(2 * PI * b) * epsilon/(abs(C12) * (Double)2 * z * (j + 1)));
            Complex y = G((a + 2 * b * K)/sqrt(2 * PI * b) + I * sqrt(2 * b) * (Double)K/sqrt(PI), (Double)1/(2 * PI), 0, r, epsilon/(abs(C12) * (Double)2 * z * (j + 1)));
            S2 = S2 + z * I_power(r) * y;
            if(verbose::IC1) {
                cout << r << ": G( " << (a + 2 * b * K)/sqrt(2 * PI * b) + I * sqrt(2 * b) * (Double)K/sqrt(PI) << ", " << (Double)1/(2 * PI) << ", 0, " << r << ") returned " << y << endl;
                cout << r << ": z * G * C11 * C12 = " << z * y * C11 * C12 << endl;
            }
        }
        S2 = S2 * C12;
        if(verbose::IC1) {
            cout << "S2 = " << S2 << endl;
        }
        S = S + S2;
        */

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
        if(verbose::IC1) {
            cout << "S2 = " << S2 << endl;
        }
        S = S + S2;
    }

    S *= C11;

    if(verbose::IC1)
        cout << "Returning IC1(" << a << ", " << b << ", " << j << ", " << K << " = " << S << endl;

    return S;
}

Complex IC1c(int K, int j, Double a, Double b, Complex C8, theta_cache * cache, Double epsilon) {
    //
    // Compute C8 * exp(-2 pi a K) int_0^K t^j exp(2 pi i a t - 4 pi b K t + 2 pi i b t^2),
    //
    // ( where we expect that C8 = -i exp(-2 pi i b K^2) )
    //
    // requires that 2bK >= 1

    if(verbose::IC1c) {
        cout << "Inside IC1c:C8 = " << C8 << endl;
        cout << "            a  = " << a << endl;
        cout << "            b  = " << b << endl;
        cout << "            K  = " << K << endl;
        cout << "            j  = " << j << endl;
        cout << "      epsilon  = " << epsilon << endl;
    }
    Complex S = (Complex)0;

    int L = min(K, max(0, to_int(ceil(-LOG(epsilon) - 2 * PI * a * (Double)K/(4 * PI * b * K)) ) ));

    if(verbose::IC1c) {
        cout << "            L = " << L << endl;
    }

    for(int l = 0; l <= j; l++) {
        Complex S1 = 0;
        Double z = K_power(l, cache)/binomial_coefficient(j, l);
        for(int n = 0; n <= L - 1; n++) {
            S1 = S1 + EXP(2.0 * PI * n * (I * a - 2.0 * b * (Double)K + I * b * (Double)n) ) 
                    * G(a + (Double)2.0 * I * b * (Double)K + (Double)2.0 * b * (Double)n, b, n, l, epsilon * exp(4 * PI * b * K * (Double)n + 2 * PI * a * K) * z, 0);
        }
        S1 = S1 * minus_I_power(l)/z;
        if(verbose::IC1c) {
            cout << "Term " << l << ": " << S1 << endl;
        }
        S = S + S1;
    }

    //for(int n = 0; n <= L - 1; n++) {
    //    S = S + EXP(2.0 * PI * n * (I * a - 2.0 * b * (Double)K + I * b * (Double)n) ) * G(a + (Double)2.0 * I * b * (Double)K + (Double)2.0 * b * (Double)n, b, epsilon * exp(4 * PI * b * K * (Double)n + 2 *
    //     PI * a * K) );
    //}

    S *= EXP(-2 * PI * a * K);
    S *= C8;

    if(verbose::IC1c) {
        cout << "Computed IC1c(";
        cout << K << ", ";
        cout << j << ", ";
        cout << a << ", ";
        cout << b << ") = ";
        cout << S << endl;
    }

    return S;
}




// IC2 not defined anywhere

// IC3  defined inline in theta_sums.h
// IC3c defined inline in theta_sums.h


Complex IC4(int K, int j, Double a, Double b, Complex C11, theta_cache * cache, Double epsilon) {
    // needs a + 2bK <= 0
    // b is always positive, so this will be satisfied if 
    
    Complex S = 0;
    for(int r = 0; r <= j; r++) {
        Double z = binomial_coefficient(j, r);
        S = S + z * minus_I_power(r) * IC9H(K, r, -(a + 2 * b * (Double)K), b, cache, epsilon/z);
    }
    return -S * C11;
}




Complex IC4c(int K, int j, Double a, Double b, Complex C11, theta_cache * cache, Double epsilon) {
    // called when a > 0 and b is small but not too small.
    
    Complex S = 0;
    for(int r = 0; r <= j; r++) {
        Double z = binomial_coefficient(j, r);
        S = S + z * I_power(r) * IC9H(K, r, (a + 2 * b * (Double)K), b, cache, epsilon/z);
    }
    return S * C11;
}

Complex IC5(int K, int j, Double a, Double b, theta_cache * cache, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi i exp(i pi /4) a - 2 pi b t^2) dt 
    // assuming that a is negative (??)
    //

    return (Double)minus_one_power(j) * I_power(j + 1) * IC7(-1, j, -a, b, cache, epsilon * K_power(j, cache)) * K_power(-j, cache);
}

Complex IC6(int K, int j, Double a, Double b, mpfr_t mp_a, theta_cache * cache, Double epsilon) {
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

    if(verbose::IC6) {
        cout << "In IC6, using L = " << L << endl;
        cout << "              x = " << x << endl;
        cout << "              y = " << y << endl;
    }

    Complex z1 = 1.0/ExpA(mp_a, K);

    //Complex z = (Double)1/sqrt(2 * PI * b) * exp(I * PI/(Double)4 + 2 * PI*(I - (Double)1) * a * (Double)K - 4 * PI * b * (Double)K * Double(K));
    //Complex z = (Double)1/sqrt(2 * PI * b) * exp(I * PI/(Double)4 - 2.0 * PI* a * (Double)K - 4 * PI * b * (Double)K * Double(K)) * z1;

    //Complex z = pow(2.0 * PI * b, -((Double)j + 1.0)/2) * exp(I * PI * ((Double)j + 1)/(Double)4 - 2.0 * PI* a * (Double)K - 4 * PI * b * (Double)K * Double(K)) * z1;
    Complex z = exp(I * PI * ((Double)j + 1)/(Double)4 - 2.0 * PI* a * (Double)K - 4 * PI * b * (Double)K * Double(K)) * z1;

    Complex alpha = (I - (Double)1)*a/(2 * sqrt(PI * b)) - 2 * sqrt(b) * (Double)K / sqrt(PI);
    Double beta = -1/(2 * PI);

    Complex S = (Complex)0;
    
    for(int r = 0; r <= j; r++) {
     
        Complex S1 = 0;

        Double z2 = K_power(r, cache) * pow(2.0 * PI * b, ((Double)r + 1)/2) / ( pow(2, ((Double)(j-r))/2) * binomial_coefficient(j, r));

        for(int n = 0; n < L; n++) {
            Complex w1 = exp(2 * PI * alpha * (Double)n + 2 * PI * beta * (Double)n * (Double)n);
            Complex w2 = G(-I * alpha - (Double)2 * I * beta * (Double)n, -I * beta, n, r, epsilon * z2/(abs(z) * abs(w1)));
            if(verbose::IC6) {
                cout << "For n = " << n << ", z * w1 * w2 = " << z * w1 * w2 << endl;
            }
            S1 = S1 + z * w1 * w2;
        }
        S1 = S1/z2;
        S = S + S1;
    }
    S = S;

    return S;
}

Complex IC7(int K, int j, Double a, Double b, theta_cache * cache, Double epsilon) {
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


    stats::IC7++;

    if(verbose::IC7) {
        cout << "Entering IC7() with " << endl;
        cout << "         K = " << K << endl;
        cout << "         j = " << j << endl;
        cout << "         a = " << a << endl;
        cout << "         b = " << b << endl;
        cout << "   epsilon = " << epsilon << endl;
    }

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

Complex IC7_method1(int K, int j, Double a, Double b, theta_cache * cache, Double epsilon, int L) {
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
        if(verbose::IC7 >= 2) {
            cout << "Computing IC7() = " << C9 << " * " << root_2pi_b_power(-(j + 1), cache) << " S " << endl;
            cout << "S == 0" << endl;
        }
        {
            Double one_over_L = 1.0/L;
            Complex alpha = C10 * a * root_2pi_b_power(-1, cache);
            Complex s1 = alpha;
            Complex y1 = exp(2 * PI * I * s1);
            Complex y = 1.0;
            Double newepsilon = epsilon * two_pi_b_power * one_over_L;
            for(int n = 0; n <= L-1; n = n + 1) {
                //Complex z3 = exp(2 * PI * I * n * s1- n * n);
                Double y3 = exp(-n * n);
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
                
                if(verbose::IC7 >= 2) {

                    cout << "S += " << z << "; (S == " << S << ")" << endl;
                }


            }
        }

        S = S * C9 * root_2pi_b_power(-(j + 1), cache);

    }
    else {
        K_to_the_j = K_power(j, cache);
        //Double two_pi_b_power = pow(2 * PI * b, ((Double)j + 1.0)/2);

        if(verbose::IC7 >= 2) {
            cout << "Computing IC7() = " << C9 << " * " << root_2pi_b_power(-(j + 1), cache) << " S " << endl;
            cout << "S == 0" << endl;
        }


        Double two_pi_b_power = root_2pi_b_power(j + 1, cache);
        {
            Double one_over_L = 1.0/L;
            Complex alpha = C10 * a * root_2pi_b_power(-1, cache);
            Complex s1 = alpha;
            Complex y1 = exp(2 * PI * I * s1);
            Complex y = 1.0;
            Double newepsilon = epsilon * two_pi_b_power * K_to_the_j * one_over_L;
 
            for(int n = 0; n <= L-1; n = n + 1) {
                //Complex z3 = exp(2 * PI * I * n * s1- n * n);
                Double y3 = exp(-n * n);


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
                if(verbose::IC7 >= 2) {
                    cout << "S += " << z << "; (S == " << S << ")" << endl;
                }

            }
        }


        S = S * C9 * K_power(-j, cache) * root_2pi_b_power(-(j + 1), cache);
    }
    if(verbose::IC7)
        cout << "IC7() returning " << S << endl;

    return S;




}

static struct {
    int a_per_unit_interval;
    int number_of_a;
    int max_j;
    Complex ** values;
    double a_spacing;
} IC7_cache;

static bool IC7_cache_initialized = false;

void build_IC7_cache(int a_per_unit_interval, Double max_a, int max_j, Double epsilon) {
    int number_of_a = max_a * a_per_unit_interval + 1;
    
    IC7_cache.number_of_a = number_of_a;
    IC7_cache.a_per_unit_interval = a_per_unit_interval;
    IC7_cache.max_j = max_j;
    IC7_cache.a_spacing = 1.0/a_per_unit_interval;


    Double memory_used = number_of_a * (max_j + 1) * 16;


    clock_t start_time = clock();

    cout << "Building IC7 cache." << endl;
    cout << "    Total memory used for this cache will be " << memory_used/1000000.0 << " million bytes." << endl;

    if(FAKE_PRECOMPUTATION) {
        IC7_cache_initialized = true;
        return;
    }





    IC7_cache.values = new Complex * [number_of_a];
    Double a = 0;
    Double a_increment = 1.0/a_per_unit_interval;


    Double percent_done = 0.0;
    Double old_percent_done = 0.0;

    for(int n = 0; n < number_of_a; n++) {
        percent_done = (Double)n/(number_of_a);
        if(old_percent_done + .005 < percent_done) {
            cout <<  "    " << percent_done * 100 << " percent done with IC7." << endl;
            old_percent_done = percent_done;
        }

        IC7_cache.values[n] = new Complex[max_j + 1];
        for(int j = 0; j <= max_j; j++) {
            IC7_cache.values[n][j] = IC7star(a, j, epsilon, false);
        }
        a += a_increment;
    }

    clock_t end_time = clock();

    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Total time to build IC7 cache: " << elapsed_time << " seconds." << endl;


    IC7_cache_initialized = true;
}

void free_IC7_cache() {
    IC7_cache_initialized = false;
}

inline Complex get_cached_IC7star_value(int a_index, int j) {
    if(a_index < IC7_cache.number_of_a && j <= IC7_cache.max_j) {
        if(FAKE_PRECOMPUTATION) {
            return 0.0;
        }
        else {
            return IC7_cache.values[a_index][j];
        }
    }
    else {
        cout << "Warning. IC7 cache called with values out of range." << endl;
        cout << "a_index = " << a_index << ", j = " << j << endl;
        return 0.0/0.0;
    }
}

Complex IC7star(Double a, int j, Double epsilon, bool use_cache) {
    if(verbose::IC7star) {
        cout << "Entering IC7star() with " << endl;
        cout << "         j = " << j << endl;
        cout << "         a = " << a << endl;
        cout << "   epsilon = " << epsilon << endl;
    }

    if(use_cache && IC7_cache_initialized) {
        int a0_index = round(IC7_cache.a_per_unit_interval * a);
        Double a0 = (Double)a0_index / (Double)IC7_cache.a_per_unit_interval;

        if(a0_index < IC7_cache.number_of_a) {
            if(stats::stats) {
                stats::IC7_taylor_expansion++;
            }

            int l = 0;
            Double error = 2 * epsilon;

            Double a_minus_a0_power = 1.0;
            
            Complex S = 0;
            int number_of_terms = 0;
            while(error > epsilon) {
                Double z = minus_one_power(l) * two_pi_over_factorial_power(l) * a_minus_a0_power;
                error = abs(z) * .5 * gamma_s_over_2(j + l + 1);

                //S = S + z * exp_minus_i_pi4(l) * IC7star(a0, j + l, epsilon, false);
                S = S + z * exp_minus_i_pi4(l) * get_cached_IC7star_value(a0_index, j + l);

                l++;
                a_minus_a0_power *= (a - a0);
                if(stats::stats) {
                    number_of_terms++;
                }
            }

            if(stats::stats) {
                stats::IC7_terms_used += number_of_terms;
            }

            //cout << l << " ";

            return S;
        }
    }

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

    if(verbose::IC7star) {
        cout << "In IC7star(): L = " << L << endl;
    }

    if(L == 0) {
        return 0.0;
    }

    Complex S = (Complex)0;

    {
        Double one_over_L = 1.0/L;
        Complex alpha = C10 * a;
        Complex s1 = alpha;
        Complex y1 = exp(2 * PI * I * s1);
        Complex y = 1.0;
        Double newepsilon = epsilon * one_over_L;
        for(int n = 0; n <= L-1; n = n + 1) {
            //Complex z3 = exp(2 * PI * I * n * s1- n * n);
            Double y3 = exp(-n * n);
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


Complex IC9E(int K, int j, Double a, Double b, theta_cache * cache, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi(a - ia + 2bK + 2ibK)t - 4 pi b t^2) dt
    // for a and b positive, and 2bK >= 1
    //
    
    if(verbose::IC9E) {
        cout << "IC9E() called with: " << endl;
        cout << "               K = " << K << endl;
        cout << "               j = " << j << endl;
        cout << "               a = " << a << endl;
        cout << "               b = " << b << endl;
        cout << "         epsilon = " << epsilon << endl;
    }
    
    int L = to_int(ceil(-LOG(epsilon)/(2 * PI *(a + 2 * b * K)) ));
    L = max(0, L);
    
    Complex S = (Complex)0;
    Complex c = a - I * a + (Double)2 * b * (Double)K + (Double)2 * I * b * (Double)K;
    Double K_to_the_j = K_power(j, cache);
    for(Double n = (Double)0; n <= L-1; n = n + 1) {
        
        Complex z2 = exp(-2 * PI * (n * c + (Double)2 * b * n * n));
        Complex z = G( I*(c + (Double)4 * b * n), (Double)2 * I * b, n, j, epsilon * K_to_the_j/((Double)L * abs(z2)));
        z = z * z2;
        S = S + z;
    }
    S = S * K_power(-j, cache);
    if(verbose::IC9E) {
        cout << "IC9E returning: " << S << endl;
    }
    return S;
}


// IC9H is an inline function in the header file

