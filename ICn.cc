#include "theta_sums.h"
#include "precomputed_tables.h"

#include <iostream>
#include <cmath>

using namespace std;


Complex IC0(int K, Double a, Double b, Complex C11, Complex C12, mpfr_t mp_a, mpfr_t mp_b, Double epsilon) {
    if(b <= (-LOG(epsilon)) * (-LOG(epsilon))/((Double)K * (Double)K)) {

        int N = to_int(ceil(std::max ((Double)1.0, -LOG(epsilon)) ));
        mpfr_t mp_a2, mp_b2, tmp;
        mpfr_init2(mp_a2, mpfr_get_prec(mp_a));
        mpfr_init2(mp_b2, mpfr_get_prec(mp_b));
        mpfr_init2(tmp, mpfr_get_prec(mp_a));

        mpfr_mul_si(mp_a2, mp_a, K, GMP_RNDN);
        mpfr_div_si(mp_a2, mp_a2, N, GMP_RNDN);

        Double a2 = mpfr_get_d(mp_a2, GMP_RNDN);
        
        mpfr_floor(tmp, mp_a2);
        mpfr_sub(mp_a2, mp_a2, tmp, GMP_RNDN);

        mpfr_mul_si(mp_b2, mp_b, K, GMP_RNDN);
        mpfr_mul_si(mp_b2, mp_b2, K, GMP_RNDN);
        mpfr_div_si(mp_b2, mp_b2, N * N, GMP_RNDN);

        Double b2 = mpfr_get_d(mp_b2, GMP_RNDN);

        mpfr_floor(tmp, mp_b2);
        mpfr_sub(mp_b2, mp_b2, tmp, GMP_RNDN);

        Double a2_mod1 = mpfr_get_d(mp_a2, GMP_RNDN);
        Double b2_mod1 = mpfr_get_d(mp_b2, GMP_RNDN);

        mpfr_clear(mp_a2);
        mpfr_clear(mp_b2);
        mpfr_clear(tmp);

        Complex S = (Complex)0;

        for(int n = 0; n < N; n++) {
            Complex C = exp((Double)2 * PI * I * (Double)n * (a2_mod1 + b2_mod1 * (Double)n));
            Complex z = G(a2 + (Double) 2 * (Double)n * b2, b2, (epsilon)/(abs(C) * K));
            
            S = S + C * z;
        }
        S = S * (Double)K;
        S = S / (Double)(N);

        return S;
    }
    if(-a/(2 * b) >= 0 && -a/(2 * b) <= K) {
        Complex A = IC8(a, b, epsilon/4);
        Complex B = IC6(K, a, b, epsilon/4);
        Complex C = IC5(a, b, epsilon/4);
        Complex D = IC1(K, a, b, C11, C12, epsilon/4);
        if(verbose::IC0) {
            std::cout << "   IC1:" << D << std::endl;
            std::cout << "   IC5:" << C << std::endl;
            std::cout << "   IC6:" << B << std::endl;
            std::cout << "   IC8:" << A << std::endl;
        }
        // Really here we are computing IC2 - IC1...
        return A - B - C - D;
    }
    else {
        if(-a/(2 * b) > K) {
            Complex A = IC3(a, b, epsilon/2);
            Complex B = IC4(K, a, b, C11, epsilon/2);
            if(verbose::IC0) {
                std::cout << "   IC3:" << A << std::endl;
                std::cout << "   IC4:" << B << std::endl;
            }
            return A - B;
        }
        else {// -a/(2b) < 0
            Complex A = IC3c(a, b, epsilon/2);
            Complex B = IC4c(K, a, b, C11, epsilon/2);
            if(verbose::IC0) {
                std::cout << "   IC3c:" << A << std::endl;
                std::cout << "   IC4c:" << B << std::endl;
            }
            return A - B;
        }
    }
}


Complex IC1(int K, Double a, Double b, Complex C11, Complex C12, Double epsilon) {
    //
    // C11 should be passed as I * exp(2 pi I a K + 2 PI i b K^2)
    // C12 should be passed as I * exp(-2 pi(a + 2bK)K - 2 pi i b K^2)/sqrt(2 * PI * b)
    // 
    // (NEED (-log(epsilon)^2/(K^2 * 2 pi)) < b <= 1/8K )
    //

    if(b < (LOG(epsilon) * LOG(epsilon)/(K * K))) {
        cout << "************Warning: b is too small in IC1" << endl;
    }

    Complex z = IC7(K, a + 2 * b * K, b, epsilon/2);
    Complex S = z;
 
    if( a + 2 * b * K <= - LOG(epsilon * sqrt(b))/((Double)2 * PI * K) + 1) {
        z = G((a + 2 * b * K)/sqrt(2 * PI * b) + I * (Double)2 * sqrt(b) * (Double)K/sqrt(2 * PI), (Double)1/(2 * PI), sqrt(2 * PI * b) * epsilon/(abs(C12) * (Double)2));
        S = S + C12 * z;
    }

    S *= C11;

    return S;
}

Complex IC1c(int K, Double a, Double b, Complex C8, Double epsilon) {
    //
    // Compute C8 * exp(-2 pi a K) int_0^K exp(2 pi i a t - 4 pi b K t + 2 pi i b t^2),
    //
    // ( where we expect that C8 = -i exp(-2 pi i b K^2) )
    //
    //
    // Requires the 2bK >= 1

    if(verbose::IC1c) {
        cout << "Inside IC1c:C8 = " << C8 << endl;
        cout << "            a  = " << a << endl;
        cout << "            b  = " << b << endl;
        cout << "            K  = " << K << endl;
        cout << "      epsilon  = " << epsilon << endl;
    }
    Complex S = (Complex)0;

    int L = min(K, max(0, to_int(ceil(-LOG(epsilon) - 2 * PI * a * (Double)K/(4 * PI * b * K)) ) ));

    if(verbose::IC1c) {
        cout << "            L = " << L << endl;
    }

    for(int n = 0; n <= L - 1; n++) {
        S = S + EXP(2.0 * PI * n * (I * a - 2.0 * b * (Double)K + I * b * (Double)n) ) * G(a + (Double)2.0 * I * b * (Double)K + (Double)2.0 * b * (Double)n, b, epsilon * exp(4 * PI * b * K * (Double)n + 2 *
         PI * a * K) );
    }

    S *= EXP(-2 * PI * a * K);
    S *= C8;

    if(verbose::IC1c) {
        cout << "Computed IC1c(";
        cout << K << ", ";
        cout << a << ", ";
        cout << b << ") = ";
        cout << S << endl;
    }

    return S;
}

Complex IC1c(int K, int j, Double a, Double b, Complex C8, Double epsilon) {
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
        cout << "      epsilon  = " << epsilon << endl;
    }
    Complex S = (Complex)0;

    int L = min(K, max(0, to_int(ceil(-LOG(epsilon) - 2 * PI * a * (Double)K/(4 * PI * b * K)) ) ));

    if(verbose::IC1c) {
        cout << "            L = " << L << endl;
    }

    Complex minus_i_power(0, 1);
    for(int l = 0; l <= j; l++) {
        minus_i_power *= -I;
        Complex S1 = 0;
        Double z = pow(K, l)/binomial_coefficient(j, l);
        for(int n = 0; n <= L - 1; n++) {
            S1 = S1 + EXP(2.0 * PI * n * (I * a - 2.0 * b * (Double)K + I * b * (Double)n) ) 
                    * G(a + (Double)2.0 * I * b * (Double)K + (Double)2.0 * b * (Double)n, b, n, l, epsilon * exp(4 * PI * b * K * (Double)n + 2 * PI * a * K) * z);
        }
        S1 = S1 * minus_i_power/z;
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
        cout << a << ", ";
        cout << b << ") = ";
        cout << S << endl;
    }

    return S;
}




// IC2 not defined anywhere

// IC3  defined inline in theta_sums.h
// IC3c defined inline in theta_sums.h
// IC4  defined inline in theta_sums.h
// IC4c defined inline in theta_sums.h

Complex IC5(Double a, Double b, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi i exp(i pi /4) a - 2 pi b t^2) dt 
    // assuming that a is negative (??)
    //

    return I * IC7(-1, -a, b, epsilon);
}

Complex IC5(int K, int j, Double a, Double b, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi i exp(i pi /4) a - 2 pi b t^2) dt 
    // assuming that a is negative (??)
    //

    return pow(-1, j) * pow(I, j + 1) * IC7(-1, j, -a, b, epsilon * pow(K, j)) * pow(K, -j);
}


Complex IC6(int K, Double a, Double b, Double epsilon) {
    //
    // After some changes of variables this integral is equal to
    //
    // (1/sqrt(2 pi b)) exp( i pi/4 + 2 pi (i - 1)a K - 4 PI b K^2) int_0^L exp(pi (i - 1) a t/sqrt(pi b) - 4 sqrt(pi b) K t - t^2
    //
    // 

    

    Double x = sqrt(PI/b) * a + 4 * sqrt(PI * b) * (Double)K;
    Double y = max(-LOG(epsilon * sqrt(2 * PI) * b) - (Double)2 * PI * (Double)K*(a + (Double)2 * b * (Double)K), (Double)0);

    int L = to_int(ceil(min(sqrt(y), y/x)));

    if(verbose::IC6) {
        cout << "In IC6, using L = " << L << endl;
        cout << "              x = " << x << endl;
        cout << "              y = " << y << endl;
    }

    Complex z = (Double)1/sqrt(2 * PI * b) * exp(I * PI/(Double)4 + 2 * PI*(I - (Double)1) * a * (Double)K - 4 * PI * b * (Double)K * Double(K));

    Complex alpha = (I - (Double)1)*a/(2 * sqrt(PI * b)) - 2 * sqrt(b) * (Double)K / sqrt(PI);
    Double beta = -1/(2 * PI);

    Complex S = (Complex)0;
    
    for(int n = 0; n < L; n++) {
        Complex w1 = exp(2 * PI * alpha * (Double)n + 2 * PI * beta * (Double)n * (Double)n);
        Complex w2 = G(-I * alpha - (Double)2 * I * beta * (Double)n, -I * beta , epsilon/(abs(z) * abs(w1)));
        if(verbose::IC6) {
            cout << "For n = " << n << ", z * w1 * w2 = " << z * w1 * w2 << endl;
        }
        S = S + z * w1 * w2;
    }

    return S;
}

Complex IC7(int K, Double a, Double b, Double epsilon) {
    //
    // We compute C9 * int_0^{sqrt(2) K} exp(-sqrt(2) PI a t + sqrt(2) PI i a t - 2 PI b t^2) dt
    // where C9 = exp(-I pi/4)
    //
    // K = -1 corresponds to infinity
    //
    // K, b, and epsilon should satisfy:
    //   
    //  (1) K > 2 * (-log(epsilon)/2pi)^2
    //  (2) 2bK >= 1
    //

    Complex C9(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
    Complex C10(sqrt(2.0)/2.0, sqrt(2.0)/2.0);

    Double x = a/(sqrt(2 * b * 4 * PI));

    //
    // Note that conditions (1) and (2) above should ensure that we
    // can make a change of variables and then truncate this
    // integral at an integer.
    //
    
    int L = 0;
    if(x > 1) {
        L = to_int(ceil( -LOG(epsilon*sqrt(b))/(2.0 * PI * x)));
    }
    else {
        L = to_int(ceil( -LOG(epsilon*sqrt(b))/(2.0 * PI)));
    }

    L = max(0, L);
    
    if(verbose::IC7) {
        cout << "In IC7(): L = " << L << endl;
    }

    //check_condition(L <= K, "Warning: In IC7, K was too small.");

    if(K != -1) {
        L = min(L, to_int(sqrt((Double)4 * PI * b) * K));
    }

    Complex S = (Complex)0;

    Complex x2 = (Double)2 * PI * Complex(-x, x);
    Complex x3 = Complex(x, x);
    for(Double n = (Double)0; n <= L-1; n = n + 1) {
        Complex z2 = exp( 2 * PI * I * C10 * a / sqrt((Double)2 * PI * b) * n - n * n);
        Complex z = G( C10 * a / sqrt((Double)2 * PI * b) + I * n/PI, I/( (Double)2 * PI), epsilon * sqrt((Double)2 * PI * b)/(abs(z2) * Double(L)));
        z = z * z2;
        S = S + z;
    }
    
    S = S * C9/sqrt((Double)2 * PI * b);

    return S;
}

Complex IC7(int K, int j, Double a, Double b, Double epsilon) {
    //
    // We compute C9 * int_0^{sqrt(2) K} exp(-sqrt(2) PI a t + sqrt(2) PI i a t - 2 PI b t^2) dt
    // where C9 = exp(-I pi/4)
    //
    // K = -1 corresponds to infinity
    //
    // K, b, and epsilon should satisfy:
    //   
    //  (1) K > 2 * (-log(epsilon)/2pi)^2
    //  (2) 2bK >= 1
    //

    Complex C9;
    switch(j % 8) {
        case 0:
            C9 = Complex(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 1:
            C9 = Complex(0, -1);
            break;
        case 2:
            C9 = Complex(-sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 3:
            C9 = Complex(-1, 0);
            break;
        case 4:
            C9 = Complex(-sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 5:
            C9 = Complex(0, 1);
            break;
        case 6:
            C9 = Complex(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 7:
            C9 = Complex(1, 0);
            break;
    }

    Complex C10(sqrt(2.0)/2.0, sqrt(2.0)/2.0);

    Double x = a/(sqrt(2 * b * 4 * PI));

    //
    // Note that conditions (1) and (2) above should ensure that we
    // can make a change of variables and then truncate this
    // integral at an integer.
    //
    
    int L = 0;
    if(x > 1) {
        L = to_int(ceil( -LOG(epsilon*sqrt(b))/(2.0 * PI * x)));
        if(K == -1)
            L = max((int)ceil(((Double)j + 1.0) * pow(log(j + (Double)1), 2) / x), L);
    }
    else {
        L = to_int(ceil( -LOG(epsilon*sqrt(b))/(2.0 * PI)));
        if(K == -1)
            L = max((int)ceil(((Double)j + 1.0) * pow(log(j + (Double)1), 2)), L);
    }
    if(K != -1)
        L = max(0, L);

    cout << L << endl;

    if(verbose::IC7) {
        cout << "In IC7(): L = " << L << endl;
    }

    //check_condition(L <= K, "Warning: In IC7, K was too small.");

    //cout << L << endl;

    if(K != -1) {
        L = min(L, to_int(sqrt((Double)4 * PI * b) * K));
    }

    Complex S = (Complex)0;

    Complex x2 = (Double)2 * PI * Complex(-x, x);
    Complex x3 = Complex(x, x);

    Double K_to_the_j;
    if(K == -1) {
        K_to_the_j = 1.0;       // This is potentially quite confusing. We treat K = -1 as infinity,
                                // so we don't want to do any normalization.
    }
    else {
        K_to_the_j = pow(K, j);
    }

    Double two_pi_b_power = pow(2 * PI * b, ((Double)j + 1.0)/2);

    for(Double n = (Double)0; n <= L-1; n = n + 1) {
        Complex z2 = exp( 2 * PI * I * C10 * a / sqrt((Double)2 * PI * b) * n - n * n);
        Complex z = G( C10 * a / sqrt((Double)2 * PI * b) + I * n/PI, I/( (Double)2 * PI), n, j, epsilon * two_pi_b_power * K_to_the_j/(abs(z2) * Double(L)));
        z = z * z2;
        S = S + z;
    }
    
    S = S * C9/(two_pi_b_power * K_to_the_j);

    return S;
}









// IC8 is an inline function in the header file

Complex IC9E(int K, Double a, Double b, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi(a - ia + 2bK + 2ibK)t - 4 pi b t^2) dt
    // for a and b positive, and 2bK >= 1
    //
    
    int L = to_int(ceil(-LOG(epsilon)/(2 * PI *(a + 2 * b * K)) ));
    L = max(0, L);
    
    Complex S = (Complex)0;
    Complex c = a - I * a + (Double)2 * b * (Double)K + (Double)2 * I * b * (Double)K;
    for(Double n = (Double)0; n <= L-1; n = n + 1) {
        
        Complex z2 = exp(-2 * PI * (n * c + (Double)2 * b * n * n));
        Complex z = G( I*(c + (Double)4 * b * n), (Double)2 * I * b, epsilon/((Double)L * abs(z2)));
        z = z * z2;
        S = S + z;
    }
    return S;
}


Complex IC9E(int K, int j, Double a, Double b, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi(a - ia + 2bK + 2ibK)t - 4 pi b t^2) dt
    // for a and b positive, and 2bK >= 1
    //
    
    int L = to_int(ceil(-LOG(epsilon)/(2 * PI *(a + 2 * b * K)) ));
    L = max(0, L);
    
    Complex S = (Complex)0;
    Complex c = a - I * a + (Double)2 * b * (Double)K + (Double)2 * I * b * (Double)K;
    Double K_to_the_j = pow(K, j);
    for(Double n = (Double)0; n <= L-1; n = n + 1) {
        
        Complex z2 = exp(-2 * PI * (n * c + (Double)2 * b * n * n));
        Complex z = G( I*(c + (Double)4 * b * n), (Double)2 * I * b, n, j, epsilon * K_to_the_j/((Double)L * abs(z2)));
        z = z * z2;
        S = S + z;
    }
    S = S/K_to_the_j;
    return S;
}



// IC9H is an inline function in the header file

