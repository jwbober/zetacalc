#include "theta_sums.h"

#include <iostream>
#include <cmath>

using namespace std;

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
    //cout << "------------------------------------------------------------------------------------" << endl;
    //cout << z << endl;
    //cout << "------------------------------------------------------------------------------------" << endl;
    Complex S = z;
 
    if( a + 2 * b * K <= - LOG(epsilon * sqrt(b))/((Double)2 * PI * K) + 1) {
        z = G((a + 2 * b * K)/sqrt(2 * PI * b) + I * (Double)2 * sqrt(b) * (Double)K/sqrt(2 * PI), (Double)1/(2 * PI), sqrt(2 * PI * b) * epsilon/(abs(C12) * (Double)2));
    //cout << "------------------------------------------------------------------------------------" << endl;
    //cout << C12 << endl;
    //cout << z << endl;
    //cout << "------------------------------------------------------------------------------------" << endl;

        S = S + C12 * z;
    }


    //cout << "------------------------------------------------------------------------------------" << endl;
    //cout << C12 << endl;
    //cout << z << endl;
    //cout << z * C12 << endl;
    //cout << C11 << endl;
    //cout << "------------------------------------------------------------------------------------" << endl;
    S *= C11;

    return S;
}

Complex IC1c(int K, Double a, Double b, Complex C8, Double epsilon) {
    //
    // Compute C8 * exp(-2 pi a K) int_0^K exp(2 pi i a t - 4 pi b K t + 2 pi i b t^2),
    //
    // ( where we expect that C8 = -i exp(-2 pi i b K^2) )
    //

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

    //cout << L << endl;

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

Complex IC5(Double a, Double b, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi i exp(i pi /4) a - 2 pi b t^2) dt 
    // assuming that a is positive
    //

    return I * IC7(-1, -a, b, epsilon);
}

//Complex IC8(Double a, Double b, Double epsilon) {
//    return exp(I * PI * (.25 - a * a/(2.0 * b)))/sqrt((Double)2 * b);
//}

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
        //Complex w1 = exp(2 * PI * (Double)n *  (alpha + 2 * beta + beta * (Double)n));
        Complex w1 = exp(2 * PI * alpha * (Double)n + 2 * PI * beta * (Double)n * (Double)n);
        Complex w2 = G(-I * alpha - (Double)2 * I * beta * (Double)n, -I * beta , epsilon/(abs(z) * abs(w1)));
        if(verbose::IC6) {
            cout << "For n = " << n << ", z * w1 * w2 = " << z * w1 * w2 << endl;
        }
        S = S + z * w1 * w2;
    }

    //S = S * z;

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
        L = to_int(ceil( -LOG(epsilon*sqrt(b))));
    }

    L = max(0, L);
    
    if(verbose::IC7) {
        cout << "In IC7(): L = " << L << endl;
    }

    //check_condition(L <= K, "Warning: In IC7, K was too small.");

    if(K != -1) {
        //L = min(L, (int)(sqrt(2.0 * b)*K));
        L = min(L, to_int(sqrt((Double)4 * PI * b) * K));
    }

    Complex S = (Complex)0;

    Complex x2 = (Double)2 * PI * Complex(-x, x);
    Complex x3 = Complex(x, x);
    for(Double n = (Double)0; n <= L-1; n = n + 1) {
        //Complex z2 = exp(x2*n/(sqrt((Double(2) * PI)) - n * n));
        Complex z2 = exp( 2 * PI * I * C10 * a / sqrt((Double)2 * PI * b) * n - n * n);
        //Complex z = G(x3 + (Double)2 * I * n, I, epsilon/(abs(z2) * Double(L)));
        Complex z = G( C10 * a / sqrt((Double)2 * PI * b) + I * n/PI, I/( (Double)2 * PI), epsilon * sqrt((Double)2 * PI * b)/(abs(z2) * Double(L)));
        //Complex z = G( (x3/sqrt( (Double)2 * PI)) + I * n/(Double)PI, I/((Double)2 * PI), epsilon/(abs(z2) * Double(L)));
        z = z * z2;
        S = S + z;
    }
    
    S = S * C9/sqrt((Double)2 * PI * b);

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

/*
Complex IC9H(Double a, Double b, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi a t - 2 pi i b t^2)
    //
    // after a change of contour, this is well approximated by IC7(K, a, b) with large K
    //

    //int K = (int)(10 * -log(epsilon * sqrt(b)));
    //if(K < 0) {
    //    K = 10;
    //}
    
    int K = (int)(10 * ceil( max(-log(epsilon * sqrt(b)), 1.0)/sqrt(b) ));
    return IC7(K, a, b, epsilon);
}
*/
