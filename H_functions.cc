#include "theta_sums.h"
#include <iostream>

#include "precomputed_tables.h"

using namespace std;

inline Complex H_method2(int j, Complex alpha, Double epsilon);
Complex H_method5(int j, Complex alpha, Double epsilon);

Complex H(int j, Complex alpha, Double epsilon) {
    //
    // Compute the integral int_0^1 t^j exp(-2 pi alpha t) dt
    // We have three different methods to compute this depending
    // on the range of input.
    //
    
    if(epsilon > 1) {
        return 0.0;
    }

    if(stats::stats) {
        const Double D = 50.0;
        if( abs(alpha) > D ) {
        //    cout << alpha << endl;
            stats::H_function_big++;
        }
        else {
            stats::H_function_small++;
        }
    }

    const Double alpha_0 = 1/(2*PI);

    if(verbose::H) {
        cout << "Function H() called with: " << endl;
        cout << "                      j = " << j << endl;
        cout << "                  alpha = " << alpha << endl;
        cout << "                epsilon = " << epsilon << endl;
    }

    if(abs(alpha) < alpha_0) {
        if(verbose::H) {
            cout << "   In function H(), using method 2" << endl;
        }
        return H_method2(j, alpha, epsilon);
    }
    else if(abs(2 * PI * alpha) > 2 * j) {
        if(verbose::H) {
            cout << "   In function H(), using method 1" << endl;
        }
        return H_method1(j,  alpha);
    }
    else {
        if(verbose::H) {
            cout << "   In function H(), using method 3" << endl;
        }
        return H_method4(j, alpha, epsilon);
    }
   
}

Complex H_method1(int j, Complex alpha) {
    stats::H_method1++;
    // In this case we compute an "exact" value using the antiderivative of the integrand.
    // 
    Complex S = (Complex)0.0;
    Double j_factorial = factorial(j);
    Complex alpha_power = (Complex)1;
    for(int v = 0; v < j + 1; v++) {
        if(v > 0) {
    //        v_factorial *= v;
            alpha_power *= alpha;
        }
        Complex z = alpha_power * two_pi_over_factorial_power(v);  //two_pi_alpha_power/factorial(v);
        S = S + z;
    }
    S = S * EXP(-2 * PI * alpha);
    S = (Double)1.0 - S;
    alpha_power *= alpha;
    S = S * j_factorial/(alpha_power * two_pi_power(j+1));

    if(verbose::H) {
        cout << "Computed H_method1(" << j << ", " << alpha << ") = " << S << endl;
    }

    return S;
}

inline Complex H_method2(int j, Complex alpha, Double epsilon) {
    // When alpha is very small, using the antiderivative can lead to bad
    // cancellation and loss of precision, so we use the taylor series
    // expansion for exp(2 pi alpha t) instead, and integrate term by term
    // using antiderivatives.

    // The taylor expansion gives
    //
    // H(j, alpha) = sum_{m=0}^\infty (-2 pi alpha)^m / (m! (j + m + 1) )
    //
    //int N = max(ceil(-log(epsilon)), ceil(2 * PI * alpha_0 * E * E));; // the number of terms we will use in the taylor expansion

    stats::H_method2++;


    Complex S = (Complex)0;
    
    Double error = epsilon + 1.0;
    int m = 0;
    Complex alpha_power = (Complex)1;
    while(error > epsilon/2.0) {
        if(m > 0) {
            alpha_power *= -alpha;
        }
        Complex z = two_pi_over_factorial_power(m) * alpha_power / Double(j + m + 1);
        S = S + z;
        error = abs(z);
        m++;
    }

    return S;
}

Complex H_method3(int j, Complex alpha, Double epsilon) {
    //
    // We use a change of variables to rewrite the integral, and compute it as
    //
    // (1/M^{j+1}) sum_{m=0}^{M-1} exp( 2 pi i alpha m / M ) sum_{r=0}^j (j choose r) m^r H(j-r, alpha/M)
    //
    // Here, in the inner sum, alpha/M will always be small enough that H(j-r, alpha/M) is computed
    // with method 2.
    //
    stats::H_method3++;

    int M = to_int(ceil(abs(2 * PI * alpha)));

    Complex S = (Complex)0;

    Complex H_table[j + 1];
    for(int m = 0; m <= j; m++) {
        H_table[m] = H_method2(m, alpha/(Double)M, j * epsilon);
    }


    Complex exponential_multiplication_factor = EXP(-2.0 * PI * alpha/ (Double)M);
    Complex exponential_factor = (Complex)1.0/exponential_multiplication_factor;

    for(int m = 0; m <= M - 1; m++) {
        Complex S2 = (Complex)0;
        Double m_power = 1;
        for(int r = 0; r <= j; r++) {
            if(r > 0) {
                m_power *= m;
            }
            Complex z = binomial_coefficient(j, r) * m_power * H_table[j-r];
            S2 = S2 + z;
        }

        exponential_factor *= exponential_multiplication_factor;
        Complex z = exponential_factor * S2; // Here we should be multiplying S2 by exp(-2 pi alpha m/M), if we are doing things correctly.
        S = S + z;


    }
    S = S/(Double)pow((Double)M,j + 1);
    return S;
}

Complex H_method4(int j, Complex alpha, Double epsilon) {
    // Compute H(j, alpha) using a continued fraction expansion
    //
    // This code is largely copied from lcalc.
   
    stats::H_method4++;

    Complex P1 = (Double)1;
    Complex P2 = (Double)(j + 1);
    Complex P3 = (Complex)0;

    Complex Q1 = (Complex)0;
    Complex Q2 = (Double)1;
    Complex Q3 = (Complex)0;
    
    Complex u = PI * alpha;
    Double z = j + 1;
    Complex w = (Double)2 * PI * alpha;
    //ttype P1=1.,P2=z,P3,Q1=0.,Q2=1.,Q3;
    //ttype u=.5*w;
    //ttype t1,t2;

    int n=0;
    Double error = epsilon + 1;
    while( (error > epsilon || n < 3) && n < 1000000) {
        n++;
        P3=(z+n)*P2-(z+(n-1)*.5)*w*P1;
        Q3=(z+n)*Q2-(z+(n-1)*.5)*w*Q1;

        //t1=z+n;
        //t2=(z+(n-1)*.5)*w;
        //P3=t1*P2-t2*P1;
        //Q3=t1*Q2-t2*Q1;

        P1=P2;P2=P3;
        Q1=Q2;Q2=Q3;

        n++;
        P3=(z+n)*P2+(Double)n*u*P1;
        Q3=(z+n)*Q2+(Double)n*u*Q1;
        //t1=t1+1; t2=n*u;
        //P3=t1*P2+t2*P1;
        //Q3=t1*Q2+t2*Q1;

        P1=P2;P2=P3;
        Q1=Q2;Q2=Q3;

        //cout << P2/Q2 << " " << norm(Q2*P1-P2*Q1) / norm(Q2*P1*tolerance) <<endl;

        //to prevent overlow
        if(n%8==0&&(real(P2)>1.e50||real(P2)<-1.e50||imag(P2)>1.e50||imag(P2)<-1.e50)){
            P1=P1*(Double)(1.e-50);
            P2=P2*(Double)(1.e-50);
            Q1=Q1*(Double)(1.e-50);
            Q2=Q2*(Double)(1.e-50);

        }

        error = abs( ((P1 * Q2) - (P2 * Q1))/(Q1 * Q2) );

        //cout << P2 << "   "  << Q2 << P2/Q2 << endl;
        //cout << P2/Q2 << endl;

    }

    Complex g=P2/Q2;

    //cout<< "using cfrac for comp inc " << t << " " << n << endl;

    if(n>999999){
         cout << "Mofu. Continued fraction for g(z,w) failed to converge. z = "
         << z << "  w = " << w << endl;
         //exit(1);
    }

    g = exp(-w)/g;
    return g;
}

Complex H_method5(int j, Complex alpha, Double epsilon) {
    const int S = 20;
    const int D = 5;
    const int N = S * D;
    const int M = S * D;

    int n = floor(real(alpha)*N);
    int m = floor(imag(alpha)*M);

//    Complex alpha_0 = Complex( n/(Double)N, m/(Double)M );
    Complex alpha_0 = alpha + 1.0;
    Complex difference = alpha - alpha_0;
    Complex difference_power = 1;

//    Double error = epsilon + 2;

    Double error = 12;

    Complex Z = 0;
    int r = 0;
    Double one_over_j_plus_1 = 1.0/(j + 1.0);
    while(error > 10) {
        Complex z1 = two_pi_over_factorial_power(r) * difference_power;
        Z = Z + (Double)minus_one_power(r) * z1; // * H_method4(j + r, alpha_0, epsilon/20);
        error = abs(z1) * one_over_j_plus_1;
        difference_power *= difference;
        r++;
    }
//        cout << r << " ";

    return Z;
}
