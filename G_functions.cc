#include "theta_sums.h"
#include "precomputed_tables.h"

#include <iostream>
#include <cmath>

using namespace std;


Complex G(Complex alpha, Complex b, Double epsilon, int method) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = exp(2 pi i alpha t + 2 pi i b t^2)
    //
    // For certain cases of input we call a different function that computes
    // this integral using Euler-Maclaurin summation. Otherwise the integrand
    // gets expanded as a power series in b and the computation is reduced
    // to a bunch of calls to H().
    //
    
    
    //if(beta < 0) {
    //    // do something
    //}
    //if(b < 0 || b > 1) {
    //    // do something
    //}

    check_condition(imag(alpha) >= 0, "In function G(), Imag(alpha) should be nonnegative, but it isn't");
    //check_condition(b >= -1 && b <= 1, "In function G(), b should be between 0 and 1, but it isn't");

    //cout << alpha << "  " << abs(alpha) << endl;

    if(verbose::G) {
        cout << "G() called with:  " << endl;
        cout << "          alpha = " << alpha << endl;
        cout << "              b = " << b << endl;
        cout << "        epsilon = " << epsilon << endl;
    }

    if(epsilon >= 1) {
        if(verbose::G) {
            cout << "in G: epsilon >= 1, so returning 0" << endl;
        }
        return (Complex)0;
    }

    //
    // TODO: these if statements need to be optimized.
    //

    if(method == 0) {
        if(abs(alpha) < 1.5) {
            method = 2;
            stats::G_method2++;
        }
        else if( (1.5 <= abs(alpha)) && (abs(alpha) <= 2.5) && abs(b) > .024) {
            //cout << "here" << endl;
            method = 2;
            stats::G_method2++;
        }
        else if( (abs(alpha) >= 2.5) && (abs(alpha) - 2.0)/10.0 < abs(b)) {
            method = 2;
            stats::G_method2++;
        }
        else {
            method = 1;
            stats::G_method1++;
        }
    }
    if(verbose::G) {
        cout << "  In G(), using method " << method << endl;
    }

    if(b == Complex(0, 0)) {
        return H(0, -I * alpha, epsilon);
    }

    if(method == 2) {
        return G_via_Euler_MacLaurin(alpha, b, epsilon);
    }

    int N = max(1, to_int(ceil( -LOG(epsilon) )));

    //cout << "in G, using N = " << N << endl;

    Complex S = (Complex)0;
//    Double j_factorial = 1;

    //for(int j = 0; j < N + 1; j++) {
    Double error = epsilon + 1;
    int j = 0;
    Complex Ib_power = (Complex)1;
    if(verbose::G >= 2) {
        cout << "  In G(), using taylor expansion." << endl;
    }
    while(error > epsilon/2) {
        if(j > 0) {
            Ib_power *= (I * b);
        }
        Complex s = Ib_power * two_pi_over_factorial_power(j);
        Complex s2 = H(2 * j, -I * alpha, epsilon/(2 * abs(s) * (Double)N));
        Complex z = s * s2;
        if(verbose::G >= 2) {
            cout << "  Term " << j << " in expansion is " << s << " * " << s2 << " = " << z << endl;
        }
        S = S + z;
        j++;
        error = abs(s/(Complex)max( PI * abs(alpha) / (j + 1.0),  (Double)(2.0 * j + 1.0)));
    }


    return S;

}

inline Complex g(Complex alpha, Complex b, Double t) {
    return EXP(2 * PI * I * t * (alpha + b * t) );
}
Complex G_via_Euler_MacLaurin(Complex alpha, Complex b, Double epsilon) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = exp(2 pi i alpha t + 2 pi i b t^2)
    // using euler maclaurin summation

    // By rough experimentation, we decide to use Euler MacLaurin summation in the cases that
    // (i) |alpha| < 1.5
    // OR (ii) 1.5 <= |alpha| <= 2.5 AND b > .024
    // OR (iii) 2.5 <= |alpha| AND (|alpha| - 2)/10 < b
    //
    // (These heuristics might not work well when b is larger, but b will be at most .25 in our applications.)

    //int n =  a * (2 * m + ceil(2 * PI * (abs(alpha) + abs((Double)2 * b))));

    if(epsilon >= 1) {
        return (Complex)0;
    }

    int n = 4 * to_int(ceil(abs(alpha) + abs((Complex)2.0 * b) - LOG(epsilon)/(2 * PI)));
    
    //cout << epsilon << endl;

    //cout << alpha << endl;
    //cout << b << endl;
    //cout << epsilon << endl;

    //cout << n << endl;

    //int n = 10000;

    //cout << n << endl;

    Complex S = (Complex)0;
    for(int j = 0; j <= n; j++) {
        S = S + g(alpha, b, (Double)j/(Double)n);
    }
    S = S - (Double)(.5) * (g(alpha, b, 0) + g(alpha, b, 1));
    S = S/(Double)n;

    //cout << S << endl;

    
    //if(2 * m - 3 > 21) {
    //    cout << "I can't compute that many derivatives. Sorry." << endl;
    //    return Complex(0, 0);
    //}

    //cout << m << endl;

//    initialize_power_arrays(21, alpha, b);
    Double n_power = 1;
//    for(int r = 1; r <= m - 1; r++) {

    Double error = epsilon + 1; 
    int r = 1;

    Complex * p;
    Complex * p_prev = new Complex[0];

    while(error > epsilon/2) {
       // cout << error << endl;
    
        p = new Complex[2 * r - 2 + 1];
        g_derivative_polynomial(2 * r - 2, p, p_prev, alpha, b);
        delete [] p_prev;
        p_prev = p;
        p = new Complex[2 * r - 1 + 1];
        g_derivative_polynomial(2 * r - 1, p, p_prev, alpha, b);
        delete [] p_prev;
        p_prev = p;

        Complex derivative_at_1 = (Complex)0;
        for(int k = 0; k <= 2 * r - 1; k++) {
            derivative_at_1 = derivative_at_1 + p[k];
        }
        derivative_at_1 *= exp(2 * PI * I * (alpha + b));
        Complex derivative_at_0 = p[0];

        n_power *= ((Double)n * (Double)n);
//        Complex z = bernoulli_table[2 * r]/(factorial(2 * r) * n_power) * (g_derivative_at_1(2 * r - 1) - g_derivative_at_0(2 * r - 1));
        Complex z = bernoulli_table[2 * r]/(factorial(2 * r) * n_power) * (derivative_at_1 - derivative_at_0);

//        cout << g_derivative_at_1(2 * r - 1) << "    " << derivative_at_1 << endl;
//        cout << g_derivative_at_0(2 * r - 1) << "    " << derivative_at_0 << endl;
//        cout << "    " << bernoulli_table[2*r] << endl;
//        cout << "    " << factorial(2 * r) << endl;
//        cout << "    " n_power << endl;

//        cout << z << "   " << g_derivative_at_1(2 * r - 1) << "   " <<  g_derivative_at_0(2 * r - 1)  << endl;
        S = S - z;
        r = r + 1;
        error = abs(z);
       // cout << error << endl;
    }

     delete [] p;
 
    if(verbose::G) {
        cout << "In G(), using Euler-Maclaurin computed G(" << alpha << ", " << b << ") = " << S << endl;
    }

    return S;

}


Complex G(Complex alpha, Complex b, int n, int j, Double epsilon, int method) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)/(2n)^j
    //
    // For certain cases of input we call a different function that computes
    // this integral using Euler-Maclaurin summation. Otherwise the integrand
    // gets expanded as a power series in b and the computation is reduced
    // to a bunch of calls to H().
    //
    
    
    //if(beta < 0) {
    //    // do something
    //}
    //if(b < 0 || b > 1) {
    //    // do something
    //}
    
    if(j == 0) {
        return G(alpha, b, epsilon, method);
    }

    check_condition(imag(alpha) >= 0, "In function G(), Imag(alpha) should be nonnegative, but it isn't");
    //check_condition(b >= -1 && b <= 1, "In function G(), b should be between 0 and 1, but it isn't");

    //cout << alpha << "  " << abs(alpha) << endl;

    if(verbose::G) {
        cout << "G() called with:  " << endl;
        cout << "          alpha = " << alpha << endl;
        cout << "              b = " << b << endl;
        cout << "        epsilon = " << epsilon << endl;
    }

    if(epsilon >= 1) {
        if(verbose::G) {
            cout << "in G: epsilon >= 1, so returning 0" << endl;
        }
        return (Complex)0;
    }

    //
    // TODO: these if statements need to be optimized.
    //

    if(method == 0) {
        if(abs(alpha) < 1.5) {
            method = 2;
            stats::G_method2++;
        }
        else if( (1.5 <= abs(alpha)) && (abs(alpha) <= 2.5) && abs(b) > .024) {
            //cout << "here" << endl;
            method = 2;
            stats::G_method2++;
        }
        else if( (abs(alpha) >= 2.5) && (abs(alpha) - 2.0)/10.0 < abs(b)) {
            method = 2;
            stats::G_method2++;
        }
        else {
            method = 1;
            stats::G_method1++;
        }
    }
    if(verbose::G) {
        cout << "  In G(), using method " << method << endl;
    }

    if(method == 2) {
        return G_via_Euler_MacLaurin(alpha, b, n, j, epsilon);
    }

    if(n > 0) {
        Complex S = 0;
        for(int s = 0; s <= j; s++) {
            Double z = binomial_coefficient(j, s) * pow(2, -j) * pow(n, -s);
            S = S + z * G(alpha, b, 0, s, epsilon/(z * j));
        }
        return S;
    }

    if(b == Complex(0, 0)) {
        return H(j, -I * alpha, epsilon);
    }

    int N = max(1, to_int(ceil( -LOG(epsilon) )));

    //cout << "in G, using N = " << N << endl;

    Complex S = (Complex)0;
//    Double j_factorial = 1;

    //for(int j = 0; j < N + 1; j++) {
    Double error = epsilon + 1;
    int r = 0;
    Complex Ib_power = (Complex)1;
    if(verbose::G >= 2) {
        cout << "  In G(), using taylor expansion." << endl;
    }
    while(error > epsilon/2) {
        if(r > 0) {
            Ib_power *= (I * b);
        }
        Complex s = Ib_power * two_pi_over_factorial_power(r);
        Complex s2 = H(j + 2 * r, -I * alpha, epsilon/(2 * abs(s) * (Double)N));
        Complex z = s * s2;
        if(verbose::G >= 2) {
            cout << "  Term " << r << " in expansion is " << s << " * " << s2 << " = " << z << endl;
        }
        S = S + z;
        r++;
        error = abs(s/(Complex)max( PI * abs(alpha) / (r + 1.0),  (Double)(2.0 * r + 1.0)));
    }


    return S;

}


inline Double POW(Double a, Double b) {
    if(a == 0 && b == 0) {
        return 1;
    }
    return pow(a, b);
}




inline Complex g(Complex alpha, Complex b, Double n, Double j, Double t) {
    return POW(t + n, j) * EXP(2 * PI * I * t * (alpha + b * t) );
}
Complex G_via_Euler_MacLaurin(Complex alpha, Complex b, int n, int j, Double epsilon) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)/(2n)^j
    // using euler maclaurin summation

    // By rough experimentation, we decide to use Euler MacLaurin summation in the cases that
    // (i) |alpha| < 1.5
    // OR (ii) 1.5 <= |alpha| <= 2.5 AND b > .024
    // OR (iii) 2.5 <= |alpha| AND (|alpha| - 2)/10 < b
    //
    // (These heuristics might not work well when b is larger, but b will be at most .25 in our applications.)

    //int n =  a * (2 * m + ceil(2 * PI * (abs(alpha) + abs((Double)2 * b))));

    if(j == 0) {
        return G_via_Euler_MacLaurin(alpha, b, epsilon);
    }

    if(epsilon >= 1) {
        return (Complex)0;
    }

    int N = 4 * to_int(ceil(abs(alpha) + abs((Complex)2.0 * b) - LOG(epsilon)/(2 * PI)));
    
    Double two_n_to_the_j;
    if(n != 0)
        two_n_to_the_j = pow(2 * n, j);
    else
        two_n_to_the_j = 1;
    Double one_over_two_n_to_the_j = 1/two_n_to_the_j;

    //cout << epsilon << endl;

    //cout << alpha << endl;
    //cout << b << endl;
    //cout << epsilon << endl;

    //cout << n << endl;

    //int n = 10000;

    //cout << n << endl;

    Complex S = (Complex)0;
    for(int s = 0; s <= N; s++) {
        S = S + g(alpha, b, n, j, (Double)s/(Double)N) * one_over_two_n_to_the_j;
    }
    S = S - (Double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1)) * one_over_two_n_to_the_j;
    S = S/(Double)N;

    //cout << S << endl;

    
    //if(2 * m - 3 > 21) {
    //    cout << "I can't compute that many derivatives. Sorry." << endl;
    //    return Complex(0, 0);
    //}

    //cout << m << endl;

//    initialize_power_arrays(21, alpha, b);
    Double N_power = 1;
//    for(int r = 1; r <= m - 1; r++) {

    Double error = epsilon + 1; 
    int r = 1;

    Complex * p;
    Complex * p_prev = new Complex[j+1];
    if(n == 0) {
        for(int s = 0; s <= j - 1; s++) {
            p_prev[s] = 0;
        }
        p_prev[j] = 1;
    }
    else {
        for(int s = 0; s <= j; s++) {
            p_prev[s] = pow(2, -j) * pow(n, -s) * binomial_coefficient(j, s);
        }
    }


    while(error > epsilon/2) {
       // cout << error << endl;
    
        if(r > 1) {
            p = new Complex[2 * r - 2 + 1 + j];
            g_derivative_polynomial(2 * r - 2 + j, p, p_prev, alpha, b);
            delete [] p_prev;
            p_prev = p;
        }
        p = new Complex[2 * r - 1 + 1 + j];
        g_derivative_polynomial(2 * r - 1 + j, p, p_prev, alpha, b);
        delete [] p_prev;
        p_prev = p;

        Complex derivative_at_1 = (Complex)0;
        for(int k = 0; k <= 2 * r - 1 + j; k++) {
            derivative_at_1 = derivative_at_1 + p[k];
        }
        derivative_at_1 *= exp(2 * PI * I * (alpha + b));
        Complex derivative_at_0 = p[0];

        N_power *= ((Double)N * (Double)N);
//        Complex z = bernoulli_table[2 * r]/(factorial(2 * r) * n_power) * (g_derivative_at_1(2 * r - 1) - g_derivative_at_0(2 * r - 1));
        Complex z = bernoulli_table[2 * r]/(factorial(2 * r) * N_power) * (derivative_at_1 - derivative_at_0);

//        cout << g_derivative_at_1(2 * r - 1) << "    " << derivative_at_1 << endl;
//        cout << g_derivative_at_0(2 * r - 1) << "    " << derivative_at_0 << endl;
//        cout << "    " << bernoulli_table[2*r] << endl;
//        cout << "    " << factorial(2 * r) << endl;
//        cout << "    " n_power << endl;

//        cout << z << "   " << g_derivative_at_1(2 * r - 1) << "   " <<  g_derivative_at_0(2 * r - 1)  << endl;
        S = S - z;
        r = r + 1;
        error = abs(z);
       // cout << error << endl;
    }

     delete [] p;
 
    if(verbose::G) {
        cout << "In G(), using Euler-Maclaurin computed G(" << alpha << ", " << b << ") = " << S << endl;
    }

    return S;

}


