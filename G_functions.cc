#include "theta_sums.h"
#include "precomputed_tables.h"
#include "log.h"

#include <iostream>
#include <cmath>

using namespace std;


Complex G(Complex alpha, Complex b, int n, int j, Double epsilon, int method) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    //
    // For certain cases of input we call a different function that computes
    // this integral using Euler-Maclaurin summation. Otherwise the integrand
    // gets expanded as a power series in b and the computation is reduced
    // to a bunch of calls to H().
    //
    
    //cout << "New code called" << endl;
    
    // Note b is _always_ (I think) either real or purely imaginary.
    // We should take advantage of this.
    //
    // alpha is typically complex, however.


    check_condition(imag(alpha) >= 0, "In function G(), Imag(alpha) should be nonnegative, but it isn't");

    if(verbose::G) {
        cout << "G() called with:  " << endl;
        cout << "          alpha = " << alpha << endl;
        cout << "              b = " << b << endl;
        cout << "              n = " << n << endl;
        cout << "              j = " << j << endl;
        cout << "        epsilon = " << epsilon << endl;
    }

    //if(epsilon >= pow((Double)n + 1, (Double)j)) {
    if(fastlog(epsilon) > j * fastlog(n+1)) {
        if(verbose::G) {
            cout << "in G: epsilon >= 1, so returning 0" << endl;
        }
        return (Complex)0;
    }

    //
    // TODO: these if statements need to be optimized.
    //

    Double norm_alpha = norm(alpha);
    Double norm_b = norm(b);

    if(method == 0) {
//        if(abs(alpha) < 5) {
//            return exp(I * alpha);
//        }
        //if(abs(alpha) < 1.5) {
        if(norm_alpha < 2.25) {
            method = 2;
            stats::G_method2++;
        }
        //else if( (1.5 <= abs(alpha)) && (abs(alpha) <= 2.5) && abs(b) > .024) {
        else if( (2.25 <= norm_alpha) && (norm_alpha <= 6.25) && norm_b > .000576) {
            method = 2;
            stats::G_method2++;
        }
        //else if( (abs(alpha) >= 2.5) && (abs(alpha) - 2.0)/10.0 < abs(b)) {
        else if( (norm_alpha >= 6.25) && (norm_alpha - 4.0) < 100 * norm_b) {
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
        Double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            Double z = binomial_coefficient(j, s) * n_power;
            S = S + z * G(alpha, b, 0, j - s, epsilon/(z * (j + 1)));
            n_power *= n;
        }
        return S;
    }

    // At this point we assume that n == 0

    if(b == Complex(0, 0)) {
        return H(j, -I * alpha, epsilon);
    }

    //int N = max(1, to_int(ceil( -LOG(epsilon) )));
    int N = max(1, fastlog(epsilon));

    Complex S = (Complex)0;

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
        error = abs(s/(Complex)max( PI * abs(alpha) / (r + 1.0 + j/2),  (Double)(2.0 * r + 1.0 + j)));
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
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    // using euler maclaurin summation

    // By rough experimentation, we decide to use Euler MacLaurin summation in the cases that
    // (i) |alpha| < 1.5
    // OR (ii) 1.5 <= |alpha| <= 2.5 AND b > .024
    // OR (iii) 2.5 <= |alpha| AND (|alpha| - 2)/10 < b
    //
    // (These heuristics might not work well when b is larger, but b will be at most .25 in our applications.)

    //if(epsilon >= 1) {
    //    return (Complex)0;
    //}

    //int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((Complex)2.0 * b) + max(-LOG(epsilon)/(2 * PI), 0.0) ) * (1 + j * log(n + 1)/4.0) )), 1);

//    int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((Complex)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) )), 1);

    int N = to_int(ceil(  ( abs(alpha) + abs((Complex)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) ));
    N = 4 * max(N, 1);

    if(verbose::G >= 2) {
        cout << "In G(), using " << N << " sample points in Euler-Maclaurin summation." << endl;
    }

    /*
    Double two_n_to_the_j = 1;
    if(n != 0)
        two_n_to_the_j = pow(2 * n, j);
    else
        two_n_to_the_j = 1;
    Double one_over_two_n_to_the_j = 1/two_n_to_the_j;
    */

    Complex S = (Complex)0;
    for(int s = 0; s <= N; s++) {
        //S = S + g(alpha, b, n, j, (Double)s/(Double)N) * one_over_two_n_to_the_j;
        S = S + g(alpha, b, n, j, (Double)s/(Double)N);
    }
    //S = S - (Double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1)) * one_over_two_n_to_the_j;
    S = S - (Double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1));
    S = S/(Double)N;

    Double N_power = 1;

    Double error = epsilon + 1; 
    int r = 1;


    Complex p1[bernoulli_range + j];            // below p will need to be as big as 2r + j
    Complex p2[bernoulli_range + j];            // for increasing r. If 2r ever gets as large
                                                // as bernoulli_range, then the computation
                                                // will fail anyway, so the biggest p
                                                // that makes sense is of size bernoulli_range + j.
                                                //
                                                // As long as bernoulli_range isn't too big, this
                                                // should be hopefully be a good choice.

    Complex * ptmp;                             // A temporary variable used for swapping
                                                // p and p_prev.

    Complex * p = p1;
    Complex * p_prev = p2;

    if(n == 0) {
        for(int s = 0; s <= j - 1; s++) {
            p_prev[s] = 0;
        }
        p_prev[j] = 1;
    }
    else {
        /*
        for(int s = 0; s <= j; s++) {
            //p_prev[s] = pow(2, -j) * pow(n, -s) * binomial_coefficient(j, s);
            p_prev[s] = pow(n, j - s) * binomial_coefficient(j, s);
        }
        */
        Double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            //p_prev[s] = pow(2, -j) * pow(n, -s) * binomial_coefficient(j, s);
            p_prev[j - s] = n_power * binomial_coefficient(j, s);
            n_power *= n;
        }
    }

    if(verbose::G >= 2) {
        cout << "In G(), starting to compute correction terms in Euler-Maclaurin summation." << endl;
    }

    while(2 * error > epsilon * epsilon) {
    
        if(r > 1) {
//            p = new Complex[2 * r - 2 + 1 + j];
            g_derivative_polynomial(2 * r - 2 + j, p, p_prev, alpha, b);
//            delete [] p_prev;
//            p_prev = p;
            ptmp = p;
            p = p_prev;
            p_prev = ptmp;
        }
//        p = new Complex[2 * r - 1 + 1 + j];
        g_derivative_polynomial(2 * r - 1 + j, p, p_prev, alpha, b);
//        delete [] p_prev;
//        p_prev = p;

        Complex derivative_at_1 = (Complex)0;
        for(int k = 0; k <= 2 * r - 1 + j; k++) {
            derivative_at_1 = derivative_at_1 + p[k];
        }
        derivative_at_1 *= exp(2 * PI * I * (alpha + b));
        Complex derivative_at_0 = p[0];

        N_power *= ((Double)N * (Double)N);
        Complex z = bernoulli_table[2 * r]/(factorial(2 * r) * N_power) * (derivative_at_1 - derivative_at_0);

        S = S - z;
        error = norm(z);
        if(verbose::G >= 2) {
            cout << "   Correction term " << r << " has size: " << error << endl;
            cout << "      derivative at 0 = " << derivative_at_0 << endl;
            cout << "      derivative at 1 = " << derivative_at_1 << endl;
            cout << "                 N^2r = " << N_power << endl;
        }
        r = r + 1;
        ptmp = p;
        p = p_prev;
        p_prev = ptmp;
    }

//    delete [] p;
 
    if(verbose::G) {
        cout << "In G(), using Euler-Maclaurin computed G(" << alpha << ", " << b << ") = " << S << endl;
    }

    return S;

}


