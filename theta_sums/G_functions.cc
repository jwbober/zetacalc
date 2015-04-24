//
// Functions to compute the integral
//
// G(alpha, b, n, j) = int_0^1 (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2) dt
//
// (We apologize for the strange naming convention of the variables. There is
// no particular reason that one is greek and one is not.)
//
// We have two different methods to compute this integral. One method
// is Euler-Maclaurin summation. The second method (which below is
// called method1) is:
//
//      1. If n != 0 and j != 0, do a binomial expansion to
//         write the integral as a linear combination of
//         integrals of the same form, but with n = 0.
//      
//      2. In the n = 0 case, do a taylor expansion of the integrand
//         in the "b" variable, which reduces the the integral to a sum
//         of integrals of the form
//
//              H(alpha, j) = \int_0^1 t^j exp(2 pi i alpha t)
//
//         These "H functions" are dealt with in a separate file.
//
//
// There are various functions in this file which duplicate functionality.
// In general, the parameters alpha and b are complex numbers, but in many
// instances one of them is real, so we have extra functions which take
// advantage of this. Also, in a certain case we compute this function many
// times with b = I/2pi, so this is also written in as a specific function.



#include "theta_sums.h"
#include "precomputed_tables.h"
#include "log.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using std::complex;
using std::pow;
using std::max;

complex<double> G_method1_I_over_twopi(complex<double> alpha, int n, int j, double epsilon);
complex<double> G_via_Euler_MacLaurin_I_over_twopi(complex<double> alpha, int n, int j, double epsilon);

complex<double> G_via_Euler_MacLaurin_deprecated(complex<double> alpha, complex<double> b, int n, int j, double epsilon);
complex<double> G_method1_deprecated(complex<double> alpha, complex<double> b, int n, int j, double epsilon);

complex<double> G_deprecated(complex<double> alpha, complex<double> b, int n, int j, double epsilon, int method) {
    //
    // OLD METHOD NO LONGER USED, SINCE b is ALWAYS EITHER REAL OR PURELY
    // IMAGINARY.
    //
    
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    //
    // For certain cases of input we call a different function that computes
    // this integral using Euler-Maclaurin summation. Otherwise the integrand
    // gets expanded as a power series in b and the computation is reduced
    // to a bunch of calls to H().
    //
    
    // Note b is _always_ (I think) either real or purely imaginary.
    // We should take advantage of this.
    //
    // alpha is typically complex, however.



    //if(imag(b) == 0) {
    //    return G_R(alpha, real(b), n, j, epsilon, method);
    //}
    //else if(real(b) == 0) {
    //    return G_I(alpha, imag(b), n, j, epsilon, method);
    //}

    //if(epsilon >= pow((double)n + 1, (double)j)) {
    if(fastlog2(epsilon) > j * fastlog2(n+1)) {
        return (complex<double>)0;
    }

    //
    // TODO: these if statements need to be optimized.
    //

    double norm_alpha = norm(alpha);
    double norm_b = norm(b);

    if(method == 0) {
//        if(abs(alpha) < 5) {
//            return exp(I * alpha);
//        }
        //if(abs(alpha) < 1.5) {
        if(norm_alpha < 2.25) {
            method = 2;
        }
        //else if( (1.5 <= abs(alpha)) && (abs(alpha) <= 2.5) && abs(b) > .024) {
        else if( (2.25 <= norm_alpha) && (norm_alpha <= 6.25) && norm_b > .000576) {
            method = 2;
        }
        //else if( (abs(alpha) >= 2.5) && (abs(alpha) - 2.0)/10.0 < abs(b)) {
        else if( (norm_alpha >= 6.25) && (norm_alpha - 4.0) < 100 * norm_b) {
            method = 2;
        }
        else {
            method = 1;
        }
    }

    if(method == 2) {
        //if(imag(b) == 0) {
        //    return G_via_Euler_MacLaurin_R(alpha, real(b), n, j, epsilon);
        //}
        //else if(real(b) == 0) {
        //    return G_via_Euler_MacLaurin_I(alpha, imag(b), n, j, epsilon);
        //}
        //else
            return G_via_Euler_MacLaurin_deprecated(alpha, b, n, j, epsilon);
    }
    else {
        //if(imag(b) == 0) {
        //    return G_method1_R(alpha, real(b), n, j, epsilon);
        //}
        //if(real(b) == 0) {
        //    return G_method1_I(alpha, imag(b), n, j, epsilon);
        //}
        return G_method1_deprecated(alpha, b, n, j, epsilon);
    }

}

complex<double> G(complex<double> alpha, double b, int n, int j, double epsilon, int method) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    //
    // Specialized for b real.
    //
    // For certain cases of input we call a different function that computes
    // this integral using Euler-Maclaurin summation. Otherwise the integrand
    // gets expanded as a power series in b and the computation is reduced
    // to a bunch of calls to H().
    //
    
    // Note b is _always_ (I think) either real or purely imaginary.
    // We should take advantage of this.
    //
    // alpha is typically complex, however.

    //if(epsilon >= pow((double)n + 1, (double)j)) {
    if(fastlog(epsilon) > j * fastlog(n+1) + j) {
        return (complex<double>)0;
    }


    //
    // TODO: these if statements need to be optimized.
    //

    double norm_alpha = norm(alpha);
    double norm_b = abs(b) * abs(b);

    if(method == 0) {
//        if(abs(alpha) < 5) {
//            return exp(I * alpha);
//        }
        //if(abs(alpha) < 1.5) {
        if(norm_alpha < 2.25) {
            method = 2;
        }
        //else if( (1.5 <= abs(alpha)) && (abs(alpha) <= 2.5) && abs(b) > .024) {
        else if( (2.25 <= norm_alpha) && (norm_alpha <= 6.25) && abs(b) > .024) {
            method = 2;
        }
        //else if( (abs(alpha) >= 2.5) && (abs(alpha) - 2.0)/10.0 < abs(b)) {
        else if( (norm_alpha >= 6.25) && (norm_alpha - 4.0) < 100 * norm_b) {
            method = 2;
        }
        else {
            method = 1;
        }
    }

    if(method == 2) {
        return G_via_Euler_MacLaurin_R(alpha, b, n, j, epsilon);
    }
    else {
        return G_method1_R(alpha, b, n, j, epsilon);
    }

}

complex<double> G_I(complex<double> alpha, double b, int n, int j, double epsilon, int method) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    //
    // Specialized for b imaginary. b has an implicit I in front.
    //
    // For certain cases of input we call a different function that computes
    // this integral using Euler-Maclaurin summation. Otherwise the integrand
    // gets expanded as a power series in b and the computation is reduced
    // to a bunch of calls to H().
    //
    
    // Note b is _always_ (I think) either real or purely imaginary.
    // We should take advantage of this.
    //
    // alpha is typically complex, however.


    //if(epsilon >= pow((double)n + 1, (double)j)) {
    if(fastlog(epsilon) > j * fastlog(n+1)) {
        return (complex<double>)0;
    }

    //
    // TODO: these if statements need to be optimized.
    //

    double norm_alpha = norm(alpha);
    double norm_b = abs(b) * abs(b);

    if(method == 0) {
//        if(abs(alpha) < 5) {
//            return exp(I * alpha);
//        }
        //if(abs(alpha) < 1.5) {
        if(norm_alpha < 2.25) {
            method = 2;
        }
        //else if( (1.5 <= abs(alpha)) && (abs(alpha) <= 2.5) && abs(b) > .024) {
        else if( (2.25 <= norm_alpha) && (norm_alpha <= 6.25) && norm_b > .000576) {
            method = 2;
        }
        //else if( (abs(alpha) >= 2.5) && (abs(alpha) - 2.0)/10.0 < abs(b)) {
        else if( (norm_alpha >= 6.25) && (norm_alpha - 4.0) < 100 * norm_b) {
            method = 2;
        }
        else {
            method = 1;
        }
    }

    if(method == 2) {
        return G_via_Euler_MacLaurin_I(alpha, b, n, j, epsilon);
    }
    else {
        return G_method1_I(alpha, b, n, j, epsilon);
    }

}

complex<double> G_I_over_twopi(complex<double> alpha, int n, int j, double epsilon, int method) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    //
    // Specialized for b = I/2pi
    //
    // For certain cases of input we call a different function that computes
    // this integral using Euler-Maclaurin summation. Otherwise the integrand
    // gets expanded as a power series in b and the computation is reduced
    // to a bunch of calls to H().
    //
    
    // Note b is _always_ (I think) either real or purely imaginary.
    // We should take advantage of this.
    //
    // alpha is typically complex, however.

    //if(epsilon >= pow((double)n + 1, (double)j)) {
    if(fastlog(epsilon) > j * fastlog(n+1) + j) {
        return (complex<double>)0;
    }

    //
    // TODO: these if statements need to be optimized.
    //

    double norm_alpha = norm(alpha);

    if(method == 0) {
//        if(abs(alpha) < 5) {
//            return exp(I * alpha);
//        }
        //if(abs(alpha) < 1.5) {
        if(norm_alpha <= 6.25) {
            method = 2;
        }
        //else if( (abs(alpha) >= 2.5) && (abs(alpha) - 2.0)/10.0 < abs(b)) {
        else if(norm_alpha - 4.0 < 25/(PI * PI)) {
            method = 2;
        }
        else {
            method = 1;
        }
    }


    if(method == 2) {
        return G_via_Euler_MacLaurin_I_over_twopi(alpha, n, j, epsilon);
        //return G_via_Euler_MacLaurin_I(alpha, 1.0/(2 * PI), n, j, epsilon);
    }
    else {
        return G_method1_I_over_twopi(alpha, n, j, epsilon);
    }

}

complex<double> G_method1_deprecated(complex<double> alpha, complex<double> b, int n, int j, double epsilon) {
    if(n > 0) {
        complex<double> S = 0;
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            double z = binomial_coefficient(j, s) * n_power;
            S = S + z * G_method1_deprecated(alpha, b, 0, j - s, epsilon/(z * (j + 1)));
            n_power *= n;
        }
        return S;
    }

    // At this point we assume that n == 0

    if(b == complex<double>(0, 0)) {
        return H(j, -I * alpha, epsilon);
    }

    //int N = max(1, to_int(ceil( -LOG(epsilon) )));
    int N = max(1, -fastlog(epsilon));

    complex<double> S = (complex<double>)0;

    double error = epsilon + 1;
    int r = 0;
    complex<double> Ib_power = (complex<double>)1;
    while(error > epsilon/2) {
        if(r > 0) {
            Ib_power *= (I * b);
        }
        complex<double> s = Ib_power * two_pi_over_factorial_power(r);
        complex<double> s2 = H(j + 2 * r, -I * alpha, epsilon/(2 * abs(s) * (double)N));
        complex<double> z = s * s2;
        S = S + z;
        r++;
        error = abs(s/(complex<double>)max( PI * abs(alpha) / (r + 1.0 + j/2),  (double)(2.0 * r + 1.0 + j)));
    }


    return S;
}

complex<double> G_method1_R(complex<double> alpha, double b, int n, int j, double epsilon) {
    // Specialized for b real.
    

    if(n > 0) {
        complex<double> S = 0;
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            double z = binomial_coefficient(j, s) * n_power;
            S = S + z * G_method1_R(alpha, b, 0, j - s, epsilon/(z * (j + 1)));
            n_power *= n;
        }
        return S;
    }

    // At this point we assume that n == 0

    if(b == 0) {
        return H(j, -I * alpha, epsilon);
    }

    //int N = max(1, to_int(ceil( -LOG(epsilon) )));
    int N = max(1, -fastlog(epsilon));

    complex<double> S = (complex<double>)0;

    double error = epsilon + 1;
    int r = 0;
    //complex<double> Ib_power = (complex<double>)1;
    double b_power = 1.0;
    while(error > epsilon/2) {
        if(r > 0) {
            //Ib_power *= (I * b);
            b_power *= b;
        }
    
        //Ib_power = (ib)^r

        double s = b_power * two_pi_over_factorial_power(r);
        complex<double> s2 = H(j + 2 * r, -I * alpha, epsilon/(2 * abs(s) * (double)N));
        complex<double> z = I_power(r) * s * s2;
        S = S + z;
        r++;
        error = abs(s/max( PI * abs(alpha) / (r + 1.0 + j/2),  (double)(2.0 * r + 1.0 + j)));
    }

    return S;
}

complex<double> G_method1_I(complex<double> alpha, double b, int n, int j, double epsilon) {
    //
    // specialized for b purely imaginary. b has an implicit I in front
    //
    if(n > 0) {
        complex<double> S = 0;
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            double z = binomial_coefficient(j, s) * n_power;
            S = S + z * G_method1_I(alpha, b, 0, j - s, epsilon/(z * (j + 1)));
            n_power *= n;
        }
        return S;
    }

    // At this point we assume that n == 0

    if(b == 0) {
        return H(j, -I * alpha, epsilon);
    }

    //int N = max(1, to_int(ceil( -LOG(epsilon) )));
    int N = max(1, -fastlog(epsilon));

    complex<double> S = (complex<double>)0;

    double error = epsilon + 1;
    int r = 0;
    double Ib_power = 1.0;
    while(error > epsilon/2) {
        if(r > 0) {
            Ib_power *= (-b);
        }
        double s = Ib_power * two_pi_over_factorial_power(r);
        complex<double> s2 = H(j + 2 * r, -I * alpha, epsilon/(2 * abs(s) * (double)N));
        complex<double> z = s * s2;
        S = S + z;
        r++;
        error = abs(s/max( PI * abs(alpha) / (r + 1.0 + j/2),  (double)(2.0 * r + 1.0 + j)));
    }


    return S;
}

complex<double> G_method1_I_over_twopi(complex<double> alpha, int n, int j, double epsilon) {
    //
    // specialized for b purely imaginary. b has an implicit I in front
    //
    if(n > 0) {
        complex<double> S = 0;
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            double z = binomial_coefficient(j, s) * n_power;
            S = S + z * G_method1_I_over_twopi(alpha, 0, j - s, epsilon/(z * (j + 1)));
            n_power *= n;
        }
        return S;
    }

    // At this point we assume that n == 0

    //int N = max(1, to_int(ceil( -LOG(epsilon) )));
    int N = max(1, fastlog(epsilon));

    complex<double> S = (complex<double>)0;

    double error = epsilon + 1;
    int r = 0;
    //double Ib_power = 1.0;
    while(error > epsilon/2) {
        //if(r > 0) {
        //    Ib_power *= (-b);
        //}
        // Ib_power = (-1/2pi)^r
        // so this is just...
        //double s = Ib_power * two_pi_over_factorial_power(r);
        double s = minus_one_power(r) / factorial(r);
        complex<double> s2 = H(j + 2 * r, -I * alpha, epsilon/(2 * abs(s) * (double)N));
        complex<double> z = s * s2;
        S = S + z;
        r++;
        error = abs(s/max( PI * abs(alpha) / (r + 1.0 + j/2),  (double)(2.0 * r + 1.0 + j)));
    }


    return S;
}

inline double POW(double a, double b) {
    if(a == 0 && b == 0) {
        return 1;
    }
    return pow(a, b);
}

inline complex<double> g(complex<double> alpha, complex<double> b, double n, double j, double t) {
    return POW(t + n, j) * EXP(2 * PI * I * t * (alpha + b * t) );
}
complex<double> G_via_Euler_MacLaurin_deprecated(complex<double> alpha, complex<double> b, int n, int j, double epsilon) {
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
    //    return (complex<double>)0;
    //}

    //int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-LOG(epsilon)/(2 * PI), 0.0) ) * (1 + j * log(n + 1)/4.0) )), 1);

//    int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) )), 1);

    int N = to_int(ceil(  ( abs(alpha) + abs((double)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) ));
    N = 2 * max(N, 1);


    /*
    double two_n_to_the_j = 1;
    if(n != 0)
        two_n_to_the_j = pow(2 * n, j);
    else
        two_n_to_the_j = 1;
    double one_over_two_n_to_the_j = 1/two_n_to_the_j;
    */

    complex<double> S = (complex<double>)0;
    {
        complex<double> alpha_term = 1.0;
        complex<double> alpha_term_multiplier = EXP(2 * PI * I * alpha/(double)N);
        double t = 0;
        double t_increment = 1.0/N;
        switch(j) {
            case 0:
                if(imag(b) == 0) {
                    for(int s = 0; s <= N; s++) {
                        S = S + complex<double>(cos(2 * PI * real(b) * t * t), sin(2 * PI * real(b) * t * t)) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                else if(real(b) == 0) {
                    for(int s = 0; s <= N; s++) {
                        S = S + EXP(-2 * PI * imag(b) * t * t) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                else {
                    for(int s = 0; s <= N; s++) {
                        S = S + EXP(2 * PI * I * b * t * t) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                break;
            case 1:
                if(imag(b) == 0) {
                    for(int s = 0; s <= N; s++) {
                        S = S + (n + t) * complex<double>(cos(2 * PI * real(b) * t * t), sin(2 * PI * real(b) * t * t)) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                else if(real(b) == 0) {
                    for(int s = 0; s <= N; s++) {
                        S = S + (n + t) * EXP(-2 * PI * imag(b) * t * t) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                else {
                    for(int s = 0; s <= N; s++) {
                        S = S + (n + t) * EXP(2 * PI * I * b * t * t) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                break;

            default:
                if(imag(b) == 0) {
                    for(int s = 0; s <= N; s++) {
                        S = S + pow(t + n, j) * complex<double>(cos(2 * PI * real(b) * t * t), sin(2 * PI * real(b) * t * t)) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                else if(real(b) == 0) {
                    for(int s = 0; s <= N; s++) {
                        S = S + pow(t + n, j) * EXP(-2 * PI * imag(b) * t * t) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
                else {
                    for(int s = 0; s <= N; s++) {
                        S = S + pow(t + n, j) * EXP(2 * PI * I * b * t * t) * alpha_term;
                        alpha_term *= alpha_term_multiplier;
                        t += t_increment;
                    }
                }
            
        }


    }
    //S = S - (double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1)) * one_over_two_n_to_the_j;
    S = S - (double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1));
    S = S/(double)N;

    double N_power = 1;

    double error = 10 * epsilon * epsilon;
    int r = 1;


    complex<double> p1[bernoulli_range + max_j];            // below p will need to be as big as 2r + j
    complex<double> p2[bernoulli_range + max_j];            // for increasing r. If 2r ever gets as large
                                                // as bernoulli_range, then the computation
                                                // will fail anyway, so the biggest p
                                                // that makes sense is of size bernoulli_range + j.
                                                //
                                                // As long as bernoulli_range isn't too big, this
                                                // should be hopefully be a good choice.

    complex<double> * ptmp;                             // A temporary variable used for swapping
                                                // p and p_prev.

    complex<double> * p = p1;
    complex<double> * p_prev = p2;

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
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            //p_prev[s] = pow(2, -j) * pow(n, -s) * binomial_coefficient(j, s);
            p_prev[j - s] = n_power * binomial_coefficient(j, s);
            n_power *= n;
        }
    }


    if(imag(b) == 0) {
        while(4 * error > epsilon * epsilon) {
            if(r > 1) {
                g_derivative_polynomial(2 * r - 2 + j, p, p_prev, alpha, real(b));
                ptmp = p;
                p = p_prev;
                p_prev = ptmp;
            }
            g_derivative_polynomial(2 * r - 1 + j, p, p_prev, alpha, real(b));

            complex<double> derivative_at_1 = (complex<double>)0;
            for(int k = 0; k <= 2 * r - 1 + j; k++) {
                derivative_at_1 = derivative_at_1 + p[k];
            }
            derivative_at_1 *= EXP(2 * PI * I * (alpha + real(b)));
            complex<double> derivative_at_0 = p[0];

            N_power *= ((double)N * (double)N);
            complex<double> z = bernoulli_table[2 * r]/(factorial(2 * r) * N_power) * (derivative_at_1 - derivative_at_0);

            S = S - z;
            error = norm(z);
            r = r + 1;
            ptmp = p;
            p = p_prev;
            p_prev = ptmp;
        }
    }
    else if(real(b) == 0) {
        while(4 * error > epsilon * epsilon) {
            if(r > 1) {
                g_derivative_polynomial_I(2 * r - 2 + j, p, p_prev, alpha, imag(b));
                ptmp = p;
                p = p_prev;
                p_prev = ptmp;
            }
            g_derivative_polynomial_I(2 * r - 1 + j, p, p_prev, alpha, imag(b));

            complex<double> derivative_at_1 = (complex<double>)0;
            for(int k = 0; k <= 2 * r - 1 + j; k++) {
                derivative_at_1 = derivative_at_1 + p[k];
            }
            derivative_at_1 *= EXP(2 * PI * I * (alpha + b));
            complex<double> derivative_at_0 = p[0];

            N_power *= ((double)N * (double)N);
            complex<double> z = bernoulli_table[2 * r]/(factorial(2 * r) * N_power) * (derivative_at_1 - derivative_at_0);

            S = S - z;
            error = norm(z);
            r = r + 1;
            ptmp = p;
            p = p_prev;
            p_prev = ptmp;
        }
    }
    else {
        while(4 * error > epsilon * epsilon) {
            if(r > 1) {
                g_derivative_polynomial(2 * r - 2 + j, p, p_prev, alpha, b);
                ptmp = p;
                p = p_prev;
                p_prev = ptmp;
            }
            g_derivative_polynomial(2 * r - 1 + j, p, p_prev, alpha, b);

            complex<double> derivative_at_1 = (complex<double>)0;
            for(int k = 0; k <= 2 * r - 1 + j; k++) {
                derivative_at_1 = derivative_at_1 + p[k];
            }
            derivative_at_1 *= EXP(2 * PI * I * (alpha + b));
            complex<double> derivative_at_0 = p[0];

            N_power *= ((double)N * (double)N);
            complex<double> z = bernoulli_table[2 * r]/(factorial(2 * r) * N_power) * (derivative_at_1 - derivative_at_0);

            S = S - z;
            error = norm(z);
            r = r + 1;
            ptmp = p;
            p = p_prev;
            p_prev = ptmp;
        }
    }


    return S;

}

complex<double> G_via_Euler_MacLaurin_R(complex<double> alpha, double b, int n, int j, double epsilon) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    // using euler maclaurin summation
    //
    // This routine is specialized for b real.

    // By rough experimentation, we decide to use Euler MacLaurin summation in the cases that
    // (i) |alpha| < 1.5
    // OR (ii) 1.5 <= |alpha| <= 2.5 AND b > .024
    // OR (iii) 2.5 <= |alpha| AND (|alpha| - 2)/10 < b
    //
    // (These heuristics might not work well when b is larger, but b will be at most .25 in our applications.)

    //if(epsilon >= 1) {
    //    return (complex<double>)0;
    //}

    //int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-LOG(epsilon)/(2 * PI), 0.0) ) * (1 + j * log(n + 1)/4.0) )), 1);

//    int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) )), 1);

    int N = to_int(ceil(  ( abs(alpha) + abs(2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) ));
    N = 2 * max(N, 1);

    /*
    double two_n_to_the_j = 1;
    if(n != 0)
        two_n_to_the_j = pow(2 * n, j);
    else
        two_n_to_the_j = 1;
    double one_over_two_n_to_the_j = 1/two_n_to_the_j;
    */
    complex<double> alpha_term_multiplier = EXP(2 * PI * I * alpha/(double)N);
    complex<double> alpha_term = alpha_term_multiplier;

    complex<double> S = (complex<double>)0;
    {
        double t_increment = 1.0/N;
        double t = t_increment;
        switch(j) {
            case 0:
                for(int s = 1; s < N; s++) {
                    S = S + complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 1:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 2:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a*a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 3:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a*a*a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 4:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a*a*a*a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 5:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a*a*a*a*a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 6:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a*a*a*a*a*a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 7:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a*a*a*a*a*a*a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 8:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + a*a*a*a*a*a*a*a * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
   
            default:
                for(int s = 1; s < N; s++) {
                    double a = (n + t);
                    S = S + pow(a, j) * complex<double>(cos(2 * PI * b * t * t), sin(2 * PI * b * t * t)) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
        }
    }
    //S = S - (double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1)) * one_over_two_n_to_the_j;
    
    complex<double> exp_factor_at_one = alpha_term * complex<double>(cos(2 * PI * b), sin(2 * PI * b));
    //S = S - .5 * (pow(n, j) + pow(n + 1, j) * EXP(2 * PI * I * (alpha + b)));
    S = S + .5 * (pow(n, j) + pow(n + 1, j) * exp_factor_at_one);
    S = S/(double)N;

    double N_power = 1;

    double error = 10 * epsilon * epsilon;
    int r = 1;


    complex<double> p1[bernoulli_range + max_j];            // below p will need to be as big as 2r + j
    complex<double> p2[bernoulli_range + max_j];            // for increasing r. If 2r ever gets as large
                                                // as bernoulli_range, then the computation
                                                // will fail anyway, so the biggest p
                                                // that makes sense is of size bernoulli_range + j.
                                                //
                                                // As long as bernoulli_range isn't too big, this
                                                // should be hopefully be a good choice.

    complex<double> * ptmp;                             // A temporary variable used for swapping
                                                // p and p_prev.

    complex<double> * p = p1;
    complex<double> * p_prev = p2;

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
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            //p_prev[s] = pow(2, -j) * pow(n, -s) * binomial_coefficient(j, s);
            p_prev[j - s] = n_power * binomial_coefficient(j, s);
            n_power *= n;
        }
    }


    while(4 * error > epsilon * epsilon) {
        if(r > 1) {
            //g_derivative_polynomial_R(2 * r - 2 + j, p, p_prev, alpha, b);
            g_derivative_polynomial(2 * r - 2 + j, p, p_prev, alpha, b);
            ptmp = p;
            p = p_prev;
            p_prev = ptmp;
        }
        //g_derivative_polynomial_R(2 * r - 1 + j, p, p_prev, alpha, b);
        g_derivative_polynomial(2 * r - 1 + j, p, p_prev, alpha, b);

        complex<double> derivative_at_1 = (complex<double>)0;
        for(int k = 0; k <= 2 * r - 1 + j; k++) {
            derivative_at_1 = derivative_at_1 + p[k];
        }
        //derivative_at_1 *= exp(2 * PI * I * (alpha + b));
        derivative_at_1 *= exp_factor_at_one;
        complex<double> derivative_at_0 = p[0];

        N_power *= ((double)N * (double)N);
        complex<double> z = bernoulli_table[2 * r]/(factorial(2 * r) * N_power) * (derivative_at_1 - derivative_at_0);

        S = S - z;
        error = norm(z);
        r = r + 1;
        ptmp = p;
        p = p_prev;
        p_prev = ptmp;
    }


    return S;

}

complex<double> G_via_Euler_MacLaurin_I(complex<double> alpha, double b, int n, int j, double epsilon) {
    // I don't think that this function is ever used.
 
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    // using euler maclaurin summation
    //
    // This routine is specialized for b imaginary, and b has an implicit I in front of it.

    // By rough experimentation, we decide to use Euler MacLaurin summation in the cases that
    // (i) |alpha| < 1.5
    // OR (ii) 1.5 <= |alpha| <= 2.5 AND b > .024
    // OR (iii) 2.5 <= |alpha| AND (|alpha| - 2)/10 < b
    //
    // (These heuristics might not work well when b is larger, but b will be at most .25 in our applications.)

    //if(epsilon >= 1) {
    //    return (complex<double>)0;
    //}

    //int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-LOG(epsilon)/(2 * PI), 0.0) ) * (1 + j * log(n + 1)/4.0) )), 1);

//    int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) )), 1);

    int N = to_int(ceil(  ( abs(alpha) + abs((double)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) ));
    N = 2 * max(N, 1);


    /*
    double two_n_to_the_j = 1;
    if(n != 0)
        two_n_to_the_j = pow(2 * n, j);
    else
        two_n_to_the_j = 1;
    double one_over_two_n_to_the_j = 1/two_n_to_the_j;
    */

    complex<double> S = (complex<double>)0;
    {
        complex<double> alpha_term = 1.0;
        complex<double> alpha_term_multiplier = EXP(2 * PI * I * alpha/(double)N);
        double t = 0;
        double t_increment = 1.0/N;
        switch(j) {
            case 0:
                for(int s = 0; s <= N; s++) {
                    S = S + EXP(-2 * PI * b * t * t) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 1:
                for(int s = 0; s <= N; s++) {
                    S = S + (n + t) * EXP(-2 * PI * b * t * t) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;

            default:
                for(int s = 0; s <= N; s++) {
                    S = S + pow(t + n, j) * EXP(-2 * PI * b * t * t) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
        }
    }
    //S = S - (double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1)) * one_over_two_n_to_the_j;


    //S = S - (double)(.5) * (g(alpha, I * b, n, j, 0) + g(alpha, I * b, n, j, 1));
    S = S - .5 * (pow(n, j) + pow(1 + n, j) * EXP(2 * PI * I * (alpha + I * b)));
    S = S/(double)N;

    double N_power = 1;

    double error = 10 * epsilon * epsilon;
    int r = 1;


    complex<double> p1[bernoulli_range + max_j];            // below p will need to be as big as 2r + j
    complex<double> p2[bernoulli_range + max_j];            // for increasing r. If 2r ever gets as large
                                                // as bernoulli_range, then the computation
                                                // will fail anyway, so the biggest p
                                                // that makes sense is of size bernoulli_range + j.
                                                //
                                                // As long as bernoulli_range isn't too big, this
                                                // should be hopefully be a good choice.

    complex<double> * ptmp;                             // A temporary variable used for swapping
                                                // p and p_prev.

    complex<double> * p = p1;
    complex<double> * p_prev = p2;

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
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            //p_prev[s] = pow(2, -j) * pow(n, -s) * binomial_coefficient(j, s);
            p_prev[j - s] = n_power * binomial_coefficient(j, s);
            n_power *= n;
        }
    }


    while(4 * error > epsilon * epsilon) {
        if(r > 1) {
            g_derivative_polynomial_I(2 * r - 2 + j, p, p_prev, alpha, b);
            ptmp = p;
            p = p_prev;
            p_prev = ptmp;
        }
        g_derivative_polynomial_I(2 * r - 1 + j, p, p_prev, alpha, b);

        complex<double> derivative_at_1 = (complex<double>)0;
        for(int k = 0; k <= 2 * r - 1 + j; k++) {
            derivative_at_1 = derivative_at_1 + p[k];
        }
        derivative_at_1 *= EXP(2 * PI * I * (alpha + I * b));
        complex<double> derivative_at_0 = p[0];

        N_power *= ((double)N * (double)N);
        complex<double> z = bernoulli_table[2 * r]/(factorial(2 * r) * N_power) * (derivative_at_1 - derivative_at_0);

        S = S - z;
        error = norm(z);
        r = r + 1;
        ptmp = p;
        p = p_prev;
        p_prev = ptmp;
    }


    return S;

}

complex<double> G_via_Euler_MacLaurin_I_over_twopi(complex<double> alpha, int n, int j, double epsilon) {
    // We compute the integral int_0^1 g(t),
    // where g(t) = (t + n)^j exp(2 pi i alpha t + 2 pi i b t^2)
    // using euler maclaurin summation
    //
    // Specialized for b = I/2pi

    // By rough experimentation, we decide to use Euler MacLaurin summation in the cases that
    // (i) |alpha| < 1.5
    // OR (ii) 1.5 <= |alpha| <= 2.5 AND b > .024
    // OR (iii) 2.5 <= |alpha| AND (|alpha| - 2)/10 < b
    //
    // (These heuristics might not work well when b is larger, but b will be at most .25 in our applications.)

    //if(epsilon >= 1) {
    //    return (complex<double>)0;
    //}

    //int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-LOG(epsilon)/(2 * PI), 0.0) ) * (1 + j * log(n + 1)/4.0) )), 1);

//    int N = 4 * max(to_int(ceil(  ( abs(alpha) + abs((complex<double>)2.0 * b) + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) )), 1);

    int N = to_int(ceil(  ( abs(alpha) + 1/PI + max(-fastlog(epsilon)/(2 * PI), 0.0) ) * (1 + j * (fastlog(n + 1) + 1)/4.0) ));
    N = 2 * max(N, 1);


    /*
    double two_n_to_the_j = 1;
    if(n != 0)
        two_n_to_the_j = pow(2 * n, j);
    else
        two_n_to_the_j = 1;
    double one_over_two_n_to_the_j = 1/two_n_to_the_j;
    */

    complex<double> alpha_term_multiplier = EXP(2 * PI * I * alpha/(double)N);
    complex<double> alpha_term = alpha_term_multiplier;
    complex<double> S = (complex<double>)0;

    {
        double t_increment = 1.0/N;
        double t = t_increment;
        switch(j) {
            case 0:
                for(int s = 1; s < N; s++) {
                    //S = S + exp(-t * t) * alpha_term;
                    S = S + exp_t_over_N_squared(s, N) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;
            case 1:
                for(int s = 1; s < N; s++) {
                    //S = S + (n + t) * exp(-t * t) * alpha_term;
                    S = S + (n + t) * exp_t_over_N_squared(s, N) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
                break;

            default:
                for(int s = 1; s < N; s++) {
                    //S = S + pow(t + n, j) * exp(-t * t) * alpha_term;
                    S = S + pow(t + n, j) * exp_t_over_N_squared(s, N) * alpha_term;
                    alpha_term *= alpha_term_multiplier;
                    t += t_increment;
                }
        }
    }
    //S = S - (double)(.5) * (g(alpha, b, n, j, 0) + g(alpha, b, n, j, 1)) * one_over_two_n_to_the_j;
    //S = S - (double)(.5) * (g(alpha, I/(2 * PI), n, j, 0) + g(alpha, I/(2 * PI), n, j, 1));
    
    //complex<double> exp_factor_at_1 = exp(2 * PI * I * alpha - 1.0);
    complex<double> exp_factor_at_1 = alpha_term/E;
    S = S + .5 * (pow(n, j) + pow(1 + n, j) * exp_factor_at_1);
    S = S/(double)N;


    double error = 10 * epsilon * epsilon;
    int r = 1;


    complex<double> p1[bernoulli_range + max_j];            // below p will need to be as big as 2r + j
    complex<double> p2[bernoulli_range + max_j];            // for increasing r. If 2r ever gets as large
                                                // as bernoulli_range, then the computation
                                                // will fail anyway, so the biggest p
                                                // that makes sense is of size bernoulli_range + j.
                                                //
                                                // As long as bernoulli_range isn't too big, this
                                                // should be hopefully be a good choice.

    complex<double> * ptmp;                             // A temporary variable used for swapping
                                                // p and p_prev.

    complex<double> * p = p1;
    complex<double> * p_prev = p2;

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
        double n_power = 1.0;
        for(int s = 0; s <= j; s++) {
            //p_prev[s] = pow(2, -j) * pow(n, -s) * binomial_coefficient(j, s);
            p_prev[j - s] = n_power * binomial_coefficient(j, s);
            n_power *= n;
        }
    }

    double N_power = 1;
    double N_power_multiplier = 1.0/(N * N);

    while(4 * error > epsilon * epsilon) {
        if(r > 1) {
            g_derivative_polynomial_I_over_twopi(2 * r - 2 + j, p, p_prev, alpha);
            ptmp = p;
            p = p_prev;
            p_prev = ptmp;
        }
        g_derivative_polynomial_I_over_twopi(2 * r - 1 + j, p, p_prev, alpha);

        complex<double> derivative_at_1 = (complex<double>)0;
        for(int k = 0; k <= 2 * r - 1 + j; k++) {
            derivative_at_1 = derivative_at_1 + p[k];
        }
        //derivative_at_1 *= exp(2 * PI * I * alpha - 1.0);
        derivative_at_1 *= exp_factor_at_1;
        complex<double> derivative_at_0 = p[0];

        N_power *= N_power_multiplier;
        //complex<double> z = N_power * bernoulli_table[2 * r]/factorial(2 * r) * (derivative_at_1 - derivative_at_0);
        complex<double> z = N_power * bernoulli_over_factorial(2*r) * (derivative_at_1 - derivative_at_0);

        S = S - z;
        error = norm(z);
        r = r + 1;
        ptmp = p;
        p = p_prev;
        p_prev = ptmp;
    }

    return S;

}

