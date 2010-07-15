#include "theta_sums.h"
#include "precomputed_tables.h"
#include "log.h"

#include <iostream>
#include <cmath>

using namespace std;

static struct{
    long number_of_a;
    long number_of_b;
    long max_j;
    long max_M;
    double a_spacing;
    double b_spacing;
    long a_per_unit_interval;
    long b_per_unit_interval;
    Complex **** values;
} F0_cache;

inline Complex get_cached_F0_value(int a_index, int b_index, int j, int M) {
    if(M == -1) {
        if(a_index <= (F0_cache.number_of_a  - 1 )* F0_cache.max_M  + 1 && b_index < F0_cache.number_of_b && j <= F0_cache.max_j) {
            if(FAKE_PRECOMPUTATION) {
                return 1.0;
            }
            else{
                return F0_cache.values[0][a_index][b_index][j];
            }
        }
        else {
            cout << "F0 called for a value out of range. Exiting." << endl;
            cout << "F0 called with: ";
            cout << "    a_index = " << a_index << endl;
            cout << "    b_index = " << b_index << endl;
            cout << "          j = " << j << endl;
            cout << "          M = " << M << endl;
            exit(1);
            //return 0.0/0.0;
            //Double a = 1.0/(F0_cache.number_of_a - 1.0) * a_index;
            //Double b = .25/(F0_cache.number_of_b - 1.0) * b_index;
            //return J_Integral_0(a, b, j, M, 1, NULL, exp(-20), false);
        }
    }

    if(a_index < F0_cache.number_of_a && b_index < F0_cache.number_of_b && j <= F0_cache.max_j && M <= F0_cache.max_M) {
        if(FAKE_PRECOMPUTATION) {
            return 0.0;
        }
        else {
            return F0_cache.values[M][a_index][b_index][j];
        }
    }
    else {
        cout << "F0 called for a value out of range. Exiting." << endl;
        cout << "F0 called with: ";
        cout << "    a_index = " << a_index << endl;
        cout << "    b_index = " << b_index << endl;
        cout << "          j = " << j << endl;
        cout << "          M = " << M << endl;
        exit(1);
        //return 0.0/0.0;
        //Double a = 1.0/(F0_cache.number_of_a - 1.0) * a_index;
        //Double b = .25/(F0_cache.number_of_b - 1.0) * b_index;
//        cout << "a = " << a << "  b = " << b << endl;
        //return J_Integral_0(a, b, j, M, 1, NULL, exp(-20), false);
    }
}

static bool F0_cache_initialized = false;

static struct{
    long number_of_a;
    long number_of_b;
    long max_j;
    double a_spacing;
    double b_spacing;
    long a_per_unit_interval;
    long b_per_unit_interval;
    Complex *** values;
} F1_cache;

static bool F1_cache_initialized = false;


static struct{
    long max_a1;
    long number_of_a1;
    long number_of_a2;
    long number_of_b;
    long max_j;
    double a1_spacing;
    double a2_spacing;
    double b_spacing;
    long a1_per_unit_interval;
    long a2_per_unit_interval;
    long b_per_unit_interval;
    Complex *** values;
} F2_cache;


static bool F2_cache_initialized = false;


inline Complex get_cached_F1_value(int a_index, int b_index, int j) {
    if(a_index < F1_cache.number_of_a && b_index < F1_cache.number_of_b && j <= F1_cache.max_j) {
        if(FAKE_PRECOMPUTATION) {
            return 0.0;
        }
        else {
            return F1_cache.values[a_index][b_index][j];
        }
    }
    else {
        cout << "Warning: requested F1 value out of range. a_index = " << a_index << ", b_index = " << b_index << ", j = " << j << endl;
        return 0.0/0.0;
//        Double a = 6.0/(F0_cache.number_of_a - 1.0) * a_index;
//        Double b = .25/(F0_cache.number_of_b - 1.0) * b_index;
//        cout << "a = " << a << "  b = " << b << endl;
//        return J_Integral_1(a, b, j, 10, 1, NULL, exp(-20), false);
    }
}

inline Complex get_cached_F2_value(int a1_index, int a2_index, int b_index, int j) {
    if(a1_index < (F2_cache.number_of_a1 - 1) * F2_cache.max_a1 + 1 && a2_index < F2_cache.number_of_a2 && b_index < F2_cache.number_of_b) {
        if(FAKE_PRECOMPUTATION) {
            return 0.0;
        }
        else {
            if(j == 0) {
                return F2_cache.values[a1_index][a2_index][b_index];
            }
            else {
                return F0_cache.values[0][a1_index][b_index][j] - F0_cache.values[0][a2_index][b_index][j];
            }
        }
    }
    else {
        cout << "Warning. F2 value not precomputed." << endl;
        cout << "Asked for a1_index = " << a1_index << " a2_index = " << a2_index << " b_index = " << b_index << " j = " << j << endl;
        return 0.0/0.0;
    }
}


Complex H_Integral_0(int j, Double a, int M, Double epsilon) {
    //
    // Compute the integral int_0^1 t^j exp(-2 pi a t) (1 - exp(-2 pi M t))/(exp(2 pi t) - 1)dt,
    // which can be written (by expanding the geometric sum) as
    //
    // sum_{m=1}^M int_0^1 t^j exp(-2 pi (a + m) t) dt
    //
    //
    // In terms of our function H(j, a, epsilon), this is
    //
    // sum_{m=1}^M H(j, a + m, epsilon/(2 j + 1)
    //
    // We don't want to compute all of these terms directly, however, because
    // M can be large. We separate this into the case of small m, where
    // we compute directly, and then large m, where the antiderivative
    // can be approximated simply. We compute "simple" approximation
    // with Euler-Maclaurin summation.
    //
    // small m is just m < C = max(j, ceil(-log(epsilon)/2pi)
    //

    stats::H_Integral_0++;

    Complex S = (Complex)0;

    //int C = max(to_int(ceil(j/(2 * PI * E))), to_int(ceil(-LOG(epsilon)/(2 * PI))) );
    int C = max(to_int(ceil(j/(2 * PI * E))), to_int(ceil(-fastlog(epsilon)/(2 * PI))) );
    //int C = ceil(max(j/E, -(Double)(fastlog(epsilon)))/2 * PI) ;

    if(M != -1) {
        C = min(M, C);
    }

    for(int m = 1; m <= C; m++) {
        Complex z = H(j, a + m, epsilon/(C + 1));
        S = S + z;
    }

    if(C == M) {
        return S;
    }
    
    Double z = factorial(j)/two_pi_power(j + 1);

    S = S + z * sum_of_offset_inverse_powers(a, C + 1, M, j+1, epsilon/z);

    return S;

}

Complex J_Integral_0(Double a, Double b, int j, int M, int K, theta_cache * cache, Double epsilon, bool use_cache) {
    //
    // Compute the integral int_0^1 exp(-2 pi a t - 2 pi i b t^2) (1 - exp(-2 pi M t))/(exp(2 pi t) - 1) dt
    //
    // We can do this by expanding the denominator in a geometric series
    // and doing a taylor expansion in b, which puts this integral in the form
    //
    // sum_{r=0}^\infty ((2 pi i b)^r)/r! \sum_{m=1}^M int_0^1 t^j exp(-2 pi(m + a)t)dt
    //
    // -1 corresponds to infinity
    stats::J_Integral_0++;
    
    return 0.0;

    Double K_pow_j;
    Double K_pow_minus_j;
    if(cache) {
        K_pow_j = K_power(j, cache);
        K_pow_minus_j = K_power(-j, cache);
    }
    else {
        K_pow_j = pow(K, j);
        K_pow_minus_j = 1.0/K_pow_j;
    }

    

    //if(2 * pow(K, -j) * log(abs((Double)M) + 10) < epsilon) {
    if( 3 + fastlog2(K_pow_minus_j) + fastlog2(fastlog(abs(M) + 10.0)) < fastlog2(epsilon)) {
        stats::J_Integral_0_zero++;
        return 0.0;
    }
    epsilon = epsilon * K_pow_j;

    if(use_cache && F0_cache_initialized) {
        if(stats::stats)
            stats::J_Integral_0_taylor_expansion++;

        int a0_index = round( (F0_cache.number_of_a - 1) * a );
        Double a0 = a0_index * 1.0/(F0_cache.number_of_a - 1);
        int b0_index = round( (F0_cache.number_of_b - 1) * b/.25);
        Double b0 = b0_index * .25/(F0_cache.number_of_b - 1);

        Double outside_error = 2 * epsilon;
        Double b_minus_b0_power = 1.0;

        int s = 0;
        Complex S = 0;

//        cout << "a = " << a << ",  b = " << b << endl;
//        cout << "a0 = " << a0 << ", b0 = " << b0 << endl;
//        cout << "M = " << M << endl;

        int number_of_terms = 0;

        while(outside_error > epsilon/2) {
            Double a_minus_a0_power = 1.0;

            Complex S2 = 0;
            int l = 0;

            Complex z1 = two_pi_over_factorial_power(s) * minus_I_power(s) * b_minus_b0_power;
            outside_error = abs(z1);
            
            Double inside_error = 2 * epsilon;
            while(l == 0 || inside_error *  outside_error > epsilon/6) {
                //Complex z = minus_one_power(l) * two_pi_over_factorial_power(l) * a_minus_a0_power * F0_cache.values[M][a0_index][b0_index][j + l + 2 * s];
                //Complex z = minus_one_power(l) * two_pi_over_factorial_power(l) * a_minus_a0_power * get_cached_F0_value(a0_index, b0_index, j + l + 2 * s, M);
                Complex z = minus_one_power(l) * two_pi_over_factorial_power(l) * a_minus_a0_power;
                inside_error = abs(z);
                z *= get_cached_F0_value(a0_index, b0_index, j + l + 2 * s, M);
                S2 = S2 + z;
                a_minus_a0_power *= (a - a0);
                l++;
                if(stats::stats)
                    number_of_terms++;
            }

            z1 = z1 * S2;
            S = S + z1;
            s++;
            b_minus_b0_power *= (b - b0);
        }

        if(stats::stats)
            stats::J_Integral_0_terms_used += number_of_terms;

        Complex answer = K_pow_minus_j * S;
        
        if(verbose::J_Integral_0)
            cout << "Computed J_Integral_0( " << a << ", " << b << ", " << j << ", " << M << ", " << K << ") = " << answer << endl;

        return answer;
    }

    Complex S = (Complex)0;

    int r = 0;
    Double error = 2 * epsilon;
    
    if(b == 0) {
        Complex answer = H_Integral_0(j, a, M, epsilon)*K_pow_minus_j;
        if(verbose::J_Integral_0)
            cout << "Computed J_Integral_0( " << a << ", " << b << ", " << j << ", " << M << ", " << K << ") = " << answer << endl;

        return answer;
    }
    
    int N = max(1, to_int(ceil(-fastlog(epsilon))));    // an approximation of the number of terms
                                                    // we compute in the taylor expansion

    Double b_power = (Double)1.0;
    while(error > epsilon) {
        Complex z = minus_I_power(r) * two_pi_over_factorial_power(r) * b_power;
        error = abs(z);
        z *= H_Integral_0(2 * r + j, a, M, epsilon/(error *N));  // TODO: We really don't need to compute this term to
                                                                                // this much precision usually. We really should figure
                                                                                // our how many terms we are going to compute
                                                                                // and then divide appropriately.
                                                                                //
        //cout << "error = " << error << "; z = " << z << endl;

        S = S + z;
        b_power *= b;
        r++;
    }

    S = S * K_pow_minus_j;

    if(verbose::J_Integral_0)
        cout << "Computed J_Integral_0( " << a << ", " << b << ", " << j << ", " << M << ", " << K << ") = " << S << endl;

    return S;
}

Complex J_Integral_1(Double a, Double b, int j, int M, int K, theta_cache * cache, Double epsilon, bool use_cache) {
    //
    // Compute the integral int_1^K exp(-2 pi a t - 2 pi i b t^2)(1 - exp(-2 pi M t))/(exp(2 pi t) - 1) dt
    //
    // If M == -1, then we treat it as positive infinity.
    //
    
    // We truncate the integral at L, where L is given by
    stats::J_Integral_1++;

    return 0.0;

    //int L = min(K, max(1, to_int(ceil(-LOG(epsilon)/(2 * PI * (1 + a))))));

    int L = max(1.0, ceil(-fastlog(epsilon)/(2 * PI * (1 + a))));
    L = max(1.0, L - j * fastlog( (Double)K/(Double)(L + 1))/(2 * PI) );
    L = min(K, L);

    if(verbose::J_Integral_1 > 1) {
        cout << "J_Integral_1 called with: " << endl;
        cout << "                           a = " << a << endl;
        cout << "                           b = " << b << endl;
        cout << "                           j = " << j << endl;
        cout << "                           M = " << M << endl;
        cout << "                           K = " << K << endl;
        cout << "   computed L = " << L << endl;
    }
    else if(verbose::J_Integral_1) {
        cout << "j, L = " << j << ", " << L << endl;
    }

    if(L <= 1) {
        if(stats::stats)
            stats::J_Integral_1_zero++;
        return 0.0;
    }



    // Now we compute the integral as a sum over unit length integrals.
    // On each unit length interval, we do a change of variables to get
    // the range of integration between 0 and 1, and reduce the computation
    // to a sum of exponentials times G()

    Complex S = (Complex)0;

    //Double one_over_K_to_the_j = pow( (Double)K, -j);
    Double one_over_K_to_the_j;
    if(cache)
        one_over_K_to_the_j = K_power(-j, cache);
    else
        one_over_K_to_the_j = pow(K, -j);

    //if(2 * pow(L, j) * one_over_K_to_the_j * log(abs((Double)M) + 10) < epsilon) {
    if( 5 + j * fastlog2(L) + fastlog2(one_over_K_to_the_j) - 2 * PI * (a + 1) * log2(E) < fastlog2(epsilon)) {
        if(stats::stats)
            stats::J_Integral_1_zero++;
        return 0.0;
    }


    if(a <= 6 && use_cache && F1_cache_initialized && (M > ceil(-fastlog(epsilon)/(PI)) || M == -1)) {

        if(stats::stats)
            stats::J_Integral_1_taylor_expansion++;

        epsilon /= one_over_K_to_the_j;
        
        int a0_index = round( (F1_cache.number_of_a - 1) * a/6.0 );
        Double a0 = a0_index * 6.0/(F1_cache.number_of_a - 1);
        int b0_index = round( (F1_cache.number_of_b - 1) * b/.25);
        Double b0 = b0_index * .25/(F1_cache.number_of_b - 1);

        Double outside_error = 2 * epsilon;
        Double b_minus_b0_power = 1.0;

        int s = 0;
        Complex S = 0;

//        cout << "a = " << a << ",  b = " << b << endl;
//        cout << "a0 = " << a0 << ", b0 = " << b0 << endl;
//        cout << "M = " << M << endl;

        Double outside_a_plus_1_power = pow(a + 1, -j);
        Double a_plus_1_inverse = 1.0/(a + 1.0);

        int number_of_terms = 0;

        //while(s == 0 || L_power1 * outside_error > epsilon/2) {
        while(s == 0 ||  outside_a_plus_1_power/(two_pi_over_factorial_power(j + 2 *  s) * 2 * PI) * outside_error > epsilon/2) {
        //while(L_power1 * outside_error > epsilon/100) {
            Double a_minus_a0_power = 1.0;


            Complex z1 = two_pi_over_factorial_power(s) * minus_I_power(s) * b_minus_b0_power;
            outside_error = abs(z1);
            
            Double inside_error = 0;
            
            int l = 0;
            Complex S2 = 0;
            Double inside_a_plus_1_power = outside_a_plus_1_power;
            while(l == 0 || inside_error *  outside_error * inside_a_plus_1_power/(two_pi_over_factorial_power(j + 2*s + l) * 2 * PI) > epsilon/6) {
            //while(inside_error *  outside_error * L_power2 > epsilon/300) {
                //Complex z = minus_one_power(l) * two_pi_over_factorial_power(l) * a_minus_a0_power * F0_cache.values[M][a0_index][b0_index][j + l + 2 * s];
                //Complex z = minus_one_power(l) * two_pi_over_factorial_power(l) * a_minus_a0_power * get_cached_F0_value(a0_index, b0_index, j + l + 2 * s, M);
                Complex z = minus_one_power(l) * two_pi_over_factorial_power(l) * a_minus_a0_power;
                inside_error = abs(z);
                z *= get_cached_F1_value(a0_index, b0_index, j + l + 2 * s);
                S2 = S2 + z;
                a_minus_a0_power *= (a - a0);
                l++;
                inside_a_plus_1_power *= a_plus_1_inverse;
                if(stats::stats)
                    number_of_terms++;
            }


            outside_a_plus_1_power *= (a_plus_1_inverse * a_plus_1_inverse);
            z1 = z1 * S2;
            S = S + z1;
            s++;
            b_minus_b0_power *= (b - b0);
        }

        if(stats::stats)
            stats::J_Integral_1_terms_used += number_of_terms;

        //cout << "number of terms computed = " << number_of_terms << endl;

        //Complex answer = S * one_over_K_to_the_j;
        //Complex direct_computed_answer = J_Integral_1(a, b, j, M, K, NULL, epsilon, false);
        //cout << "For J_Integral_1(" << a << ", " << b << ", " << j << ", " << M << ", " << K << "):" << endl;
        //cout << "Using taylor expansion returning " << answer << endl;
        //cout << "Directly computing, get " << direct_computed_answer << endl;
        
        Complex answer = S * one_over_K_to_the_j;

        if(verbose::J_Integral_1 > 1)
            cout << "Computed J_Integral_1( " << a << ", " << b << ", " << j << ", " << M << ", " << K << ", " << epsilon << ") = " << answer << endl;
        
        return answer;

    }

    Double exp_minus_twopi = exp(-2.0 * PI);
    Double exp_minus_twopi_n = 1.0;

    //Double a_factor = exp(-2.0 * PI * (1 + a));
    //Double a_multiplier = a_factor;

    for(Double n = (Double)1; n <= L - 1; n = n + 1) {
        exp_minus_twopi_n *= exp_minus_twopi;
        int end_point;
        if(M == -1) {
            end_point = to_int(ceil(-fastlog(epsilon/(2 * L))/(2 * PI * n) -  j * fastlog((Double)K/(n + 1))/(2 * PI) ));
        }
        else {
            end_point = min(M, to_int(ceil(-fastlog(epsilon/(2 * L))/(2 * PI * n) - j * fastlog((Double)K/(n + 1))/(2 * PI))  )  );
            //end_point = min(M, to_int(ceil(-log(epsilon/(2 * L))/(2 * PI * n)    )    )      );
        }

        end_point = max(1, end_point);

        //end_point = max(1, (int)ceil(end_point - j * log((Double)K/(n + 1))/(2 * PI)));
        //end_point = max(1, (int)(end_point - j * fastlog((Double)K/(n + 1))/(2 * PI)));

        Complex S1 = 0;

        Complex x = EXP(-2.0 * PI * n * (1.0 + a + I * b * n));
        //Double d = -2 * PI * b * n * n;
        //Complex x = a_factor * Complex(cos(d), sin(d));

        for(Double m = (Double)1; m <= end_point; m = m + 1) {
            if(m > 1)
                x = x * exp_minus_twopi_n;
            //Complex z =  G(I*(m + a + (Complex)2.0 * I * b * n), -b, n, j, epsilon/(abs(x) * end_point * Double(L - 1) * one_over_K_to_the_j));
            Complex z =  G_R(I*(m + a + (Complex)2.0 * I * b * n), -b, n, j, epsilon/(abs(x) * end_point * Double(L - 1) * one_over_K_to_the_j));
            z *= x;
            S1 = S1 + z;
        }
        S1 = S1 * one_over_K_to_the_j;
        S = S + S1;
        //a_factor *= a_multiplier;
    }


    if(verbose::J_Integral_1 > 1)
        cout << "Computed J_Integral_1( " << a << ", " << b << ", " << j << ", " << M << ", " << K << ", " << epsilon << ") = " << S << endl;


    return S;

}


Complex J_Integral_1_precomputation(Double a, Double b, int j, Double epsilon) {
    //
    // Compute the integral int_1^K exp(-2 pi a t - 2 pi i b t^2)(1 - exp(-2 pi M t))/(exp(2 pi t) - 1) dt
    //
    // If M == -1, then we treat it as positive infinity.
    //
    
    int M = -1;

    // We truncate the integral at L, where L is given by
    //int L = min(K, max(1, to_int(ceil(-LOG(epsilon)/(2 * PI * (1 + a))))));

    int L = max(1.0, ceil(-fastlog(epsilon)/(2 * PI * (1 + a))));

    L = (j + 1) * L;

    if(L <= 1) {
        return 0.0;
    }



    // Now we compute the integral as a sum over unit length integrals.
    // On each unit length interval, we do a change of variables to get
    // the range of integration between 0 and 1, and reduce the computation
    // to a sum of exponentials times G()

    Complex S = (Complex)0;

    if( 5 + j * fastlog2(L) - 2 * PI * (a + 1) * log2(E) < fastlog2(epsilon)) {
        return 0.0;
    }

    Double exp_minus_twopi = exp(-2.0 * PI);
    Double exp_minus_twopi_n = 1.0;

    //Double a_factor = exp(-2.0 * PI * (1 + a));
    //Double a_multiplier = a_factor;

    //bool a_is_small;
    //if(2 * PI * (1 + a) * E > j) {
    //    a_is_small = false;
    //}
    //else {
    //    a_is_small = true;
    //}

    for(Double n = (Double)1; n <= L - 1; n = n + 1) {
        exp_minus_twopi_n *= exp_minus_twopi;
        int end_point;
        if(M == -1) {
            end_point = (j + 1) * to_int(ceil(-fastlog(epsilon)/(2 * PI * n)));
        }
        else {
            end_point = min(M, to_int(ceil(-fastlog(epsilon)/(2 * PI * n) ) ));
        }

        end_point = max(end_point, 1);

        Complex S1 = 0;

        Complex x = EXP(-2.0 * PI * n * (1.0 + a + I * b * n));

        for(Double m = (Double)1; m <= end_point; m = m + 1) {
            if(m > 1)
                x = x * exp_minus_twopi_n;
            //Complex z =  G(I*(m + a + (Complex)2.0 * I * b * n), -b, n, j, epsilon/(abs(x) * end_point * Double(L - 1) * one_over_K_to_the_j));
            Complex z =  G_R(I*(m + a + (Complex)2.0 * I * b * n), -b, n, j, epsilon/(abs(x) * end_point * Double(L - 1)));
            z *= x;
            S1 = S1 + z;
        }
        S = S + S1;
        //a_factor *= a_multiplier;
        
        //if(a_is_small) {
        //    if( j * fastlog2(j * n) - 2 * PI * (1 + a) * L * log2(E) < fastlog2(epsilon) )
        //        break;
        //}
        //else {
        //    if( j * fastlog2(n) - 2 * PI * (1 + a) * L * log2(E) < fastlog2(epsilon) )
        //        break;
        //}
    }

    return S;

}






Complex H_Integral_2(int j, Double a1, Double a2, Double epsilon) {
    stats::H_Integral_2++;

    Complex S = (Complex)0;

    int C = max(to_int(ceil(j/(2*PI*E))), to_int(ceil(-fastlog(epsilon)/(2 * PI) ) ));

    for(int m = 1; m <= C; m++) {
        Complex z = H(j, m + a1, epsilon/(C + 1)) - H(j, m + a2, epsilon/(C + 1));
        S = S + z;
    }

    S = S + factorial(j)/two_pi_power(j + 1) * infinite_sum_of_differenced_inverse_powers(a1, a2, C + 1, j+1, epsilon);

    return S;

}

Complex J_Integral_2(Double a1, Double a2, Double b, theta_cache * cache, Double epsilon, bool use_cache) {
    //
    // Compute the integral int_0^1 exp(-2 pi i b t^2) (exp(-2 pi a1 t) - exp(-2 pi a2 t))(exp(2 pi t) - 1) dt,
    //
    // which is equal to 
    //
    // lim_{M --> oo} J_Integral_0(a1, b, M) - J_Integral_0(a2, b, M)
    //
    stats::J_Integral_2++;

    return 0.0;

    if(a2 > a1) {
        //if (2.0 * ( (a2 - a1)/a2 * log(1 + a2) ) < epsilon) {
        if(1 + fastlog2( (a2 - a1)/a2 * (fastlog(1 + a2) + 1.0)) < fastlog2(epsilon)) {
            stats::J_Integral_2_zero++;
            return 0.0;
        }
    }
    else {
        //if (2.0 * ( (a1 - a2)/a1 * log(1 + a1) ) < epsilon) {
        if(1 + fastlog2( (a1 - a2)/a1 * (fastlog(1 + a1) + 1.0)) < fastlog2(epsilon)) {
            stats::J_Integral_2_zero++;
            return 0.0;
        }
    }

    if(use_cache && F2_cache_initialized) {

        if(stats::stats)
            stats::J_Integral_2_taylor_expansion++;

        int sign = 1;
        if(a2 > a1) {
            Double tmp = a1;
            a1 = a2;
            a2 = tmp;
            sign = -1;
        }


        int a1_0_index = round( (F2_cache.number_of_a1 - 1) * a1);
        Double a1_0 = (Double)a1_0_index /(Double)(F2_cache.number_of_a1 - 1);

        int a2_0_index = round( (F2_cache.number_of_a2 - 1) * a2);
        Double a2_0 = a2_0_index * 1.0/(F2_cache.number_of_a2 - 1);

        int b0_index = round( (F0_cache.number_of_b - 1) * b/.25);
        Double b0 = b0_index * .25/(F2_cache.number_of_b - 1);

        if(verbose::J_Integral_2) {
            cout << "With a1 = " << a1 << ", a2 = " << a2 << ", b = " << b << endl;
            cout << "    using a1_0 = " << a1_0 << endl;
            cout << "    using a2_0 = " << a2_0 << endl;
            cout << "    using b0 = " << b0 << endl;
            cout << "  a1_0 index = " << a1_0_index << endl;
            cout << "  a2_0 index = " << a2_0_index << endl;
            cout << "    b0 index = " << b0_index << endl;
        }

        Complex S = 0;

        int number_of_terms = 1;

        S = S + get_cached_F2_value(a1_0_index, a2_0_index, b0_index, 0);
        //S = S + J_Integral_2(a1_0, a2_0, b0, NULL, epsilon, 0);

        Complex S1 = 0;
        Double error = 2 * epsilon;

        {
            int r = 1;
            Double a2_minus_a2_0_power = (a2 - a2_0);
            while(error > epsilon/3) {
                Double z = minus_one_power(r) * two_pi_over_factorial_power(r) * a2_minus_a2_0_power;
                error = abs(z);
                S1 = S1 + z * get_cached_F0_value(a2_0_index, b0_index, r, -1);
                a2_minus_a2_0_power *= (a2 - a2_0);
                r++;
                if(stats::stats)
                    number_of_terms++;
            }
            if(verbose::J_Integral_2)
                cout << r << " terms in first loop." << endl;
        }

        S = S - S1;

        S1 = 0.0;
        {
            Double a1_minus_a1_0_power = (a1 - a1_0);
            int l = 1;
            error = 2 * epsilon;
            while(error > epsilon/3) {
                Double z = minus_one_power(l) * two_pi_over_factorial_power(l) * a1_minus_a1_0_power;
                error = abs(z);
                S1 = S1 + z * get_cached_F0_value(a1_0_index, b0_index, l, -1);
                a1_minus_a1_0_power *= (a1 - a1_0);
                l++;
                if(stats::stats)
                    number_of_terms++;
            }
            if(verbose::J_Integral_2)
                cout << l << " terms in second loop." << endl;
        }

        S = S + S1;

        S1 = 0;
        {
            int s = 1;
            Double b_minus_b0_power = b - b0;
            error = 2 * epsilon;
            while(error > epsilon/3) {
                Double z = two_pi_over_factorial_power(s) * b_minus_b0_power;
                error = abs(z);
                S1 = S1 + z * minus_I_power(s) * (J_Integral_0(a1, b0, 2 * s, -1, 1, NULL, epsilon/error) - J_Integral_0(a2, b0, 2 * s, -1, 1, NULL, epsilon/error));
                b_minus_b0_power *= b - b0;
                s++;
                if(stats::stats)
                    number_of_terms++;
            }
            if(verbose::J_Integral_2)
                cout << s << " terms in third loop." << endl;
        }

        if(stats::stats)
            stats::J_Integral_2_terms_used += number_of_terms;

        S = S + S1;

        Complex answer = (Double)sign * S;

        if(verbose::J_Integral_2)
            cout << "Computed J_Integral_2( " << a1 << ", " << a2 << ", " << b << ") = " << answer << endl;

        //cout << "                   Computed J_Integral_2( " << a1 << ", " << a2 << ", " << b << ") = " << answer << endl;
        //Complex answer2 = J_Integral_2(a1, a2, b, NULL, epsilon, 0);
        //cout << "should have computed " << answer2 << endl;


        return answer;
    }

    Complex S = (Complex)0;

    int r = 0;
    Double error = 2 * epsilon;
    int sign = -1;
    
    if(b == 0) {
        Complex answer = H_Integral_2(0, a1, a2, epsilon);
        if(verbose::J_Integral_2)
            cout << "Computed J_Integral_2( " << a1 << ", " << a2 << ", " << b << ") = " << answer << endl;
        return answer;
    }
   
    int N = max(1, to_int(-fastlog(epsilon)));

    Complex Ib_power = (Complex)1.0/(I * b);
    while(error > epsilon) {
        sign = -sign;
        Ib_power *= (I * b);
        Complex z = sign * two_pi_over_factorial_power(r) * Ib_power;
        error = abs(z);
        
        z *= H_Integral_2(2 * r, a1, a2, epsilon/(error * N));  // TODO: We really don't need to compute this term to
                                                                                // this much precision usually. We really should figure
                                                                                // our how many terms we are going to compute
                                                                                // and then divide appropriately.
        S = S + z;
        r++;
    }


    if(verbose::J_Integral_2)
        cout << "Computed J_Integral_2( " << a1 << ", " << a2 << ", " << b << ") = " << S << endl;

    return S;

}

/*
void build_F0_cache(long number_of_a, long number_of_b, long max_j, long max_M, Double epsilon) {

    F0_cache.number_of_a = number_of_a;
    F0_cache.number_of_b = number_of_b;
    F0_cache.max_j = max_j;
    F0_cache.max_M = max_M;

    F0_cache.values = new Complex *** [number_of_a];

    Double a = 0;
    Double a_increment = 1.0/(number_of_a - 1);
    Double b = 0;
    Double b_increment = .25/(number_of_b - 1);

    for(long k = 0; k < number_of_a; k++) {
        cout << "a index = " << k << ", a = " << a << endl;
        F0_cache.values[k] = new Complex ** [number_of_b];
        for(long l = 0; l < number_of_b; l++) {
            cout << "b index = " << l << ", b = " << b << endl;
            F0_cache.values[k][l] = new Complex * [max_M + 1];
            F0_cache.values[k][l][0] = new Complex[max_j + 1];
            for(long j = 1; j <= max_j; j++) {
                    F0_cache.values[k][l][0][j] = J_Integral_0(a, b, j, -1, 1, NULL, epsilon);
            }
            for(long M = 1; M <= max_M; M++) {
//                cout << "M = " << M << endl;
                F0_cache.values[k][l][M] = new Complex [max_j + 1];
                for(long j = 0; j <= max_j; j++) {
                    F0_cache.values[k][l][M][j] = J_Integral_0(a, b, j, M, 1, NULL, epsilon);
                }
            }
            b += b_increment;
        }
        a += a_increment;
        b = 0.0;
    }

}
*/

//void build_F0_cache(long number_of_a, long number_of_b, long max_j, long max_M, Double epsilon) {
void build_F0_cache(long a_per_unit_interval, long b_per_unit_interval, long max_j, long max_M, Double epsilon) {
    
    long number_of_a = a_per_unit_interval + 1;
    long number_of_b = b_per_unit_interval/4 + 1;

    F0_cache.number_of_a = number_of_a;
    F0_cache.number_of_b = number_of_b;
    F0_cache.a_spacing = 1.0/a_per_unit_interval;
    F0_cache.b_spacing = 1.0/b_per_unit_interval;
    F0_cache.a_per_unit_interval = a_per_unit_interval;
    F0_cache.b_per_unit_interval = b_per_unit_interval;
    F0_cache.max_j = max_j;
    F0_cache.max_M = max_M;

    Double stage_1_memory = (((Double)max_j + 1) * (Double)(b_per_unit_interval/4 + 1) * (Double)a_per_unit_interval * max_M + 1.0) * 16.0;
    Double stage_2_memory = (((Double)max_j + 1) * (Double)(b_per_unit_interval/4 + 1) * (Double)(a_per_unit_interval + 1) * max_M ) * 16.0;

    clock_t start_time = clock();

    cout << "Building F0 cache." << endl;
    cout << "    Total used memory for this will be " << (stage_1_memory + stage_2_memory) / 1000000.0 << " million bytes." << endl;

    if(FAKE_PRECOMPUTATION) {
        F0_cache_initialized = true;
        return;
    }

    F0_cache.values = new Complex *** [max_M + 1];

    Double a = 0;
    Double a_increment = 1.0/(number_of_a - 1);
    Double b = 0;
    Double b_increment = .25/(number_of_b - 1);

    cout << "Building F0 cache stage 1 (for M == 0)" << endl;

    Double percent_done = 0.0;
    Double old_percent_done = 0.0;

    {
        // M == 0 case
        F0_cache.values[0] = new Complex ** [(number_of_a - 1) * max_M + 2];
        for(int k = 0; k < (number_of_a - 1) * max_M + 2; k++) {
            percent_done = (Double)k/( (number_of_a - 1) * max_M + 1);
            if(old_percent_done + .005 < percent_done) {
                cout <<  "    " << percent_done * 100 << " percent done with F0 stage 1." << endl;
                old_percent_done = percent_done;
            }
            //if(k % 100 == 0)
            //    cout << "a index = " << k << ", a = " << a << endl;
            F0_cache.values[0][k] = new Complex * [number_of_b];
            for(int l = 0; l < number_of_b; l++) {
//                cout << "b index = " << l << ", b = " << b << endl;
//                cout << a << " " << b << endl;
                F0_cache.values[0][k][l] = new Complex [max_j + 1];
                for(int j = 1; j <= max_j; j++) {
                    F0_cache.values[0][k][l][j] = J_Integral_0(a, b, j, -1, 1, NULL, epsilon, false);
                }
                b += b_increment;
            }
            a += a_increment;
            b = 0;
        }
    }

    a = 0;
    b = 0;

    cout << "Building F0 cache stage 2." << endl;


    percent_done = 0.0;
    old_percent_done = 0.0;
    for(int M = 1; M <= max_M; M++) {
        percent_done = (Double)M/(Double)max_M;
        if(old_percent_done + .005 < percent_done) {
            cout <<  "    " << percent_done * 100 << " percent done with F0 stage 2." << endl;
            old_percent_done = percent_done;
        }
        F0_cache.values[M] = new Complex ** [number_of_a];
        //cout << "M = " << M << endl;
        for(int k = 0; k < number_of_a; k++) {
            //cout << "a index = " << k << ", a = " << a << endl;
            F0_cache.values[M][k] = new Complex * [number_of_b];
            for(int l = 0; l < number_of_b; l++) {
            //    cout << a << " " << b << endl;
                F0_cache.values[M][k][l] = new Complex [max_j + 1];
                for(int j = 0; j <= max_j; j++) {
                    F0_cache.values[M][k][l][j] = J_Integral_0(a, b, j, M, 1, NULL, epsilon, false);
                }
                b += b_increment;
            }
            a += a_increment;
            b = 0;
        }
        a = 0;
    }

    clock_t end_time = clock();

    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Total time to build F0 cache: " << elapsed_time << " seconds." << endl;

    F0_cache_initialized = true;
}


void build_F1_cache(long a_per_unit_interval, long b_per_unit_interval, long max_j, Double epsilon) {


    clock_t start_time = clock();

    int number_of_a = a_per_unit_interval * 6 + 1;
    int number_of_b = b_per_unit_interval/4 + 1;

    F1_cache.number_of_a = number_of_a;
    F1_cache.number_of_b = number_of_b;
    F1_cache.a_spacing = 1.0/a_per_unit_interval;
    F1_cache.b_spacing = 1.0/b_per_unit_interval;
    F1_cache.max_j = max_j;
    F1_cache.a_per_unit_interval = a_per_unit_interval;
    F1_cache.b_per_unit_interval = b_per_unit_interval;

    Double memory_used = (Double)(a_per_unit_interval + 1) * (Double)(b_per_unit_interval/4 + 1) * (Double)(max_j + 1) * 16;

    cout << "Building F1 cache." << endl;
    cout << "    Total memory used for this cache will be " << memory_used/1000000.0 << " million bytes." << endl;

    if(FAKE_PRECOMPUTATION) {
        F1_cache_initialized = true;
        return;
    }



    F1_cache.values = new Complex ** [number_of_a];

    Double a = 0;
    Double a_increment = 6.0/(number_of_a - 1);
    Double b = 0;
    Double b_increment = .25/(number_of_b - 1);


    Double percent_done = 0.0;
    Double old_percent_done = 0.0;

    for(long k = 0; k < number_of_a; k++) {
        percent_done = (Double)k/(Double)number_of_a;
        if(old_percent_done + .005 < percent_done) {
            cout <<  "    " << percent_done * 100 << " percent done with F1." << endl;
            old_percent_done = percent_done;
        }

        //cout << "a index = " << k << ", a = " << a << endl;
        F1_cache.values[k] = new Complex * [number_of_b];
        for(long l = 0; l < number_of_b; l++) {
//            cout << a << " " << b << endl;
//            cout << "b index = " << l << ", b = " << b << endl;
            F1_cache.values[k][l] = new Complex[max_j + 1];
            for(long j = 0; j <= max_j; j++) {
                //F1_cache.values[k][l][j] = pow(100, j) * J_Integral_1(a, b, j, -1, 100, NULL, epsilon/pow(100, j), false);
                F1_cache.values[k][l][j] = J_Integral_1_precomputation(a, b, j, epsilon);
            }
            b += b_increment;
        }
        a += a_increment;
        b = 0.0;
    }


    clock_t end_time = clock();

    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Total time to build F1 cache: " << elapsed_time << " seconds." << endl;

    F1_cache_initialized = true;
}

void build_F2_cache(long max_a1, long a1_per_unit_interval, long a2_per_unit_interval, long b_per_unit_interval, Double epsilon) {
    int number_of_a1 = a1_per_unit_interval + 1;
    int number_of_a2 = a2_per_unit_interval + 1;
    int number_of_b = b_per_unit_interval/4 + 1;


    F2_cache.max_a1 = max_a1;
    F2_cache.number_of_a1 = number_of_a1;
    F2_cache.number_of_a2 = number_of_a2;
    F2_cache.a1_spacing = 1.0/a1_per_unit_interval;
    F2_cache.a2_spacing = 1.0/a2_per_unit_interval;
    F2_cache.b_spacing = 1.0/b_per_unit_interval;
    F2_cache.number_of_b = number_of_b;
    F2_cache.a1_per_unit_interval = a1_per_unit_interval;
    F2_cache.a2_per_unit_interval = a2_per_unit_interval;
    F2_cache.b_per_unit_interval = b_per_unit_interval;

    Double memory_used = (Double)(a1_per_unit_interval * max_a1 + 1) * (Double)(a2_per_unit_interval + 1) * (Double)(b_per_unit_interval/4 + 1) * 16;


    clock_t start_time = clock();

    cout << "Building F2 cache." << endl;
    cout << "    Total memory used for this cache will be " << memory_used/1000000.0 << " million bytes." << endl;

    if(FAKE_PRECOMPUTATION) {
        F2_cache_initialized = true;
        return;
    }


    Double a1_increment = 1.0/(number_of_a1 - 1.0);
    Double a2_increment = 1.0/(number_of_a2 - 1.0);
    Double b_increment = .25/(number_of_b - 1);

    Double a1 = 0;
    Double a2 = 0;
    Double b = 0;


    Double percent_done = 0.0;
    Double old_percent_done = 0.0;

    F2_cache.values = new Complex ** [(number_of_a1 - 1) * max_a1 + 2];
    for(long k = 0; k <= a1_per_unit_interval * max_a1 + 1; k++) {

        percent_done = (Double)k/(  a1_per_unit_interval * max_a1 + 1 );
        if(old_percent_done + .005 < percent_done) {
            cout <<  "    " << percent_done * 100 << " percent done with F2." << endl;
            old_percent_done = percent_done;
        }


        F2_cache.values[k] = new Complex * [number_of_a2];
        for(long j = 0; j < number_of_a2; j++) {
            F2_cache.values[k][j] = new Complex[number_of_b];
            for(long l = 0; l < number_of_b; l++) {
                //cout << a1 << " " << a2 << " " << b << endl;
                F2_cache.values[k][j][l] = J_Integral_2(a1, a2, b, NULL, epsilon, false);
                b += b_increment;
            }
            a2 += a2_increment;
            b = 0.0;
        }
        a1 += a1_increment;
        a2 = 0.0;
    }


    clock_t end_time = clock();

    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Total time to build F2 cache: " << elapsed_time << " seconds." << endl;

    F2_cache_initialized = true;
}


void free_F0_cache() {
    F0_cache_initialized = false;

    for(long k = 0; k < F0_cache.number_of_a; k++) {
        for(long l = 0; l < F0_cache.number_of_b; l++) {
            for(long M = 0; M <= F0_cache.max_M; M++) {
                delete [] F0_cache.values[k][l][M];
            }
            delete [] F0_cache.values[k][l];
        }
        delete [] F0_cache.values[k];
    }
    delete [] F0_cache.values;
}

void free_F1_cache() {
    F1_cache_initialized = false;

    for(long k = 0; k < F1_cache.number_of_a; k++) {
        for(long l = 0; l < F1_cache.number_of_b; l++) {
            delete [] F1_cache.values[k][l];
        }
        delete [] F1_cache.values[k];
    }
    delete [] F1_cache.values;
}

void free_F2_cache() {
    F2_cache_initialized = false;

    for(long k = 0; k < F2_cache.number_of_a1; k++) {
        for(long j = 0; j < F2_cache.number_of_a2; j++) {
            delete F2_cache.values[k][j];
        }
        delete [] F2_cache.values[k];
    }
    delete [] F2_cache.values;
}
