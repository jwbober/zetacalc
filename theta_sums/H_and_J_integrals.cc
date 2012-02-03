#include "theta_sums.h"
#include "precomputed_tables.h"
#include "log.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include <cstdlib>

using namespace std;

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


void J_Integral_0(complex<double> * J, double a, double b, int j, int M, theta_cache * cache, double * epsilon) {
    //
    // Compute and add to J[l] the integral
    //
    //   /1
    //   |  j                          2    1 - exp(-2 pi M t)
    //   | t  exp( -2pi a t - 2pi i b t  ) -------------------- dt
    //   |                                    exp(2 pi t) - 1
    //   /0
    //
    // to precision at least epsilon[j].
    //
    // Note (since elsewhere the documentation might misleading) that we
    // DO NOT divide by K^j. We let the caller handle that, since it will
    // have to be done a few times. (Hence we do not even take K as an
    // argument.)
    //
    //
    //
    //
    // To compute the integral for a given value of j, we do a taylor
    // expansion on the term exp(-2pi i b t^2), and rewrite the integral as
    //
    //  oo           
    // ---         r /1
    // \   (-2pi b)  |   j + 2r                  1 - exp(-2 pi M t)  
    //  )   ------   |  t       exp( - 2pi a t) ------------------- dt
    // /      r!     |                             exp(2 pi t) - 1
    // ---           /0
    // r=0
    //
    // These subintegrals are computed by the function H_integral_0(j, a, M).
    // When j > 0, H_integral_0(j, a, M) < 1, which means that the summands
    // quickly decay and we can truncate the sum for a given j as soon as
    //
    //         r
    //  (2pi b)
    //  -------   < epsilon[j].
    //    r!
    //
    // For different values of j, many terms are repeated, so we first compute
    // and save a bunch of the subintegrals H_integral_0(j, a, M). It may be
    // a little tricky to decide just how much precision we need when
    // computing each one, so to be conservative we
    //
    // -- figure out the smallest epsilon

    double small_epsilon = epsilon[0];
    for(int n = 1; n <= j; n++) {
        small_epsilon = min(small_epsilon, epsilon[n]);
    }

    // -- then compute a list of values (2 pi b)^r / r! until they are
    //    smaller than small_epsilon
    
    vector<double> Z; Z.reserve(30); //  I expect that it should be rare for r to be this
                                     //  large. (I don't know what is optimal here. We
                                     //  could just allocate space for an array of size
                                     //  200 or something on the stack, since eventually
                                     //  we run out of precision and things are 0.)
    double largest_Z = 1;

    int R; // this is going to be the actual largest r that we need
           // but for artificial reasons (which are revealed below)
           // we are going to inflate Z so that its size is a multiple
           // of 4.
    {
        double b_power = 1.0;
        double this_term = 1.0;
        int r = 0;
        Z.push_back(this_term);
        while(this_term > small_epsilon) {
            r++;
            b_power *= b;
            this_term = b_power * two_pi_over_factorial_power(r);
            Z.push_back(this_term);
            largest_Z = max(largest_Z, this_term);
        }
        R = r;
        int N = Z.size() % 4;
        if(N != 0) {
            for(int n = 0; n < 4 - N; n++) {
                Z.push_back(0);
            }
        }
    }

    // -- Now we compute a bunch of values of H_integral_0(l, a, M). The
    //    largest j we could need is l = j + 2 * (Z.size() - 1). When l is
    //    large be probably don't need as much precision as when l is
    //    small, if we take into account the various epsilon[l], and the
    //    various Z[l]. Instead we just always make epsilon the smallest
    //    it could ever need to be.
    //

    double H_integrals[j + 2 * (Z.size() - 1) + 1];
    double new_epsilon[j + 2 * (Z.size() - 1) + 1];
    for(unsigned int n = 0; n < j + 2 * (Z.size() - 1) + 1; n++) {
        new_epsilon[n] = 1000;
    }
    for(int l = 0; l <= j; l++) {
        for(int r = 0; r <= R; r++) {
            new_epsilon[l + 2 * r] = min(epsilon[l]/Z[r], new_epsilon[l + 2*r]);
        }
    }
    //double smaller_epsilon = small_epsilon/(largest_Z + Z.size());
    for(int l = 0; l <= j + 2 * R; l++) {
        H_integrals[l] = real(H_Integral_0(l, a, M, new_epsilon[l]/Z.size())); // REMARK: H_Integral_0
                                                                               // is always real, so
                                                                               // it is stupid that it
                                                                               // does complex arithmetic.
    }

    // -- Now we finally compute the integrals that we are looking for. We have
    //
    //      J[l] = sum_{r} (-i)^r Z[r] H_integrals[l + 2*r]
    //

    for(int l = 0; l <= j; l++) {
        for(int r = 0; r <= R; r += 4) {
            real(J[l]) += Z[r]     * H_integrals[l + 2*r];          // (-i)^0 ==  1
            imag(J[l]) -= Z[r + 1] * H_integrals[l + 2*(r + 1)];    // (-i)^1 == -i
            real(J[l]) -= Z[r + 2] * H_integrals[l + 2*(r + 2)];    // (-i)^2 == -1
            imag(J[l]) += Z[r + 3] * H_integrals[l + 2*(r + 3)];    // (-i)^3 == i
            if(Z[r] < epsilon[l])
                break;
        }
    }

    //for(int l = 0; l <= j; l++) {
    //    if(!isinf(epsilon[l])) {
    //        J[l] += J_Integral_0(a, b, l, M, K, cache, epsilon[l]);
    //    }
    //}
}




void J_Integral_2(complex<double> * J, double a1, double a2, double b, int j, theta_cache * cache, double * epsilon) {
    //
    // Compute and add to J[l] the integral
    //
    //   /1
    //   |  j                2    (exp(-2 pi a1 t) - exp(-2 pi a2 t))
    //   | t  exp(- 2pi i b t  ) ------------------------------------- dt
    //   |                                   exp(2 pi t) - 1
    //   /0
    //
    // to precision at least epsilon[j].
    //
    // Note (since elsewhere the documentation might misleading) that we
    // DO NOT divide by K^j. We let the caller handle that, since it will
    // have to be done a few times. (Hence we do not even take K as an
    // argument.)
    //
    // This integral is in fact the same as
    //
    //   lim_{M ---> oo} J_Integral_0(a1, b, j, M) - J_Integral_0(a2, b, j, M)
    //
    // and we could compute it as a different of 2 such integrals, but this
    // will not work in the j = 0 case, and in any event it may be more efficient
    // to just compute the difference directly.
    //
    // The method of computation is _exactly_ the same as the method
    // for J_Integral_0 (so it is not fully explained again) except that the
    // "subintegral" is different, and we compute it in the function H_Integral_2().
    // Basically all of the code is duplicated, though.
    //
    // -- figure out the smallest epsilon

    double small_epsilon = epsilon[0];
    for(int n = 1; n <= j; n++) {
        small_epsilon = min(small_epsilon, epsilon[n]);
    }

    // -- then compute a list of values (2 pi b)^r / r! until they are
    //    smaller than small_epsilon
    
    vector<double> Z; Z.reserve(30); //  I expect that it should be rare for r to be this
                                     //  large. (I don't know what is optimal here. We
                                     //  could just allocate space for an array of size
                                     //  200 or something on the stack, since eventually
                                     //  we run out of precision and things are 0.)
    double largest_Z = 1;

    int R; // this is going to be the actual largest r that we need
           // but for artificial reasons (which are revealed below)
           // we are going to inflate Z so that its size is a multiple
           // of 4.
    {
        double b_power = 1.0;
        double this_term = 1.0;
        int r = 0;
        Z.push_back(this_term);
        while(this_term > small_epsilon) {
            r++;
            b_power *= b;
            this_term = b_power * two_pi_over_factorial_power(r);
            Z.push_back(this_term);
            largest_Z = max(largest_Z, this_term);
        }
        R = r;
        int N = Z.size() % 4;
        if(N != 0) {
            for(int n = 0; n < 4 - N; n++) {
                Z.push_back(0);
            }
        }
    }

    // -- Now we compute a bunch of values of H_integral_0(l, a, M). The
    //    largest j we could need is l = j + 2 * (Z.size() - 1). When l is
    //    large be probably don't need as much precision as when l is
    //    small, if we take into account the various epsilon[l], and the
    //    various Z[l]. Instead we just always make epsilon the smallest
    //    it could ever need to be.
    //

    double H_integrals[j + 2 * (Z.size() - 1) + 1];
    double new_epsilon[j + 2 * (Z.size() - 1) + 1];
    for(unsigned int n = 0; n < j + 2 * (Z.size() - 1) + 1; n++) {
        new_epsilon[n] = 1000;
    }
    for(int l = 0; l <= j; l++) {
        for(int r = 0; r <= R; r++) {
            new_epsilon[l + 2 * r] = min(epsilon[l]/Z[r], new_epsilon[l + 2*r]);
        }
    }
    //double smaller_epsilon = small_epsilon/(largest_Z + Z.size());
    for(int l = 0; l <= j + 2 * R; l++) {
        H_integrals[l] = real(H_Integral_2(l, a1, a2, new_epsilon[l]/Z.size())); // REMARK: H_Integral_0
                                                                               // is always real, so
                                                                               // it is stupid that it
                                                                               // does complex arithmetic.
    }

    // -- Now we finally compute the integrals that we are looking for. We have
    //
    //      J[l] = sum_{r} (-i)^r Z[r] H_integrals[l + 2*r]
    //

    for(int l = 0; l <= j; l++) {
        for(int r = 0; r <= R; r += 4) {
            real(J[l]) += Z[r]     * H_integrals[l + 2*r];          // (-i)^0 ==  1
            imag(J[l]) -= Z[r + 1] * H_integrals[l + 2*(r + 1)];    // (-i)^1 == -i
            real(J[l]) -= Z[r + 2] * H_integrals[l + 2*(r + 2)];    // (-i)^2 == -1
            imag(J[l]) += Z[r + 3] * H_integrals[l + 2*(r + 3)];    // (-i)^3 == i
            if(Z[r] < epsilon[l])
                break;
        }
    }

    //for(int l = 0; l <= j; l++) {
    //    if(!isinf(epsilon[l])) {
    //        J[l] += J_Integral_0(a, b, l, M, K, cache, epsilon[l]);
    //    }
    //}
}

Complex J_Integral_0(Double a, Double b, int j, int M, int K, theta_cache * cache, Double epsilon) {
    //
    // Compute the integral (1/K^j) int_0^1 t^j exp(-2 pi a t - 2 pi i b t^2) (1 - exp(-2 pi M t))/(exp(2 pi t) - 1) dt
    //
    // We can do this by expanding the denominator in a geometric series
    // and doing a taylor expansion in b, which puts this integral in the form
    //
    // sum_{r=0}^\infty ((2 pi i b)^r)/r! \sum_{m=1}^M int_0^1 t^j exp(-2 pi(m + a)t)dt
    //
    // -1 corresponds to infinity

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
        return 0.0;
    }
    epsilon = epsilon * K_pow_j;

    Complex S = (Complex)0;

    int r = 0;
    Double error = 2 * epsilon;
    
    if(b == 0) {
        Complex answer = H_Integral_0(j, a, M, epsilon)*K_pow_minus_j;
        return answer;
    }
    
    int N = max(1, to_int(ceil(-fastlog(epsilon))));    // an approximation of the number of terms
                                                        // we compute in the taylor expansion

    Double b_power = (Double)1.0;
    while(error > epsilon) {
        Complex z = minus_I_power(r) * two_pi_over_factorial_power(r) * b_power;
        error = abs(z);
        z *= H_Integral_0(2 * r + j, a, M, epsilon/(error *N));     // TODO: We really don't need to compute this term to
                                                                    // this much precision usually. We really should figure
                                                                    // our how many terms we are going to compute
                                                                    // and then divide appropriately.
                                                                    //
        S = S + z;
        b_power *= b;
        r++;
    }

    S = S * K_pow_minus_j;

    return S;
}


Complex J_Integral_1(Double a, Double b, int j, int M, int K, theta_cache * cache, Double epsilon) {
    //
    // Compute the integral int_1^K t^j exp(-2 pi a t - 2 pi i b t^2)(1 - exp(-2 pi M t))/(exp(2 pi t) - 1) dt
    //
    // If M == -1, then we treat it as positive infinity.
    //
    
    // We truncate the integral at L, where L is given by

    //int L = min(K, max(1, to_int(ceil(-LOG(epsilon)/(2 * PI * (1 + a))))));

    int L = ceil(-fastlog(epsilon) * (1/(2 * PI)) * inverse((int)(1 + a)));
    if(L <= 1) {
        return 0.0;
    }
    L = max(1.0, L - j * fastlog( (Double)K * inverse(L + 1)) * (1.0/(2 * PI)) );
    L = min(K, L);


    if(L <= 1) {
        return 0.0;
    }

    // Now we compute the integral as a sum over unit length integrals.
    // On each unit length interval, we do a change of variables to get
    // the range of integration between 0 and 1, and reduce the computation
    // to a sum of exponentials times G()

    Complex S = (Complex)0;

    Double one_over_K_to_the_j;
    if(cache)
        one_over_K_to_the_j = K_power(-j, cache);
    else
        one_over_K_to_the_j = pow(K, -j);

    if( 5 + j * fastlog2(L) + fastlog2(one_over_K_to_the_j) - 2 * PI * (a + 1) * log2(E) < fastlog2(epsilon)) {
        return 0.0;
    }

    Double exp_minus_twopi = EXP(-2.0 * PI);
    Double exp_minus_twopi_n = 1.0;

    for(Double n = (Double)1; n <= L - 1; n = n + 1) {
        exp_minus_twopi_n *= exp_minus_twopi;
        int end_point;
        if(M == -1) {
            end_point = to_int(ceil(-fastlog(epsilon/(2 * L))/(2 * PI * n) -  j * fastlog((Double)K/(n + 1))/(2 * PI) ));
        }
        else {
            end_point = min(M, to_int(ceil(-fastlog(epsilon/(2 * L))/(2 * PI * n) - j * fastlog((Double)K/(n + 1))/(2 * PI))  )  );
        }

        end_point = max(1, end_point);

        Complex S1 = 0;

        Complex x = EXP(-2.0 * PI * n * (1.0 + a + I * b * n));

        for(Double m = (Double)1; m <= end_point; m = m + 1) {
            
            if(m > 1)
                x = x * exp_minus_twopi_n;

            Complex z =  G(I*(m + a + (Complex)2.0 * I * b * n), -b, n, j, epsilon/(abs(x) * end_point * Double(L - 1) * one_over_K_to_the_j));
            z *= x;
            S1 = S1 + z;
        }
        S1 = S1 * one_over_K_to_the_j;
        S = S + S1;
    }


    return S;

}

Complex H_Integral_2(int j, Double a1, Double a2, Double epsilon) {

    Complex S = (Complex)0;

    int C = max(to_int(ceil(j/(2*PI*E))), to_int(ceil(-fastlog(epsilon)/(2 * PI) ) ));

    if(j % 2 == 0) {
        for(int m = 1; m <= C; m++) {
            Complex z = H(j, m + a1, epsilon/(C + 1)) - H(j, m + a2, epsilon/(C + 1));
            S = S + z;
        }

        S = S + factorial(j)/two_pi_power(j + 1) * infinite_sum_of_differenced_inverse_powers(a1, a2, C + 1, j + 1, epsilon);
    }
    else {
        for(int m = 1; m <= C; m++) {
            Complex z = H(j, m + a1, epsilon/(C + 1)) + H(j, m + a2, epsilon/(C + 1));
            S = S + z;
        }
        S = S + factorial(j)/two_pi_power(j + 1) * (sum_of_offset_inverse_powers(a1, C + 1, -1, j + 1, epsilon) + sum_of_offset_inverse_powers(a2, C + 1, -1, j + 1, epsilon));
    }

    return S;

}

Complex J_Integral_2(Double a1, Double a2, Double b, theta_cache * cache, Double epsilon) {
    //
    // Compute the integral int_0^1 exp(-2 pi i b t^2) (exp(-2 pi a1 t) - exp(-2 pi a2 t))(exp(2 pi t) - 1) dt,
    //
    // which is equal to 
    //
    // lim_{M --> oo} J_Integral_0(a1, b, M) - J_Integral_0(a2, b, M)
    //

    if(a2 > a1) {
        //if (2.0 * ( (a2 - a1)/a2 * log(1 + a2) ) < epsilon) {
        if(1 + fastlog2( (a2 - a1)/a2 * (fastlog(1 + a2) + 1.0)) < fastlog2(epsilon)) {
            return 0.0;
        }
    }
    else {
        //if (2.0 * ( (a1 - a2)/a1 * log(1 + a1) ) < epsilon) {
        if(1 + fastlog2( (a1 - a2)/a1 * (fastlog(1 + a1) + 1.0)) < fastlog2(epsilon)) {
            return 0.0;
        }
    }

    Complex S = (Complex)0;

    int r = 0;
    Double error = 2 * epsilon;
    int sign = -1;
    
    if(b == 0) {
        Complex answer = H_Integral_2(0, a1, a2, epsilon);
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



    return S;

}

void JBulk(complex<double> * J, Double a, Double b,
           int j, int M, int K, theta_cache * cache, double * epsilon) {

    for(int l = 0; l <= j; l++) {
        J[l] = 0.0;
    }

    double Kpower = 1.0;
    double new_epsilon[j + 1]; new_epsilon[0] = epsilon[0];

    for(int l = 1; l <= j; l++) {
        Kpower *= K;
        new_epsilon[l] = epsilon[l] * Kpower;
    }

    J_Integral_0(J, a, b, j, M, cache, new_epsilon);
    Kpower = 1.0;
    double Kinv = 1.0/K;
    for(int l = 1; l <= j; l++) {
        Kpower *= Kinv;
        J[l] *= Kpower;
        //
        // FOR TESTING PURPOSES. TO BE REMOVED.
        //
        //complex<double> z = J_Integral_0(a, b, l, M, K, cache, epsilon[l]/2);
        //double a = abs(J[l] - z);
        //if(a > epsilon[l]) {
        //    cout << a << " ";
        //    cout << b << " ";
        //    cout << l << " ";
        //    cout << M << " ";
        //    cout << K << " ";
        //    cout << epsilon[l] << " " << endl;
        //    cout << "ohno" << endl;
        //    exit(-1);
        //}
    }

    for(int l = 0; l <= j; l++) {
        double x = epsilon[l]/2;
        complex<double> B = 0.0;
        if(!isinf(x)) {
            B = J_Integral_1(a, b, l, M, K, cache, x);
        }
        J[l] += B;
    }
}                                                                                       


void JBoundary(complex<double> * J, double a1, double a2, double b, int j, int K, theta_cache * cache, double * epsilon){ 
    for(int l = 0; l <= j; l++) {
        J[l] = 0.0;
    }

    double Kpower = 1.0;
    double new_epsilon[j + 1]; new_epsilon[0] = epsilon[0];

    for(int l = 1; l <= j; l++) {
        Kpower *= K;
        new_epsilon[l] = epsilon[l] * Kpower;
    }

    J_Integral_2(J, a1, a2, b, j, cache, new_epsilon);
    Kpower = 1.0;
    double Kinv = 1.0/K;
    for(int l = 1; l <= j; l++) {
        Kpower *= Kinv;
        J[l] *= Kpower;
        //
        // FOR TESTING PURPOSES. TO BE REMOVED.
        //
        //complex<double> z;
        //if(l == 0)
        //    z = J_Integral_2(a1, a2, b, cache, epsilon[l]/2);
        //else
        //    z = J_Integral_0(a1, b, l, -1, K, cache, epsilon[l]/2) + (Double)minus_one_power(l+1) * J_Integral_0(a2, b, l, -1, K, cache, epsilon[l]/2);
        //double a = abs(J[l] - z);
        //if(a > epsilon[l]/2) {
        //    cout << z << " " << J[l] << " " << a << endl;
        //    cout << a1 << " ";
        //    cout << a2 << " ";
        //    cout << b << " ";
        //    cout << l << " ";
        //    cout << K << " ";
        //    cout << epsilon[l] << " " << endl;
        //    cout << "ohno" << endl;
        //}
    }

    for(int l = 0; l <= j; l++) {
        double x = epsilon[l]/2;
        complex<double> B = 0.0;
        if(!isinf(x)) {
            B = J_Integral_1(a1, b, l, -1, K, cache, x) + (Double)minus_one_power(l+1) * J_Integral_1(a2, b, l, -1, K, cache, x);
        }
        J[l] += B;
    }

}                                                                                     
