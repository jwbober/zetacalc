#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include "gmpfrxx/gmpfrxx.h"

using namespace std;

mpfr_class sinc(mpfr_class u, mpfr_class Delta) {
    mpfr_class x = Delta * u * const_pi();
    if(abs(x) < .00001) {
        return -pow(x, 6)/5040 + pow(x, 4)/120 - pow(x, 2)/6 + 1;
    }
    mpfr_class a = sin(x);
    a = a/x;

    return a * a;
}

mpfr_class tri(mpfr_class u, mpfr_class Delta) {
    mpfr_class x = abs(u/Delta);
    x = 1 - x;
    x = x/Delta;
    return x;
}

mpfr_class phi(mpfr_class u, mpfr_class Delta) {
    return exp(-Delta * u * u);
}

mpfr_class phi_hat(mpfr_class u, mpfr_class Delta) {
    return sqrt(const_pi()/Delta) * exp(-const_pi() * const_pi() * u * u/Delta);
}

mpfr_class prime_sum_term(mpfr_class t, mpfr_class Delta, mpfr_class (*fhat)(mpfr_class, mpfr_class), int num_terms = 0) {
    mpfr_class endpoint;
    if(num_terms == 0) 
        endpoint = exp(2 * const_pi() * Delta);
    else
        endpoint = num_terms;
    mpz_class p;
    p = 1;
    mpfr_class S(0);
    while(p < endpoint) {
        if(mpz_probab_prime_p(p.get_mpz_t(), 10)) {
            int n = 1;
            mpfr_class logp = log(mpfr_class(p));
            mpfr_class p_power(p);
            while(p_power <= endpoint) {
                S = S - logp * cos(t * log(p_power))/sqrt(p_power) * fhat(log(p_power)/(mpfr_class(2) * const_pi()), Delta);

                n += 1;
                p_power *= p;
            }
        }
        p = p + 1;
        if(p % 10000 == 0) {
            cout << p.get_d()/endpoint << " done with prime sum. Answer so far: " << S/const_pi() << endl;
        }
    }
    return S/const_pi();
}

mpfr_class zeros_term(mpfr_class t0, double * zero_list, int number_of_zeros, mpfr_class t, mpfr_class Delta, mpfr_class (*f)(mpfr_class, mpfr_class)) {
    mpfr_class S(0);

    for(int n = 0; n < number_of_zeros; n++) {
        S = S + f(t0 + zero_list[n] - t, Delta);
    }
    return S;
}

complex<mpfr_class> digamma(complex<mpfr_class> z) {
    return log(z) - mpfr_class(1)/(mpfr_class(2) * z) + mpfr_class(1)/(mpfr_class(12) * z * z);
}

mpfr_class integral_term(mpfr_class t, mpfr_class Delta, mpfr_class (*f)(mpfr_class, mpfr_class)) {
    // return (an approximation to) the integral
    // 1/(2pi) * int_{-oo}^{oo} f(u - t) Re(digamma(.25 + iu/2)du)

    // For now we just compute a Riemann sum over a relatively
    // small interval
    
    double length = 50;
    mpfr_class u = t - length/2;
    mpfr_class end = t + length/2;
    double delta = .002;

    unsigned long num_terms = length/delta;

    mpfr_class S = 0;
    complex<mpfr_class> z(.25, u/2);
    mpfr_class x;

    unsigned long count = 0;
    while(u < end) {
        x = delta * f(u - t, Delta) * digamma(z).real();
        S = S + x;
        u = u + delta;
        z.imag() += delta/2;
        count++;
        if(count % 10000 == 0) {
            cout << (double)count/(double)num_terms << " finished with integral." << endl;
        }
    }

    return S/(2 * const_pi());

}

double test( mpfr_class t0, double * zeros, int number_of_zeros, mpfr_class t, mpfr_class Delta, mpfr_class (*f)(mpfr_class, mpfr_class), mpfr_class (*fhat)(mpfr_class, mpfr_class), int num_terms = 0) {
    mpfr_class fhat_term, integral, prime_sum, LHS;
    
    fhat_term = -fhat(0, Delta) * log(const_pi()) / (mpfr_class(2) * const_pi());
    integral = integral_term(t, Delta, f);
    prime_sum = prime_sum_term(t, Delta, fhat, num_terms);
    
    LHS = zeros_term(t0, zeros, number_of_zeros, t, Delta, f);

    mpfr_class error = LHS - (fhat_term + integral + prime_sum);

    cout << LHS << " - (" << fhat_term << " + " << integral << " + " << prime_sum << ") = " << LHS - (fhat_term + integral + prime_sum) << endl;
    return error.get_d();
}


int main(int argc, char ** argv) {
    if(argc < 2) {
        cout << "Missing input file name." << endl;
        return 1;
    }

    mpfr_class::set_dprec(300);

    ifstream infile(argv[1]);
    mpfr_class t0;
    infile >> t0;

    int N;
    infile >> N;

    double zeros[N];
    for(int n = 0; n < N; n++) {
        infile >> zeros[n];
    }

    cout << setprecision(17);

    mpfr_class t, Delta;
    //t = t0 + 29.5925;
    //t = t0 + 19.592;
    t = t0 + 15.5;
    Delta = 2;

    //mpfr_class x("13066434408794883841332.2315014023");
    //cout << prime_sum_term(x, Delta, phi_hat, 100000) << endl;

    //test(t0, zeros, N, t, Delta, sinc, tri);

    int M = 100;
    double spacing = .1;
    double errors[M];

    for(int m = 0; m < M; m++) {
        t = t0 + 15.0 + spacing * m;
        errors[m] = test(t0, zeros, N, t, Delta, phi, phi_hat, 100000);
    }

    for(int m = 0; m < M; m++) {
        cout << errors[m] << endl;
    }

    //test(t0, zeros, N, t, Delta, phi, phi_hat, 10000);
    //test(t0, zeros, N, t, Delta, sinc, tri);
    
    //cout << integral_term(t, 1, sinc) << endl;

    return 0;
}
