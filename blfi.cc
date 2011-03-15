/*
 * Program to use bandlimited interpolation to compute values of the zeta
 * function from a precomputed list of values of the riemann siegel sum
 */



#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <mpfr.h>
#include <gmp.h>
#include <cmath>
#include <complex>
#include <string>
#include <cstdlib>
#include <vector>
#include "gmpfrxx/gmpfrxx.h"

typedef double Double;
typedef std::complex<double> Complex;
using namespace std;

#include "rs_sum.h"


double sinc(double x);
double blfi_kernel(double u, double c, double epsilon1);

mpfr_class siegel_theta(mpfr_class t) {
    return (t/2) * log(t/(2 * const_pi())) - t/2 - const_pi()/8;

    /*
    mpfr_t temp;
    mpfr_init2(temp, mpfr_get_prec(theta));
    
    
    mpfr_const_pi(theta, GMP_RNDN);             // theta = pi
    mpfr_set_ui(temp, 1u, GMP_RNDN);
    mpfr_exp(temp, temp, GMP_RNDN);           // temp = e



    mpfr_mul(theta, theta, temp, GMP_RNDN);     // theta = e * pi
    mpfr_mul_2ui(theta, theta, 1, GMP_RNDN);    // theta = 2 e pi




    mpfr_div(theta, t, theta, GMP_RNDN);        // theta = t/2epi
    mpfr_log(theta, theta, GMP_RNDN);           // theta = log(t/2epi)
    mpfr_mul(theta, theta, t, GMP_RNDN);        // theta = t log (t/2epi)
    mpfr_div_2ui(theta, theta, 1, GMP_RNDN);    // theta = .../2

    mpfr_const_pi(temp, GMP_RNDN);              //
    mpfr_div_2ui(temp, temp, 3, GMP_RNDN);      //
    mpfr_sub(theta, theta, temp, GMP_RNDN);     // theta -= pi/8

    mpfr_clear(temp);
    */
}

mpz_class N_approx(mpfr_class t) {
    return siegel_theta(t)/const_pi();
    //mpfr_t theta;
    //mpfr_init2(theta, 300);
    //siegel_theta(theta, t);
    //mpfr_div(theta, theta, mp_pi, GMP_RNDN);
    //mpfr_add_ui(theta, theta, 1u, GMP_RNDN);
    //mpfr_get_z(N, theta, GMP_RNDD);
    //mpfr_clear(theta);
}

mpfr_class grampoint(mpz_class n) {
    // return the nth grampoint
    // (to precision roughly 1/t)

    // Let's be lazy and just start with
    // test points n/log n and n, and then use
    // a binary search.

    mpfr_class t1;
    mpfr_class t2;
    mpfr_class t3;
    mpfr_class theta;
    mpfr_class target;
    mpfr_class d;

    target = const_pi() * n;
    t1.set_prec(300);
    t2.set_prec(300);
    t3.set_prec(300);
    theta.set_prec(300);
    target.set_prec(300);
    d.set_prec(300);

    t1 = n;
    t1 = t1/log(t1);
    t2 = n;

    target = n * const_pi();


    //mpfr_set_z(t2, n, GMP_RNDN);
    //mpfr_log(t1, t2, GMP_RNDN);
    //mpfr_div(t1, t2, t2, GMP_RNDN);

    double epsilon = .00000000001;

    // we now maintain the invariant that t1 is too
    // small and t2 is too large

    //mpfr_sub(d, t2, t1, GMP_RNDN);
    d = t2 - t1;
    while(d > epsilon) {
        //cout << t1 << endl;
        //cout << t2 << endl;
        t3 = (t1 + t2)/2;
        //mpfr_add(t3, t1, t2, GMP_RNDN);
        //mpfr_div_2ui(t3, t3, 1, GMP_RNDN);

        theta = siegel_theta(t3);
        //siegel_theta(theta, t3);
        //if(mpfr_cmp(theta, target) < 0) {
        if(theta < target) {
            t1 = t3;
            //mpfr_set(t1, t3, GMP_RNDN);
        }
        else {
            //mpfr_set(t2, t3, GMP_RNDN);
            t2 = t3;
        }
        d = t2 - t1;
        //mpfr_sub(d, t2, t1, GMP_RNDN);
    }

    return (t1 + t2)/2;
    //mpfr_add(g, t1, t2, GMP_RNDN);
    //mpfr_div_2ui(g, g, 1, GMP_RNDN);

    //mpfr_clear(t1);
    //mpfr_clear(t2);
    //mpfr_clear(t3);
    //mpfr_clear(theta);
    //mpfr_clear(d);
    //mpfr_clear(target);
}

Complex I(0, 1);



class ZetaComputation {

public:

    mpfr_class t0;
    double delta;
    
    double tau;
    double beta;
    double lambda;
    double epsilon1;
    double alpha;

    int N;
    Complex * rs_sum;

    ZetaComputation(mpfr_class _t0, double _delta, int _N, Complex * _rs_sum) {
        //mpfr_init2(t0, mpfr_get_prec(_t0));
        //mpfr_set(t0, _t0, GMP_RNDN);
        
        t0.set_prec(_t0.get_prec());
        t0 = _t0;

        delta = _delta;
        N = _N;

        rs_sum = new Complex[N];

        for(int n = 0; n < N; n++) {
            rs_sum[n] = _rs_sum[n];
        }

        //mpfr_t temp;
        //mpfr_init2(temp, 300);

        //mpfr_const_pi(temp, GMP_RNDN);
        //mpfr_mul_2ui(temp, temp, 1, GMP_RNDN);

        //mpfr_div(temp, t0, temp, GMP_RNDN);
        //mpfr_log(temp, temp, GMP_RNDN);

        //tau = .5 * mpfr_get_d(temp, GMP_RNDN);
        mpfr_class temp;
        temp.set_prec(300);
        temp = log(t0/(2 * const_pi()));
        tau = .5 * temp.get_d();
        beta = M_PI/delta;
        lambda = (beta + tau)/2.0;
        epsilon1 = (beta - tau)/2.0;
        alpha = tau;

        //mpfr_clear(temp);

    }

    double Z(double t) {
        // return the value of Z(t0 + _t)
        

        Complex S = 0;

        double c = 20 * M_PI * epsilon1/beta;

        for(int n = ceil(t/delta) - 19; n <= ceil(t/delta) + 18; n++) {
            double u = n * delta - t;
            
            Complex z = rs_sum[n] * exp(-I * alpha * u) * sinc(lambda * u) * blfi_kernel(u, c, epsilon1);
            S = S + z;

        }

        S = S * lambda/beta;

        mpfr_t temp;
        mpfr_init2(temp, 300);

        mpfr_add_d(temp, t0.get_mpfr_t(), t, GMP_RNDN);

        S = 2.0 * S * rs_rotation(temp);
        
        mpfr_clear(temp);

        return S.real();
    }

    double find_zero(double t1, double t2, double epsilon = .00000001) {
        // given points t1 and t2 such that Z(t1) * Z(t2) < 0,
        // return an approximation of a zero in this interval

        if(Z(t1) * Z(t2) > 0) {
            return -1;
        }
        
        double t3;

        if(Z(t1) < 0) {
            t3 = t1;
            t1 = t2;
            t2 = t3;
        }

        while( abs(t2 - t1) > epsilon) {
            t3 = (t1 + t2)/2.0;
            if(Z(t3) > 0)
                t1 = t3;
            else
                t2 = t3;
        }

        return (t1 + t2)/2.0;
    }

    vector<double> find_zeros(double start, double end, double delta, bool verbose=false) {
        double t1 = start;
        double t2 = start + delta;

        vector<double> zeros(0);
        int count = 0;

        while(t2 < end) {
            if(Z(t1) * Z(t2) < 0) {
                double zero = find_zero(t1, t2);
                if(verbose) {
                    count++;
                    cout << "Found zero number " << count << " at " << zero << endl;
                }
                zeros.push_back(find_zero(t1, t2));
            }
            t1 += delta;
            t2 += delta;
        }

        return zeros;
    }

    //void calculate_N(mpz_class N, mpfr_class t) {
    //    mpfr_class g;
    //    mpfr_class T;

    //    g.set_prec(300);
    //    T.set_prec(300);

    //    T = t + t0;
    //    N_approx(N, T);



    //}

};


double sinc(double x){
    //Taylor coefficients in sin expansions
    const double sinc_coeff[6]={1, -1.0 / 6.0, 1.0 / 120.0, -1.0 / 5040.0, 1 / 362880.0, -1.0 / 39916800.0};

    //Parameters needed in sinc function; help overcome instability when sinc(x) is called with very small x 
    double sinc_tol=1.e-5;
    int sinc_terms = 3;

    if(abs(x) < sinc_tol) {
        //double x_squared = x*x;
        //double x_power = x_squared;
        //double ans = 1;

        //for(int j = 1; j < sinc_terms; j++){
        //    ans = ans + sinc_coeff[j] * x_power;
        //    x_power = x_squared * x_power;
        //}
        
        //return ans;
        
        return 1/120.0 * pow(x, 4) - 1.0/6.0 * pow(x, 2) + 1;

    }
    double ans = sin(x) / x;
    return ans;
}


//computes the kernel function in BLFI
double blfi_kernel(double u, double c, double epsilon_1){
    //Taylor coefficients in sinh and sin expansions
    const double sinh_coeff[6]={1, 1 / 6.0, 1 / 120.0, 1 / 5040.0, 1 / 362880.0, 1 / 39916800.0};
    //Parameters needed in kernel function; help overcome instability when kernel(x) is called with very small x 
    double sinh_tol=1.e-5;
    int sinh_terms = 3;

    double x = c*c - pow( epsilon_1*u, 2);

    if( x < 0 ) {
        cout << endl << "kernel: u is too large!" << endl;
        exit(1);
    }

    double w;

    if(abs(x) < sinh_tol) {

        w = 1/39916800.0*pow(x, 5) + 1/362880.0*pow(x, 4) + 1/5040.0*pow(x, 3) + 1/120.0 *pow(x,2) + 1/6 * x + 1;

        //double x_power = x_squared;
        //double ans = 1;
        
        //for(int j = 1; j < sinh_terms; j++) {
        //    ans = ans + sinh_coeff[j] * x_power;
        //    x_power = x_squared * x_power;
        //}
        //ans = ans * 2 * c/(exp(-c) + exp(c));
        //return ans;
    }
    else {
        w = sinh(sqrt(x))/sqrt(x);
    }

    //double sinh_x = 0.5 * (exp(sqrt(x)) + exp(- sqrt(x)));
    //double ans = sinh_x / sqrt(x);
    //ans = ans * 2 * c/(exp(-c) + exp(c));
    return w * c/cosh(c);
}


int main(int argc, char * argv[]) {

    cout << setprecision(200);

    //mpfr_t g;
    //mpfr_t theta;
    //mpfr_init2(g, 300);
    //mpfr_init2(theta, 300);
    //mpz_t n;
    //mpz_init(n);

    //mpfr_init2(mp_pi, 300);
    //mpfr_const_pi(mp_pi, GMP_RNDN);

    mpfr_class g;
    mpz_class n;

    g.set_prec(300);

    //while(1) {
    //    cin >> n;
        //mpfr_set_z(g, n, GMP_RNDN);
        //siegel_theta(theta, g);
    //    g = grampoint(n);
    //    cout << g << endl;
        //grampoint(g, n);
        ///mpfr_get_z(n, g, GMP_RNDD);
        //mpfr_sub_z(g, g, n, GMP_RNDN);
        //cout << n << " + " << mpfr_get_d(g, GMP_RNDN) << endl;
        //mpfr_get_z(n, theta, GMP_RNDD);
        //mpfr_sub_z(theta, theta, n, GMP_RNDN);
        //cout << n << " + " << mpfr_get_d(theta, GMP_RNDN) << endl;

    //}


    if(argc < 2) {
        cout << "missing required input file." << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    
    mpfr_class t;
    t.set_prec(300);
    double delta;
    int N = 1000;
    complex<double> rs_sum[N];

    infile >> t;
    infile >> delta;

    for(int n = 0; n < N; n++) {
        infile >> rs_sum[n];
    }

    infile.close();

    ZetaComputation Z(t, delta, N, rs_sum);

    //cout << mpfr_get_d(t, GMP_RNDN) << endl;
    
    //for(int n = 0; n < 10000; n++) {
    //    cout << Z.Z(1 + 38/10000.0 * n) << endl;
    //}

    //for(int n = 21; n < N - 21; n++) {
    //    cout << Z.Z(n * delta) << endl;
    //}
    
    Z.find_zeros(1.0, 39.0, .001, true);

    return 0;
}
