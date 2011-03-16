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
#include <algorithm>

typedef double Double;
typedef std::complex<double> Complex;
using namespace std;

#include "rs_sum.h"


double sinc(double x);
double blfi_kernel(double u, double c, double epsilon1);

bool reverse_cmp(double x, double y) {
    return x > y;
}

int minus_one_power(mpz_class n) {
    if(n % 2 == 0) {
        return 1;
    }
    else
        return -1;
}

mpfr_class siegel_theta(mpfr_class t) {
    return (t/2) * log(t/(2 * const_pi())) - t/2 - const_pi()/8;
}

mpfr_class N_approx(mpfr_class t) {
    mpfr_class z;
    z = siegel_theta(t)/const_pi();
    z = z + 1;
    return z;
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

    t1 = n;
    t1 = t1/log(t1);
    t2 = n;

    target = n * const_pi();


    double epsilon = .00000000001;

    // we now maintain the invariant that t1 is too
    // small and t2 is too large

    d = t2 - t1;
    while(d > epsilon) {
        t3 = (t1 + t2)/2;
        theta = siegel_theta(t3);
        if(theta < target) {
            t1 = t3;
        }
        else {
            t2 = t3;
        }
        d = t2 - t1;
    }

    return (t1 + t2)/2;
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

    vector<double> zeros;
    mpz_class zero_index;
    double g0;

    ZetaComputation(mpfr_class _t0, double _delta, int _N, Complex * _rs_sum) {
        
        t0 = _t0;

        delta = _delta;
        N = _N;

        rs_sum = new Complex[N];

        for(int n = 0; n < N; n++) {
            rs_sum[n] = _rs_sum[n];
        }

        mpfr_class temp;
        temp = log(t0/(2 * const_pi()));
        tau = .5 * temp.get_d();
        beta = M_PI/delta;
        lambda = (beta + tau)/2.0;
        epsilon1 = (beta - tau)/2.0;
        alpha = tau;


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

    double Z(mpfr_class t) {
        mpfr_class T;
        T = t - t0;
        return Z(T.get_d());
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

        int count = 0;

        zeros.clear();

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


    mpz_class calculate_N(mpfr_class t) {
        mpfr_class starting_g;
        mpfr_class g;
        mpfr_class g2;
        mpfr_class T;

        mpz_class N;
        mpz_class starting_N;

        T = t + t0;
        N = floor(N_approx(T) - 1);
        g = grampoint(N);
        double delta = .001;

        while( minus_one_power(N) * Z(g) < .001 ) {
            N = N + 1;
            g = grampoint(N);
        }

        starting_N = N;

        mpfr_class S_upper_bound;

        mpfr_class a = 2.067;
        mpfr_class b = .059;

        starting_g = g;
        mpfr_class h;
        mpfr_class h_sum;

        h = 0;
        h_sum = 0;
        S_upper_bound = 100;
        while(S_upper_bound > 1.95) {
            N = N + 1;
            g2 = grampoint(N);
            if(g + h < g2 && minus_one_power(N) * Z(g2) > 0) {
                h = 0;
                g = g2;
                S_upper_bound = (a + b * log(g) + h_sum)/(g - starting_g);
                S_upper_bound += 1;
                cout << "Found upper bound of " << S_upper_bound << endl;
            }
            else {
                h = g + h - g2 + delta;
                while(minus_one_power(N) * Z(g2 + h) > 0) {
                    h = h + delta;
                }
                g = g2;
                h_sum += h;
            }
        }

        //cout << "Successfully proved that S(g) <= 0." << endl;

        mpfr_class S_lower_bound;
        h = 0;
        h_sum = 0;
        S_lower_bound = -100;
        N = starting_N;

        while(S_lower_bound < -1.95) {
            N = N - 1;
            g2 = grampoint(N);

            if(g + h > g2 && minus_one_power(N) * Z(g2) > 0) {
                h = 0;
                g = g2;
                S_lower_bound = (a + b * log(g) - h_sum)/(starting_g - g);
                S_lower_bound += 1;
                S_lower_bound = -S_lower_bound;
                cout << "Found lower bound of " << S_lower_bound << endl;
            }
            else {
                h = g + h - g2 - delta;
                while(minus_one_power(N) * Z(g2 + h) > 0) {
                    h = h - delta;
                }
                g = g2;
                h_sum += h;
            }
        }

        //cout << "Successfully proved that S(g) >= 0." << endl;

        zero_index = starting_N + 1;
        g = starting_g - t0;
        g0 = g.get_d();

        return starting_N + 1;

    }

    double calculate_S_values() {
        double current_max = 0;
        cout << setprecision(10);
        mpz_class current_zero_index;
        current_zero_index = zero_index - (lower_bound(zeros.begin(), zeros.end(), g0) - zeros.begin());
        mpfr_class x;

        for(vector<double>::iterator i = zeros.begin(); i < zeros.end(); i++) {
            x = current_zero_index - N_approx(t0 + *i);
            double S_value = x.get_d();
            current_max = max(current_max, abs(S_value));
            cout << *i << " " << S_value << endl;
            S_value += 1;
            current_max = max(current_max, abs(S_value));
            cout << *i << " "  << S_value << endl;
            current_zero_index += 1;
        }

        return current_max;
    }

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

    mpfr_class::set_dprec(300);

    cout << setprecision(200);


    bool filename_set = false;
    const char * filename = "";
    int list_values = 0;
    int list_zeros = 0;
    int list_S_values = 0;
    int list_midpoint_values = 0;
    double spacing = .001;
    double start = 1;
    double end = 39;

    while(1) {
        enum {FILENAME = 2, ZEROS, S_OPTION, START, END, SPACING};
        static struct option options[] = 
            {
                {"filename", required_argument, 0, FILENAME},
                {"values", no_argument, &list_values, 1},
                {"zeros", no_argument, &list_zeros, 1},
                {"S", no_argument, &list_S_values, 1},
                {"start", required_argument, 0, START},
                {"end", required_argument, 0, END},
                {"spacing", required_argument, 0, SPACING},
                {"list_midpoint_values", no_argument, &list_midpoint_values, 1},
                {0, 0, 0, 0}
            }; 

        int option_index = 0;
        int c = getopt_long(argc, argv, "", options, &option_index);
        if (c == -1)
            break;

        switch(options[option_index].val) {
            case FILENAME:
                filename_set = true;
                filename = optarg;
                break;
            case SPACING:
                spacing = atof(optarg);
                break;
            case START:
                start = atof(optarg);
                break;
            case END:
                end = atof(optarg);
                break;
                
        }
    }

    if(!filename_set) {
        cout << "Missing input file name." << endl;
        return 0;
    }

    ifstream infile(filename);
    
    mpfr_class t0;
    double delta;
    int N = 1000;
    complex<double> rs_sum[N];

    infile >> t0;
    infile >> delta;

    for(int n = 0; n < N; n++) {
        infile >> rs_sum[n];
    }

    infile.close();

    ZetaComputation Z(t0, delta, N, rs_sum);

    //cout << mpfr_get_d(t, GMP_RNDN) << endl;
    
    //for(int n = 0; n < 10000; n++) {
    //    cout << Z.Z(1 + 38/10000.0 * n) << endl;
    //}

    //for(int n = 21; n < N - 21; n++) {
    //    cout << Z.Z(n * delta) << endl;
    //}
    
    if(list_zeros)
        Z.find_zeros(start, end, spacing, true);

    if(list_values) {
        cout << setprecision(15);
        double t = start;
        while(t < end) {
            cout << t << " " << Z.Z(t) << endl;
            t = t + spacing;
        }
    }

    if(list_S_values) {
        Z.find_zeros(start, end, spacing, true);
        Z.calculate_N(7);
        double S = Z.calculate_S_values();
        cout << "maximal value of S was " << S << endl;
    }

    if(list_midpoint_values) {
        vector<double> zeros = Z.find_zeros(start, end, spacing, true);
        vector<double> values(0);
        double previous = zeros[0];
        for(vector<double>::iterator i = zeros.begin() + 1; i < zeros.end(); i++) {
            values.push_back( abs(Z.Z( (*i + previous)/2 ) ));
            previous = *i;
        }
        
        sort(values.begin(), values.end(), reverse_cmp);
        for(vector<double>::iterator i = values.begin(); i < values.end(); i++) {
            cout << *i << endl;
        }
    }

    return 0;
}
