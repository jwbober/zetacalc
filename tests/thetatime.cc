#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <thread>
#include <random>

#include "theta_sums.h"

typedef std::mt19937 RNG;
std::uniform_real_distribution<double> random01(0, 1);

using std::max;
using std::stringstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::setprecision;
//using std::log2;
using std::complex;

inline double random_double(RNG &rng, const double start = 0, const double end = 1) {
    return random01(rng) * (end - start) + start;
}

inline complex<double> random_complex(RNG &rng, const double start = 0, const double end = 1) {
    return complex<double>(random_double(rng, start, end), random_double(rng, start, end));
}

inline int random_int(RNG &rng, const int start, const int end) {
    //
    // includes endpoint
    //
    return (int)random_double(rng, start, end + 1);
}

int run_theta_algorithm(int Kmin, int number_of_runs,
                         int approx_K,
                         double epsilon = pow(2, -29),
                         bool EM = false,
                         unsigned int seed = 0) {

    RNG rng;
    rng.seed(seed);

    const int j_max = 18;
    complex<double> * v = new complex<double>[j_max + 1];

    complex<double> S1;

    double a,b;
    
    for(int n = 0; n < number_of_runs; n++) {

        int K = random_int(rng, approx_K - 500, approx_K + 501);
        int j = random_int(rng, 0, 18);
    
        if(EM) {
            a = random_double(rng);
            b = random_double(rng)/(2 * K); 
        }
        else {
            a = random_double(rng, -10, 10);
            b = random_double(rng, -10, 10);
        }


        for(int k = 0; k <= j; k++) {
            //v[k] = (random_complex(-1, 1) * 2.0 - complex<double>(1.0, 1.0))/(k*k*k*k + 1.0);
            //v[k] = (random_complex(-1, 1) * 2.0 - complex<double>(1.0, 1.0));
            v[k] = 1.0;
        }

        S1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
    }

    delete [] v;
    return 0;
}


int main(int argc, char ** _argv) {
    cout << setprecision(17);
    vector<string> argv(_argv, _argv + argc);
    
    // expect input in the form seed, nthreads, Kmin, number_of_tests, Kapprox, epsilon, EM
    //
    // a seed of 0 will mean we use the choice seed = time().
    //
    // everything is optional and will be set to default values if not present.
    
    stringstream ss;
    for(unsigned int n = 1; n < argv.size(); n++) {
        ss << argv[n];
        if(n != argv.size() - 1) {
            ss << " ";
        }
    }

    unsigned int seed = 0, nthreads = 1, number_of_tests = 1000, Kapprox = 1000, Kmin = 100;
    bool EM = false;
    double epsilon = pow(2, -30);

    ss >> seed;
    ss >> nthreads;
    ss >> Kmin;
    ss >> number_of_tests;
    ss >> Kapprox;
    ss >> epsilon;
    ss >> EM;

    if(seed == 0) {
        seed = time(NULL);
    }
    
    if(nthreads == 1) {
        run_theta_algorithm(Kmin, number_of_tests, Kapprox, epsilon, EM, seed);
    }
    else {
        std::thread threads[nthreads];
        for(unsigned int k = 0; k < nthreads; k++) {
            threads[k] = std::thread(run_theta_algorithm, Kmin, number_of_tests, Kapprox, epsilon, EM, seed + k);
        }
        for(unsigned int k = 0; k < nthreads; k++) {
            threads[k].join();
        }
    }
    return 0;
}
