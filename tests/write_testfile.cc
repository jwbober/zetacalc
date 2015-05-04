#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <thread>
#include <mutex>

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <complex>
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

const int histogram_size = 60;
const int histogram_start = 0;
int error_histogram[histogram_size] = {0};
complex<double> maxerror = 0.0;

string largest_error_string;

std::mutex report_mutex;

void report(double a, double b, int j, int K, complex<double> * v, complex<double> S) {
    report_mutex.lock();
    cout.write((char*)&a, sizeof(a));
    cout.write((char*)&b, sizeof(b));
    cout.write((char*)&j, sizeof(j));
    cout.write((char*)&K, sizeof(K));
    cout.write((char*)v, sizeof(complex<double>)*(j+1));
    cout.write((char*)&S, sizeof(S));
    cout.flush();
    report_mutex.unlock();
}

int compute_some_sums(int number_of_tests,
                      int approx_K,
                      double epsilon = pow(2, -29),
                      bool EM = false,
                      unsigned int seed = 0) {

    RNG rng;
    rng.seed(seed);

    const int j_max = 18;
    complex<double> * v = new complex<double>[j_max + 1];

    complex<double> S1;
    complex<double> S2;
    complex<double> error;

    double a,b;

    for(int n = 0; n < number_of_tests; n++) {

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

        complex<double> S = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 4);
        report(a, b, j, K, v, S);
    }
    delete [] v;
    return 0;
}


int main(int argc, char ** _argv) {
    vector<string> argv(_argv, _argv + argc);

    if(argc < 2) {
        cout << "select a test." << endl;
        return 0;
    }

    if(argv[1] == "many") {
        //
        // expect input in the form seed, nthreads, number_of_tests, Kapprox, epsilon, EM
        //
        // a seed of 0 will mean we use the choice seed = time().
        //
        // everything is optional and will be set to default values if not present.
        
        stringstream ss;
        for(unsigned int n = 2; n < argv.size(); n++) {
            ss << argv[n];
            if(n != argv.size() - 1) {
                ss << " ";
            }
        }

        unsigned int seed = 0, nthreads = 1, number_of_tests = 1000, Kapprox = 1000;
        bool EM = false;
        double epsilon = pow(2, -30);

        ss >> seed;
        ss >> nthreads;
        ss >> number_of_tests;
        ss >> Kapprox;
        ss >> epsilon;
        ss >> EM;

        if(seed == 0) {
            seed = time(NULL);
        }
        
        if(nthreads == 1) {
            compute_some_sums(number_of_tests, Kapprox, epsilon, EM, seed);
        }
        else {
            std::thread threads[nthreads];
            for(unsigned int k = 0; k < nthreads; k++) {
                threads[k] = std::thread(compute_some_sums, number_of_tests/nthreads, Kapprox, epsilon, EM, seed + k);
            }
            for(unsigned int k = 0; k < nthreads; k++) {
                threads[k].join();
            }
        }

        return 0;
    }
}
