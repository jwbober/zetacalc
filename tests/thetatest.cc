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

void report(complex<double> error, string test_string) {
    report_mutex.lock();
    double logerror = -log2(abs(error));
    if(logerror < histogram_start) {
        error_histogram[0]++;
    }
    else if(logerror >= histogram_size - histogram_start){
        error_histogram[histogram_size - 1]++;
    }
    else {
        error_histogram[int(floor(logerror)) + histogram_start]++;
    }

    if(abs(error) > abs(maxerror)) {
        maxerror = error;
        largest_error_string = test_string;
        cout << "Worst so far:" << endl;
        cout << test_string;
    }
    report_mutex.unlock();
}

void print_histogram() {

    int hist_max_entry = 0;
    for(int n = 0; n < histogram_size; n++) {
        hist_max_entry = max(error_histogram[n], hist_max_entry);
    }

    double hist_normalization = 1.0;

    if(hist_max_entry > 70) {
        hist_normalization = hist_max_entry/70.0;
    }

    cout << "Error histogram:" << endl;
    for(int n = 0; n < histogram_size; n++) {
        cout << n + histogram_start << "\t";
        for(int m = 0; m < ceil(error_histogram[n]/hist_normalization); m++) {
            cout << "+";
        }
        cout << endl;
    }
}


complex<double> test_one(complex<double> * S1,
                         complex<double> * S2,
                         double a,
                         double b,
                         int j,
                         int K,
                         complex<double> * v,
                         double epsilon,
                         string &teststring) {

    *S1 = compute_exponential_sums(a, b, j, K, v, epsilon, 100, 0);
    *S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 4);
 
    stringstream ss;
    ss << std::setprecision(17);
    ss << "./thetatest one ";
    ss << a << " ";
    ss << b << " ";
    ss << j << " ";
    ss << K << " ";
    ss << epsilon << " ";
    for(int k = 0; k < j; k++) {
        ss << "\\(" << v[k].real() << "," << v[k].imag() << "\\)" << " ";
    }
    ss << "\\(" << v[j].real() << "," << v[j].imag() << "\\)" << endl;

    complex<double> error = *S1 - *S2;
    ss << "               ERROR: " << error << endl;
    ss << "         log2(ERROR): " << log2(abs(error)) << endl;
    ss << "    algorithm answer: " << *S1 << endl;
    ss << "       direct answer: " << *S2 << endl;
    ss << endl;

    teststring = ss.str();

    return error;
}

int test_theta_algorithm(int number_of_tests,
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

        string teststring;
        error = test_one(&S1, &S2, a, b, j, K, v, epsilon, teststring);
        error = S1 - S2;
        report(error, teststring);
    }
    delete [] v;
    return 0;
}


int main(int argc, char ** _argv) {
    vector<string> argv(_argv, _argv + argc);
    cout << std::setprecision(17);

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
        
        cout << "seeding random() with " << seed << "." << endl;
        if(nthreads == 1) {
            test_theta_algorithm(number_of_tests, Kapprox, epsilon, EM, seed);
        }
        else {
            std::thread threads[nthreads];
            for(unsigned int k = 0; k < nthreads; k++) {
                threads[k] = std::thread(test_theta_algorithm, number_of_tests/nthreads, Kapprox, epsilon, EM, seed + k);
            }
            for(unsigned int k = 0; k < nthreads; k++) {
                threads[k].join();
            }
        }

        print_histogram();
        return 0;
    }

    if(argv[1] == "one") {
        //
        // expect input in the form 
        // a, b, j, K, epsilon, v_0, v_1, v_2, ..., v_j
        //
        // Everything is optional, however. As soon as something is missing,
        // it is set to a default value.
        //

        stringstream ss;
        for(unsigned int n = 2; n < argv.size(); n++) {
            ss << argv[n];
            if(n != argv.size() - 1) {
                ss << " ";
            }
        }
        double a = .1234, b = .2345, epsilon = pow(2,-30);
        int j = 18, K = 501;
        ss >> a;
        ss >> b;
        ss >> j;
        ss >> K;
        ss >> epsilon;
        complex<double> * v = new complex<double>[j + 1];
        for(int n = 0; n <= j; n++) {
           v[n] = 1.0;
        }
        for(int n = 0; n <= j; n++) {
            if(ss.eof())
                break;
            ss >> v[n];
        }
        complex<double> S1;
        complex<double> S2;
        string teststring;
        test_one(&S1, &S2, a, b, j, K, v, epsilon, teststring);
        cout << teststring;

        delete [] v;
        return 0;
    }

    //if(argc > 1) {
    //    seed = atoi(argv[1].c_str());
    //}


    return 0;

    // The following had a bad problem at some point.
    //complex<double> v[19];
    //for(int j = 0; j < 19; j++) {
    //    v[j] = 1;
    //}
    //complex<double> S1;
    //complex<double> S2;
    //complex<double> error = test_one(&S1, &S2, -1.1759932897985, -2.12365078739822, 18, 10414, v, pow(2, -30));
    //complex<double> error = test_one(&S1, &S2, -1.1759932897985, -2.12365078739822, 18, 10352, v, pow(2, -30));
    //cout << error << " ";
    //cout << S1 << " ";
    //cout << S2 << endl;

    return 0;
}
