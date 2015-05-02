#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <complex>

#include "theta_sums.h"

using std::max;
using std::stringstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;
//using std::log2;
using std::complex;

inline double random_double(const double start = 0, const double end = 1) {
    return ((double)std::rand()/((double)RAND_MAX + 1.0)) * (end - start) + start;
}

inline complex<double> random_complex(const double start = 0, const double end = 1) {
    return complex<double>(random_double(start, end), random_double(start, end));
}

inline int random_int(const int start, const int end) {
    //
    // includes endpoint
    //
    return (int)random_double(start, end + 1);
}


complex<double> test_one(complex<double> * S1,
                         complex<double> * S2,
                         double a,
                         double b,
                         int j,
                         int K,
                         complex<double> * v,
                         double epsilon ) {

    *S1 = compute_exponential_sums(a, b, j, K, v, epsilon, 100, 0);
    *S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 4);
    
    cout << "./thetatest one ";
    cout << a << " ";
    cout << b << " ";
    cout << j << " ";
    cout << K << " ";
    cout << epsilon << " ";
    for(int k = 0; k < j; k++) {
        cout << "\\(" << v[k].real() << "," << v[k].imag() << "\\)" << " ";
    }
    cout << "\\(" << v[j].real() << "," << v[j].imag() << "\\)" << endl;

    complex<double> error = *S1 - *S2;
    cout << "               ERROR: " << error << endl;
    cout << "         log2(ERROR): " << log2(abs(error)) << endl;
    cout << "    algorithm answer: " << *S1 << endl;
    cout << "       direct answer: " << *S2 << endl;
    cout << endl;

    return error;
}

int test_theta_algorithm(int number_of_tests,
                         int approx_K,
                         double epsilon = pow(2, -29),
                         bool EM = false,
                         int run_only = -1) {

    const int j_max = 18;
    complex<double> * v = new complex<double>[j_max + 1];

    cout << endl;
    cout << endl;
    cout << endl;
    cout << "Testing theta algorithm with various random parameters ";
    cout <<         number_of_tests << " times, with log2(epsilon) = " << log2(epsilon) << endl;

    complex<double> maxerror = 0.0;
    complex<double> S1, S1error;
    complex<double> S2, S2error;
    complex<double> error;

    double a,b;
    
    const int histogram_size = 60;
    const int histogram_start = 0;
    int error_histogram[histogram_size];
    for(int n = 0; n < histogram_size; n++) {
        error_histogram[n] = 0;
    }

    int bad_test_number = -1;
    for(int n = 0; n < number_of_tests; n++) {

        int K = random_int(approx_K - 500, approx_K + 501);
        int j = random_int(0, 18);
    
        if(EM) {
            a = random_double();
            b = random_double()/(2 * K); 
        }
        else {
            a = random_double(-10, 10);
            b = random_double(-10, 10);
        }


        for(int k = 0; k <= j; k++) {
            //v[k] = (random_complex(-1, 1) * 2.0 - complex<double>(1.0, 1.0))/(k*k*k*k + 1.0);
            //v[k] = (random_complex(-1, 1) * 2.0 - complex<double>(1.0, 1.0));
            v[k] = 1.0;
        }

        if(run_only == -1 || run_only == n) {
            cout << "Test number " << n << endl;
            error = test_one(&S1, &S2, a, b, j, K, v, epsilon);
            error = S1 - S2;

            double logerror = -log2(abs(S1 - S2));
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
                bad_test_number = n;
                maxerror = error;
                S1error = S1;
                S2error = S2;
            }

            
            //cout    << "Test "  << n
            //        << ": a = " << a
            //        << ", b = " << b
            //        << ", j = " << j
            //        << ", K = " << K
            //        << ":                     log2(error) = " << log2(abs(error)) << endl
            //        << "          error = " << error                              << endl
            //        << "          answer = " << S2                                << endl;
        }
    }

    cout << "Largest error was "    << maxerror                 << endl
         << "    log2(maxerror) = "   << log2(abs(maxerror))    << endl
         << "    requested was "  << log2(epsilon)              << endl
         << "This happened on test number " << bad_test_number  << endl;

    cout << endl;
    cout << endl;
 
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
        // expect input in the form seed, number_of_tests, Kapprox, epsilon, EM
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

        unsigned int seed = 0, number_of_tests = 1000, Kapprox = 1000;
        bool EM = false;
        double epsilon = pow(2, -30);

        ss >> seed;
        ss >> number_of_tests;
        ss >> Kapprox;
        ss >> epsilon;
        ss >> EM;

        if(seed == 0) {
            seed = time(NULL);
        }
        
        srand(seed);
        cout << "seeding rand() with " << seed << "." << endl;
        test_theta_algorithm(number_of_tests, Kapprox, epsilon, EM);
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
        test_one(&S1, &S2, a, b, j, K, v, epsilon);

        delete [] v;
        return 0;
    }

    //if(argc > 1) {
    //    seed = atoi(argv[1].c_str());
    //}


    return 0;

    // The following had a bad problem at some point.
    complex<double> v[19];
    for(int j = 0; j < 19; j++) {
        v[j] = 1;
    }
    complex<double> S1;
    complex<double> S2;
    //complex<double> error = test_one(&S1, &S2, -1.1759932897985, -2.12365078739822, 18, 10414, v, pow(2, -30));
    complex<double> error = test_one(&S1, &S2, -1.1759932897985, -2.12365078739822, 18, 10352, v, pow(2, -30));
    cout << error << " ";
    cout << S1 << " ";
    cout << S2 << endl;

    return 0;
}
