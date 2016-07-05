#include <iostream>
#include <fstream>
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

using std::max;
using std::stringstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::cerr;
//using std::log2;
using std::complex;

const int j_max = 20;

const int histogram_size = 60;
const int histogram_start = 0;
int error_histogram[histogram_size] = {0};
int relerror_histogram[histogram_size] = {0};
double maxerror = 0.0;
unsigned long testcount = 0;
std::mutex report_mutex;
void report(double a, double b, int j, int K, complex<double> * v, double epsilon, complex<double> S, complex<double> S1) {
    std::lock_guard<std::mutex> lock(report_mutex);
    testcount++;
    complex<double> diff = S - S1;
    double abserror = abs(diff);
    double relerror = abserror/abs(S);
    //double error = std::min(abserror, relerror);
    double error = abserror;
    double logerror = -log2(error);
    double logrelerror = -log2(relerror);
    if(logerror < histogram_start) {
        error_histogram[0]++;
    }
    else if(logerror >= histogram_size - histogram_start){
        error_histogram[histogram_size - 1]++;
    }
    else {
        int z = int(floor(logerror)) + histogram_start;
        if(z < 0) {
            z = 0;
        }
        if(z >= histogram_size) z = histogram_size - 1;
        error_histogram[z]++;
    }

    if(logrelerror < histogram_start) {
        relerror_histogram[0]++;
    }
    else if(logrelerror >= histogram_size - histogram_start){
        relerror_histogram[histogram_size - 1]++;
    }
    else {
        int z = int(floor(logrelerror)) + histogram_start;
        if(z < 0) {
            z = 0;
        }
        if(z >= histogram_size) z = histogram_size - 1;
        relerror_histogram[z]++;
    }

    if(error > maxerror || isnan(error) || maxerror == 0) {
        maxerror = error;
        cout << "Worst so far:" << endl;
        cout << std::setprecision(17);
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

        cout << "      RELATIVE ERROR: " << relerror << endl;
        cout << "               ERROR: " << abserror << endl;
        cout << "          DIFFERENCE: " << diff << endl;
        cout << "         log2(ERROR): " << log2(abserror) << endl;
        cout << "log2(RELATIVE ERROR): " << log2(relerror) << endl;
        cout << "    algorithm answer: " << S1 << endl;
        cout << "       direct answer: " << S << endl;
        cout << endl;
    }
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

    int first_nonzero = 0;
    int last_nonzero = 0;
    for(int n = 0; n < histogram_size; n++) {
        if(error_histogram[n] > 0) last_nonzero = n;
    }
    for(int n = histogram_size - 1; n >= 0; n--) {
        if(error_histogram[n] > 0) first_nonzero = n;
    }

    for(int n = first_nonzero; n <= last_nonzero; n++) {
        cout << n + histogram_start << " " << error_histogram[n] << endl;
    }
    cout << endl;
    cout << endl;

    for(int n = first_nonzero; n <= last_nonzero; n++) {
        cout << n + histogram_start << "\t";
        for(int m = 0; m < ceil(error_histogram[n]/hist_normalization); m++) {
            cout << "+";
        }
        cout << endl;
    }

    cout << "Relative error histogram:" << endl;

    first_nonzero = 0;
    last_nonzero = 0;
    for(int n = 0; n < histogram_size; n++) {
        if(relerror_histogram[n] > 0) last_nonzero = n;
    }
    for(int n = histogram_size - 1; n >= 0; n--) {
        if(relerror_histogram[n] > 0) first_nonzero = n;
    }
    for(int n = first_nonzero; n <= last_nonzero; n++) {
        cout << n + histogram_start << "\t";
        for(int m = 0; m < ceil(relerror_histogram[n]/hist_normalization); m++) {
            cout << "+";
        }
        cout << endl;
    }

    cout << "Completed " << testcount << " tests." << endl;
}

ifstream infile;
std::mutex next_mutex;
void get_next(double * a, double * b, int * j, int * K, complex<double> * v, complex<double> * S) {
    std::lock_guard<std::mutex> lock(next_mutex);
    if(!infile.read((char*)a, sizeof(double))) {*j = -1; return;}
    if(!infile.read((char*)b, sizeof(double))) {*j = -1; return;}
    if(!infile.read((char*)j, sizeof(int))) {*j = -1; return;}
    if(*j < 0 || *j > j_max) {
        cout << "possible error reading input." << endl;
        cout << "got j = " << *j << endl;
        cout << "largest allowed j is " << j_max << endl;
        *j = -1; return;
    }
    if(!infile.read((char*)K, sizeof(int))) {*j = -1; return;}
    if(!infile.read((char*)v, sizeof(complex<double>) * (*j + 1))) {*j = -1; return;}
    if(!infile.read((char*)S, sizeof(complex<double>))) {*j = -1; return;}
}

void run_tests(double epsilon, const int Kmin) {
    double a, b;
    int j, K;
    complex<double> v[j_max + 1];
    for(int j = 0; j <= j_max; j++) v[j] = 0.0;
    complex<double> S;

    get_next(&a, &b, &j, &K, v, &S);
    while(j != -1) {
        complex<double> S1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
        report(a, b, j, K, v, epsilon, S, S1);
        get_next(&a, &b, &j, &K, v, &S);
    }
}


int main(int argc, char ** argv) {
    if(argc < 2) {
        cerr << "Missing input." << endl;
        return -1;
    }
    infile.open(argv[1]);
    double epsilon = pow(2.0, -30);
    int nthreads = 1;
    if(argc > 2) {
        epsilon = atof(argv[2]);
    }
    if(argc > 3) {
        nthreads = atoi(argv[3]);
    }
    int Kmin = 100;
    if(argc > 4) {
        Kmin = atoi(argv[4]);
    }
    if(nthreads == 1) {
        run_tests(epsilon, Kmin);
    }
    else {
        std::thread threads[nthreads];
        for(int k = 0; k < nthreads; k++) {
            threads[k] = std::thread(run_tests, epsilon, Kmin);
        }
        for(int k = 0; k < nthreads; k++) {
            threads[k].join();
        }
    }
    print_histogram();
    return 0;
}
