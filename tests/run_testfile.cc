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

const int j_max = 18;

const int histogram_size = 60;
const int histogram_start = 0;
int error_histogram[histogram_size] = {0};
complex<double> maxerror = 0.0;
std::mutex report_mutex;
void report(double a, double b, int j, int K, complex<double> * v, double epsilon, complex<double> S, complex<double> S1) {
    std::lock_guard<std::mutex> lock(report_mutex);
    complex<double> error = S - S1;
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

        cout << "               ERROR: " << error << endl;
        cout << "         log2(ERROR): " << log2(abs(error)) << endl;
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
    for(int n = 0; n < histogram_size; n++) {
        cout << n + histogram_start << "\t";
        for(int m = 0; m < ceil(error_histogram[n]/hist_normalization); m++) {
            cout << "+";
        }
        cout << endl;
    }
}

ifstream infile;
std::mutex next_mutex;
void get_next(double * a, double * b, int * j, int * K, complex<double> * v, complex<double> * S) {
    std::lock_guard<std::mutex> lock(next_mutex);
    if(!infile.read((char*)a, sizeof(double))) {*j = -1; return;}
    if(!infile.read((char*)b, sizeof(double))) {*j = -1; return;}
    if(!infile.read((char*)j, sizeof(int))) {*j = -1; return;}
    if(*j < 0 || *j > 20) {
        cout << "possible error reading input." << endl;
        cout << "got j = " << *j << endl;
        cout << "largest allowed j is " << j_max << endl;
        *j = -1; return;
    }
    if(!infile.read((char*)K, sizeof(int))) {*j = -1; return;}
    if(!infile.read((char*)v, sizeof(complex<double>) * (*j + 1))) {*j = -1; return;}
    if(!infile.read((char*)S, sizeof(complex<double>))) {*j = -1; return;}
}

void run_tests(double epsilon) {
    double a, b;
    int j, K;
    complex<double> v[j_max];
    complex<double> S;

    get_next(&a, &b, &j, &K, v, &S);
    while(j != -1) {
        complex<double> S1 = compute_exponential_sums(a, b, j, K, v, epsilon, 100, 0);
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
    if(nthreads == 1) {
        run_tests(epsilon);
    }
    else {
        std::thread threads[nthreads];
        for(int k = 0; k < nthreads; k++) {
            threads[k] = std::thread(run_tests, epsilon);
        }
        for(int k = 0; k < nthreads; k++) {
            threads[k].join();
        }
    }
    print_histogram();
    return 0;
}
