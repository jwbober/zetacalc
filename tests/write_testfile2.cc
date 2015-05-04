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
using std::cerr;
using std::endl;
using std::complex;


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

double last_a = 0;
double last_b = 0;
double a_step = .01;
double b_step = .0001;
long total_tests = 1/a_step + .25/b_step;
long count = 0;
int last_percent_printed = 0;

std::mutex next_mutex;
void next_ab(double * a, double * b) {
    std::lock_guard<std::mutex> lock(next_mutex);
    count++;
    last_b += b_step;
    if(last_b > .25) {
        last_b = 0;
        last_a += a_step;
        if(last_a > 1) {*a = -1; *b = -1; return;}
    }
    *a = last_a;
    *b = last_b;
    int percent_finished = (int)((double)count/(double)total_tests * 100);
    if(percent_finished > last_percent_printed) {
        cerr << percent_finished << endl;
        last_percent_printed = percent_finished;
    }
}

int compute_some_sums(int K, int j) {
    complex<double> * v = new complex<double>[j + 1];
    for(int k = 0; k <= j; k++) {
        v[k] = 1.0;
    }

    double a = 1;
    double b = 1;

    next_ab(&a, &b);
    while(a >= 0) {
        complex<double> S = compute_exponential_sums(a, b, j, K, v, pow(2, -50), 0, 4);
        report(a, b, j, K, v, S);
        next_ab(&a, &b);
    }
    delete [] v;
    return 0;
}

int main(int argc, char ** argv) {
    if(argc < 3) {
        cerr << "missing input" << endl;
    }
    int j = atoi(argv[1]);
    int K = atoi(argv[2]);
    int nthreads = 1;
    if(argc > 3) {
        nthreads = atoi(argv[3]);
    }
        
    if(nthreads == 1) {
        compute_some_sums(K, j);
    }
    else {
        std::thread threads[nthreads];
        for(int k = 0; k < nthreads; k++) {
            threads[k] = std::thread(compute_some_sums, K, j);
        }
        for(int k = 0; k < nthreads; k++) {
            threads[k].join();
        }
    }
}
