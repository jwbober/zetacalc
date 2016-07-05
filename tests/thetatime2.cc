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

#include "theta_sums.h"

using std::max;
using std::stringstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::setprecision;
//using std::log2;
using std::complex;

int main(int argc, char ** _argv) {
    cout << setprecision(17);
    vector<string> argv(_argv, _argv + argc);

    // expect input in the form a_subdivisions, b_subdivisions, J, Kstart, Kend, multiplier, Kmin, epsilon

    stringstream ss;
    for(unsigned int n = 1; n < argv.size(); n++) {
        ss << argv[n];
        if(n != argv.size() - 1) {
            ss << " ";
        }
    }

    unsigned int a_subdivisions = 100;
    unsigned int b_subdivisions = 100;
    unsigned int J = 6;
    double epsilon = 1e-15;
    unsigned long Kstart = 10000;
    unsigned long Kend = 10001;
    double multiplier = 1.1;
    unsigned int Kmin = 100;

    ss >> a_subdivisions;
    ss >> b_subdivisions;
    ss >> J;
    ss >> Kstart;
    ss >> Kend;
    ss >> multiplier;
    ss >> Kmin;
    ss >> epsilon;

    complex<double> v[J + 1];
    complex<double> S1 = 0.0;

    unsigned long K = Kstart;
    while(K < Kend) {
        clock_t starttime = clock();
        double astep = .25/a_subdivisions;
        double bstep = .25/b_subdivisions;
        long count = 0;
        for(double a = astep; a < .25; a += astep) {
            for(double b = bstep; b < .25; b += bstep) {
                count++;
                //cout << a << " " << b << " " << J << " " << K << " " << epsilon << " " << Kmin << endl;
                for(int j = 0; j <= J; j++) v[j] = complex<double>(1);
                S1 += compute_exponential_sums(a, b, J, K, v, epsilon, Kmin, 0);
            }
        }
        clock_t endtime = clock();
        cout << K << " " << ((endtime - starttime)/(double)CLOCKS_PER_SEC)/count << endl;
        K = K * multiplier;
    }
    return 0;
}
