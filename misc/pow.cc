#include "misc.h"

//#define __pow_stats__

#include <iostream>
#ifdef __pow_stats__
class Stats {
public:
    long counter = 0;
    int max_exponent = 0;
    long exponent_counts[41];
    long out_of_range = 0;
    void report(int exponent) {
        counter++;
        if(exponent > max_exponent) max_exponent = exponent;
        if(exponent < 41) {
            exponent_counts[exponent]++;
        }
        else {
            out_of_range++;
        }
    }
    Stats() {
        for(int k = 0; k < 41; k++) exponent_counts[k] = 0;
    }
    ~Stats() {
        std::cout << "pow called " << counter << " times." << std::endl;
        std::cout << "Max exponent was " << max_exponent << "." << std::endl;
        std::cout << "Out of range exponents: " << out_of_range << "." << std::endl;
        std::cout << "Exponent histogram:" << std::endl;
        long max_count = 0;
        for(int k = 0; k < 41; k++) {
            if(exponent_counts[k] > max_count) max_count = exponent_counts[k];
        }
        for(int k = 0; k < 41; k++) {
            std::cout << k << " ";
            for(int j = 0; j < exponent_counts[k]/(double)max_count * 20; j++) {
                std::cout << '*';
            }
            std::cout << std::endl;
        }

    }
};
Stats stats;
#endif // __pow_stats__

double misc::pow(double x, int exponent) {
    //std::cout << x << " " << exponent << std::endl;
    if(exponent < 0) {
        exponent = -exponent;
        x = 1.0/x;
    }
#ifdef __pow_stats__
    stats.report(exponent);
#endif
    switch(exponent) {
        case 0:
            return 1;
        case 1:
            return x;
        case 2:
            return x*x;
        case 3:
            return x*x*x;
        case 4:
            {double x2 = x*x; return x2*x2;}
        case 5:
            {double x2 = x*x; return x2*x2*x;}
        case 6:
            {double x3 = x*x*x; return x3*x3;}
        case 7:
            {double x3 = x*x*x; return x3*x3*x;}
        case 8:
            {double x2 = x*x; double x4 = x2*x2; return x4*x4;}
        case 9:
            {double x2 = x*x; double x4 = x2*x2; return x4*x4*x;}
        case 10:
            {double x2 = x*x; double x4 = x2*x2; return x4*x4*x2;}
        case 11:
            {double x2 = x*x; double x4 = x2*x2; return x4*x4*x2*x;}
        case 12:
            {double x2 = x*x; double x4 = x2*x2; return x4*x4*x4;}
        case 13:
            {double x2 = x*x; double x4 = x2*x2; return x4*x4*x4*x;}
        case 14:
            {double x2 = x*x; double x4 = x2*x2; return x4*x4*x4*x2;}
        case 15:
            {double x3 = x*x*x; double x6 = x3*x3; return x6*x6*x3;}
        case 16:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4; return x8*x8;}
        case 17:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4; return x8*x8*x;}
        case 18:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4; return x8*x8*x2;}
        case 19:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4; return x8*x8*x2*x;}
        case 20:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4; return x8*x8*x4;}
        case 21:
            {double x3 = x*x*x; double x6 = x3*x3; return x6*x6*x6*x3;}
        case 22:
            {double x2 = x*x; double x4 = x2*x2; double x11 = x4*x4*x2*x; return x11*x11;}
        case 23:
            {double x2 = x*x; double x3 = x*x2; double x5 = x2*x3;
             double x10 = x5*x5; return x10*x10*x3;}
        case 24:
            {double x2 = x*x; double x4 = x2*x2; double x12 = x4*x4*x4; return x12*x12;}
        case 25:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4; return x8*x8*x8*x;}
        case 26:
            {double x2 = x*x; double x4 = x2*x2; double x13 = x4*x4*x4*x; return x13*x13;}
        case 27:
            {double x2 = x*x; double x4 = x2*x2; double x9 = x4*x4*x; return x9*x9*x9;}
        case 28:
            {double x2 = x*x; double x4 = x2*x2; double x14 = x4*x4*x4*x2; return x14*x14;}
        case 29:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4; return x4*x8*x8*x8*x;}
        case 30:
            {double x3 = x*x*x; double x6 = x3*x3; double x15 = x6*x6*x3; return x15*x15;}
        case 31:
            {double x3 = x*x*x; double x6 = x3*x3; double x15 = x6*x6*x3; return x15*x15*x;}
        case 32:
            {double x2 = x*x; double x4 = x2*x2;
             double x8 = x4*x4; double x16 = x8*x8; return x16*x16;}
        case 33:
            {double x2 = x*x; double x4 = x2*x2;
             double x8 = x4*x4; double x16 = x8*x8; return x*x16*x16;}
        case 34:
            {double x2 = x*x; double x4 = x2*x2;
             double x8 = x4*x4; double x17 = x8*x8*x; return x17*x17;}
        case 35:
            {double x2 = x*x; double x4 = x2*x2; double x9 = x4*x4*x;
             double x13 = x4*x9; return x13*x13*x9;}
        case 36:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4;
             double x18 = x8*x8*x2; return x18*x18;}
        case 37:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4;
             double x18 = x8*x8*x2; return x*x18*x18;}
        case 38:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4;
             double x19 = x8*x8*x2*x; return x19*x19;}
        case 39:
            {double x2 = x*x; double x3 = x*x2; double x6 = x3*x3;
             double x12 = x6*x6; return x12*x12*x12*x3;}
        case 40:
            {double x2 = x*x; double x4 = x2*x2; double x8 = x4*x4;
             double x20 = x8*x8*x4; return x20*x20;}

        default:
            return std::pow(x, exponent);
    }
}
