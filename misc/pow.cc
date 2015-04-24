#include "misc.h"

double misc::pow(double x, int exponent) {
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

        default:
            return std::pow(x, exponent);
    }
}
