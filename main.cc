#include "theta_sums.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

using namespace std;


namespace stats {
    int H_method1 = 0;
    int H_method2 = 0;
    int H_method3 = 0;
    int H_method4 = 0;

    int G_method1 = 0;
    int G_method2 = 0;

    int exp = 0;

    int exponential_sum_called = 0;
    int exponential_sum_euler_maclaurin = 0;
    int exponential_sum_taylor_expansion = 0;

    int H_Integral_0;
    int H_Integral_2;
    int J_Integral_0;
    int J_Integral_1;
    int J_Integral_2;
}





void test(int K, Double epsilon, Double error_allowance) {
    for(Double a = 0; a < 1.0; a += .01) {
        for(Double b = 1.0/((Double)2 * K); b <= .25; b += .01) {
            Complex z1, z2;
            z1 = compute_exponential_sum(a, b, K, epsilon, 2);
            z2 = compute_exponential_sum(a, b, K, epsilon, 2);
            if( abs(z1 - z2) > error_allowance) {
                cout << "For a = " << a << " b = " << b << " K = " << K << " epsilon = " << epsilon << ": " << endl;
                cout << "                                       Error was   " << abs(z1 - z2)  << endl;
                cout << "                                       log(error): " << log(abs(z1 - z2)) << endl;
                cout << "  method 1 gives: " << z1 << endl;
                cout << "  method 2 gives: " << z2 << endl;
            }
        }
    }
}

void test2(Double epsilon, Double error_allowance) {
    for(int n = 0; n < 100; n++) {
        Double a = (Double)rand()/(Double)RAND_MAX;
        Double b = (Double)rand()/(Double)RAND_MAX * .25;

        int K = to_int((Double)rand()/(Double)RAND_MAX * 100000);

        Complex z1 = compute_exponential_sum(a, b, K, epsilon);
        Complex z2 = direct_exponential_sum_evaluation(a, b, 0, K, 300);

        cout << z1 << endl;
        cout << z2 << endl;

        cout << a << "  "  << b << "  " << K << "  " << abs(z1 - z2) << "   " << log(abs(z1 - z2)) << endl;
    }
}

void timing_test(Double epsilon, int number_of_tests) {
    
}

/*
int main() {
    cout << setprecision(15);

    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);
    mpfr_set_str(mp_a, ".13", 10, GMP_RNDN);
    int K = 10003;
    mpfr_set_str(mp_b, ".001", 10, GMP_RNDN);
    //mpfr_div_si(mp_b, mp_b, K, GMP_RNDN);
    //mpfr_div_ui(mp_b, mp_b, 4, GMP_RNDN);

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);

    Complex C11 = compute_C11(mp_a, mp_b, K);
    Complex C12 = compute_C12(mp_a, mp_b, K);

    //cout << C12 << endl;

    Double epsilon = exp(-20);

    Complex S = 0;

    //for(int n = 0; n < 10; n++) {
    //    S = S + compute_exponential_sum(mp_a, mp_b, K, epsilon);
        //S = S + direct_exponential_sum_evaluation(a, b, 0, K);
    //}
    S = compute_exponential_sum(mp_a, mp_b, K, epsilon);
    cout << S << endl;
    cout << direct_exponential_sum_evaluation(mp_a, mp_b, 0, K) << endl;
    cout << direct_exponential_sum_evaluation(a, b, 0, K) << endl;

    print_stats();


}
*/

int main() {


    cout << setprecision(15);

//    cout << G(.2, .5, 0, 10, exp(-20), 1) << endl;
//    cout << G(.2, .5, 0, 10, exp(-20), 2) << endl;

    //cout << J_Integral_1(.41, .1, 1, 20, 400, exp(-60)) << endl;

    //cout << J_Integral_2(.2, .3, .2, 1, 400, exp(-20)) << endl;
    
    //cout << G(.1, .1, 5, 20, 100) << endl;
    //cout << G(.1, .1, 5, 20, exp(-20)) << endl;
    
    mpfr_t mp_a;
    mpfr_t mp_b;

    mpfr_init2(mp_a, 100);
    mpfr_set_str(mp_a, ".5134513", 10, GMP_RNDN);

    mpfr_init2(mp_b, 100);
    mpfr_set_str(mp_b, "0", 10, GMP_RNDN);
    mpfr_set_str(mp_b, ".124131435123", 10, GMP_RNDN);


    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
       
    int K = 10000;
    int j = 1;

    
    Complex C11 = compute_C11(mp_a, mp_b, K);
    Complex C12 = compute_C12(mp_a, mp_b, K);

    //cout << IC0(K, 5, a, b, C11, C12, mp_a, mp_b, exp(-20)) << endl;
    //cout << IC0(K, 6, a, b, C11, C12, mp_a, mp_b, exp(-20)) << endl;
    //cout << IC0(K, 3, a, b, C11, C12, mp_a, mp_b, exp(-20)) << endl;
    //cout << IC0(K, 11, a, b, C11, C12, mp_a, mp_b, exp(-20)) << endl;
    //cout << IC0(K, a, b, C11, C12, mp_a, mp_b, exp(-20)) << endl;

    int m = 0;
    int M = -1;

    //cout << H_Integral_0(0, .1, -1, exp(-20)) << endl;

    //cout << JBoundary(.2, .1, .2, 2, 100, exp(-20)) << endl;

    //Complex coeffs[] = {0, 0, 0, 1};
    
    Complex coeffs[101];
    for(int k = 0; k <= 100; k++) {
        coeffs[k] = 0;
    }
    coeffs[0] = 1.0;
    coeffs[1] = 1.0;
    coeffs[2] = 1.0;
    coeffs[3] = 2.0 * I;
    coeffs[4] = 1.0;
    coeffs[5] = 1.0;
    coeffs[6] = 4.1;
    coeffs[7] = 1.0 + I;
    coeffs[8] = 1.2;
    coeffs[9] = -5;
    coeffs[10] = I;


    Complex coeffs2[] = {1};

    Complex z1, z2, z3, z4;

    
    z1 = compute_exponential_sums(mp_a, mp_b, 10, 100000, coeffs, exp(-20));
//    z2 = compute_exponential_sums(mp_a, mp_b, 0, 100000, coeffs2, exp(-20));
//    z3 = direct_exponential_sum_evaluation2(mp_a, mp_b, 1, 0, 100000);
//    z4 = direct_exponential_sum_evaluation2(mp_a, mp_b, 1, 0, 100000);
    z4 = compute_exponential_sums_directly(mp_a, mp_b, 10, 100000, coeffs, exp(-40));
    //z4 = direct_exponential_sum_evaluation(a, b, 0, 1000);
    cout << z1 << endl;
//    cout << z2 << endl;
//    cout << z3 << endl;
    cout << z4 << endl;

    cout << z1 - z4 << endl;
    cout << abs(z1 - z4) << endl;

//    cout << IC0(K, 1, a, b, C11, C12, mp_a, mp_b, exp(-20)) << endl;
//    cout << IC0(K, 1, a, b, C11, C12, mp_a, mp_b, exp(-40)) << endl;
//    cout << IC0(K, 1, a, b, C11, C12, mp_a, mp_b, exp(-80)) << endl;
//    cout << IC0(K, a, b, C11, C12, mp_a, mp_b, exp(-20)) << endl;

//    for(int s = 0; s <= 5; s++) {
//        for(int j = s; j <= 5; j++) {
//            cout << s << ", " << j << ": " << w_coefficient(mp_a, mp_b, 10000, s, j) << endl;
//        }
//    }


//    cout << JBulk(.99, .12, 0, 23999, 100000, exp(-60)) << endl;
//    cout << IC7(100000, 1, .99, .12, exp(-60)) << endl;

    //cout << JBoundary(.99, .21, .12, 1, 100000, exp(-20)) << endl;
    //cout << JBoundary(.99, .21, .12, 1, 1000, exp(-20)) << endl;
    //cout << JBoundary(.99, .21, .12, 1, 100, exp(-20)) << endl;
    //cout << JBoundary(.99, .21, .12, 1, 10, exp(-20)) << endl;

    //cout << J_Integral_0(.1, .2, 2, -1, 100, exp(-20)) << endl;

    //cout << sum_of_offset_inverse_powers(a, m, M, j, 1) << endl;
    //cout << sum_of_offset_inverse_powers(a, m, M, j, exp(-1)) << endl;
    //cout << sum_of_offset_inverse_powers(a, m, M, j, exp(-30)) << endl;
    //cout << sum_of_offset_inverse_powers(a, m, M, j, exp(-40)) << endl;

/* 
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 300);
    mpfr_init2(mp_b, 300);
    mpfr_set_str(mp_a, ".24288677062973695886", 10, GMP_RNDN);
    mpfr_set_str(mp_b, ".034307894196504679085", 10, GMP_RNDN);
    int K = 80417;
    Double epsilon = exp(-30);

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);


    Complex z1 = direct_exponential_sum_evaluation(mp_a, mp_b, 0, K);
    Complex z2 = compute_exponential_sum(mp_a, mp_b, K, epsilon);
    Double z3 = abs(z2 - z1);
    cout << z1 << endl;
    cout << z2 << endl;
    cout << endl;
    cout << z3 << endl;    
    cout << -log(z3) << endl;    
    cout << endl;
*/

 //   srand(time(NULL));

    
    //test2(exp(-20), exp(-20));
//    return 0;
/*
    mpfr_t mp_a, mp_b;
    mpfr_init2(mp_a, 100);
    mpfr_init2(mp_b, 100);
    
    int K = 10000013;

    mpfr_set_str(mp_a, "1", 10, GMP_RNDN);
    mpfr_div_si(mp_a, mp_a, 99, GMP_RNDN);
    mpfr_set_str(mp_b, "1", 10, GMP_RNDN);
    mpfr_div_si(mp_b, mp_b, K, GMP_RNDN);
    mpfr_div_si(mp_b, mp_b, K, GMP_RNDN);

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
    Double epsilon = exp(-20);

    Complex S = (Complex)0;

    int n = 0;
    
    while(a < 1) {
        while(b < .25) {
            S = S + compute_exponential_sum(a, b, K, epsilon);
            //S = S + direct_exponential_sum_evaluation(a, b, 0, K);
            b = b + 1000.1/(Double)K;
            n++;
        }
        a = a + .09;
        b = 1/(Double)K + .000002;
    }
*/

//    cout << H(13, Complex(2, 1), exp(-25)) << endl;
//    cout << H_method3(13, Complex(2, 1), exp(-25)) << endl;
    
//    cout << "---"  << IC0(K, a, b, 0.0, 0.0, mp_a, mp_b, epsilon) << endl;

   // cout << compute_exponential_sum(a, b, K, epsilon) << endl;
   // cout << direct_exponential_sum_evaluation(a, b, 0, K) << endl;

    print_stats();

//    cout << n << endl;
//    cout << S << endl;
    
    return 0;
}
