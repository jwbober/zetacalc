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

    
    test2(exp(-20), exp(-20));
    return 0;
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

    cout << H(13, Complex(2, 1), exp(-25)) << endl;
    cout << H_method3(13, Complex(2, 1), exp(-25)) << endl;
    
//    cout << "---"  << IC0(K, a, b, 0.0, 0.0, mp_a, mp_b, epsilon) << endl;

   // cout << compute_exponential_sum(a, b, K, epsilon) << endl;
   // cout << direct_exponential_sum_evaluation(a, b, 0, K) << endl;

    print_stats();

//    cout << n << endl;
//    cout << S << endl;
    
    return 0;
}
