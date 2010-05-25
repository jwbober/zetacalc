#include "theta_sums.h"
#include "zeta.h"

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

void test3(Double epsilon, Double error_allowance, int j, int K) {
    Complex coeffs[j + 1];
    for(int l = 0; l < j; l++) {
        coeffs[l] = 0;
    }
    coeffs[j] = 1;

    int n = 0;

    for(Double a = 0; a < 1.0; a += .0097) {
        for(Double b = 1.0/((Double)2 * K * K); b <= .5; b += .0087) {
            n++;
            if(n % 100 == 0) {
                cout << "Running test number " << n << " with a = " << a << " b = " << b << endl;
            }
            Complex z1, z2;
            z1 = compute_exponential_sums(a, b, j, K, coeffs, epsilon, 0);
            z2 = z1;
            //z2 = compute_exponential_sums(a, b, j, K, coeffs, epsilon, 1);
            if( abs(z1 - z2) > error_allowance) {
                cout << "For a = " << a << " b = " << b << " K = " << K << " epsilon = " << epsilon << ": " << endl;
                cout << "                                       Error was   " << abs(z1 - z2)  << endl;
                cout << "                                       log(error): " << log(abs(z1 - z2)) << endl;
                cout << "  method 0 gives: " << z1 << endl;
                cout << "  method 1 gives: " << z2 << endl;
            }
        }
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

    mpfr_t v;
    mpfr_t t;
    mpfr_init2(v, 150);
    mpfr_init2(t, 150);
    mpfr_set_str(v, "1000000000000000", 10, GMP_RNDN);
    mpfr_set_str(t, "1000000000000000000000000000000", 10, GMP_RNDN);
    int K = 30012;

    

    Complex Z[13];
    compute_taylor_coefficients(t, Z);

    Complex z1 = 0.0;
    Complex z2, z3;

    int total = 10000;
    mpfr_t increment;
    mpfr_init2(increment, 150);
    mpfr_div_si(increment, v, total, GMP_RNDN);
    mpfr_floor(increment, increment);

    for(int l = 0; l < total; l++) {
        if(l % 10 == 0) {
            cout << l << endl;
        }
        //z3 = zeta_block(v, K, t, Z, 2);
        z2 = zeta_block(v, K, t, Z);
        z3 = z2;
        z1 = z1 + z2;
        //cout << z2 << " " << z3 << " " << abs(z2 - z3) << endl;
        mpfr_add(v, v, increment, GMP_RNDN);
    }

    cout << z1 << endl;

//    Complex z1 = zeta_block(v, K, t, Z);
//    cout << z1 << " " << abs(z1) << endl;
//    Complex z2 = zeta_block(v, K, t, Z, 2);
//    cout << z2 << " " << abs(z2) << endl;

//    cout << z1 - z2 << endl;


//    zeta_sum(t);

//    test3(exp(-20), exp(-19), 9, 40001);

/*    
    mpfr_t mp_a;
    mpfr_t mp_b;

    mpfr_init2(mp_a, 100);
    mpfr_set_str(mp_a, ".28", 10, GMP_RNDN);

    mpfr_init2(mp_b, 100);
    mpfr_set_str(mp_b, "0", 10, GMP_RNDN);
    mpfr_set_str(mp_b, ".21005", 10, GMP_RNDN);


    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
       
    int K = 10000;
    int j = 1;

    
    Complex C11 = compute_C11(mp_a, mp_b, K);
    Complex C12 = compute_C12(mp_a, mp_b, K);

    Complex coeffs[101];
    for(int k = 0; k <= 100; k++) {
        coeffs[k] = 0;
    }
    coeffs[0] = 1.0;

    Complex z1, z2;
    
    z1 = compute_exponential_sums(mp_a, mp_b, 0, 10000, coeffs, exp(-20));
    z2 = compute_exponential_sums_directly(mp_a, mp_b, 0, 10000, coeffs, exp(-20));
    cout << z1 << endl;
    cout << z2 << endl;

    cout << z1 - z2 << endl;
    cout << abs(z1 - z2) << endl;
*/
    
    print_stats();

    return 0;
}
