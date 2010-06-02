#include "theta_sums.h"
#include "zeta.h"
#include "log.h"

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

int run_exp_itlogn_test() {
    mpfr_t v;
    mpfr_t t;
    mpfr_t M;
    mpz_t MM;
    mpfr_init2(v, 150);
    mpfr_init2(t, 158);
    mpfr_init2(M, 150);
    mpz_init(MM);
    mpfr_set_str(M, "4e6", 10, GMP_RNDN);
    mpfr_get_z(MM, M, GMP_RNDN);
    mpfr_set_str(v, "100000000", 10, GMP_RNDN);
    mpfr_set_str(t, "1e30", 10, GMP_RNDN);


    Complex z1, z2, z3;
    Double x1, x2, x3;
    
    mpz_t n;
    mpz_init(n);

    mpfr_t twopi;
    mpfr_init2(twopi, mpfr_get_prec(t));
    mpfr_const_pi(twopi, GMP_RNDN);

    mpfr_t w1, w2;
    mpfr_init2(w1, mpfr_get_prec(t));
    
    create_exp_itlogn_table(t);

    mpz_set_str(n, "3410873413481", 10);

//    cout << exp_itlogn(n) << endl;
    
    z1 = 0;
    z2 = 0;

    Double tt = mpfr_get_d(t, GMP_RNDN);

    for(int k = 1; k <= 1000000; k++) {
        mpz_add_ui(n, n, 1);
//        mpfr_set_z(w1, n, GMP_RNDN);

//        mpfr_log(w1, w1, GMP_RNDN);
//        mpfr_mul(w1, w1, t, GMP_RNDN);
//        mpfr_fmod(w1, w1, twopi, GMP_RNDN);
//        z1 = exp(I * mpfr_get_d(w1, GMP_RNDN));

        z2 += exp_itlogn3(n);

//        z3 += exp(I * tt * log(k));
//        Double error = abs(z1 - z2);
//        if(error > 1e-14)
//            cout << n << "   " << abs(z1 - z2) << endl;
    }

    cout << z1 << endl;
    cout << z2 << endl;
    cout << z3 << endl;

    cout << exp_itlogn_stats::bigger_than_one << endl;
    cout << exp_itlogn_stats::smaller_than_one << endl;




    return 0;
}


int run_zeta_test() {
    
//    z1 = hardy_Z(t, R);
//    z1 = zeta_sum_basic(t);
//    z2 = zeta_sum(t);
//    z3 = zeta_sum_mpfr(t);

    print_zeta_stats();
    return 0;

}

int run_theta_sums_test() {

    int j = 9;
    Double epsilon = exp(-20);
    Double error_allowance = exp(-18);
    int K = 10000;


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

    return 0;
}

int main() {


    cout << setprecision(15);

    return run_theta_sums_test();
}
