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
    mpfr_t M;
    mpz_t MM;
    mpfr_init2(v, 150);
    mpfr_init2(t, 110);
    mpfr_init2(M, 150);
    mpz_init(MM);
    mpfr_set_str(M, "4e6", 10, GMP_RNDN);
    mpfr_get_z(MM, M, GMP_RNDN);
    mpfr_set_str(v, "100000000", 10, GMP_RNDN);
    mpfr_set_str(t, "1e15", 10, GMP_RNDN);
//    mpfr_set_str(t, "1111010010045", 10, GMP_RNDN);

    mpz_set_str(MM, "1", 10);
//    Complex zz = zeta_block_d_stupid(MM, 12615662, t);
    
//    cout << zz << endl;
//    return 0;

    Complex z1, z2, z3;
    Double x1, x2, x3;
    
    mpz_t n;
    mpz_init(n);

    mpfr_t twopi;
    mpfr_init2(twopi, mpfr_get_prec(t));
    mpfr_const_pi(twopi, GMP_RNDN);


//    make_tlog_table(t, 50000);

    mpz_set_str(n, "1", 10);

    x1 = 0;

    for(int l = 1; l <= 300000; l++) {
        mpz_add_ui(n, n, 1);
//        x1 = fmod(tlog(t, n), 2 * PI);

//        x1 += tlog(t, n);

        mpfr_set_z(M, n, GMP_RNDN);
        mpfr_log(M, M, GMP_RNDN);
        mpfr_mul(M, M, t, GMP_RNDN);

        
        mpfr_fmod(M, M, twopi, GMP_RNDN);

        x2 += mpfr_get_d(M, GMP_RNDN);

//        cout << n << "   " << x1 - x2 << endl;
    }
    
    cout << x2 << endl;

    return 0;
    

    Double vv = mpfr_get_d(v, GMP_RNDN);
    Double tt = mpfr_get_d(t, GMP_RNDN);

    Complex R;

    z1 = hardy_Z(t, R);
//    z1 = zeta_sum_basic(t);
//    z2 = zeta_sum(t);
//    z3 = zeta_sum_mpfr(t);
    cout << z1 << endl;
    cout << z1 * R << endl;
//    cout << z2 << endl;
//    cout << z3 << endl;

//    cout << abs(z2 - z3) << endl;

    print_zeta_stats();
    return 0;
}
