#include "theta_sums.h"
#include "log.h"
#include "zeta.h"

#include <ctime>
#include <iostream>
#include <iomanip>
#include "precomputed_tables.h"

using namespace std;

inline double random_double() {
    return (double)rand()/(double)RAND_MAX;   
}
inline complex<double> random_complex() {
    return complex<double>(random_double(), random_double());
}

int test_fastlog() {
    const int number_of_tests = 10000;

    cout << "Testing fastlog() on " << number_of_tests << " uniformly random numbers between 1 and 1000000...";

    int failures1 = 0;
    int smaller1 = 0;
    int exact1 = 0;

    for(int n = 0; n < number_of_tests && failures1 <= 10; n++) {
        double d = (double)rand()/(double)RAND_MAX * 999999.0 + 1.0;
        int a = fastlog(d);
        int b = (int)floor(log(d));
        if(a != b && a+1 != b) {
            failures1++;
            if(failures1 < 11) {
                cout << endl;
                cout << "fastlog(" << d << ") gives wrong result. Got: " << a << ". Expected: " << b << ".";
            }
        }
        if(a == b)
            exact1++;
        else
            smaller1++;
    }
    if(failures1) {
        cout << "done." << endl;
        cout << "**************************************************" << endl;
        cout << "TEST FAILURES FOUND IN fastlog()." << endl;
        cout << "**************************************************" << endl;
    }
    else {
        cout << "OK." << endl;
        cout << "       Answer was exact " << exact1 << " times." << endl;
        cout << "  Answer was off by one " << smaller1 << " times." << endl;
    }

    cout << "Testing fastlog() on " << number_of_tests << " uniformly random numbers between 0 and 1...";

    int failures2 = 0;
    int smaller2 = 0;
    int exact2 = 0;

    for(int n = 0; n < number_of_tests && failures2 <= 10; n++) {
        double d = (rand() + 1.0)/(RAND_MAX + 1.0);
        int a = fastlog(d);
        int b = (int)floor(log(d));
        if(a != b && a + 1 != b) {
            failures2++;
            if(failures2 < 11) {
                cout << endl;
                cout << "fastlog(" << d << ") gives wrong result. Got: " << a << ". Expected: " << b << ".";
            }
        }
        if(a == b)
            exact2++;
        else
            smaller2++;
    }
    if(failures2) {
        cout << "done." << endl;
        cout << "**************************************************" << endl;
        cout << "TEST FAILURES FOUND IN fastlog()." << endl;
        cout << "**************************************************" << endl;
    }
    else {
        cout << "OK." << endl;
        cout << "       Answer was exact " << exact2 << " times." << endl;
        cout << "  Answer was off by one " << smaller2 << " times." << endl;
    }

    cout << "Testing fastlog() on " << number_of_tests << " uniformly random numbers between 0 and .0000001...";

    int failures3 = 0;
    int smaller3 = 0;
    int exact3 = 0;

    for(int n = 0; n < number_of_tests && failures3 <= 10; n++) {
        double d = (rand() + 1.0)/(RAND_MAX + 1.0) * .0000001;
        int a = fastlog(d);
        int b = (int)floor(log(d));
        if(a != b && a + 1 != b) {
            failures3++;
            if(failures3 < 11) {
                cout << endl;
                cout << "fastlog(" << d << ") gives wrong result. Got: " << a << ". Expected: " << b << ".";
            }
        }
        if(a == b)
            exact3++;
        else
            smaller3++;
    }
    if(failures3) {
        cout << "done." << endl;
        cout << "**************************************************" << endl;
        cout << "TEST FAILURES FOUND IN fastlog()." << endl;
        cout << "**************************************************" << endl;
    }
    else {
        cout << "OK." << endl;
        cout << "       Answer was exact " << exact3 << " times." << endl;
        cout << "  Answer was off by one " << smaller3 << " times." << endl;
    }
    

    return failures1 + failures2 + failures3;
}




int test_fastlog2() {
    const int number_of_tests = 10000;

    cout << "Testing fastlog2() on " << number_of_tests << " uniformly random numbers between 1 and 1000000...";

    int failures1 = 0;

    for(int n = 0; n < number_of_tests && failures1 <= 10; n++) {
        double d = (double)rand()/(double)RAND_MAX * 999999.0 + 1.0;
        int a = fastlog2(d);
        int b = (int)floor(log2(d));
        if(a != b) {
            failures1++;
            if(failures1 < 11) {
                cout << endl;
                cout << "fastlog2(" << d << ") gives wrong result. Got: " << a << ". Expected: " << b << ".";
            }
        }
    }
    if(failures1) {
        cout << "done." << endl;
        cout << "**************************************************" << endl;
        cout << "TEST FAILURES FOUND IN fastlog2()." << endl;
        cout << "**************************************************" << endl;
    }
    else {
        cout << "OK." << endl;
    }

    cout << "Testing fastlog2() on " << number_of_tests << " uniformly random numbers between 0 and 1...";

    int failures2 = 0;

    for(int n = 0; n < number_of_tests && failures2 <= 10; n++) {
        double d = (rand() + 1.0)/(RAND_MAX + 1.0);
        int a = fastlog2(d);
        int b = (int)floor(log2(d));
        if(a != b) {
            failures2++;
            if(failures2 < 11) {
                cout << endl;
                cout << "fastlog2(" << d << ") gives wrong result. Got: " << a << ". Expected: " << b << ".";
            }
        }
    }
    if(failures2) {
        cout << "done." << endl;
        cout << "**************************************************" << endl;
        cout << "TEST FAILURES FOUND IN fastlog2()." << endl;
        cout << "**************************************************" << endl;
    }
    else {
        cout << "OK." << endl;
    }

    cout << "Testing fastlog2() on " << number_of_tests << " uniformly random numbers between 0 and .00001...";

    int failures3 = 0;

    for(int n = 0; n < number_of_tests && failures3 <= 10; n++) {
        double d = (rand() + 1.0)/(RAND_MAX + 1.0) * .00001;
        int a = fastlog2(d);
        int b = (int)floor(log2(d));
        if(a != b) {
            failures3++;
            if(failures3 < 11) {
                cout << endl;
                cout << "fastlog2(" << d << ") gives wrong result. Got: " << a << ". Expected: " << b << ".";
            }
        }
    }
    if(failures3) {
        cout << "done." << endl;
        cout << "**************************************************" << endl;
        cout << "TEST FAILURES FOUND IN fastlog2()." << endl;
        cout << "**************************************************" << endl;
    }
    else {
        cout << "OK." << endl;
    }


    return failures1 + failures2 + failures3;
}

int test_theta_algorithm(int number_of_tests, int approx_K, Double epsilon = pow(2, -29), int run_only = -1) {
    const int j_max = 18;

    Complex v[j_max + 1];

    cout << endl;
    cout << endl;
    cout << endl;
    cout << "Testing theta algorithm with various random parameters " << number_of_tests << " times, with log2(epsilon) = " << log2(epsilon) << endl;
    Complex maxerror = 0.0;
    int bad_test_number = -1;
    for(int n = 0; n < number_of_tests; n++) {
        Double a = (double)rand()/(double)RAND_MAX * 20.0 - 10.0;
        Double b = (double)rand()/(double)RAND_MAX * 20.0 - 10.0;

        int K = (int)((double)rand()/(double)RAND_MAX * 500.0 + approx_K);
        int j = (int)((double)rand()/(double)RAND_MAX * j_max);
    
        for(int k = 0; k <= j; k++) {
            //v[k] = random_complex() * 2.0 - complex<double>(1.0, 1.0);
            v[k] = (random_complex() * 2.0 - complex<double>(1.0, 1.0))/(k * k * k * k + 1.0);
        }

        if(run_only == -1) {
            Complex S1 = compute_exponential_sums(a, b, j, K, v, epsilon, 100, 0);
            Complex S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 1);

            Complex error = S1 - S2;

            if(abs(error) > abs(maxerror)) {
                bad_test_number = n;
                maxerror = error;
            }

            cout << "Test " << n << ": a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ":                     log2(error) = " << log2(abs(error)) 
                << endl << "          error = " << error << endl << "          answer = " << S2 << endl;
        }
        else if(run_only == n) {
            Complex S1 = compute_exponential_sums(a, b, j, K, v, epsilon, 100, 0);
            Complex S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 1);

            Complex error = S1 - S2;
            
            cout << "Test " << n << ": a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ": log2(error) = " << log2(abs(error)) << "; error = " << error << ", answer = " << S2 << endl;
        }
    }
    cout << "Largest error was " << maxerror << "; log2(maxerror) = " << log2(abs(maxerror)) << " when requested was " << log2(epsilon) << endl;
    cout << "This happened on test number " << bad_test_number << endl;

    return 0;
}

int test_theta_algorithm_EM_case(int number_of_tests, int approx_K, Double epsilon = pow(2, -29), int run_only = -1) {
    const int j_max = 18;

    Complex v[j_max + 1];

    cout << endl;
    cout << endl;
    cout << endl;
    cout << "Testing theta algorithm with various random parameters in the (mostly) Euler-Maclaurin case, with K ~= " << approx_K << " and log2(epsilon) = " << log2(epsilon) << endl;
    cout << "       Running " << number_of_tests << " tests." << endl;
    Double maxerror = 0.0;
    int bad_test_number = -1;
    for(int n = 0; n < number_of_tests; n++) {
        Double a = (double)rand()/(double)RAND_MAX;

        int K = (int)((double)rand()/(double)RAND_MAX * 500.0 + approx_K);
        int j = (int)((double)rand()/(double)RAND_MAX * j_max);
 
        Double b = ((double)rand()/(double)RAND_MAX) / (2.0 * K);

        for(int k = 0; k <= j; k++) {
            v[k] = (random_complex() * 2.0 - complex<double>(1.0, 1.0))/(k * k * k * k + 1.0);
        }

        if(run_only == -1) {
            Complex S1 = compute_exponential_sums(a, b, j, K, v, epsilon, 100, 0);
            Complex S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 1);

            Double error = abs(S1 - S2);
            if(error > maxerror) {
                bad_test_number = n;
                maxerror = error;
            }
            //cout << "Test " << n << ": a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ": log2(error) = " << log2(error) << "; error = " << error << ", answer = " << S2 << endl;

            cout << "Test " << n << ": a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ":                     log2(error) = " << log2(abs(error)) 
                << endl << "          error = " << error << endl << "          answer = " << S2 << endl;

        }
        else if(run_only == n) {
            Complex S1 = compute_exponential_sums(a, b, j, K, v, epsilon, 100, 0);
            Complex S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 1);

            Double error = abs(S1 - S2);
            maxerror = max(error, maxerror);
            
            cout << "Test " << n << ": a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ": log2(error) = " << log2(error) << "; error = " << error << ", answer = " << S2 << endl;
        }
    }
    cout << "Largest error was " << maxerror << "; log2(maxerror) = " << log2(maxerror) << " when requested was " << log2(epsilon) << endl;
    cout << "This happened on test number " << bad_test_number << endl;

    return 0;
}




void test2() {
    Complex v[10];

    Double a = -9.25347071106242;
    Double b = 7.26188676304272;

    int K = 2462;
    int j = 1;

    K = 120;
    a = .0001;
    //a = 0;
    b = 1.1/K;
    j = 1;


    for(int k = 0; k <= j; k++) {
        v[k] = 1;
    }
    v[0] = 1;
    v[1] = 1;

    Complex S1 = compute_exponential_sums(a, b, j, K, v, pow(2.0, -29), 100, 0);
    Complex S2 = compute_exponential_sums(a, b, j, K, v, pow(2.0, -29), 0, 1);

    Double error = abs(S1 - S2);
    cout << "a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ": log2(error) = " << log2(error) << endl;

}


void test3() {
    Complex v[12];

    Double a = -9.25347071106242;
    Double b = 7.26188676304272;

    int K;
    int j = 10;

    K = 20010;
    a = .817876690448204;
    b = 0.999958985019456;

    for(int k = 0; k <= j; k++) {
        v[k] = 1.0/(k * k + 1);
    }
    v[0] = 1;
    v[1] = 1;

    Complex S1 = compute_exponential_sums(a, b, j, K, v, exp(-20), 700, 0);
    Complex S2 = compute_exponential_sums(a, b, j, K, v, exp(-20), 0, 1);

    Double error = abs(S1 - S2);
    cout << "a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ": log2(error) = " << log2(error) << endl;

}

void test_specific_inputs(Double a, Double b, int K, int j, Double epsilon, int Kmin = 0) {
    Complex v[j + 1];
    for(int l = 0; l <= j; l++) {
        v[l] = 1.0/(l * l * l * l + 1);
    }

    Complex S1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
    Complex S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 1);

    Complex error = S1 - S2;
    cout << "a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ": log2(error) = " << log2(abs(error)) << "; error = " << error << ", answer = " << S2 << endl;

}


double time_theta_algorithm(int j, int K) {
    Complex v[j + 1];
    for(int k = 0; k <= j; k++) {
        v[j] = 1.0;
    }

    clock_t start_time = clock();

    Double epsilon = exp(-20);

    int n = 0;
    Complex z1 = 0.0;

    const int number_of_tests = 1000;

    cout << "Timing theta_algorithm with K = " << K << " and j = " << j << endl;
    cout << "   Running approximately 1000 iterations total." << endl;

    for(Double a = 0; a < .5; a += .5/(number_of_tests/1000.0) ) {
        for(Double b = 1.1/((Double)K); b <= 1.0; b += 1.0/1000.0) {
            n++;
            if(n % 100 == 0) {
                cout << "   Running iteration number " << n << " with a = " << a << " b = " << b << endl;
            }
            z1 += compute_exponential_sums(a, b, j, K, v, epsilon);
        }
    }
    cout << "Sum was " << z1 << endl;

    clock_t end_time = clock();

    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Number of seconds for this run was " << elapsed_time << endl;

    return elapsed_time;
    
}

double time_theta_algorithm_varying_Kmin(int j, int approx_K, int number_of_tests, Double epsilon) {
    Complex v[j + 1];
    for(int k = 0; k <= j; k++) {
        v[k] = 1.0/(k*k*k*k*k*k*k*k + 1);
    }


    Complex z1 = 0.0;

    cout << "Timing theta_algorithm with K ~= " << approx_K << " and j = " << j << endl;

    int Kmin_start = 100;
    int Kmin_end = 1500;
    int Kmin_increment = 100;
    cout << "   Running " << number_of_tests << " iterations total for various Kmin from " << Kmin_start << " to " << Kmin_end << "." << endl;
    for(int Kmin = Kmin_start; Kmin <= Kmin_end; Kmin+=Kmin_increment) {
        clock_t start_time = clock();
        z1 = 0;
        for(int k = 0; k < number_of_tests; k++) {
            int K = rand()/(Double)RAND_MAX * 2500.0 + approx_K;
            Double a = rand()/(Double)RAND_MAX;
            Double b = rand()/(Double)RAND_MAX;
            //z1 += compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
            z1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);


        }
        cout << "Sum was " << z1 << endl;
        clock_t end_time = clock();
        double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
        cout << "Number of seconds with Kmin = " << Kmin << " was " << elapsed_time << endl;
    }

    return 1.0;
    
}


double time_theta_algorithm_zeta_case_varying_Kmin(int j, int approx_K, int number_of_tests, Double epsilon) {
    Complex v[j + 1];
    for(int k = 0; k <= j; k++) {
        v[k] = pow(2 * PI/3.0, k/3) * pow(2, -k)/factorial(k/3);
    }


    Complex z1 = 0.0;

    cout << "Timing theta_algorithm with K ~= " << approx_K << " and j = " << j << endl;

    int Kmin_start = 100;
    int Kmin_end = 1500;
    int Kmin_increment = 100;
    cout << "   Running " << number_of_tests << " iterations total for various Kmin from " << Kmin_start << " to " << Kmin_end << "." << endl;
    for(int Kmin = Kmin_start; Kmin <= Kmin_end; Kmin+=Kmin_increment) {
        clock_t start_time = clock();
        z1 = 0;
        for(int k = 0; k < number_of_tests; k++) {
            int K = rand()/(Double)RAND_MAX * 2500.0 + approx_K;
            Double a = rand()/(Double)RAND_MAX;
            Double b = rand()/(Double)RAND_MAX;
            //z1 += compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
            z1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);


        }
        cout << "Sum was " << z1 << endl;
        clock_t end_time = clock();
        double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
        cout << "Number of seconds with Kmin = " << Kmin << " was " << elapsed_time << endl;
    }

    return 1.0;
    
}

double time_theta_algorithm_zeta_case_fixed_Kmin(int j, int approx_K, int number_of_tests, int Kmin, Double epsilon) {
    Complex v[j + 1];
    for(int k = 0; k <= j; k++) {
        v[k] = pow(2 * PI/3.0, k/3) * pow(2, -k)/factorial(k/3);
    }


    Complex z1 = 0.0;

    cout << "Timing theta_algorithm with K ~= " << approx_K << " and j = " << j << endl;

    cout << "   Running " << number_of_tests << " iterations total with Kmin = " << Kmin << endl;
    clock_t start_time = clock();
    z1 = 0;
    for(int k = 0; k < number_of_tests; k++) {
        int K = rand()/(Double)RAND_MAX * 2500.0 + approx_K;
        Double a = rand()/(Double)RAND_MAX;
        Double b = rand()/(Double)RAND_MAX;
        //z1 += compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
        z1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);


    }
    cout << "Sum was " << z1 << endl;
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Number of seconds with Kmin = " << Kmin << " was " << elapsed_time << endl;

    return 1.0;
    
}






double time_theta_algorithm_EM_case(int j, int K, int number_of_tests, Double epsilon) {
    Complex v[j + 1];
    for(int k = 0; k <= j; k++) {
        v[k] = 1.0/(k*k + 1);
    }


    Complex z1 = 0.0;

    cout << "Timing theta_algorithm euler-maclaurin case with K = " << K << " and j = " << j << endl;
    clock_t start_time = clock();
    z1 = 0;
    for(int k = 0; k < number_of_tests; k++) {
        Double a = (rand()/(Double)RAND_MAX);
        Double b = (rand()/(Double)RAND_MAX)/K;
        //z1 += compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
        z1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);

    }
    cout << "Sum was " << z1 << endl;
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Number of seconds was " << elapsed_time << endl;

    return elapsed_time;
    
}







int test_exp_itlogn(gmp_randstate_t state) {
    mpfr_t t;
    mpfr_init2(t, 158);
    mpfr_set_str(t, "1e30", 10, GMP_RNDN);

    Complex z1, z2;
    
    mpz_t n, m;
    mpz_init(n);
    mpz_init(m);
    mpz_set_str(m, "100000000000000", 10);

    mpfr_t twopi;
    mpfr_init2(twopi, mpfr_get_prec(t));
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_t w1;
    mpfr_init2(w1, mpfr_get_prec(t));
    
    create_exp_itlogn_table(t);

    mpz_set_str(n, "100000000", 10);

//    cout << exp_itlogn(n) << endl;
    
    z1 = 0;
    z2 = 0;

    int failures = 0;

    cout << "Testing exp_itlogn() against mpfr 200000 times ...";
    cout.flush();

    for(int k = 1; k <= 200000; k++) {
        mpz_urandomm(n, state, m);
        mpfr_set_z(w1, n, GMP_RNDN);

        mpfr_log(w1, w1, GMP_RNDN);
        mpfr_mul(w1, w1, t, GMP_RNDN);
        mpfr_fmod(w1, w1, twopi, GMP_RNDN);
        z1 = exp(I * mpfr_get_d(w1, GMP_RNDN));

        z2 = exp_itlogn(n);

        Double error = abs(z1 - z2);

        if(error > 1e-14) {
            cout << endl;
            cout << "Large error found with n = " << n << ". Error was " << abs(z1 - z2);
            failures++;
        }
        if(failures == 11) {
            break;
        }
    }
    if(failures) {
        cout << "done." << endl;
    }
    else {
        cout << "OK." << endl;
    }

    return failures;
}

int time_exp_itlogn() {
    mpfr_t t;
    mpfr_init2(t, 158);
    mpfr_set_str(t, "1e30", 10, GMP_RNDN);

    Complex z2;
    
    mpz_t n;
    mpz_init(n);

    create_exp_itlogn_table(t);

    mpz_set_str(n, "100000000000", 10);

//    cout << exp_itlogn(n) << endl;
    
    z2 = 0;

    const int number_of_iterations = 2000000;
    cout << "Timing exp_itlogn() over " << number_of_iterations << " iterations ...";
    cout.flush();
    
    clock_t start_time = clock();

    for(int k = 0; k < number_of_iterations; k++) {
        mpz_add_ui(n, n, 1u);
        z2 += exp_itlogn(n);
    }
    cout << z2 << "...";

    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << elapsed_time << " seconds." << endl;

    return elapsed_time;
}

int time_exp_itlogn_mpfr() {
    mpfr_t t;
    mpfr_init2(t, 158);
    mpfr_set_str(t, "1e30", 10, GMP_RNDN);

    Complex z1;
    
    mpz_t n;
    mpz_init(n);

    mpfr_t twopi;
    mpfr_init2(twopi, mpfr_get_prec(t));
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_t w1;
    mpfr_init2(w1, mpfr_get_prec(t));
    
    create_exp_itlogn_table(t);

    mpz_set_str(n, "100000000000", 10);

//    cout << exp_itlogn(n) << endl;
    
    z1 = 0;

    const int number_of_iterations = 2000000;
    cout << "Timing exp_itlogn using mpfr over " << number_of_iterations << " iterations ...";
    cout.flush();

    clock_t start_time = clock();
    for(int k = 0; k <= number_of_iterations; k++) {
        mpz_add_ui(n, n, 1u);
        mpfr_set_z(w1, n, GMP_RNDN);

        mpfr_log(w1, w1, GMP_RNDN);
        mpfr_mul(w1, w1, t, GMP_RNDN);
        mpfr_fmod(w1, w1, twopi, GMP_RNDN);
        z1 += exp(I * mpfr_get_d(w1, GMP_RNDN));
    }

    cout << z1 << "...";

    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << elapsed_time << " seconds." << endl;

    return elapsed_time;

}

int test_zeta_sum_stage1(gmp_randstate_t rand_state) {
    mpfr_t t, big_number;
    mpfr_init2(t, 200);
    mpfr_init2(big_number, 200);
    mpfr_set_str(big_number, "1e30", 10, GMP_RNDN);

    mpz_t n, one;
    mpz_init(n);
    mpz_init(one);
    mpz_set_ui(one, 1u);

    unsigned int K = 101331;
    mpz_set_ui(n, K);

    Complex S1 = 0.0;
    Complex S2 = 0.0;


    for(int k = 0; k < 10; k++) {
        mpfr_urandomb(t, rand_state);
        mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
        mpfr_add(t, t, big_number, GMP_RNDN);

        S1 = zeta_block_mpfr(one, K - 1, t);
        zeta_sum_stage1(n, t, 1, 1, &S2);
        //S2 = zeta_block_stage1(one, K, t);

        cout << S1 << " " << S1 - S2 << endl;
    }

    mpz_clear(n);
    mpz_clear(one);
    mpfr_clear(t);
    mpfr_clear(big_number);

    return 0;
}

int test_zeta_sum_stage2(gmp_randstate_t rand_state) {
    mpfr_t t, big_number;  
    mpfr_init2(t, 200);
    mpfr_init2(big_number, 200);
    mpfr_set_str(big_number, "1e30", 10, GMP_RNDN);

    mpz_t v;
    mpz_init(v);

    int length = 100000;
    
    mpz_t mp_length;
    mpz_init(mp_length);
    mpz_set_si(mp_length, length);

    Complex S1;
    Complex S2;

    cout << "Testing stage 2 sum on block sizes of 100000 ten times just past stage1_bound." << endl;

    for(int k = 0; k < 10; k++) {
        mpfr_urandomb(t, rand_state);
        mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
        mpfr_add(t, t, big_number, GMP_RNDN);

        stage_1_bound(v, t);

        create_exp_itlogn_table(t);

        S1 = zeta_block_mpfr(v, length, t);
        zeta_sum_stage2(v, mp_length, t, 1, 1, &S2);

        cout << S1 << "   " << S2 << "   " << S1 - S2 << endl;
    }

    length = 102013;
    mpz_set_si(mp_length, length);

    cout << "Now testing stage 2 sum on block sizes of " << length << " ten times with v = 2 * stage1_bound." << endl;

    for(int k = 0; k < 10; k++) {
        mpfr_urandomb(t, rand_state);
        mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
        mpfr_add(t, t, big_number, GMP_RNDN);

        stage_1_bound(v, t);
        mpz_mul_ui(v, v, 2u);

        create_exp_itlogn_table(t);

        S1 = zeta_block_mpfr(v, length, t);
        zeta_sum_stage2(v, mp_length, t, 1, 1, &S2);

        cout << S1 << "   " << S2 << "   " << S1 - S2 << endl;
    }




    length = 2020317;
    mpz_set_si(mp_length, length);

    cout << "Now testing once on a block size of length " << length << " with v = 3 * stage1_bound. This may take a few minutes." << endl;
    mpfr_urandomb(t, rand_state);
    mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
    mpfr_add(t, t, big_number, GMP_RNDN);

    stage_1_bound(v, t);
    mpz_mul_ui(v, v, 3u);

    create_exp_itlogn_table(t);

    S1 = zeta_block_mpfr(v, length, t);
    zeta_sum_stage2(v, mp_length, t, 1, 1, &S2);

    cout << S1 << "   " << S2 << "   " << S1 - S2 << endl;

    mpz_clear(v);
    mpz_clear(mp_length);

    return 0;

}

int test_zeta_sum_stage3(gmp_randstate_t rand_state) {
    mpfr_t t, big_number;  
    mpfr_init2(t, 200);
    mpfr_init2(big_number, 200);
    mpfr_set_str(big_number, "1e30", 10, GMP_RNDN);

    mpz_t v;
    mpz_init(v);

    int length = 100000;
    
    mpz_t mp_length;
    mpz_init(mp_length);
    mpz_set_si(mp_length, length);

    Complex S1;
    Complex S2;

    cout << "Testing stage 3 sum on block sizes of 100000 ten times just past stage2_bound." << endl;

    for(int k = 0; k < 10; k++) {
        mpfr_urandomb(t, rand_state);
        mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
        mpfr_add(t, t, big_number, GMP_RNDN);

        stage_2_bound(v, t);

        create_exp_itlogn_table(t);

        S1 = zeta_block_mpfr(v, length, t);
        zeta_sum_stage3(v, mp_length, t, 1, 1, &S2);

        cout << S1 << "   " << S2 << "   " << S1 - S2 << endl;
    }

    length = 102013;
    mpz_set_si(mp_length, length);

    cout << "Now testing stage 3 sum on block sizes of " << length << " ten times with v = 2 * stage2_bound." << endl;

    for(int k = 0; k < 10; k++) {
        mpfr_urandomb(t, rand_state);
        mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
        mpfr_add(t, t, big_number, GMP_RNDN);

        stage_2_bound(v, t);
        mpz_mul_ui(v, v, 2u);

        create_exp_itlogn_table(t);

        S1 = zeta_block_mpfr(v, length, t);
        zeta_sum_stage3(v, mp_length, t, 1, 1, &S2);

        cout << S1 << "   " << S2 << "   " << S1 - S2 << endl;
    }




    length = 2020317;
    mpz_set_si(mp_length, length);

    cout << "Now testing once on a block size of length " << length << " with v = (t/2pi)^1/2. This may take a few minutes." << endl;
    mpfr_urandomb(t, rand_state);
    mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
    mpfr_add(t, t, big_number, GMP_RNDN);

    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));
    mpfr_div_d(x, t, 2 * PI, GMP_RNDN);
    mpfr_sqrt(x, x, GMP_RNDN);
    mpfr_get_z(v, x, GMP_RNDN);

    create_exp_itlogn_table(t);

    S1 = zeta_block_mpfr(v, length, t);
    zeta_sum_stage3(v, mp_length, t, 1, 1, &S2);

    cout << S1 << "   " << S2 << "   " << S1 - S2 << endl;

    mpz_clear(v);
    mpz_clear(mp_length);

    return 0;

}


int test_zeta_sum() {
    mpfr_t t;
    mpfr_init2(t, 155);
    mpfr_set_str(t, "1e15", 10, GMP_RNDN);
    //mpfr_set_str(t, "2.38137487412446e12", 10, GMP_RNDN);
    //mpfr_set_str(t, "144176897509546973538.2912188",  10, GMP_RNDN); // the 10^21st zero

    Complex rotation_factor;
    Complex S = hardy_Z(t, rotation_factor);

    cout << endl << S << endl;
    cout << S * rotation_factor << endl;

    return 0;
}



double time_zeta_sum_stage1() {
    mpfr_t t;
    mpfr_init2(t, 158);
    mpfr_set_str(t, "1e30", 10, GMP_RNDN);

    mpz_t n;
    mpz_init(n);
    mpz_set_ui(n, 1000000u);
    
    cout << "Timing zeta_sum_stage_1() on a sum of length " << n << " with t = 1e30....";

    clock_t start_time = clock();

    Complex S1;
    zeta_sum_stage1(n, t, 1, 1, &S1);

    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << elapsed_time << " seconds." << endl;
    
    mpz_clear(n);
    mpfr_clear(t);

    return elapsed_time;
}


int time_zeta_sum_stage3(gmp_randstate_t rand_state) {
    mpfr_t t, big_number;  
    mpfr_init2(t, 200);
    mpfr_init2(big_number, 200);
    mpfr_set_str(big_number, "1e30", 10, GMP_RNDN);

    mpz_t v, v_increment;
    mpz_init(v);
    mpz_init(v_increment);

    int length = 20000000;
    
    mpz_t mp_length;
    mpz_init(mp_length);
    mpz_set_si(mp_length, length);

    Complex S1;
    Complex S2;

    cout << "Timing stage 3 sum ten times on large block sizes of length " << length << " starting at v = 20 stage2_bound for random large t...";
    cout.flush();
    int Kmin_start = 500;
    int Kmin_end = 1000;
    int Kmin_increment = 100;
    Complex Z[30];

    int number_of_gridpoints = 200;

    for(int Kmin = Kmin_start; Kmin <= Kmin_end; Kmin+=Kmin_increment) {
        clock_t start_time = clock();

        for(int k = 0; k < 10; k++) {
            cout << k << " ";
            cout.flush();
            mpfr_urandomb(t, rand_state);
            mpfr_mul_ui(t, t, 1000000000, GMP_RNDN);
            mpfr_add(t, t, big_number, GMP_RNDN);

            compute_taylor_coefficients(t, Z);

            stage_2_bound(v, t);
            mpz_cdiv_q_ui(v_increment, v, number_of_gridpoints);

//            cout << v << endl;

            mpz_mul_ui(v, v, 20u);
            create_exp_itlogn_table(t);

            for(int n = 0; n < number_of_gridpoints; n++) {
                Complex z;
                S2 += zeta_block_stage3(v, 200000, t, Z, 1, 1, &z, Kmin);
                //S2 += zeta_block_stage2(v, 200000, t);
                mpz_add(v, v, v_increment);
            }

            //S2 = zeta_sum_stage3(v, mp_length, t, Kmin);
        }

        clock_t end_time = clock();
        double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
        cout << "With Kmin = " << Kmin << ": Done in " << elapsed_time << " seconds." << endl;
    }

    mpz_clear(v);
    mpz_clear(v_increment);
    mpz_clear(mp_length);
    mpfr_clear(t);
    mpfr_clear(big_number);

    print_stats();

    return 0;

}

void H_test() {

    cout << H_method1(58, 72.0 * I) << endl;
    cout << H_method1(58, 80.0 * I) << endl;
    cout << H_method1(65, 80.0 * I) << endl;
    cout << H_method1(70, 80.0 * I) << endl;
    cout << H_method1(80, 80.0 * I) << endl;
    cout << H_method1(85, 80.0 * I) << endl;
}


int main() {
    unsigned int seed = time(NULL);
    //seed = 1276487827;
    //seed = 1276414014;
    //seed = 1277852897;
    //seed = 1277923991;

    //seed = 1278127602;
    //seed = 1278182770; // this seed makes the first test in test_theta_algorithm(20, 5000)
    //seed = 1278268854;
    //seed = 1278723182;
    //seed = 1279180855;
    cout << "Seeding rand() and gmp with " << seed << "." << endl;
    srand(seed);
    
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, seed);

    cout << setprecision(17);
    test_zeta_sum_stage1(rand_state);

    return 0;

    //build_F0_cache(11, 6, 25, 11000, exp(-20));
    //build_F1_cache(181, 51, 25, exp(-30));
    //build_F1_cache(181, 51, 25, exp(-30));
    //build_F1_cache(361, 101, 25, exp(-30));
    
    //build_F0_cache(11, 6, 25, 2800, exp(-25));
    //build_F0_cache(11, 6, 25, 10500, exp(-25));
    //build_F1_cache(1801, 501, 30, exp(-30));
    //build_F2_cache(10500, 11, 11, 6, exp(-30));
    //build_IC7_cache(600, 200, 25, exp(-30));

    //build_F0_cache(11, 6, 25, 3100, exp(-25));
    //build_F1_cache(181, 51, 30, exp(-30));
    //build_F2_cache(3100, 11, 11, 6, exp(-30));
    //build_IC7_cache(600, 200, 25, exp(-30));

    //build_F0_cache(5, 4, 25, 50, exp(-30));
    //build_F0_cache(10, 20, 25, 20000, exp(-30));
    //build_F0_cache(10, 100, 25, 30000, exp(-30));
    //build_F1_cache(30, 2000, 30, exp(-30));
    //build_F2_cache(30000, 10, 10, 100, exp(-30));
    //build_IC7_cache(600, 200, 25, exp(-30));

    //time_zeta_sum_stage3(rand_state);

    //test3();

    //test_specific_inputs(-3.79326471769869, 0.000369601256708863, 1726, 16, exp(-30));
    //test_specific_inputs(-8.01233758126029, 0.00151760844655872, 339, 15, exp(-36));
    int K = 30000;
    Double b = (36 * 36 + 1000)/ (30000.0 * 30000.0);
    Double a = -2 * b * K;

    
    //K = 61;
    //b = 1/(8 * K);
    //a = .49;

    a = .80263398159418;
    b = 6.01723246090917;
    K = 5813;

    a = 0.2115681902677;
    a = 0.7515681902677;
    a = 0.22;
    b = 0.00750425709434452;
    b = 0.0075;
    K = 201;
    K = 186;
    int j = 19;
//compute_exponential_sums() called with a = 0.462984559807454, b = 8.62172660502867e-06, K = 10300, j = 13, epsilon = 8.88178419700125e-16

    a = 0.462984559807454;
    b = 8.62172660502867e-06;
    K = 10300;
    j = 13;


    a = 0.896067520089479;
    b = 0.00054575118798739;
    j = 15;
    K = 489;

    a = 0.81720095212441;
    b = 1.61006127354104e-05; j = 9; K = 10094;

    a = 6.4375082293699997; b = -6.9832637240100954; j = 16; K = 8383;

    a = 2.8761590052750705; b = 9.0024296282801899; j = 19; K = 234; // in this example the coefficients of the subsum get really big, and we
                                                                     // lose some precision if the coefficients of the original sum
                                                                     // were all the same size.
    
    b = .00245;

    a = -7.7368794371080023, b = -3.1673489944857307, j = 2, K = 474;
    
    a = 0.052528253780923899, b = 2.9356077194412081e-06, j = 0, K = 100235;

    //test_specific_inputs(a, b, K, j, pow(2, -50));
    //test_specific_inputs(a, b, K, j, pow(2, -60));
    //test_specific_inputs(a, b, K, j, pow(2, -70));
    //test_specific_inputs(a, b, K, j, pow(2, -80));

    //print_stats();
    //exit(0);

    //test_specific_inputs(.0000001, .0000000001, 10000, 0, pow(2, -30));


    //test_theta_algorithm_EM_case(20, 10000, pow(2, -30));
    //test_theta_algorithm_EM_case(100, 1000, pow(2, -50));
    //test_theta_algorithm_EM_case(100, 1000, pow(2, -30));
    //test_theta_algorithm_EM_case(20, 100000, pow(2, -50));
    //test_theta_algorithm_EM_case(100, 517, pow(2, -30));

    //cout << H_method1(50, 1000.0 * I) << endl;
    //cout << H(50, 1000.0 * I, exp(-30)) << endl;

    
    //time_theta_algorithm_varying_Kmin(15, 40000, 30000, pow(2, -22));
    //time_theta_algorithm_zeta_case_varying_Kmin(15, 50000, 30000, pow(2, -20));
    //time_theta_algorithm_zeta_case_fixed_Kmin(15, 50000, 40000, 900, pow(2, -20));


    //test_specific_inputs(0.809496893458765, 9.61035310379055e-06, 10311, 0, pow(2, -50));

    //test_theta_algorithm(10, 5000);
    //test_theta_algorithm(10, 3000);
    //test_theta_algorithm(10, 1000);
    //test_theta_algorithm(20, 1023);
    //test_theta_algorithm(5, 10231);
    //test_theta_algorithm(5, 20231);
    //test_theta_algorithm(50, 5013);
    

    //test_theta_algorithm(1, 1000000, pow(2, -30));
    //test_theta_algorithm(1, 1000000, pow(2, -50));

    //test_theta_algorithm(100, 517, pow(2, -30));
    //test_theta_algorithm(100, 517, pow(2, -40));;
    //test_theta_algorithm(50, 5432, pow(2, -30));
    //test_theta_algorithm(50, 5432, pow(2, -40));
    //test_theta_algorithm(50, 1432, pow(2, -30));
    //test_theta_algorithm(50, 1432, pow(2, -40));
    //test_theta_algorithm(500, 123, pow(2, -30));
    //test_theta_algorithm(500, 123, pow(2, -40));
    //test_theta_algorithm(10, 10000, pow(2, -40));
    
    /*
    test_theta_algorithm(10, 10000, pow(2, -30));
    test_theta_algorithm(10, 10000, pow(2, -50));
    test_theta_algorithm(10, 8000, pow(2, -30));
    test_theta_algorithm(10, 8000, pow(2, -50));
    test_theta_algorithm(15, 20000, pow(2, -30));
    test_theta_algorithm(15, 20000, pow(2, -50));
    test_theta_algorithm(5, 40000, pow(2, -30));
    test_theta_algorithm(5, 40000, pow(2, -50));
    test_theta_algorithm(2, 100000, pow(2, -30));
    test_theta_algorithm(2, 100000, pow(2, -50));
    test_theta_algorithm(100, 1000, pow(2, -30));
    test_theta_algorithm(100, 1000, pow(2, -50));
 
    test_theta_algorithm_EM_case(10, 10000, pow(2, -30));
    test_theta_algorithm_EM_case(10, 10000, pow(2, -50));
    
    test_theta_algorithm_EM_case(3, 100000, pow(2, -30));
    test_theta_algorithm_EM_case(3, 100000, pow(2, -50));
    
    test_theta_algorithm_EM_case(20, 1000, pow(2, -30));
    test_theta_algorithm_EM_case(20, 1000, pow(2, -50));
    test_theta_algorithm_EM_case(100, 517, pow(2, -30));
    test_theta_algorithm_EM_case(100, 517, pow(2, -50));;
    test_theta_algorithm_EM_case(20, 5432, pow(2, -30));
    test_theta_algorithm_EM_case(20, 5432, pow(2, -50));
    test_theta_algorithm_EM_case(50, 1432, pow(2, -30));
    test_theta_algorithm_EM_case(50, 1432, pow(2, -50));
    test_theta_algorithm_EM_case(500, 123, pow(2, -30));
    test_theta_algorithm_EM_case(500, 123, pow(2, -50));
    */



    //time_theta_algorithm_varying_Kmin(10, 20010, 10000, exp(-14));
    

    //int new_seed = 1278263575;
    //srand(new_seed);
    //time_theta_algorithm_EM_case(10, 20010, 1000, exp(-14));

    //time_theta_algorithm(15, 10000);
    //test2();

    print_stats();

    //free_F1_cache();

    //test_fastlog2();
    //test_fastlog();
    //time_theta_algorithm(18, 10010);
    //test_exp_itlogn(rand_state);
    //time_exp_itlogn();
    //time_exp_itlogn_mpfr();
    //time_zeta_sum_stage1();
    //test_zeta_sum_stage2(rand_state);
    test_zeta_sum_stage3(rand_state);
    //test_zeta_sum();
    //time_theta_algorithm_varying_Kmin(10, 10010, 10000);
    //time_theta_algorithm_varying_Kmin(18, 10010, 10000);
    //time_theta_algorithm_varying_Kmin(18, 2010, 10000);
    //time_zeta_sum_stage3(rand_state);

    cout << "Used random seed " << seed << endl;
    gmp_randclear(rand_state);
    return 0;
}
