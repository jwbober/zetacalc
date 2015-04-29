#include "theta_sums.h"
#include "log.h"
#include "main_sum.h"

#include <ctime>
#include <iostream>
#include <iomanip>
#include "theta_sums/precomputed_tables.h"

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
    Complex * v = new Complex[j + 1];
    for(int l = 0; l <= j; l++) {
        v[l] = 1.0/(l * l * l * l + 1);
    }

    Complex S1 = compute_exponential_sums(a, b, j, K, v, epsilon, Kmin, 0);
    Complex S2 = compute_exponential_sums(a, b, j, K, v, epsilon, 0, 1);

    Complex error = S1 - S2;
    cout << "a = " << a << ", b = " << b << ", j = " << j << ", K = " << K << ": log2(error) = " << log2(abs(error)) << "; error = " << error << ", answer = " << S2 << endl;

    delete [] v;
}


double time_theta_algorithm(int j, int K, int number_of_tests = 1000) {
    Complex * v = new Complex[j + 1];
    for(int k = 0; k <= j; k++) {
        v[j] = 1.0/(k * k * k * k * k * k + 1);
        cout << v[j] << endl;
        //v[j] = 1.0;
    }

    clock_t start_time = clock();

    Double epsilon = exp(-20);

    int n = 0;
    Complex z1 = 0.0;

    cout << "Timing theta_algorithm with K = " << K << " and j = " << j << endl;
    cout << "   Running approximately 1000 iterations total." << endl;

    for(Double a = 0; a < .5; a += .5/(number_of_tests/1000.0) ) {
        for(Double b = 1.1/((Double)K); b <= 1.0; b += 1.0/1000.0) {
            n++;
       //     if(n % 100 == 0) {
       //         cout << "   Running iteration number " << n << " with a = " << a << " b = " << b << endl;
       //     }
            z1 += compute_exponential_sums(a, b, j, K, v, epsilon);
        }
    }
    cout << "Sum was " << z1 << endl;

    clock_t end_time = clock();

    double elapsed_time = (double)(end_time - start_time)/(double)CLOCKS_PER_SEC;
    cout << "Number of seconds for this run was " << elapsed_time << endl;

    delete [] v;

    return elapsed_time;
    
}

double time_theta_algorithm_varying_Kmin(int j, int approx_K, int number_of_tests, Double epsilon) {
    Complex * v = new Complex[j + 1];
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

    delete [] v;

    return 1.0;
    
}


double time_theta_algorithm_zeta_case_varying_Kmin(int j, int approx_K, int number_of_tests, Double epsilon) {
    Complex * v = new Complex[j + 1];
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

    delete [] v;

    return 1.0;
    
}

double time_theta_algorithm_zeta_case_fixed_Kmin(int j, int approx_K, int number_of_tests, int Kmin, Double epsilon) {
    Complex * v = new Complex[j + 1];
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

    delete [] v;

    return 1.0;
    
}






double time_theta_algorithm_EM_case(int j, int K, int number_of_tests, Double epsilon) {
    Complex * v = new Complex[j + 1];
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

    delete [] v;

    return elapsed_time;
    
}







int test_exp_itlogn(gmp_randstate_t state) {
    mpfr_t t;
    mpfr_init2(t, 158);
    mpfr_set_str(t, "1374481970215211977756", 10, GMP_RNDN);
    mpfr_set_str(t, "1633275126614341648036", 10, GMP_RNDN);
    mpfr_set_str(t, "133275126614341648036", 10, GMP_RNDN);
    mpfr_set_str(t, "233275126614341648036", 10, GMP_RNDN);
    mpfr_set_str(t, "433275126614341648036", 10, GMP_RNDN);

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
            cout << "On test number " << k << ": Large error found with n = " << n << ". Error was " << abs(z1 - z2);
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

void test1() {
    mpfr_t a;
    mpfr_init2(a, 200);
    mpfr_set_str(a, "0.46666552511643095", 10, GMP_RNDN);

    mpfr_t b;
    mpfr_init2(b, 200);
    mpfr_set_str(b, "0.18333333337945587", 10, GMP_RNDN);
    
    int j = 18;
    unsigned int block_size = 17431;
    Complex Z[19];
    Z[0] = Complex(1,0);
    Z[1] = Complex(-4.1202789467000038e-11,0);
    Z[2] = Complex(2.546504789792894e-21,0);
    Z[3] = Complex(0,0.24294252905079758);
    Z[4] = Complex(-0,-2.502477469265138e-11);
    Z[5] = Complex(0,2.2271555299401193e-21);
    Z[6] = Complex(-0.029510536210798802,0);
    Z[7] = Complex(4.8636656422072966e-12,0);
    Z[8] = Complex(-5.7864515793819967e-22,0);
    Z[9] = Complex(-0,-0.0023897881002322013);
    Z[10] = Complex(0,5.4156264780535148e-13);
    Z[11] = Complex(-0,-8.0938571022957095e-23);
    Z[12] = Complex(0.00014514529124147794,0);
    Z[13] = Complex(0,0);
    Z[14] = Complex(0,0);
    Z[15] = Complex(7.0523928268038457e-06,0);
    Z[16] = Complex(0,0);
    Z[17] = Complex(0,0);
    Z[18] = Complex(2.8555435820057152e-07,0);
    Z[18] = exp(-15);

    //Complex S = compute_exponential_sums(a, b, j, block_size - 1, Z, exp(-20), 800, 0);
    Complex S1 = compute_exponential_sums(a, b, j, block_size - 1, Z, exp(-20), 500, 0);
    Complex S2 = compute_exponential_sums(a, b, j, block_size - 1, Z, exp(-20), 500, 1);

    cout << S1 - S2 << endl;


    //cout << "S = " << S << endl;
    cout << "S1 = " << S1 << endl;

}

int main(int argc, char ** argv) {

    unsigned int seed = time(NULL);

    if(argc > 1) {
        seed = strtoul(argv[1], NULL, 10);
    }

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


    test_theta_algorithm(10, 5000);
    test_theta_algorithm(10, 3000);

    test_theta_algorithm(100, 517, pow(2, -30));
    test_theta_algorithm(100, 517, pow(2, -40));;
    test_theta_algorithm(50, 5432, pow(2, -30));
    test_theta_algorithm(50, 5432, pow(2, -40));
    test_theta_algorithm(50, 1432, pow(2, -30));
    test_theta_algorithm(50, 1432, pow(2, -40));
    test_theta_algorithm(500, 123, pow(2, -30));
    test_theta_algorithm(500, 123, pow(2, -40));
    test_theta_algorithm(5, 40000, pow(2, -30));
    test_theta_algorithm(5, 40000, pow(2, -50));
    test_theta_algorithm(2, 100000, pow(2, -30));
    test_theta_algorithm(2, 100000, pow(2, -50));
 
    test_theta_algorithm_EM_case(20, 10000, pow(2, -30));
    test_theta_algorithm_EM_case(100, 1000, pow(2, -50));
    test_theta_algorithm_EM_case(100, 1000, pow(2, -30));

    test_theta_algorithm_EM_case(10, 10000, pow(2, -30));
    test_theta_algorithm_EM_case(10, 10000, pow(2, -50));
    
    test_theta_algorithm_EM_case(3, 100000, pow(2, -30));
    test_theta_algorithm_EM_case(3, 100000, pow(2, -50));

    cout << "Used random seed " << seed << endl;
    gmp_randclear(rand_state);
    return 0;
}
