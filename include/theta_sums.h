#include <iostream>
#include <complex>
#include <string>
#include <cmath>
#include <cstdlib>

#include "mpfr.h"
#include "gmp.h"


#define BUILTIN_EXPECT(a, b) __builtin_expect( (a), (b) )

#include "misc.h"

// For timing purposes, we can set FAKE_PRECOMPUTATION = true
// in this case, no precomputation will be done, and the answers
// will be nonsense, but we can quickly get an idea of the
// performance of the code
const bool FAKE_PRECOMPUTATION = false;
//const bool FAKE_PRECOMPUTATION = true;

const bool FAKE_J_INTEGRALS = false;
const bool FAKE_IC7 = false;
const bool FAKE_EULER_MACLAURIN = false;

const int Kmin = 800;
const int mpfr_Kmin = 2000;
//const int mpfr_Kmin = 2000;


inline Complex I_power(int n) {
    Complex S = 0;
    switch(n % 4) {
        case 0:
            S = 1;
            break;
        case 1:
            S = I;
            break;
        case 2:
            S = -1;
            break;
        case 3:
            S = -I;
            break;
    }
    return S;
}

inline Complex minus_I_power(int n) {
    Complex S = 0.0;
    switch(n % 4) {
        case 0:
            S = Complex(1.0, 0.0);
            break;
        case 1:
            S = Complex(0.0, -1.0);
            break;
        case 2:
            S = Complex(-1.0, 0.0);
            break;
        case 3:
            S = Complex(0.0, 1.0);
            break;
    }
    return S;
}

inline int minus_one_power(int n) {
    int S = 0;
    switch(n % 2) {
        case 0:
            S = 1;
            break;
        case 1:
            S = -1;
            break;
    }
    return S;
}

inline Complex exp_i_pi4(int n) {
    //
    // Return exp( n * I * PI / 4)
    //

    Complex S = 0;
    switch(n % 8) {
        case 0:
            S = Complex(1, 0);
            break;
        case 1:
            S = Complex(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 2:
            S = Complex(0, 1);
            break;
        case 3:
            S = Complex(-sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 4:
            S = Complex(-1, 0);
            break;
        case 5:
            S = Complex(-sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 6:
            S = Complex(0, -1);
            break;
        case 7:
            S = Complex(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
    }
    return S;
}


inline Complex exp_minus_i_pi4(int n) {
    //
    // Return exp(-n * I * PI /4)
    //

    Complex S = 0;
    switch(n % 8) {
        case 0:
            S = Complex(1, 0);
            break;
        case 1:
            S = Complex(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 2:
            S = Complex(0, -1);
            break;
        case 3:
            S = Complex(-sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 4:
            S = Complex(-1, 0);
            break;
        case 5:
            S = Complex(-sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 6:
            S = Complex(0, 1);
            break;
        case 7:
            S = Complex(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
    }
    return S;
}


typedef struct{
    Double a;
    Double b;
    int K;
    int j;
    int q;
    Complex ExpAK;
    Complex ExpBK;
    Complex ExpAK_inverse;
    Complex ExpBK_inverse;
    Complex ExpAB;
    Complex ExpABK;
    Complex C1;
    Complex C5;
    Complex C7;
    Complex C8;
} theta_cache;



inline Double K_power(int l, const theta_cache * cache) {
    return * ( (Double *)((intptr_t)(cache) + sizeof(theta_cache) + (cache->j + l) * sizeof(Double) ));
}

inline Double root_2pi_b_power(int l, const theta_cache * cache) {
    return * ( (Double *)((intptr_t)(cache) + sizeof(theta_cache) + (3 * cache->j + 2 + l) * sizeof(Double) ));
}

theta_cache * build_theta_cache(mpfr_t mp_a, mpfr_t mp_b, int j, int K);
void free_theta_cache(theta_cache * cache);

namespace verbose {
    const int IC0 = 0;
    const int IC1 = 0;
    const int IC1c = 0;
    const int IC6 = 0;
    const int IC7 = 0;
    const int IC7star = 0;
    const int IC9E = 0;
    const int compute_exponential_sum = 0;
    const int S1 = 0;
    const int S2 = 0;
    const int G = 0;
    const int H = 0;

    const int direct_evaluation = 0;

    const int J_Integral_0 = 0;
    const int J_Integral_1 = 0;
    const int J_Integral_2 = 0;

    const int JBulk = 0;
}

namespace stats {
    const bool stats = false;
    extern int H_method1;
    extern int H_method2;
    extern int H_method3;
    extern int H_method4;

    extern int H_function_big;
    extern int H_function_small;

    extern int G_method1;
    extern int G_method2;
    
    extern int G_itwopi_method1;
    extern int G_itwopi_method2;

    extern int G_R_method1;
    extern int G_R_method2;

    extern int exp;

    extern int exponential_sum_called;
    extern int exponential_sum_euler_maclaurin;
    extern int exponential_sum_taylor_expansion;

    extern int H_Integral_0;
    extern int H_Integral_2;
    extern int J_Integral_0;
    extern int J_Integral_1;
    extern int J_Integral_2;

    extern int J_Integral_0_taylor_expansion;
    extern int J_Integral_1_taylor_expansion;
    extern int J_Integral_2_taylor_expansion;

    extern int J_Integral_0_zero;
    extern int J_Integral_1_zero;
    extern int J_Integral_2_zero;

    extern int J_Integral_0_terms_used;
    extern int J_Integral_1_terms_used;
    extern int J_Integral_2_terms_used;

    extern int IC7;
    extern int IC7zero;
    
    extern int IC7_taylor_expansion;
    extern int IC7_terms_used;

    extern int IC0;
    extern int IC0_method1;
    extern int IC0_method2;
    extern int IC0_method3;
    extern int IC0_method4;
}

void print_stats();

Complex ExpA(mpfr_t A, int K);
Complex ExpAK(mpfr_t A, int K);
Complex ExpB(mpfr_t B, int K);
Complex ExpBK(mpfr_t B, int K);
Complex ExpAB(mpfr_t A, mpfr_t B);
Complex ExpABK(mpfr_t A, mpfr_t B, int K);



inline Complex EXP(Complex z) {                                                 //--------------------------------------------
    stats::exp++;                                                               // This is just here for counting purposes.
    return std::exp(z);                                                         // There is pretty much no overhead in using it,
}                                                                               // and commenting out the first line of this
                                                                                // function will get rid of any overhead at all.
                                                                                // ---------------------------------------------

inline bool check_condition(bool condition, char * message) {              //----------------------------------------------
    if(!condition) {                                                            //
        std::cout << "WARNING: " << message << std::endl;                       //
    }                                                                           //
    return condition;                                                           //
}                                                                               //----------------------------------------------


                                                                                //----------------------------------------------
                                                                                //  Functions to compute the integral
Complex H(int j, Complex alpha, Double epsilon);                                //    / 1 
Complex H_method1(int j, Complex alpha);                                        //    |   j
//inline Complex H_method2(int j, Complex alpha, Double epsilon);               //    |  t  exp(- 2 pi alpha) dt
Complex H_method3(int j, Complex alpha, Double epsilon);                        //    |
Complex H_method4(int j, Complex alpha, Double epsilon);                        //    / 0
                                                                                //  (Defined in H_functions.cc)
                                                                                //----------------------------------------------


                                                                                //----------------------------------------------
                                                                                //  Functions to compute the integral         
                                                                                //    / 1
                                                                                //    |
Complex G(Complex alpha, Complex b, int n, int j, Double epsilon, int method = 0);  //
Complex G_R(Complex alpha, Double b, int n, int j, Double epsilon, int method = 0);  //
Complex G_I(Complex alpha, Double b, int n, int j, Double epsilon, int method = 0);  //
Complex G_I_over_twopi(Complex alpha, int n, int j, Double epsilon, int method = 0);  //
Complex G_method1(Complex alpha, Complex b, int n, int j, Double epsilon);  //
Complex G_method1_R(Complex alpha, Double b, int n, int j, Double epsilon);  //
Complex G_method1_I(Complex alpha, Double b, int n, int j, Double epsilon);  //
Complex G_via_Euler_MacLaurin(Complex alpha, Complex b, int n, int j, Double epsilon); //    / 0
Complex G_via_Euler_MacLaurin_R(Complex alpha, Double b, int n, int j, Double epsilon); //    / 0
Complex G_via_Euler_MacLaurin_I(Complex alpha, Double b, int n, int j, Double epsilon); //    / 0
                                                                                //  For complex parameters alpha and b.
                                                                                //  (Defined in G_functions.cc)
                                                                                //----------------------------------------------

Complex H_Integral_0(int j, Double a, int M, Double epsilon);                   //----------------------------------------------
Complex J_Integral_0(Double a, Double b, int j, int M, int K, theta_cache * cache, Double epsilon, bool use_cache = true);  //
Complex J_Integral_1(Double a, Double b, int j, int M, int K, theta_cache * cache, Double epsilon, bool use_cache = true);  //
Complex H_Integral_2(int j, Double a1, Double a2, Double epsilon);              //
Complex J_Integral_2(Double a1, Double a2, Double b, theta_cache * cache, Double epsilon, bool use_cache = true);// Various integrals.

void build_F0_cache(long number_of_a, long number_of_b, long max_j, long max_M, Double epsilon, string cache_directory);
void build_F1_cache(long a_per_unit_interval, long b_per_unit_interval, long max_j, Double epsilon, string cache_directory);
void build_F1_cache(string cache_directory);
void build_F2_cache(long max_a1, long number_of_a1, long number_of_a2, long number_of_b, Double epsilon, string cache_directory);
void free_F0_cache();
void free_F1_cache();
void free_F2_cache();
                
Complex JBulk2(Double a, Double b, int j, int M, int K, theta_cache * cache, Complex * Z, Double * epsilon);                

inline Complex JBulk(Double a, Double b, int j, int M, int K, theta_cache * cache, Double epsilon) {         //                         
    Double x = epsilon/2;
    Complex A = J_Integral_0(a, b, j, M, K, cache, x);
    Complex B = J_Integral_1(a, b, j, M, K, cache, x);
    if(verbose::JBulk) {
        cout << "JBulk returning " << A << " + " << B << " = " << A + B << endl;
    }
    return A + B;
    //return J_Integral_0(a, b, j, M, K, cache, epsilon/2)                                          // See H_and_J_integrals.cc
    //                                    + J_Integral_1(a, b, j, M, K, cache, epsilon/2);       //
}                                                                                       //
                                                                                        //
inline Complex JBoundary(Double a1, Double a2, Double b, int j, int K, theta_cache * cache, Double epsilon){ 
    if(j == 0) {
        Double x = epsilon * .33333333333333333333;
        return J_Integral_2(a1, a2, b, cache, x)                                    
                                        + J_Integral_1(a1, b, j, -1, K, cache, x)     
                                        - J_Integral_1(a2, b, j, -1, K, cache, x);
    }
    else {
        Double x = epsilon * .25;
        return J_Integral_0(a1, b, j, -1, K, cache, x) + J_Integral_1(a1, b, j, -1, K, cache, x)
                + (Double)minus_one_power(j + 1) * (  J_Integral_0(a2, b, j, -1, K, cache, x) + J_Integral_1(a2, b, j, -1, K, cache, x) );
    }
}                                                                                     





Complex IC0(int j, mpfr_t mp_a, mpfr_t mp_b, const theta_cache * cache, Double epsilon);
Complex IC1(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon);//----------------------------------------------
Complex IC1c(int K, int j, Double a, Double b, Complex C8, const theta_cache * cache, Double epsilon);     //
inline Complex IC3(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon);           //
inline Complex IC3c(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon);          //
Complex IC4(int K, int j, Double a, Double b, Complex C11, const theta_cache * cache, Double epsilon);     //
Complex IC4c(int K, int j, Double a, Double b, Complex C11, const theta_cache * cache, Double epsilon);    //
Complex IC5(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon);                  //
Complex IC6(int K, int j, Double a, Double b, mpfr_t mp_a, const theta_cache * cache, Double epsilon);     //      
Complex IC7(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon);                  //
Complex IC7_method1(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon, int L);                  //
Complex IC7star(Double a, int j, Double epsilon, bool use_cache = true);                  //
Complex IC8(int K, int j, mpfr_t mp_a, mpfr_t mp_b, const theta_cache * cache);                            //
Complex IC9E(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon);                 //
inline Complex IC9H(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon);          //
                                                                                //  (Defined in ICn.cc, unless defined inline below)
                                                                                //
                                                                                //

void build_IC7_cache(int a_per_unit_interval, Double max_a, int max_j, Double epsilon, string cache_directory); //max_a should be passed as about sqrt(K) to compute sums of length K

inline Complex IC3(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon) {
    // needs a <= 0
    //       b >= 0
    //return pow(-I, j+1) * IC9H(K, j, -a, b, cache, epsilon);
    return minus_I_power(j+1) * IC9H(K, j, -a, b, cache, epsilon);
}

inline Complex IC3c(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon) {
    // needs a >= 0
    //       b >= 0
    //return pow(I, j+1) * IC9H(K, j, a, b, cache, epsilon);
    return I_power(j+1) * IC9H(K, j, a, b, cache, epsilon);
}

inline Complex IC9H(int K, int j, Double a, Double b, const theta_cache * cache, Double epsilon) {
    //
    // Compute the integral (1/K^j)int_0^\infty t^j exp(-2 pi a t - 2 pi i b t^2)
    //
    // after a change of contour, this is well approximated by IC7(K, a, b) with large K
    //

    //int endpoint = to_int(10 * ceil( max(-LOG(epsilon * sqrt(b)), 1.0)/sqrt(b) ));
    //Double z = pow((Double)endpoint/(Double)K, j);
    //Double z = pow(K, j);
    Complex S = IC7(-1, j, a, b, cache, epsilon * K_power(j, cache));
    S = S * K_power(-j, cache);
    return S;
}


                                                                                //-----------------------------------------------------------
                                                                                //
                                                                                // Derivatives of the function
                                                                                //
                                                                                // g(t) = exp(2 pi i alpha t + 2 pi i b t^2)
                                                                                //
//void initialize_power_arrays(int n, Complex alpha, Complex b);                // at 0, 1, and K. Only supports derivatives up
//Complex g_derivative_at_1(int n);                                             // to the 21st, and initialize_power_arrays() must
//Complex g_derivative_at_0(int n);                                             // be called before the other functions, with n
//Complex g_derivative_at_K_without_exponential_factor(int n, int K);           // at least as large as any derivatives to be
                                                                                // computed. For the function ...without_exponential_factor()
                                                                                // the answer needs to be multiplied by
                                                                                //
                                                                                //  exp(2 pi i a K + 2 pi i b K^2)
                                                                                //
                                                                                // to get the actual derivative.
                                                                                //
                                                                                // Aside from the stats:: stuff, this might be the only
                                                                                // place that this code is not thread-safe
                                                                                // ----------------------------------------------------------
//Complex g_derivative_at_1(int n, Complex a, Complex b);
void g_derivative_polynomial(int n, Complex * p, Complex * q, const Complex a, const Complex b);
void g_derivative_polynomial_R(int n, Complex * p, Complex * q, const Complex a, const Double b);
void g_derivative_polynomial_I(int n, Complex * p, Complex * q, const Complex a, const Double b);
void g_derivative_polynomial_I_over_twopi(int n, Complex * p, Complex * q, const Complex a);


//------------------------------------------------------------------------------------------------------------
//
// All of the functions below are currently defined in theta_sums.cc


Complex S1(int K, mpfr_t mp_a, mpfr_t mp_b, Double epsilon);
Complex S2(int K, mpfr_t mp_a, mpfr_t mp_b, Double epsilon);

int normalize(Double &a, Double &b);
int normalize(mpfr_t a, mpfr_t b);


Complex direct_exponential_sum_evaluation2(mpfr_t a, mpfr_t b, int j, int m, int M);


inline Complex compute_C11(mpfr_t a, mpfr_t b, int K) {
    //
    // Compute C11 = I exp(2 pi i K a + 2 pi i b K^2)
    //
    mpfr_t tmp1;
    mpfr_t tmp2;
    mpfr_init2(tmp1, mpfr_get_prec(a));
    mpfr_init2(tmp2, mpfr_get_prec(a));
    
    mpfr_const_pi(tmp1, GMP_RNDN);
    mpfr_mul_si(tmp1, tmp1, 2, GMP_RNDN);
    mpfr_mul(tmp1, tmp1, a, GMP_RNDN);
    mpfr_mul_si(tmp1, tmp1, K, GMP_RNDN);  // now tmp1 = 2 pi a K

    mpfr_const_pi(tmp2, GMP_RNDN);
    mpfr_mul_si(tmp2, tmp2, 2, GMP_RNDN);
    mpfr_mul(tmp2, tmp2, b, GMP_RNDN);
    mpfr_mul_si(tmp2, tmp2, K, GMP_RNDN);
    mpfr_mul_si(tmp2, tmp2, K, GMP_RNDN);  // now tmp2 = 2 pi b K^2
    
    mpfr_add(tmp1, tmp1, tmp2, GMP_RNDN);

    mpfr_sin_cos(tmp1, tmp2, tmp1, GMP_RNDN);
    
    Complex S(mpfr_get_d(tmp2, GMP_RNDN), mpfr_get_d(tmp1, GMP_RNDN));
    
    mpfr_clear(tmp1);
    mpfr_clear(tmp2);
    return I * S;
}

inline Complex compute_C12(mpfr_t mp_a, mpfr_t mp_b, int K) {
    Complex z = ExpBK(mp_b, K);
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
    
    Complex w = I * exp(-2 * PI * (a + 2 * b * K) * K) / (sqrt(2 * PI * b));
 
    return w/z;
}

Complex direct_exponential_sum_evaluation2(Double a, Double b, int j, int m, int M, int working_precision = 53);
Complex compute_exponential_sums_using_theta_algorithm(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int _Kmin);
Complex compute_exponential_sums_using_theta_algorithm2(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int _Kmin);
Complex compute_exponential_sums_directly(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon);
Complex compute_exponential_sums_for_small_b(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon);
Complex compute_exponential_sums(mpfr_t mp_a, mpfr_t mp_b, int j, int K, Complex * v, Double epsilon, int _Kmin = 0, int method=0);
Complex compute_exponential_sums(Double a, Double b, int j, int K, Complex * v, Double epsilon, int _Kmin = 0, int method=0);



//Complex w_coefficient(Double * a_powers, Double * b_powers, Double * q_powers, int s, int j, theta_cache * cache, Complex * inner_sums);
void compute_subsum_coefficients(Complex * v2, Complex * v, const theta_cache * cache);