#include <iostream>
#include <complex>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdint>

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
const int max_j = 30;
//const int mpfr_Kmin = 2000;


inline std::complex<double> I_power(int n) {
    std::complex<double> S = 0;
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

inline std::complex<double> minus_I_power(int n) {
    std::complex<double> S = 0.0;
    switch(n % 4) {
        case 0:
            S = std::complex<double>(1.0, 0.0);
            break;
        case 1:
            S = std::complex<double>(0.0, -1.0);
            break;
        case 2:
            S = std::complex<double>(-1.0, 0.0);
            break;
        case 3:
            S = std::complex<double>(0.0, 1.0);
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

inline std::complex<double> exp_i_pi4(int n) {
    //
    // Return exp( n * I * PI / 4)
    //

    std::complex<double> S = 0;
    switch(n % 8) {
        case 0:
            S = std::complex<double>(1, 0);
            break;
        case 1:
            S = std::complex<double>(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 2:
            S = std::complex<double>(0, 1);
            break;
        case 3:
            S = std::complex<double>(-sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 4:
            S = std::complex<double>(-1, 0);
            break;
        case 5:
            S = std::complex<double>(-sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 6:
            S = std::complex<double>(0, -1);
            break;
        case 7:
            S = std::complex<double>(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
    }
    return S;
}


inline std::complex<double> exp_minus_i_pi4(int n) {
    //
    // Return exp(-n * I * PI /4)
    //

    std::complex<double> S = 0;
    switch(n % 8) {
        case 0:
            S = std::complex<double>(1, 0);
            break;
        case 1:
            S = std::complex<double>(sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 2:
            S = std::complex<double>(0, -1);
            break;
        case 3:
            S = std::complex<double>(-sqrt(2.0)/2.0, -sqrt(2.0)/2.0);
            break;
        case 4:
            S = std::complex<double>(-1, 0);
            break;
        case 5:
            S = std::complex<double>(-sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
        case 6:
            S = std::complex<double>(0, 1);
            break;
        case 7:
            S = std::complex<double>(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
            break;
    }
    return S;
}


typedef struct{
    double a;
    double b;
    int K;
    int j;
    int q;
    std::complex<double> ExpAK;
    std::complex<double> ExpBK;
    std::complex<double> ExpAK_inverse;
    std::complex<double> ExpBK_inverse;
    std::complex<double> ExpAB;
    std::complex<double> ExpABK;
    std::complex<double> C1;
    std::complex<double> C5;
    std::complex<double> C7;
    std::complex<double> C8;
} theta_cache;



inline double K_power(int l, const theta_cache * cache) {
    return * ( (double *)((intptr_t)(cache) + sizeof(theta_cache) + (cache->j + l) * sizeof(double) ));
}

inline double root_2pi_b_power(int l, const theta_cache * cache) {
    return * ( (double *)((intptr_t)(cache) + sizeof(theta_cache) + (3 * cache->j + 2 + l) * sizeof(double) ));
}

theta_cache * build_theta_cache(mpfr_t mp_a, mpfr_t mp_b, int j, int K);
void free_theta_cache(theta_cache * cache);

std::complex<double> ExpA(mpfr_t A, int K);
std::complex<double> ExpAK(mpfr_t A, int K);
std::complex<double> ExpB(mpfr_t B, int K);
std::complex<double> ExpBK(mpfr_t B, int K);
std::complex<double> ExpAB(mpfr_t A, mpfr_t B);
std::complex<double> ExpABK(mpfr_t A, mpfr_t B, int K);



inline std::complex<double> EXP(std::complex<double> z) {                                                 //--------------------------------------------
                                                                                // This is just here for counting purposes.
    /*
    static __thread double rlast = 0.0;
    static __thread double ilast = 0.0;
    static __thread unsigned long total = 0;
    static __thread unsigned long same_as_last = 0;
    total++;
    if(real(z) == rlast && imag(z) == ilast)
        same_as_last++;
    rlast = real(z);
    ilast = imag(z);
    if(total % 1000000 == 0) {
        cout << "complex: " << same_as_last/(double)total << endl;
    }*/
    return std::exp(z);                                                         // There is pretty much no overhead in using it,
}                                                                               // and commenting out the first line of this
                                                                                // function will get rid of any overhead at all.
                                                                                // ---------------------------------------------
inline double EXP(double z) {
    /*
    static __thread double last = 0.0;
    static __thread unsigned long total = 0;
    static __thread unsigned long same_as_last = 0;
    cout << z << endl;
    total++;
    if(z == last)
        same_as_last++;
    last = z;
    if(total % 1000000 == 0) {
        cout << "   real: " << same_as_last/(double)total << endl;
    }
    */
    return std::exp(z);
}


inline bool check_condition(bool condition, char * message) {              //----------------------------------------------
    if(!condition) {                                                            //
        std::cout << "WARNING: " << message << std::endl;                       //
    }                                                                           //
    return condition;                                                           //
}                                                                               //----------------------------------------------


                                                                                //----------------------------------------------
                                                                                //  Functions to compute the integral
template<typename T> std::complex<double> H(int j, T alpha, double epsilon);                 //    / 1 
template<typename T> std::complex<double> H_method1(int j, T alpha);                         //    |   j
//inline std::complex<double> H_method2(int j, std::complex<double> alpha, double epsilon);               //    |  t  exp(- 2 pi alpha) dt
//std::complex<double> H_method3(int j, std::complex<double> alpha, double epsilon);                      //    |
template<typename T> T H_method4(int j, T alpha, double epsilon);               //    / 0
                                                                                //  (Defined in H_functions.cc)
                                                                                //----------------------------------------------


                                                                                        //----------------------------------------------
                                                                                        //  Functions to compute the integral         
                                                                                        //    
                                                                                        //      / 
std::complex<double> G(std::complex<double> alpha, double b, int n, int j, double epsilon, int method = 0);       //      |        j                                   2
std::complex<double> G_I(std::complex<double> alpha, double b, int n, int j, double epsilon, int method = 0);     //      | (t + n)  exp(2 pi i alpha t + 2 pi i beta t ) dt
std::complex<double> G_I_over_twopi(std::complex<double> alpha, int n, int j, double epsilon, int method = 0);    //      |
std::complex<double> G_method1_R(std::complex<double> alpha, double b, int n, int j, double epsilon);             //      /
std::complex<double> G_method1_I(std::complex<double> alpha, double b, int n, int j, double epsilon);             //
std::complex<double> G_via_Euler_MacLaurin_R(std::complex<double> alpha, double b, int n, int j, double epsilon); //    
std::complex<double> G_via_Euler_MacLaurin_I(std::complex<double> alpha, double b, int n, int j, double epsilon); //    
                                                                                        //  For complex parameters alpha and b.
                                                                                        //  (Defined in G_functions.cc)
                                                                                        //----------------------------------------------

std::complex<double> H_Integral_0(int j, double a, int M, double epsilon);                                        //--------------------
std::complex<double> J_Integral_0(double a, double b, int j, int M, int K, theta_cache * cache, double epsilon);  //
std::complex<double> J_Integral_1(double a, double b, int j, int M, int K, theta_cache * cache, double epsilon);  // Various integrals
std::complex<double> H_Integral_2(int j, double a1, double a2, double epsilon);                                   //
std::complex<double> J_Integral_2(double a1, double a2, double b, theta_cache * cache, double epsilon);           //--------------------

inline std::complex<double> JBulk(double a, double b, int j, int M, int K, theta_cache * cache, double epsilon) {
    double x = epsilon/2;
    std::complex<double> A = J_Integral_0(a, b, j, M, K, cache, x);
    std::complex<double> B = J_Integral_1(a, b, j, M, K, cache, x);
    return A + B;
}                                                                                       

void JBulk(std::complex<double> * J, double a, double b, int j, int M, int K, theta_cache * cache, double * epsilon);
void JBoundary(std::complex<double> * J, double a1, double a2, double b, int j, int K, theta_cache * cache, double * epsilon);

inline std::complex<double> JBoundary(double a1, double a2, double b, int j, int K, theta_cache * cache, double epsilon){ 
    if(j == 0) {
        double x = epsilon * .33333333333333333333;
        return J_Integral_2(a1, a2, b, cache, x)                                    
                                        + J_Integral_1(a1, b, j, -1, K, cache, x)     
                                        - J_Integral_1(a2, b, j, -1, K, cache, x);
    }
    else {
        double x = epsilon * .25;
        return J_Integral_0(a1, b, j, -1, K, cache, x) + J_Integral_1(a1, b, j, -1, K, cache, x)
                + (double)minus_one_power(j + 1) * (  J_Integral_0(a2, b, j, -1, K, cache, x) + J_Integral_1(a2, b, j, -1, K, cache, x) );
    }
}                                                                                     





std::complex<double> IC0(int j, mpfr_t mp_a, mpfr_t mp_b, const theta_cache * cache, double epsilon);
std::complex<double> IC1(int K, int j, double a, double b, const theta_cache * cache, double epsilon);//----------------------------------------------
std::complex<double> IC1c(int K, int j, double a, double b, std::complex<double> C8, const theta_cache * cache, double epsilon);     //
inline std::complex<double> IC3(int K, int j, double a, double b, const theta_cache * cache, double epsilon);           //
inline std::complex<double> IC3c(int K, int j, double a, double b, const theta_cache * cache, double epsilon);          //
std::complex<double> IC4(int K, int j, double a, double b, std::complex<double> C11, const theta_cache * cache, double epsilon);     //
std::complex<double> IC4c(int K, int j, double a, double b, std::complex<double> C11, const theta_cache * cache, double epsilon);    //
std::complex<double> IC5(int K, int j, double a, double b, const theta_cache * cache, double epsilon);                  //
std::complex<double> IC6(int K, int j, double a, double b, mpfr_t mp_a, const theta_cache * cache, double epsilon);     //      
std::complex<double> IC7(int K, int j, double a, double b, const theta_cache * cache, double epsilon);                  //
std::complex<double> IC7_method1(int K, int j, double a, double b, const theta_cache * cache, double epsilon, int L);                  //
std::complex<double> IC7star(double a, int j, double epsilon);                  //
std::complex<double> IC8(int K, int j, mpfr_t mp_a, mpfr_t mp_b, const theta_cache * cache);                            //
std::complex<double> IC9E(int K, int j, double a, double b, const theta_cache * cache, double epsilon);                 //
inline std::complex<double> IC9H(int K, int j, double a, double b, const theta_cache * cache, double epsilon);          //
                                                                                //  (Defined in ICn.cc, unless defined inline below)
                                                                                //
                                                                                //

inline std::complex<double> IC3(int K, int j, double a, double b, const theta_cache * cache, double epsilon) {
    // needs a <= 0
    //       b >= 0
    //return pow(-I, j+1) * IC9H(K, j, -a, b, cache, epsilon);
    return minus_I_power(j+1) * IC9H(K, j, -a, b, cache, epsilon);
}

inline std::complex<double> IC3c(int K, int j, double a, double b, const theta_cache * cache, double epsilon) {
    // needs a >= 0
    //       b >= 0
    //return pow(I, j+1) * IC9H(K, j, a, b, cache, epsilon);
    return I_power(j+1) * IC9H(K, j, a, b, cache, epsilon);
}

inline std::complex<double> IC9H(int K, int j, double a, double b, const theta_cache * cache, double epsilon) {
    //
    // Compute the integral (1/K^j)int_0^\infty t^j exp(-2 pi a t - 2 pi i b t^2)
    //
    // after a change of contour, this is well approximated by IC7(K, a, b) with large K
    //

    //int endpoint = to_int(10 * ceil( max(-LOG(epsilon * sqrt(b)), 1.0)/sqrt(b) ));
    //double z = pow((double)endpoint/(double)K, j);
    //double z = pow(K, j);
    std::complex<double> S = IC7(-1, j, a, b, cache, epsilon * K_power(j, cache));
    S = S * K_power(-j, cache);
    return S;
}


                                                                                //-----------------------------------------------------------
                                                                                //
                                                                                // Derivatives of the function
                                                                                //
                                                                                // g(t) = exp(2 pi i alpha t + 2 pi i b t^2)
                                                                                //
//void initialize_power_arrays(int n, std::complex<double> alpha, std::complex<double> b);                // at 0, 1, and K. Only supports derivatives up
//std::complex<double> g_derivative_at_1(int n);                                             // to the 21st, and initialize_power_arrays() must
//std::complex<double> g_derivative_at_0(int n);                                             // be called before the other functions, with n
//std::complex<double> g_derivative_at_K_without_exponential_factor(int n, int K);           // at least as large as any derivatives to be
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
//std::complex<double> g_derivative_at_1(int n, std::complex<double> a, std::complex<double> b);
void g_derivative_polynomial(int n, std::complex<double> * p, std::complex<double> * q, const std::complex<double> a, const std::complex<double> b);
void g_derivative_polynomial_R(int n, std::complex<double> * p, std::complex<double> * q, const std::complex<double> a, const double b);
void g_derivative_polynomial_I(int n, std::complex<double> * p, std::complex<double> * q, const std::complex<double> a, const double b);
void g_derivative_polynomial_I_over_twopi(int n, std::complex<double> * p, std::complex<double> * q, const std::complex<double> a);


//------------------------------------------------------------------------------------------------------------
//
// All of the functions below are currently defined in theta_sums.cc


std::complex<double> S1(int K, mpfr_t mp_a, mpfr_t mp_b, double epsilon);
std::complex<double> S2(int K, mpfr_t mp_a, mpfr_t mp_b, double epsilon);

int normalize(double &a, double &b);
int normalize(mpfr_t a, mpfr_t b);


std::complex<double> direct_exponential_sum_evaluation2(mpfr_t a, mpfr_t b, int j, int m, int M);


inline std::complex<double> compute_C11(mpfr_t a, mpfr_t b, int K) {
    //
    // Compute C11 = I exp(2 pi i K a + 2 pi i b K^2)
    //

    MPFR_DECL_INIT(tmp1, mpfr_get_prec(a));
    MPFR_DECL_INIT(tmp2, mpfr_get_prec(a));

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
    
    std::complex<double> S(mpfr_get_d(tmp2, GMP_RNDN), mpfr_get_d(tmp1, GMP_RNDN));
    return I * S;
}

inline std::complex<double> compute_C12(mpfr_t mp_a, mpfr_t mp_b, int K) {
    std::complex<double> z = ExpBK(mp_b, K);
    double a = mpfr_get_d(mp_a, GMP_RNDN);
    double b = mpfr_get_d(mp_b, GMP_RNDN);
    
    std::complex<double> w = I * exp(-2 * PI * (a + 2 * b * K) * K) / (sqrt(2 * PI * b));
 
    return w/z;
}

std::complex<double> direct_exponential_sum_evaluation2(double a, double b, int j, int m, int M, int working_precision = 53);
std::complex<double> compute_exponential_sums_using_theta_algorithm(mpfr_t mp_a, mpfr_t mp_b, int j, int K, std::complex<double> * v, double epsilon, int _Kmin);
std::complex<double> compute_exponential_sums_using_theta_algorithm2(mpfr_t mp_a, mpfr_t mp_b, int j, int K, std::complex<double> * v, double epsilon, int _Kmin);
std::complex<double> compute_exponential_sums_directly(mpfr_t mp_a, mpfr_t mp_b, int j, int K, std::complex<double> * v, double epsilon, int method);
std::complex<double> compute_exponential_sums_for_small_b(mpfr_t mp_a, mpfr_t mp_b, int j, int K, std::complex<double> * v, double epsilon);
std::complex<double> compute_exponential_sums(mpfr_t mp_a, mpfr_t mp_b, int j, int K, std::complex<double> * v, double epsilon, int _Kmin = 0, int method=0);
std::complex<double> compute_exponential_sums(double a, double b, int j, int K, std::complex<double> * v, double epsilon, int _Kmin = 0, int method=0);



//std::complex<double> w_coefficient(double * a_powers, double * b_powers, double * q_powers, int s, int j, theta_cache * cache, std::complex<double> * inner_sums);
void compute_subsum_coefficients(std::complex<double> * v2, std::complex<double> * v, const theta_cache * cache);
