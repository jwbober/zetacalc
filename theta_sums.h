#include <iostream>
#include <complex>
#include <string>
#include <cmath>

#include "mpfr.h"

typedef double Double;
typedef std::complex<double> Complex;

const Double PI = 3.14159265358979323846264338327950288;
const Double E = 2.7182818284590452353602874713526624977572470936999595749670;
const Complex I = Complex(0, 1);

const int Kmin = 100;

using namespace std;

inline int to_int(int x) {
    return x;
}

inline double to_double(double x) {
    return x;
}

inline double LOG(double x) {
    return log(x);
}


namespace verbose {
    const int IC0 = 0;
    const int IC1c = 0;
    const int IC6 = 0;
    const int IC7 = 0;
    const int compute_exponential_sum = 0;
    const int S1 = 0;
    const int S2 = 0;
    const int G = 0;
    const int H = 0;
}

namespace stats {
    extern int H_method1;
    extern int H_method2;
    extern int H_method3;
    extern int H_method4;

    extern int G_method1;
    extern int G_method2;

    extern int exp;

    extern int exponential_sum_called;
    extern int exponential_sum_euler_maclaurin;
    extern int exponential_sum_taylor_expansion;

    extern int H_Integral_0;
    extern int H_Integral_2;
    extern int J_Integral_0;
    extern int J_Integral_1;
    extern int J_Integral_2;
}

inline void print_stats() {
    std::cout << "In function H():" << std::endl;
    std::cout << "    Method 1 used " << stats::H_method1 << " times." << std::endl;
    std::cout << "    Method 2 used " << stats::H_method2 << " times." << std::endl;
    std::cout << "    Method 3 used " << stats::H_method3 << " times." << std::endl;
    std::cout << "    Method 4 used " << stats::H_method4 << " times." << std::endl;

    std::cout << "In function G():" << std::endl;
    std::cout << "    Method 1 used " << stats::G_method1 << " times." << std::endl;
    std::cout << "    Method 2 used " << stats::G_method2 << " times." << std::endl;

    std::cout << "EXP() called " << stats::exp << " times." << std::endl;

    std::cout << "compute_exponential_sum() called " << stats::exponential_sum_called << " times." << std::endl;
    std::cout << "[should have] used Euler-Maclaurin " << stats::exponential_sum_euler_maclaurin << " times." << std::endl;
    std::cout << "[should have] used Taylor expansions " << stats::exponential_sum_taylor_expansion << " times." << std::endl;

    std::cout << "H_Integral_0 called " << stats::H_Integral_0 << " times." << std::endl;
    std::cout << "H_Integral_2 called " << stats::H_Integral_2 << " times." << std::endl;
    std::cout << "J_Integral_0 called " << stats::J_Integral_0 << " times." << std::endl;
    std::cout << "J_Integral_1 called " << stats::J_Integral_1 << " times." << std::endl;
    std::cout << "J_Integral_2 called " << stats::J_Integral_2 << " times." << std::endl;

}

inline Complex ExpA(mpfr_t A, int K) {
    // return exp(-2 pi i A K)
    //
    // We use mpfr here because even if we only want the answer to 53 bits of
    // precision we might need to specify A and A*K to high precision.
    // 
    
    mpfr_t real_part;
    mpfr_t imag_part;
    mpfr_t tmp;

    mpfr_init2(real_part, 53);
    mpfr_init2(imag_part, 53);
    mpfr_init2(tmp, mpfr_get_prec(A));
    
    mpfr_const_pi(tmp, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, -2, GMP_RNDN);
    mpfr_mul(tmp, tmp, A, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, K, GMP_RNDN);

    mpfr_sin_cos(imag_part, real_part, tmp, GMP_RNDN);

    Complex S(mpfr_get_d(real_part, GMP_RNDN), mpfr_get_d(imag_part, GMP_RNDN));

    mpfr_clear(real_part);
    mpfr_clear(imag_part);
    mpfr_clear(tmp);

    return S;
}

inline Complex ExpB(mpfr_t B, int K) {
    // return exp(-2 pi i B K^2)
    //
    // We use mpfr here because even if we only want the answer to 53 bits of
    // precision we might need to specify B and BK^2 to high precision.
    // 
    
    mpfr_t real_part;
    mpfr_t imag_part;
    mpfr_t tmp;

    mpfr_init2(real_part, 53);
    mpfr_init2(imag_part, 53);
    mpfr_init2(tmp, mpfr_get_prec(B));
    
    mpfr_const_pi(tmp, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, -2, GMP_RNDN);
    mpfr_mul(tmp, tmp, B, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, K, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, K, GMP_RNDN);

    mpfr_sin_cos(imag_part, real_part, tmp, GMP_RNDN);

    Complex S(mpfr_get_d(real_part, GMP_RNDN), mpfr_get_d(imag_part, GMP_RNDN));

    mpfr_clear(real_part);
    mpfr_clear(imag_part);
    mpfr_clear(tmp);

    return S;
}


inline Complex ExpAB(mpfr_t A, mpfr_t B){
	// return exp(- pi i A^2 / (2 B) )
	//
	// We use mpfr to avoid loss of precision. Loss of 
	// precision can occur because B can be as small 
	// as 1/K, so A^2/B can  be as large as K.
	//
	
    mpfr_t real_part;
    mpfr_t imag_part;
    mpfr_t tmp;

    mpfr_init2(real_part, 53);
    mpfr_init2(imag_part, 53);
    mpfr_init2(tmp, mpfr_get_prec(A));
    
    mpfr_const_pi(tmp, GMP_RNDN);
    mpfr_mul_si(tmp, tmp, -1, GMP_RNDN);
    mpfr_mul(tmp, tmp, A, GMP_RNDN);
    mpfr_mul(tmp, tmp, A, GMP_RNDN);
    mpfr_div(tmp, tmp, B, GMP_RNDN);
    mpfr_div_ui(tmp, tmp, 2, GMP_RNDN);

    mpfr_sin_cos(imag_part, real_part, tmp, GMP_RNDN);

    Complex S(mpfr_get_d(real_part, GMP_RNDN), mpfr_get_d(imag_part, GMP_RNDN));

    mpfr_clear(real_part);
    mpfr_clear(imag_part);
    mpfr_clear(tmp);

    return S;
}


inline Complex EXP(Complex z) {                                                 //--------------------------------------------
    stats::exp++;                                                               // This is just here for counting purposes.
    return std::exp(z);                                                         // There is pretty much no overhead in using it,
}                                                                               // and commenting out the first line of this
                                                                                // function will get rid of any overhead at all.
                                                                                // ---------------------------------------------

inline bool check_condition(bool condition, std::string message) {              //----------------------------------------------
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
Complex G(Complex alpha, Complex b, Double epsilon, int method = 0);            //    |  exp(2 pi i alpha t + 2 pi i b t^2) dt
Complex G(Complex alpha, Complex b, int n, int j, Double epsilon, int method = 0);  //
Complex G_via_Euler_MacLaurin(Complex alpha, Complex b, Double epsilon);        //    |
Complex G_via_Euler_MacLaurin(Complex alpha, Complex b, int n, int j, Double epsilon); //    / 0
                                                                                //  For complex parameters alpha and b.
                                                                                //  (Defined in G_functions.cc)
                                                                                //----------------------------------------------

Complex H_Integral_0(int j, Double a, int M, Double epsilon);                   //----------------------------------------------
Complex J_Integral_0(Double a, Double b, int M, Double epsilon);                //
Complex J_Integral_0(Double a, Double b, int j, int M, int K, Double epsilon);  //
Complex J_Integral_1(Double a, Double b, int M, int K, Double epsilon);         // 
Complex J_Integral_1(Double a, Double b, int j, int M, int K, Double epsilon);  //
Complex H_Integral_2(int j, Double a1, Double a2, Double epsilon);              //
Complex J_Integral_2(Double a1, Double a2, Double b, Double epsilon);           //
Complex J_Integral_2(Double a1, Double a2, Double b, int j, int K, Double epsilon);// Various integrals.
                                                                                //
inline Complex JBulk(Double a, Double b, int M, int K, Double epsilon) {        //                         
    return J_Integral_0(a, b, M, epsilon/2)                                     // See H_and_J_integrals.cc
                                        + J_Integral_1(a, b, M, K, epsilon/2);  //
}                                                                               //
                                                                                //
inline Complex JBoundary(Double a1, Double a2, Double b, int K, Double epsilon){//
    return J_Integral_2(a1, a2, b, epsilon/3)                                   //
                                        + J_Integral_1(a1, b, -1, K, epsilon/3) //
                                        - J_Integral_1(a2, b, -1, K, epsilon/3);//
}                                                                               //----------------------------------------------

inline Complex JBulk(Double a, Double b, int j, int M, int K, Double epsilon) {         //                         
    return J_Integral_0(a, b, j, M, K, epsilon/2)                                          // See H_and_J_integrals.cc
                                        + J_Integral_1(a, b, j, M, K, epsilon/2);       //
}                                                                                       //
                                                                                        //
inline Complex JBoundary(Double a1, Double a2, Double b, int j, int K, Double epsilon){ //
    return J_Integral_2(a1, a2, b, j, K, epsilon/3)                                     //
                                        + J_Integral_1(a1, b, j, -1, K, epsilon/3)      //
                                        - J_Integral_1(a2, b, j, -1, K, epsilon/3);     //
}                                                                                       //----------------------------------------------





Complex IC0(int K, Double a, Double b, Complex C11, Complex C12, mpfr_t mp_a, mpfr_t mp_b, Double epsilon);
Complex IC1(int K, Double a, Double b, Complex C11, Complex C12, Double epsilon);//----------------------------------------------
Complex IC1c(int K, Double a, Double b, Complex C8, Double epsilon);            //
Complex IC1c(int K, int j, Double a, Double b, Complex C8, Double epsilon);     //
inline Complex IC3(Double a, Double b, Double epsilon);                         //
inline Complex IC3c(Double a, Double b, Double epsilon);                        //
inline Complex IC4(int K, Double a, Double b, Complex C11, Double epsilon);     //
inline Complex IC4c(int K, Double a, Double b, Complex C11, Double epsilon);    //
Complex IC5(Double a, Double b, Double epsilon);                                //  Computations of the integral of the function
Complex IC5(int K, int j, Double a, Double b, Double epsilon);                  //
Complex IC6(int K, Double a, Double b, Double epsilon);                         //      
Complex IC7(int K, Double a, Double b, Double epsilon);                         //       exp(2 pi i a t + 2 pi i b t^2)
Complex IC7(int K, int j, Double a, Double b, Double epsilon);                  //
inline Complex IC8(Double a, Double b, Double epsilon);                         //
Complex IC9E(int K, Double a, Double b, Double epsilon);                        //
Complex IC9E(int K, int j, Double a, Double b, Double epsilon);                 //
inline Complex IC9H(Double a, Double b, Double epsilon);                        //  along various contours in the complex plane.
                                                                                //  (Defined in ICn.cc, unless defined inline below)
                                                                                //
                                                                                //
inline Complex IC8(Double a, Double b, Double epsilon) {                        //
    return exp(I * PI * (.25 - a * a/(2.0 * b)))/sqrt((Double)2 * b);           //----------------------------------------------
}

inline Complex IC3(Double a, Double b, Double epsilon) {
    // needs a <= 0
    //       b >= 0
    return -I * IC9H(-a, b, epsilon);
}

inline Complex IC3c(Double a, Double b, Double epsilon) {
    // needs a >= 0
    //       b >= 0
    return I * IC9H(a, b, epsilon);
}

inline Complex IC4(int K, Double a, Double b, Complex C11, Double epsilon) {
    // needs a + 2bK <= 0
    // b is always positive, so this will be satisfied if 
    return -C11 * IC9H(-(a + 2 * b * (Double)K), b, epsilon);
}

inline Complex IC4c(int K, Double a, Double b, Complex C11, Double epsilon) {
    // needs a + 2bK <= 0
    // b is always positive, so this will be satisfied if 
    return C11 * IC9H((a + 2 * b * (Double)K), b, epsilon);
}

inline Complex IC9H(Double a, Double b, Double epsilon) {
    //
    // Compute the integral int_0^\infty exp(-2 pi a t - 2 pi i b t^2)
    //
    // after a change of contour, this is well approximated by IC7(K, a, b) with large K
    //

    int K = to_int(10 * ceil( max(-LOG(epsilon * sqrt(b)), 1.0)/sqrt(b) ));
    return IC7(K, a, b, epsilon);
}

inline Complex IC9H(int K, int j, Double a, Double b, Double epsilon) {
    //
    // Compute the integral (1/K^j)int_0^\infty t^j exp(-2 pi a t - 2 pi i b t^2)
    //
    // after a change of contour, this is well approximated by IC7(K, a, b) with large K
    //

    int endpoint = to_int(10 * ceil( max(-LOG(epsilon * sqrt(b)), 1.0)/sqrt(b) ));
    Double z = pow((Double)endpoint/(Double)K, j);
    Complex S = IC7(endpoint, a, b, epsilon/z);
    S = S * z;
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


//------------------------------------------------------------------------------------------------------------
//
// All of the functions below are currently defined in theta_sums.cc


Complex S1(int K, mpfr_t mp_a, mpfr_t mp_b, Double epsilon);
Complex S2(int K, mpfr_t mp_a, mpfr_t mp_b, Double epsilon);
Complex compute_exponential_sum_via_Euler_Maclaurin(Double a, Double b, int K, Double epsilon);
Complex compute_exponential_sum_via_Euler_Maclaurin(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon);
Complex compute_exponential_sum_for_small_b(Double a, Double b, int K, Double epsilon);
Complex compute_exponential_sum_for_small_b(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon);
Complex compute_exponential_sum(mpfr_t mp_a, mpfr_t mp_b, int K, Double epsilon, int method = 0);
Complex compute_exponential_sum(Double a, Double b, int K, Double epsilon, int method = 0);

int normalize(Double &a, Double &b);
int normalize(mpfr_t a, mpfr_t b);


Complex direct_exponential_sum_evaluation(Double alpha, Double beta, int m, int M, int working_precision = 53);
Complex direct_exponential_sum_evaluation(mpfr_t a, mpfr_t b, int m, int M);
Double sum_of_offset_inverse_powers(Double a, int m, int M, int j, Double epsilon, int method = 0);
Double infinite_sum_of_differenced_inverse_powers(Double a1, Double a2, int m, int j, Double epsilon);



inline Complex compute_C11(mpfr_t a, mpfr_t b, int K) {
    //
    // Compute C11 = I exp(2 pi i a + 2 pi i b K2)
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
    Complex z = ExpB(mp_b, K);
    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    Double b = mpfr_get_d(mp_b, GMP_RNDN);
    
    Complex w = I * exp(-2 * PI * (a + 2 * b * K) * K) / (sqrt(2 * PI * b));
 
    return w * z;
}

Complex theta_sum2(Double a, Double b, int K, Double epsilon);
