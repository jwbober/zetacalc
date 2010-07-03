#include "theta_sums.h"
#include "zeta.h"
#include "log.h"

namespace stats {
    int H_method1 = 0;
    int H_method2 = 0;
    int H_method3 = 0;
    int H_method4 = 0;

    int H_function_big = 0;
    int H_function_small = 0;

    int G_method1 = 0;
    int G_method2 = 0;

    int exp = 0;

    int exponential_sum_called = 0;
    int exponential_sum_euler_maclaurin = 0;
    int exponential_sum_taylor_expansion = 0;

    int H_Integral_0 = 0;
    int H_Integral_2 = 0;
    int J_Integral_0 = 0;
    int J_Integral_1 = 0;
    int J_Integral_2 = 0;

    int J_Integral_0_taylor_expansion = 0;
    int J_Integral_1_taylor_expansion = 0;
    int J_Integral_2_taylor_expansion = 0;

    int J_Integral_0_zero = 0;
    int J_Integral_1_zero = 0;
    int J_Integral_2_zero = 0;

    int J_Integral_0_terms_used = 0;
    int J_Integral_1_terms_used = 0;
    int J_Integral_2_terms_used = 0;

    int IC7 = 0;
    int IC7zero = 0;
    int IC7_taylor_expansion = 0;
    int IC7_terms_used = 0;
}


void print_stats() {
    std::cout << "In function H():" << std::endl;
    std::cout << "    Method 1 used " << stats::H_method1 << " times." << std::endl;
    std::cout << "    Method 2 used " << stats::H_method2 << " times." << std::endl;
    std::cout << "    Method 3 used " << stats::H_method3 << " times." << std::endl;
    std::cout << "    Method 4 used " << stats::H_method4 << " times." << std::endl;

    std::cout << "  abs(alpha) was big " << stats::H_function_big << " times." << std::endl;
    std::cout << "  abs(alpha) was small " << stats::H_function_small << " times." << std::endl;

    std::cout << "In function G():" << std::endl;
    std::cout << "    Method 1 used " << stats::G_method1 << " times." << std::endl;
    std::cout << "    Method 2 used " << stats::G_method2 << " times." << std::endl;

    std::cout << "EXP() called " << stats::exp << " times." << std::endl;

    std::cout << "compute_exponential_sum() called " << stats::exponential_sum_called << " times." << std::endl;
    std::cout << "[should have] used Euler-Maclaurin " << stats::exponential_sum_euler_maclaurin << " times." << std::endl;
    std::cout << "[should have] used Taylor expansions " << stats::exponential_sum_taylor_expansion << " times." << std::endl;

    std::cout << "H_Integral_0 called " << stats::H_Integral_0 << " times." << std::endl;
    std::cout << "H_Integral_2 called " << stats::H_Integral_2 << " times." << std::endl;

    std::cout << "J_Integral_0: " << std::endl;
    std::cout << "      called " << stats::J_Integral_0 << " times." << std::endl;
    std::cout << "      used taylor expansions " << stats::J_Integral_0_taylor_expansion << " times." << endl;
    std::cout << "          average number of terms used was " << (Double)stats::J_Integral_0_terms_used/(Double)stats::J_Integral_0_taylor_expansion << endl;
    std::cout << "      quickly returned 0 " << stats::J_Integral_0_zero << " times. " << std::endl;

    std::cout << "J_Integral_1: " << std::endl;
    std::cout << "      called " << stats::J_Integral_1 << " times." << std::endl;
    std::cout << "      used taylor expansions " << stats::J_Integral_1_taylor_expansion << " times." << endl;
    std::cout << "          average number of terms used was " << (Double)stats::J_Integral_1_terms_used/(Double)stats::J_Integral_1_taylor_expansion << endl;
    std::cout << "      quickly returned 0 " << stats::J_Integral_1_zero << " times. " << std::endl;

    std::cout << "J_Integral_2: " << std::endl;
    std::cout << "      called " << stats::J_Integral_2 << " times." << std::endl;
    std::cout << "      used taylor expansions " << stats::J_Integral_2_taylor_expansion << " times." << endl;
    std::cout << "          average number of terms used was " << (Double)stats::J_Integral_2_terms_used/(Double)stats::J_Integral_2_taylor_expansion << endl;
    std::cout << "      quickly returned 0 " << stats::J_Integral_2_zero << " times. " << std::endl;

    std::cout << "IC7: " << endl;
    std::cout << "      called " << stats::IC7 << " times." << endl;
    std::cout << "      used taylor expansions " << stats::IC7_taylor_expansion << " times." << endl;
    std::cout << "          average number of terms used was " << (Double)stats::IC7_terms_used/(Double)stats::IC7_taylor_expansion << " times." << endl;
    std::cout << "      quickly returned 0 " << stats::IC7zero << " times." << endl;


}

