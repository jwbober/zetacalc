
#define BUILTIN_EXPECT(a, b) __builtin_expect( (a), (b) )

typedef double Double;
typedef std::complex<double> Complex;

std::complex<double> w_coefficient_slow(mpfr_t mp_a, mpfr_t mp_b, int K, int s, int j, std::complex<double> CF);
void print_w_coefficient_stats();

double sum_of_offset_inverse_powers(double a, int m, int M, int j, double epsilon, int method = 0);
double infinite_sum_of_differenced_inverse_powers(double a1, double a2, int m, int j, double epsilon);

const double PI = 3.14159265358979323846264338327950288;
const double E = 2.7182818284590452353602874713526624977572470936999595749670;
const std::complex<double> I = std::complex<double>(0, 1);


inline double pow(int a, int b) {
    return pow((double)a, (double)b);
}

inline double pow(int a, double b) {
    return pow((double)a, (double)b);
}

inline double pow(int a, unsigned int b) {
    return pow((double)a, (double)b);
}

inline int to_int(int x) {
    return x;
}

inline double to_double(double x) {
    return x;
}

inline double LOG(double x) {
    return log(x);
}


