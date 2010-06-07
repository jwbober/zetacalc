
#define BUILTIN_EXPECT(a, b) __builtin_expect( (a), (b) )

typedef double Double;
typedef std::complex<double> Complex;

Complex w_coefficient(Double * a_powers, Double * b_powers, Double * q_powers, Double * K_powers, int s, int j, Complex CF);
void print_w_coefficient_stats();

Double sum_of_offset_inverse_powers(Double a, int m, int M, int j, Double epsilon, int method = 0);
Double infinite_sum_of_differenced_inverse_powers(Double a1, Double a2, int m, int j, Double epsilon);

const Double PI = 3.14159265358979323846264338327950288;
const Double E = 2.7182818284590452353602874713526624977572470936999595749670;
const Complex I = Complex(0, 1);



using namespace std;

inline Double pow(int a, int b) {
    return pow((Double)a, (Double)b);
}

inline Double pow(int a, Double b) {
    return pow((Double)a, (Double)b);
}

inline Double pow(int a, unsigned int b) {
    return pow((Double)a, (Double)b);
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


