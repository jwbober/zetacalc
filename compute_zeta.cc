#include "theta_sums.h"
#include "zeta.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//computes sinc function
Double sinc(Double x){
    //Taylor coefficients in sin expansions
    const Double sinc_coeff[6]={1, -1.0 / 6.0, 1.0 / 120.0, -1.0 / 5040.0, 1 / 362880.0, -1.0 / 39916800.0};

    //Parameters needed in sinc function; help overcome instability when sinc(x) is called with very small x 
    Double sinc_tol=1.e-5;
    int sinc_terms = 3;

    if(abs(x) < sinc_tol) {
        Double x_squared = x*x;
        Double x_power = x_squared;
        Double ans = 1;

        for(int j = 1; j < sinc_terms; j++){
            ans = ans + sinc_coeff[j] * x_power;
            x_power = x_squared * x_power;
        }
        
        return ans;
    }
    Double ans = sin(x) / x;
    return ans;
}


//computes the kernel function in BLFI
Double kernel(Double u, Double c, Double epsilon_1){
    //Taylor coefficients in sinh and sin expansions
    const Double sinh_coeff[6]={1, 1 / 6.0, 1 / 120.0, 1 / 5040.0, 1 / 362880.0, 1 / 39916800.0};
    //Parameters needed in kernel function; help overcome instability when kernel(x) is called with very small x 
    Double sinh_tol=1.e-5;
    int sinh_terms = 3;

    Double x = c*c - pow( epsilon_1*u, 2);

    if( x < 0 ) {
        cout << endl << "kernel: u is too large!" << endl;
        exit(1);
    }

    if(abs(x) < sinh_tol) {

        Double x_squared = x*x;
        Double x_power = x_squared;
        Double ans = 1;
        
        for(int j = 1; j < sinh_terms; j++) {
            ans = ans + sinh_coeff[j] * x_power;
            x_power = x_squared * x_power;
        }
        ans = ans * 2 * c/(exp(-c) + exp(c));
        return ans;
    }
 	
    Double sinh_x = 0.5 * (exp(sqrt(x)) + exp(- sqrt(x)));
    Double ans = sinh_x / sqrt(x);
    ans = ans * 2 * c/(exp(-c) + exp(c));
    return  ans;
}


//We compute the main sum using the band-limited function interpolation formula.
// We interpolate main_sum(1/2 + it) in  approximately the interval (t0 + c , t0 + max_n*PI/beta - c) to within t^0.25 * exp( - c )

Complex blfi_inter(Double t_minus_t0, Double max_n, Double c, Double start, Double end, Complex* main_sum_values){
 
    //we assume the blfi parameters are chosen as in (2.13) in "amortized complexity method...."
    //Double tau = log ( sqrt (t / (2 * PI)) ); // this is the heighest frequency in the main sum
    Double tau = (log(end) - log(start))/2;
    //Double beta = 1.5 * tau;  // note we assume main_sum_values are computed at integer multiples of PI / beta
    Double beta = 100 * tau;
    Double lambda = ( beta + tau ) / 2; // this is an extraneous blfi parameter
    Double epsilon_1 = ( beta - tau ) / 2; // an extraneous blfi parameter
    Double alpha = (log(end) + log(start))/2;

    int range = beta * c / ( PI * epsilon_1);  // 2*range is the number of terms that will be used in the bfli formula
    int n0 =  round( (t_minus_t0) * beta / PI ); // this is the closest grid point to t - t0

    //if t - t0 is too close to 0 or max_n*PI/beta then blfi doesn't work, so we exit
    if( n0  <  range  ||  n0  >  max_n  - range || 2 * range > max_n || range < 1) {
        cout << endl << "Error: blfi formula called out of range!" <<endl;
        exit(1);
    }

    // we interpolate the main sum here using precomputed values main_sum_values[n]
    Complex sum = 0;

    for(int n = n0 - range + 1; n < n0 +  range;  n++){
        Double u = n * PI /beta - (t_minus_t0);
        Complex z = main_sum_values[n] * exp(-I * alpha * u) * sinc(lambda * u) * kernel (u , c , epsilon_1);
        sum = sum + z;
    }

    sum = sum * lambda / beta;
    return sum;
}

void compute_hardy_on_unit_interval(mpfr_t t) {
    Double tt = mpfr_get_d(t, GMP_RNDN);
 
    Double end =  floor(sqrt(tt/(2 * PI)));
    Double start = 1.0;

    Double tau = (log(end) - log(start))/2;
    Double beta = 100 * tau;  // note we assume main_sum_values are computed at integer multiples of PI / beta
 
    Double delta = PI/beta;

    cout << "using delta = " << delta << endl;

    mpfr_t t0;
    mpfr_t tmp;
    mpfr_init2(t0, mpfr_get_prec(t));
    mpfr_init2(tmp, mpfr_get_prec(t));
    mpfr_sub_d(t0, t, 250 * delta, GMP_RNDN);

    mpfr_set_d(tmp, delta, GMP_RNDN);
    mpfr_fmod(tmp, t0, tmp, GMP_RNDN);
    mpfr_sub(t0, t0, tmp, GMP_RNDN);

    int N = (int)(1.0/delta) + 500;
    Complex main_sum[N];

    cout << "initially computing zeta at " << N << " points." << endl;

    zeta_sum(t0, delta, N, main_sum);
 
    cout << "Output from BLFI formula:" << endl;
    Double t_minus_t0;
    mpfr_sub(tmp, t, t0, GMP_RNDN);
    t_minus_t0 = mpfr_get_d(tmp, GMP_RNDN);

    ofstream outfile;
    outfile.open("zeta_data3.sage");

    outfile << "zeta_data2 = [";

    for(Double d = 0; d < 1; d += .001) {
        if(d != 0)
            outfile << ", ";
        Complex z = blfi_inter(t_minus_t0, N - 1, 30, start, end, main_sum);
        Double Z_value = 2 * real(rs_rotation(t) * z) + rs_remainder(t);
        Complex zeta_value = Z_value * rs_rotation(t);
        //cout << d << ", " << Z_value << "  " << zeta_value << endl;
        outfile << "(" << d << ", " << Z_value << ")";
        t_minus_t0 += .001;
        mpfr_add_d(t, t, .001, GMP_RNDN);
    }
    outfile << "]" << endl;

    outfile.close();
}


int main() {
    cout << setprecision(10) << endl;

    mpfr_t t;
    mpfr_init2(t, 150);
    mpfr_set_str(t, "1e18", 10, GMP_RNDN);

    //build_F0_cache(10, 100, 25, 500, exp(-30));
    //build_F1_cache(30, 200, 30, exp(-30));
    //build_F2_cache(500, 10, 10, 100, exp(-30));
    //build_IC7_cache(600, 200, 25, exp(-30));



    compute_hardy_on_unit_interval(t);

    print_stats();

    int N = 10;
    double delta = .01; 
    Complex S[N];
    return 0;


    hardy_Z(t, delta, N, S);
    for(int l = 0; l < N; l++) {
        cout << S[l] << "    " << S[l] * rs_rotation(t) << endl;
        mpfr_add_d(t, t, delta, GMP_RNDN);
    }

    return 0;
}
