#include "theta_sums.h"
#include "zeta.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

const char * PRECOMPUTATION_LOCATION = "/home/bober/math/experiments/theta_sums/caches/";

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

    outfile << "zeta_data3 = [";

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

void * F0_thread(void * bound ) {
    build_F0_cache(10, 100, 30, (long)(bound), exp(-30), PRECOMPUTATION_LOCATION);
    pthread_exit(NULL);
}
void * F1_thread(void * unused) {
    build_F1_cache(30, 200, 30, exp(-30), PRECOMPUTATION_LOCATION);
    pthread_exit(NULL);
}
void * F2_thread(void * bound) {
    build_F2_cache((long)(bound), 10, 10, 100, exp(-30), PRECOMPUTATION_LOCATION);
    pthread_exit(NULL);
}
void * IC7_thread(void * bound) {
    build_IC7_cache(200, (long)bound, 35, exp(-30), PRECOMPUTATION_LOCATION);
    pthread_exit(NULL);
}

void do_precomputation(mpfr_t t) {
    mpz_t v;
    mpz_init(v);

    stage_3_bound(v, t);
    int largest_block_size = stage_3_block_size(mpz_get_d(v), mpfr_get_d(t, GMP_RNDN)) + 500;

    cout << "Doing a precomputation for block sizes up to " << largest_block_size << endl;

    pthread_t thread0, thread1, thread2, thread3;

    pthread_create(&thread0, NULL, F0_thread, (void *)(largest_block_size/2));
    pthread_create(&thread1, NULL, F1_thread, (void *)(largest_block_size/2));
    pthread_create(&thread2, NULL, F2_thread, (void *)(largest_block_size/2));
    pthread_create(&thread3, NULL, IC7_thread, (void *)( (int)(sqrt(largest_block_size)) ));

    void * status;
    pthread_join(thread0, &status);
    pthread_join(thread1, &status);
    pthread_join(thread2, &status);
    pthread_join(thread3, &status);
    
    //build_F0_cache(10, 100, 30, largest_block_size/2, exp(-30));
    //build_F1_cache(30, 200, 30, exp(-30));
    //build_F2_cache(largest_block_size/2, 10, 10, 100, exp(-30));
    //build_IC7_cache((int)sqrt(largest_block_size), 200, 35, exp(-30));

    mpz_clear(v);
}

namespace zeta_config {
    extern int stage2_number_of_threads;
    extern int stage3_number_of_threads;
}

void do_computation() {
    mpfr_t t;
    mpfr_init2(t, 150);

    int N = 500;
    Double delta = .01;
    Complex main_sum_values[N];

    ofstream logfile;
    logfile.open("zeta_logfile");
    logfile << setprecision(17);

    for(int n = 16; n <= 20; n++) {
        time_t start_time = time(NULL);
        mpfr_set_si(t, 10, GMP_RNDN);
        mpfr_pow_si(t, t, n, GMP_RNDN);
        mpfr_sub_si(t, t, 5, GMP_RNDN);

        mpz_t tz;
        mpz_init(tz);
        mpfr_get_z(tz, t, GMP_RNDN);

        logfile << "Starting computation at t = " << tz << endl;
        logfile << "    Using delta = " << delta << endl;
        logfile << "    Computing at " << N << " points." << endl;
        logfile.flush();

        zeta_sum(t, delta, N, main_sum_values);

        ofstream datafile;
        stringstream filename_stream;
        filename_stream << "zeta_1e" << n << ".data";
        
        time_t end_time = time(NULL);

        logfile << "    Finished computation in " << end_time - start_time << " seconds." << endl;
        logfile << "    Writing data to file " << filename_stream.str() << endl;
        logfile.flush();

        datafile.open(filename_stream.str().c_str());
        datafile << setprecision(17);

        //mpfr_exp_t exponent;
        //char * t_string = mpfr_get_str(NULL, &exponent, 10, 0, GMP_RNDN);
        
        datafile << tz << endl;
        datafile << delta << endl;
        
        for(int k = 0; k < N; k++) {
            datafile << main_sum_values[k] << endl;
        }

        logfile << "    Finished writing data." << endl;
        logfile.flush();

        datafile.close();
        mpz_clear(tz);
    }
}

void usage() {
    const char * usage_text =
"Usage: zeta FILENAME\n\
\n\
FILENAME should be a plain text file, the first line\n\
of which is the largest t that zeta(1/2 + it) will be\n\
computed for, if the program is to using precomputation.\n\
Otherwise, the first line should be 0, indicating that\n\
no precomputation should be done.\n\
\n\
After the first line, each line should contain a single\n\
integer. Currently the program computes the main sum\n\
in the Riemann-Siegel formula at a grid of 500 points\n\
over an interval of length 5 starting at this integer\n\
and writes this data to a file zeta_N.data, where\n\
N is replaced by this integer.\n\
";
    cout << usage_text;

}

int main(int argc, char * argv[]) {
    if(argc != 2 && argc != 7) {
        usage();
        return 0;
    }

    if(argc == 7) {
        // we assume that the program was called with
        // ./zeta t delta N start length number_of_threads
        mpfr_t t;
        mpz_t start;
        mpz_t length;
        mpz_init(start);
        mpz_init(length);
        mpfr_init2(t, 200);
        mpfr_set_str(t, argv[1], 10, GMP_RNDN);
        double delta = strtod(argv[2], NULL);
        int N = atoi(argv[3]);
        mpz_set_str(start, argv[4], 10);
        mpz_set_str(length, argv[5], 10);
        
        int number_of_threads = atoi(argv[6]);

        zeta_config::stage2_number_of_threads = number_of_threads;
        zeta_config::stage3_number_of_threads = number_of_threads;

        Complex S[N];
        
        cout << setprecision(17) << endl;

        do_precomputation(t);
        partial_zeta_sum(start, length, t, delta, N, S);

        cout << "BEGIN RESULT" << endl;
        for(int k = 0; k < N; k++) {
            cout << S[k] << endl;
        }

        mpfr_clear(t);
        mpz_clear(start);
        mpz_clear(length);
        return 0;
    }

    zeta_config::stage2_number_of_threads = 4;
    zeta_config::stage3_number_of_threads = 4;

    ifstream input_file;
    input_file.open(argv[1]);
    if(!input_file.is_open()) {
        cout << "Error: Could not open file " << argv[1] << endl << endl;
        usage();
        return 0;
    }

    string line;
    getline(input_file, line);

    mpfr_t t;
    mpfr_init2(t, 200);
    mpfr_set_str(t, line.c_str(), 10, GMP_RNDN);
    if(mpfr_cmp_d(t, 0.0) != 0) {
        do_precomputation(t);
        // TODO: Write something to a logfile.
    }
    else {
        // TODO: Write something to a logfile.
    }

    ofstream logfile;
    logfile.open("zeta_logfile");
    logfile << setprecision(17);

    int N = 500;
    Double delta = .01;
    Complex main_sum_values[N];

    mpz_t tz;
    mpz_init(tz);

    while(getline(input_file, line)) {
        time_t start_time = time(NULL);
        
        mpz_set_str(tz, line.c_str(), 10);
        mpfr_set_z(t, tz, GMP_RNDN);

        logfile << "Starting computation at t = " << tz << endl;
        logfile << "    Using delta = " << delta << endl;
        logfile << "    Computing at " << N << " points." << endl;
        logfile.flush();

        zeta_sum2(t, delta, N, main_sum_values);

        ofstream datafile;
        stringstream filename_stream;
        filename_stream << "zeta_" << tz << ".data";
        
        time_t end_time = time(NULL);

        logfile << "    Finished computation in " << end_time - start_time << " seconds." << endl;
        logfile << "    Writing data to file " << filename_stream.str() << endl;
        logfile.flush();

        datafile.open(filename_stream.str().c_str());
        datafile << setprecision(17);

        //mpfr_exp_t exponent;
        //char * t_string = mpfr_get_str(NULL, &exponent, 10, 0, GMP_RNDN);
        
        datafile << tz << endl;
        datafile << delta << endl;
        
        for(int k = 0; k < N; k++) {
            datafile << main_sum_values[k] << endl;
        }

        logfile << "    Finished writing data." << endl;
        logfile.flush();

        datafile.close();

    }
    
    mpz_clear(tz);
    logfile.close();
    input_file.close();
    return 0;

    cout << setprecision(17) << endl;

    //mpfr_t t;
    //mpfr_init2(t, 150);
    //mpfr_set_str(t, "1e20", 10, GMP_RNDN);

    //do_precomputation(t);

    //do_computation();

    //print_stats();

    //return 0;

    //build_F0_cache(10, 100, 25, 500, exp(-30));
    //build_F1_cache(30, 200, 30, exp(-30));
    //build_F2_cache(500, 10, 10, 100, exp(-30));
    //build_IC7_cache(600, 200, 25, exp(-30));



    //compute_hardy_on_unit_interval(t);


    //int N = 1;
    //double delta = .01; 
    //Complex S[N];


    //mpfr_set_str(t, "1e18", 10, GMP_RNDN);
    //hardy_Z(t, delta, N, S);
    //cout << "1e18" << "  " << S[0] << "   " << S[0] * rs_rotation(t) << endl;

    //mpfr_set_str(t, "1e20", 10, GMP_RNDN);
    //hardy_Z(t, delta, N, S);
    //cout << "1e20" << "  " << S[0] << "   " << S[0] * rs_rotation(t) << endl;
    
    
    //mpfr_set_str(t, "1e20", 10, GMP_RNDN);
    //hardy_Z(t, delta, N, S);
    //cout << "1e20" << "  " << S[0] << "   " << S[0] * rs_rotation(t) << endl;
    

    //mpfr_set_str(t, "1e21", 10, GMP_RNDN);
    //hardy_Z(t, delta, N, S);
    //cout << "1e21" << "  " << S[0] << "   " << S[0] * rs_rotation(t) << endl;

    //mpfr_set_str(t, "1e22", 10, GMP_RNDN);
    //hardy_Z(t, delta, N, S);
    //cout << "1e22" << "  " << S[0] << "   " << S[0] * rs_rotation(t) << endl;

    //mpfr_set_str(t, "1e23", 10, GMP_RNDN);
    //hardy_Z(t, delta, N, S);
    //cout << "1e23" << "  " << S[0] << "   " << S[0] * rs_rotation(t) << endl;

    //mpfr_set_str(t, "1e24", 10, GMP_RNDN);
    //hardy_Z(t, delta, N, S);
    //cout << "1e24" << "  " << S[0] << "   " << S[0] * rs_rotation(t) << endl;

    //for(int l = 0; l < N; l++) {
    //    cout << S[l] << "    " << S[l] * rs_rotation(t) << endl;
    //    mpfr_add_d(t, t, delta, GMP_RNDN);
    //}


    print_stats();
    return 0;
}
