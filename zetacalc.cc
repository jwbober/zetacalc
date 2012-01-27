#include "theta_sums.h"
#include "main_sum.h"

#include <getopt.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace std;

extern int stage3_start;

string output_filename = "";

void usage() {
    const char * usage_text =
"Usage: zeta FILENAME\n\
\n\
FILENAME should be a plain text file, which contains\n\
a list of numbers and a string. The numbers will be interpreted, in\n\
order, as t, start, length, N, delta, output_filename.\n\
\n\
The program will then compute the main sum in the Riemann-Siegel\n\
formula for zeta(.5 + it), from start to start + length - 1, at N points\n\
spaced delta apart, and write the output to output_filename.\n\
If length is 0 no work will be done. If it is negative, I don't know what\n\
will happen. In the future it will be setup so that the whole sum is computed.\n\
";
    cout << usage_text;

}

int process_input_file(char const* filename, int Kmin) {

    ifstream input_file;
    input_file.open(filename);
    if(!input_file.is_open()) {
        cout << "Error: Could not open file " << filename << endl << endl;
        usage();
        return 1;
    }

    clock_t start_cpu_time = clock();
    time_t start_wall_time = time(NULL);

    mpfr_t t;
    mpz_t start, length;

    mpfr_init2(t, 200);
    mpz_init(start);
    mpz_init(length);

    double delta;
    int N;
    string t_string;


    input_file >> t_string;
    mpfr_set_str(t, t_string.c_str(), 10, GMP_RNDN);
    input_file >> start;
    input_file >> length;
    input_file >> N;
    input_file >> delta;
    if(output_filename == "") {
        input_file >> output_filename;
    }
    else {
        string blah;
        input_file >> blah;
    }

    cout << start << endl;
    cout << length << endl;
    cout << N << endl;
    cout << delta << endl;
    cout << output_filename << endl;

    Complex * S = new Complex[N];
    partial_zeta_sum(start, length, t, delta, N, S, Kmin);

    ofstream output_file;
    output_file.open(output_filename.c_str());
    
    output_file << setprecision(17);

    output_file << t_string << " ";
    output_file << start << " ";
    output_file << length << " ";
    output_file << N << " ";
    output_file << delta << " ";
    output_file << output_filename << " ";

    clock_t end_cpu_time = clock();
    time_t end_wall_time = time(NULL);

    output_file << ((double)end_cpu_time - (double)start_cpu_time)/CLOCKS_PER_SEC << " ";
    output_file << end_wall_time - start_wall_time << endl;
    
    for(int k = 0; k < N; k++) {
         output_file << S[k] << endl;
    }
    delete [] S;
    return 0;

}



int main(int argc, char * argv[]) {
    
    // initial setup of default arguments

    int Kmin = 0;
    char const* filename = NULL;

    bool length_set = false;
    bool use_input_file = false;

    double delta = .05;
    int N = 20;
    int repeat = 1;

    int Z_flag = 0;
    int zeta_flag = 0;
    int verbose = 1;
    int number_of_threads = 2;
    int fraction = -1;

    mpfr_t t;
    mpfr_init2(t, 250);
    mpfr_set_str(t, "1000", 10, GMP_RNDN);

    mpz_t start;
    mpz_t length;

    mpz_init(start);
    mpz_init(length);

    mpz_set_str(start, "1", 10);

    while (1) {
        enum {KMIN = 2, T, START, LENGTH, DELTA, DOPRECOMPUTATION, USEPRECOMPUTATION, FILENAME, NUMTHREADS, N_OPTION, STAGE3_START, OUTPUT, REPEAT, FRACTION};

        static struct option options[] = 
            {
                {"Kmin", required_argument, 0, KMIN},
                {"t", required_argument, 0, T},
                {"start", required_argument, 0, START},
                {"length", required_argument, 0, LENGTH},
                {"delta", required_argument, 0, DELTA},
                {"Z", no_argument, &Z_flag, 1},
                {"zeta", no_argument, &zeta_flag, 1},
                {"filename", required_argument, 0, FILENAME},
                {"number_of_threads", required_argument, 0, NUMTHREADS},
                {"fraction", required_argument, 0, FRACTION},
                {"N", required_argument, 0, N_OPTION},
                {"stage3_start", required_argument, 0, STAGE3_START},
                {"output", required_argument, 0, OUTPUT},
                {"repeat", required_argument, 0, REPEAT},
                {"terse", no_argument, &verbose, 0},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        int c = getopt_long(argc, argv, "", options, &option_index);
        if (c == -1)
            break;

        switch(options[option_index].val) {
            case KMIN:
                Kmin = atoi(optarg);
                break;
            case T:
                mpfr_set_str(t, optarg, 10, GMP_RNDN);
                break;
            case START:
                mpz_set_str(start, optarg, 10);
                break;
            case LENGTH:
                length_set = true;
                mpz_set_str(length, optarg, 10);
                break;
            case DELTA:
                delta = atof(optarg);
                break;
            case FILENAME:
                filename = optarg;
                use_input_file = true;
                break;
            case NUMTHREADS:
                number_of_threads = atoi(optarg);
                break;
            case FRACTION:
                fraction = atoi(optarg);
                break;
            case N_OPTION:
                N = atoi(optarg);
                break;
            case STAGE3_START:
                stage3_start = strtoul(optarg, NULL, 10);
                break;
            case OUTPUT:
                output_filename = optarg;
                break;
            case REPEAT:
                repeat = atoi(optarg);
                break;
        }
    }

    if(use_input_file) {
        return process_input_file(filename, Kmin);
    }

    complex<double> * S = new Complex[N];

    for(int k = 0; k < repeat; k++) {
        if(!length_set) {
            mpfr_t z;
            mpfr_init2(z, 250);
            mpfr_const_pi(z, GMP_RNDN);
            mpfr_mul_2ui(z, z, 1, GMP_RNDN);
            mpfr_div(z, t, z, GMP_RNDN);
            mpfr_sqrt(z, z, GMP_RNDN);
            mpfr_get_z(length, z, GMP_RNDD);
            // that is going to be the endpoint, really, or the endpoint + 1,
            // so we adjust the length if we are not starting at the endpoint.

            mpz_sub(length, length, start);
            mpz_add_ui(length, length, 1u);

            mpfr_clear(z);
        }


        partial_zeta_sum(start, length, t, delta, N, S, Kmin, number_of_threads, fraction, verbose);

        if(Z_flag) {
            compute_Z_from_rs_sum(t, delta, N, S, S);
        }

        cout << setprecision(17);

        for(int n = 0; n < N; n++) {
            cout << delta * n + k * delta * N << " " << S[n].real() << " " << S[n].imag() << endl;
        }
        mpfr_add_d(t, t, delta * N, GMP_RNDN);
    }

    delete [] S;

    return 0;


    if(argc != 2) {
        usage();
        return 1;
    }

}
