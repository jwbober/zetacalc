namespace verbose {
    const int zeta_sum = 1;
    const int zeta_block = 0;
    const int zeta_block_d = 0;
    const int initial_zeta_sum_mpfr = 1;
    const int initial_zeta_sum = 1;
    const int hardy_Z = 1;
};

Complex zeta_block_d_stupid(mpz_t v, int K, mpfr_t t);
Complex zeta_block_mpfr(mpz_t v, unsigned int K, mpfr_t t);

Complex initial_zeta_sum(mpz_t M, mpfr_t t, Double epsilon);
Complex initial_zeta_sum_mpfr(mpz_t M, mpfr_t t);

void print_zeta_stats();

Complex zeta_sum2(mpfr_t, Double delta, int N, Complex * S);
Complex zeta_sum_mpfr(mpfr_t);
Complex zeta_sum_basic(mpfr_t);

Complex zeta(mpfr_t t);

Complex rs_rotation(mpfr_t t);
Double rs_remainder(mpfr_t t);

void compute_Z_from_rs_sum(mpfr_t t0, double delta, int N, Complex * S, Complex * Z);

/************************************************************************/

// It's reasonable to remove all the functions above the mark from this file since we're 
// only concerned with the _main sum_ here, and not yet with computing zeta. 


Complex partial_zeta_sum(mpz_t start, mpz_t length, mpfr_t t, double delta, int M, Complex * S, int Kmin, int number_of_threads = 2, int fraction = -1, int verbose=1);

Complex zeta_block_stage1(mpz_t v, unsigned int K, mpfr_t t, Double delta, int M, Complex * S);
Complex zeta_block_stage2(mpz_t n, unsigned int N, mpfr_t t, Double delta, int M, Complex * S);
Complex zeta_block_stage3(mpz_t n, unsigned int N, mpfr_t t, Complex Z[30], Double delta, int M, Complex * S, int Kmin = 0);

Complex zeta_block_stage2_basic(mpz_t v, unsigned int *K, mpfr_t t, Double epsilon);
Complex zeta_block_stage3_basic(mpz_t v, unsigned int *K, mpfr_t t, Complex ZZ[30], Double epsilon, int Kmin = 0);

void stage_2_start(mpz_t v, mpfr_t t);
void stage_3_start(mpz_t v, mpfr_t t);

unsigned int stage_2_block_size(Double v, Double t);
unsigned int stage_3_block_size(Double v, Double t);

void compute_taylor_coefficients(mpfr_t t, Complex Z[30]);
