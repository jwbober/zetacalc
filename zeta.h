namespace verbose {
    const int zeta_block = 0;
    const int zeta_block_d = 0;
    const int initial_zeta_sum_mpfr = 1;
};


Complex zeta_block(mpfr_t v, int K, mpfr_t t, Complex Z[13], int method = 0);
Complex zeta_sum(mpfr_t t);
void compute_taylor_coefficients(mpfr_t t, Complex Z[13]);


Complex zeta_block_d(mpfr_t v, int K, mpfr_t t, Double epsilon);
Complex zeta_block_d_stupid(mpfr_t v, int K, mpfr_t t);
Complex zeta_block_mpfr(mpfr_t v, int K, mpfr_t t);

Complex initial_zeta_sum(mpfr_t M, mpfr_t t, Double epsilon);
Complex initial_zeta_sum_mpfr(mpz_t M, mpfr_t t);

void print_zeta_stats();

Complex zeta_sum(mpfr_t);
Complex zeta_sum_mpfr(mpfr_t);
Complex zeta_sum_basic(mpfr_t);
