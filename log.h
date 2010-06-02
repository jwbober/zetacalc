typedef double Double;
void create_exp_itlogn_table(mpfr_t t);
Complex exp_itlogn(mpz_t n);
Complex exp_itlogn2(mpz_t n);
Complex exp_itlogn3(mpz_t n);

namespace exp_itlogn_stats {
    extern int bigger_than_one;
    extern int smaller_than_one;
};
