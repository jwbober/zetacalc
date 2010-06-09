#include "theta_sums.h"
#include "zeta.h"

#include "log.h"

using namespace std;

namespace zeta_stats {
    int zeta_block_d = 0;
    int zeta_block_d_using_mpfr = 0;
    int zeta_block_d_using_mpfr_x_large = 0;
};

void print_zeta_stats() {
    cout << "zeta_block_d() called " << zeta_stats::zeta_block_d << " times." << endl;
    cout << "zeta_block_d() used zeta_block_mpfr() " << zeta_stats::zeta_block_d_using_mpfr << " times." << endl;
    cout << "       " << zeta_stats::zeta_block_d_using_mpfr_x_large << " times because K/v was too big." << endl;
}

void stage_1_bound(mpz_t v, mpfr_t t) {
    //
    // Compute the endpoint for the stage 1 sum for computing zeta(.5 + it)
    // and put it into v.
    //

    // Right now this is 3 * t^{1/4}. With this choice, stage 2 will
    // start up when it can use a block size of about 3.

    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));

    mpfr_root(x, t, 4, GMP_RNDN);
    mpfr_mul_ui(x, x, 3u, GMP_RNDN);
    mpfr_get_z(v, x, GMP_RNDN);

    mpfr_clear(x);
}


void stage_2_bound(mpz_t v, mpfr_t t) {
    //
    // Compute the endpoint for the stage 2 sum for computing zeta(.5 + it)
    // and put it into v.
    //

    // Right now this is 112 t^{1/3}. With this choice, we won't
    // start stage 3 until it can use a block size of approximately 1000.
    // (It is possible that we can to make this even larger right now.)

    mpfr_t x;
    mpfr_init2(x, mpfr_get_prec(t));

    mpfr_cbrt(x, t, GMP_RNDN);
    mpfr_mul_ui(x, x, 1120u, GMP_RNDN);

    mpfr_get_z(v, x, GMP_RNDN);

    mpfr_clear(x);
}

void stage_3_bound(mpz_t v, mpfr_t t) {
    //
    // Compute the endpoint for stage 3.
    //
    // This is just floor( sqrt(t/2pi) ) + 1;
    // The "+ 1" is because we stop BEFORE we reach stage_3_bound.
    
    mpfr_t twopi, x;
    mpfr_init2(twopi, mpfr_get_prec(t));
    mpfr_init2(x, mpfr_get_prec(t));
   
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2u, GMP_RNDN);
    
    mpfr_div(x, t, twopi, GMP_RNDN);
    mpfr_sqrt(x, x, GMP_RNDN);
    mpfr_floor(x, x);

    mpfr_get_z(v, x, GMP_RNDN);
    mpz_add_ui(v, v, 1u);

    mpfr_clear(twopi);
    mpfr_clear(x);
}

inline unsigned int stage_2_block_size(Double v, Double t) {
    //
    // For a given v and t determine a good block size to use
    // in stage 2.
    //
    
    // For now we set the block size to be v/t^{1/4}

    unsigned int block_size = (unsigned int)( v * pow(t, -.25) );

    return block_size;
}

inline unsigned int stage_3_block_size(Double v, Double t) {
    // With this choice, the block size will hit N when
    // v = N/.9 * t^(1/3)
    //
    // So if we want to start stage 3 with a block size of 100, then we should
    // set stage_2_bound to 112 * t^(1/3)
    unsigned int block_size = (unsigned int)(.9 * v * pow(t, -.3333333333333333333333));
    return block_size;
}

Complex zeta_block_stage1(mpz_t v, unsigned int K, mpfr_t t) {
    if(K == 0) {
        return 0.0;
    }
    
    mpz_t n;
    mpz_init(n);
        

    if(0) {
        mpz_add_ui(n, v, K - 1);
        cout << v << " " << n << endl;
    }

    Complex S = 0.0;

    mpz_set(n, v);


    for(unsigned int k = 0; k < K; k++) {
        S = S + exp_itlogn4(n)/sqrt(mpz_get_d(n));
        mpz_add_ui(n, n, 1u);
    }

    mpz_clear(n);

    return S;
}

Complex zeta_block_stage2_basic(mpz_t v, unsigned int *K, mpfr_t t, Double epsilon) {
    //
    // This routine calculates the sum
    //
    // sum_{n=v}^{v + K - 1} n^{.5 + it) = sum_{n=v}^{v + K} exp(it log n)/sqrt(n)
    //
    // to a nominal precision of epsilon*K/sqrt(v)
    //
    // Unless *K = 0 or 1, we always compute a good value of *K, and store it
    // in the passed variable. The *K that is passed is the maximum block
    // size that we will use, however.
    //
    // To deal with precision issues in the calculation of the exponential
    // and to get significant speedup, we do a change of variables to
    // write this sum as 
    //
    // exp(it log v)/sqrt(v) sum_{k=0}^{K-1} exp(it log(1 + k/v))/sqrt(1 + k/v).
    //
    // Then we write 1/sqrt(v + k) as exp(log -.5(1 + k/v)),
    // so the sum is
    //
    // exp(i t log(1 + k/v) + .5 log(1 + k/v))
    //
    // We Taylor expand the logarithms with a few terms, and instead
    // of computing each term using mpfr, we just calculate the first
    // terms mod 2pi using mpfr, and then multiply by powers of k
    // in the inner sum using double arithmetic.
    //
    // Let x = K/v.
    // The number of terms we compute in the taylor expansion of the log is
    //
    // -(log(t) - log(epsilon))/log(x)
    //
    // We compute the initial terms using mpfr up to -log(t)/log(x).
    //
    // The number of terms of the in the taylor expansion for the square root
    // term is log(epsilon)/log(x).

    Complex S = 0;

    if(*K == 0) {
        return 0.0;
    }
    if(*K == 1) {
        if(0)
            cout << v << "  " << v << endl;
        return exp_itlogn4(v)/sqrt(mpz_get_d(v));
    }  
    
     Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    
    unsigned int block_size = min(stage_2_block_size(vv, tt), *K);
    if(block_size < 2) {
        cout << "Error: in stage 2, computed a block size of " << block_size << ", which is smaller than 2. Shouldn't reach stage two until we can use a larger block size." << endl;
        cout << "Refusing to continue, and returning NAN, without setting K." << endl;
        return 0.0/0.0;
    }

    *K = block_size;

    if(0) {
        mpz_t n;
        mpz_init(n);
        mpz_add_ui(n, v, block_size - 1);
        cout << v << " " << n << endl;
        mpz_clear(n);
    }

    Double x = block_size/vv;

    // It seems like the following code to estimate the number
    // of terms that we need in the taylor expansion might be
    // useful (and fast) because it doesn't take any logarithms.
    // Instead we just grab the exponent of the relevant numbers
    // to estimate the log base 2. In practice it seems like
    // this slows things down a little bit, though, perhaps
    // because we sometimes (often?) take one extra term
    // in the taylor expansion, so we don't use it right now.

    // UPDATE: After writing the above comment, I changed the
    // way that we decide whether or not to use mpfr to calculate
    // the terms, and now this this method seems about as good,
    // and maybe better, so I am using it again.
    //

    int logepsilon = fastlog2(epsilon);
    int logtt = fastlog2(tt);
    int logx = fastlog2(x);

    int number_of_log_terms = (logepsilon - logtt)/logx;
    int number_of_sqrt_terms = logepsilon/logx;

    Double a[number_of_log_terms + 1];
    Double b[number_of_sqrt_terms + 1];

    int vsize = mpz_sizeinbase(v, 2);
    int precision = mpfr_get_exp(t) - vsize + 53;       
                                                // We want to accurately compute
                                                // the quantities t/v mod pi. To do so
                                                // it should be enough
                                                // to use log2(t) - log2(v) + 53 bits.
    

    mpfr_t mp_v_power, z, twopi, one_over_v, twopi_l, z1;
    
    mpfr_init2(mp_v_power, precision);
    mpfr_init2(one_over_v, precision);
    mpfr_init2(z, precision);
    mpfr_init2(twopi, precision);
    mpfr_init2(twopi_l, precision);
    mpfr_init2(z1, 53);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_set_z(one_over_v, v, GMP_RNDN);
    mpfr_ui_div(one_over_v, 1, one_over_v, GMP_RNDN);

    int sign = 1;
    Double one_over_vv = 1.0/vv;
    Double v_power = one_over_vv;
    mpfr_set(mp_v_power, one_over_v, GMP_RNDN);
    mpfr_set(twopi_l, twopi, GMP_RNDN);
    for(int l = 1; l <= number_of_log_terms; l++) {
        //if(l <= number_of_log_terms_mpfr) {
        if(precision >= 53) {
            //
            // The following calls to reduce the precision on each iteration
            // seem to slow things down just a very little bit in the tests I have run, and maybe
            // they shouldn't be there. But I suspect that in cases
            // where v is very large there should be some gain
            // in using them, so I am leaving them for now.
            //
            mpfr_set_prec(z, precision);
            mpfr_prec_round(mp_v_power, precision, GMP_RNDN);
            mpfr_prec_round(twopi_l, precision, GMP_RNDN);
            mpfr_mul(z, t, mp_v_power, GMP_RNDN);
            
            //mpfr_div_si(z, z, l, GMP_RNDN);
            //mpfr_frac(z, z, GMP_RNDN);
            mpfr_fmod(z1, z, twopi_l, GMP_RNDN);
            a[l] = sign * mpfr_get_d(z1, GMP_RNDN)/l;
            mpfr_mul(mp_v_power, mp_v_power, one_over_v, GMP_RNDN);
            v_power = v_power * one_over_vv;
            mpfr_add(twopi_l, twopi_l, twopi, GMP_RNDN);
            precision = precision - vsize;
        }
        else {
            a[l] = tt * sign * v_power/l;
            v_power = v_power * one_over_vv;
        }
        sign = -sign;
    }


    Double s = 1;
    for(int l = 1; l <= number_of_sqrt_terms; l++) {
        s = s * (-.5/vv);
        b[l] = s/l;
    }

    for(unsigned int k = 0; k < block_size; k++) {
        Double k_power = 1;
        Double x = 0;
        Double y = 0;
        for(int l = 1; l <= number_of_log_terms; l++) {
            k_power = k_power * k;
            y = y + a[l] * k_power;
            if(l <= number_of_sqrt_terms)
                x = x + b[l] * k_power;
        }
        S = S + exp(x + I * y);
    }
    S = S / sqrt(vv);

    S = S * exp_itlogn4(v);

    mpfr_clear(mp_v_power);
    mpfr_clear(z);
    mpfr_clear(twopi);
    mpfr_clear(twopi_l);
    mpfr_clear(one_over_v);
    mpfr_clear(z1);

    return S;
}


Complex zeta_block_stage3_basic(mpz_t v, unsigned int *K, mpfr_t t, Complex ZZ[30], Double epsilon) {
    if(*K == 0) {
        return 0.0;
    }
    if(*K == 1) {
        if(0) {
            cout << v << " " << v << endl;
        }
        return exp_itlogn4(v)/sqrt(mpz_get_d(v));
    }
    
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    
    unsigned int block_size = min(stage_3_block_size(vv, tt), *K);
    if(block_size < 50 && *K >= 50) {
        cout << "Error: in stage 3, computed a block size of " << block_size << ", which is smaller than 50. Shouldn't reach stage two until we can use a larger block size." << endl;
        cout << "Refusing to continue, and returning NAN, without setting K." << endl;
        return 0.0/0.0;
    }

    *K = block_size;
    if(0) {
        mpz_t n;
        mpz_init(n);
        mpz_add_ui(n, v, block_size - 1);
        cout << v << " " << n << endl;
        mpz_clear(n);
    }




    Double w = (block_size-1)/vv;
    Double w_power = 1;

    Complex Z[30];

    for(int l = 0; l < 19; l++) {
        Z[l] = ZZ[l] * w_power;
        w_power *= w;
    }

    int j = 18;

    // Compute Z[l]
 
    mpfr_t a, b, x;

    int precision = mpfr_get_prec(t);

    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(x, precision);

    mpfr_const_pi(x, GMP_RNDN);             // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);         // x = 2 pi
    mpfr_mul_z(x, x, v, GMP_RNDN);            // x = 2 pi v
    mpfr_div(a, t, x, GMP_RNDN);            // a = t / (2 pi v)

//    mpfr_mul_si(x, x, -2, GMP_RNDN);        // x = -4 pi v
//    mpfr_mul(x, x, v, GMP_RNDN);            // x = -4 pi v^2
//    mpfr_div(b, t, x, GMP_RNDN);            // b = -t/ (4 pi v^2)

//    mpfr_mul_si(b, v, -2, GMP_RNDN);
    mpfr_div_z(b, a, v, GMP_RNDN);
    mpfr_div_si(b, b, -2, GMP_RNDN);

    Complex S = compute_exponential_sums(a, b, j, block_size-1, Z, epsilon);

    // we don't need the previous values of a and b anymore, so
    // we can erase them.

    Complex z = exp_itlogn4(v);
    z = z / sqrt(mpz_get_d(v));
    S = S * z;

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(x);


    return S;
}


Complex zeta_block_stage2(mpz_t n, unsigned int N, mpfr_t t) {
    if(N == 0) {
        return 0;
    }

    Complex S = 0;

    mpz_t v;
    mpz_init(v);

    mpz_set(v, n);
    unsigned int K = N;
    while(N > 0) {
        S = S + zeta_block_stage2_basic(v, &K, t, exp(-20));
        N = N - K;
        mpz_add_ui(v, v, K);
        K = N;
    }

    mpz_clear(v);

    return S;
}

Complex zeta_block_stage3(mpz_t n, unsigned int N, mpfr_t t, Complex Z[30]) {
    if(N == 0) {
        return 0;
    }

    Complex S = 0;

    mpz_t v;
    mpz_init(v);

    mpz_set(v, n);
    unsigned int K = N;
    while(N > 0) {
        S = S + zeta_block_stage3_basic(v, &K, t, Z, exp(-20));
        N = N - K;
        mpz_add_ui(v, v, K);
        K = N;
    }

    mpz_clear(v);

    return S;
}




Complex zeta_sum_stage1(mpz_t N, mpfr_t t) {
    //
    // Compute and return the sum
    //
    // \sum_{n=1}^{N - 1} n^{-.5 + it)
    //
    // NOTE: Unlike the other zeta_sum_* functions, N is NOT the
    // number of terms that will be computed. This will compute N - 1
    // terms. (Or, it can be thought of as computing N terms,
    // but with the 0th term being 0).
    //
    // We repeatedly call zeta_block_stage1 with a block size of 10000
    // and add up all of the terms. We could just call zeta_block_stage1
    // directly, but we anticipate shortly changing this routine
    // to use multiple threads.

    create_exp_itlogn_table(t);

    Complex S = 0;

    const unsigned int block_size = 10000;

    mpz_t number_of_terms;
    mpz_init(number_of_terms);
    mpz_sub_ui(number_of_terms, N, 1);

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, number_of_terms, block_size);

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_set_si(v, 1 - block_size);
    for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
        mpz_add_ui(v, v, block_size);
        S = S + zeta_block_stage1(v, block_size, t);
        //S = S + zeta_block_mpfr(v, block_size, t);
    }
    mpz_add_ui(v, v, block_size);

    S = S + zeta_block_stage1(v, remainder, t);

    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);
    mpz_clear(number_of_terms);

    return S;
}

Complex zeta_sum_stage2(mpz_t n, mpz_t N, mpfr_t t) {
    //
    // Compute and return the sum
    //
    // \sum_{k=n}^{n + N - 1} n^{.5 + it}
    //
    // We do this by repeatedly calling zeta_block_stage2 with a block size of 1000000
    // and add up all the terms. We could just call it directly, but we anticipate
    // shortly changing this routine to use multiple threads.
    //
    Complex S = 0;

    const unsigned int block_size = 1000000;

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, N, block_size);

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_set(v, n);
    mpz_sub_ui(v, v, block_size);
    for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
        mpz_add_ui(v, v, block_size);
        S = S + zeta_block_stage2(v, block_size, t);
        //S = S + zeta_block_mpfr(v, block_size, t);
        if(mpz_divisible_ui_p(k, 5u)) {
            cout << "In stage2, completed " << k << " large blocks out of " << number_of_blocks << "." << endl;
        }
    }
    mpz_add_ui(v, v, block_size);

    S = S + zeta_block_stage2(v, remainder, t);
    //S = S + zeta_block_mpfr(v, remainder, t);

    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);

    return S;

}

Complex zeta_sum_stage3(mpz_t n, mpz_t N, mpfr_t t) {
    //
    // Compute and return the sum
    //
    // \sum_{k=n}^{n + N - 1} n^{.5 + it}
    //
    // We do this by repeatedly calling zeta_block_stage3 with a block size of 1000000
    // and add up all the terms. We could just call it directly, but we anticipate
    // shortly changing this routine to use multiple threads.
    //
    Complex S = 0;

    const unsigned int block_size = 1000000;

    mpz_t number_of_blocks;
    mpz_init(number_of_blocks);
    unsigned int remainder = mpz_fdiv_q_ui(number_of_blocks, N, block_size);

    Complex Z[30];
    compute_taylor_coefficients(t, Z);

    mpz_t k, v;
    mpz_init(k);
    mpz_init(v);

    mpz_set(v, n);
    mpz_sub_ui(v, v, block_size);
    for(mpz_set_ui(k, 0u); mpz_cmp(k, number_of_blocks) < 0; mpz_add_ui(k, k, 1u)) {
        mpz_add_ui(v, v, block_size);
        S = S + zeta_block_stage3(v, block_size, t, Z);
        //S = S + zeta_block_mpfr(v, block_size, t);
        if(mpz_divisible_ui_p(k, 5u)) {
            cout << "In stage3, completed " << k << " large blocks out of " << number_of_blocks << "." << endl;
        }
    }
    mpz_add_ui(v, v, block_size);

    S = S + zeta_block_stage3(v, remainder, t, Z);
    //S = S + zeta_block_mpfr(v, remainder, t);

    mpz_clear(v);
    mpz_clear(k);
    mpz_clear(number_of_blocks);

    return S;

}






void compute_taylor_coefficients(mpfr_t t, Complex Z[30]) {
    Double tt1 = mpfr_get_d(t, GMP_RNDN);
    Double tt2 = tt1 * tt1;
    Double tt3 = tt2 * tt1;
    Double tt4 = tt3 * tt1;
    Double tt5 = tt4 * tt1;
    Double tt6 = tt5 * tt1;
    Double tt7 = tt6 * tt1;

    Z[0] = 1.0;
    Z[1] = -.5;
    Z[2] =  .375;
    Z[3] = I * tt1 / 3.0;
    Z[4] = -I * tt1 * 5.0/12.0;
    Z[5] = I * tt1 * 9.0/20.0;
    Z[6] = - tt2 / 18.0;
    Z[7] =  tt2 / 9.0;
    Z[8] = - tt2 * 77.0/480.0;
    Z[9] = -I * tt3 / 162.0;
    Z[10] = I  * tt3 * 11.0/648.0;
    Z[11] = -I * tt3 * 133.0/4320.0;
    Z[12] = tt4 / 1944.0;
    Z[13] = 0.0;
    Z[14] = 0.0;
    Z[15] = tt5 / 29160.0;
    Z[16] = 0.0;
    Z[17] = 0.0;
    Z[18] = tt6 / 524880.0;
    Z[19] = 0;
    Z[20] = 0;
    Z[21] = tt7 / 11022480.0;
    Z[22] = 0;

}


Complex zeta_block(mpz_t v, int K, mpfr_t t, Complex ZZ[30], int method) {
    //
    // Compute sum_{n=v}^{v+K} exp(it log n)/sqrt(n)
    //
    // We assume that K is significantly smaller than v -- something
    // like K/V = t^{1/3}. We do a taylor expansion on exp(i t log n)
    // using just two terms, which is where this assumption is used.
    //

    if(method == 2) {
        Complex S = zeta_block_mpfr(v, K, t);
        return S;
    }

    Double vv = mpz_get_d(v);

    Double w = (K-1)/vv;
    Double w_power = 1;

    Complex Z[30];

    for(int l = 0; l < 30; l++) {
        Z[l] = ZZ[l] * w_power;
        w_power *= w;
    }

    int j = 18;

    // Compute Z[l]
 
    mpfr_t a, b, x;

    int precision = mpfr_get_prec(t);

    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(x, precision);

    mpfr_const_pi(x, GMP_RNDN);             // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);         // x = 2 pi
    mpfr_mul_z(x, x, v, GMP_RNDN);            // x = 2 pi v
    mpfr_div(a, t, x, GMP_RNDN);            // a = t / (2 pi v)

//    mpfr_mul_si(x, x, -2, GMP_RNDN);        // x = -4 pi v
//    mpfr_mul(x, x, v, GMP_RNDN);            // x = -4 pi v^2
//    mpfr_div(b, t, x, GMP_RNDN);            // b = -t/ (4 pi v^2)

//    mpfr_mul_si(b, v, -2, GMP_RNDN);
    mpfr_div_z(b, a, v, GMP_RNDN);
    mpfr_div_si(b, b, -2, GMP_RNDN);

    Complex S = compute_exponential_sums(a, b, j, K-1, Z, exp(-14));

    // we don't need the previous values of a and b anymore, so
    // we can erase them.

    mpfr_set_z(x, v, GMP_RNDN);
    mpfr_log(x, x, GMP_RNDN);               // x = log v
    mpfr_mul(x, x, t, GMP_RNDN);            // x = t log v

    mpfr_const_pi(a, GMP_RNDN);
    mpfr_mul_si(a, a, 2, GMP_RNDN);
    mpfr_fmod(x, x, a, GMP_RNDN);
    
    Complex z = exp(I * mpfr_get_d(x, GMP_RNDN));
    z = z / sqrt(mpz_get_d(v));
    S = S * z;

//    mpfr_sin_cos(b, a, x, GMP_RNDN);        // a + ib = exp(i t log v)
//    mpfr_sqrt(x, v, GMP_RNDN);              // x = sqrt(v)
    
//    mpfr_div(a, a, x, GMP_RNDN);
//    mpfr_div(b, b, x, GMP_RNDN);            // a + ib is now (exp(i t log v)/sqrt(v))

//    Complex z(mpfr_get_d(a, GMP_RNDN), mpfr_get_d(b, GMP_RNDN));
//    S = S * z;

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(x);

    if(verbose::zeta_block) {
        cout << "zeta block returning  " << S << endl;
        //cout << "using mpfr, answer is " << zeta_block_mpfr(v, K, t) << endl;
        cout << "   v = " << v << endl;
        cout << "   K = " << K << endl;
    }

    return S;
}

Complex zeta_block_mpfr(mpfr_t v, unsigned int K, mpfr_t t) {
    mpz_t vv;
    mpz_init(vv);

    mpfr_get_z(vv, v, GMP_RNDN);

    Complex S = zeta_block_mpfr(vv, K, t);

    mpz_clear(vv);
    return S;
}

Complex zeta_block_mpfr(mpz_t v, unsigned int K, mpfr_t t) {
    //
    // Compute the sum_{k=v}^{v + K - 1) n^{-.5 + it}
    //
    // We use machine doubles for the terms in the sum, but we use
    // mpfr to accurately calculate the quantity t log n mod 2 pi
    // for each term in the sum.
 
    //
    // First we figure out how much precision we will need from mpfr.
    //
    // We want to accurately calculate t log n mod 2pi to 53 bits, which
    // means that we need to compute t log n to 53 + log_2(t log n) bits.
    //
    // For safety we add an extra 2 bits.
    //

    mpfr_t x, y, n, twopi;
    mpfr_init2(x, 53);
    mpfr_init2(y, 53);

    mpfr_log2(x, t, GMP_RNDN);                          // Now x = log2(t)
    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log(y, y, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // And y = log2(log n) (We calculate these quantities
                                                        // to low precision because we are only interested
                                                        // in their size, really. There is probably 
                                                        // a clever way to do faster using less precision.

    mpfr_add(x, x, y, GMP_RNDN);
    int mod_precision = mpfr_get_ui(x, GMP_RNDN) + 55;  // This is the precision that we need when
                                                        // need to calculate the quantity t log n
                                                        // when we mod by 2 pi

    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // y = log2(v + K) now.
    int n_precision = mpfr_get_ui(y, GMP_RNDN) + 2;     // This is the precision that we need to exactly
                                                        // represent the largest integer that will occur
                                                        // in the summation.
 
    mpfr_init2(n, n_precision);

    mpfr_init2(twopi, mod_precision);
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_clear(x);
    mpfr_init2(x, mod_precision);

    Complex S = 0.0;
    
    mpfr_set_z(n, v, GMP_RNDN);             // The summation starts with n = v
    for(unsigned int k = 0; k <= K-1; k++) {
        mpfr_log(x, n, GMP_RNDN);           // x = log(n)
        mpfr_mul(x, x, t, GMP_RNDN);        // a = t log n
        mpfr_fmod(x, x, twopi, GMP_RNDN);
        Complex z = exp(I * mpfr_get_d(x, GMP_RNDN));
        z = z/sqrt(mpfr_get_d(n, GMP_RNDN));
        S = S + z;
    
        mpfr_add_ui(n, n, 1, GMP_RNDN);     // now on the next iteration, n will be v + k + 1
    }

    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(n);
    mpfr_clear(twopi);

    return S;
}

Complex zeta_block_d(mpz_t v, int K, mpfr_t t, Double epsilon) {
    //
    // This routine calculates the sum
    //
    // sum_{n=v}^{v + K - 1} n^{.5 + it) = sum_{n=v}^{v + K} exp(it log n)/sqrt(n)
    //
    // to a nominal precision of epsilon*K/sqrt(v)
    // 
    // To deal with precision issues in the calculation of the exponential
    // and to get significant speedup, we do a change of variables to
    // write this sum as 
    //
    // exp(it log v)/sqrt(v) sum_{k=0}^{K-1} exp(it log(1 + k/v))/sqrt(1 + k/v).
    //
    // Then we write 1/sqrt(v + k) as exp(log -.5(1 + k/v)),
    // so the sum is
    //
    // exp(i t log(1 + k/v) + .5 log(1 + k/v))
    //
    // We Taylor expand the logarithms with a few terms, and instead
    // of computing each term using mpfr, we just calculate the first
    // terms mod 2pi using mpfr, and then multiply by powers of k
    // in the inner sum using double arithmetic.
    //
    // Let x = K/v.
    // The number of terms we compute in the taylor expansion of the log is
    //
    // -(log(t) - log(epsilon))/log(x)
    //
    // We compute the initial terms using mpfr up to -log(t)/log(x).
    //
    // The number of terms of the in the taylor expansion for the square root
    // term is log(epsilon)/log(x).

    zeta_stats::zeta_block_d++;

    Complex S = 0;

    if(K == 0) {
        return S;
    }
    if(K <= 1) {
        zeta_stats::zeta_block_d_using_mpfr++;
        return zeta_block_mpfr(v, K, t);
    }

    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);
    Double x = K/vv;

    if(x > pow(tt, -1.0/6.0)) {     // TODO: maybe think of a better way to approximate this
                                    // condition so that we don't call pow.
        zeta_stats::zeta_block_d_using_mpfr_x_large++;
        zeta_stats::zeta_block_d_using_mpfr++;
        return zeta_block_mpfr(v, K, t);
    }


    // It seems like the following code to estimate the number
    // of terms that we need in the taylor expansion might be
    // useful (and fast) because it doesn't take any logarithms.
    // Instead we just grab the exponent of the relevant numbers
    // to estimate the log base 2. In practice it seems like
    // this slows things down a little bit, though, perhaps
    // because we sometimes (often?) take one extra term
    // in the taylor expansion, so we don't use it right now.

    // UPDATE: After writing the above comment, I changed the
    // way that we decide whether or not to use mpfr to calculate
    // the terms, and now this this method seems about as good,
    // and maybe better, so I am using it again.
    //

     int logepsilon;
     int logtt;
     int logx;

     frexp(epsilon, &logepsilon);
     frexp(tt, &logtt);
     frexp(x, &logx);

     int number_of_log_terms = (logepsilon - logtt)/logx;
//     int number_of_log_terms_mpfr = -logtt/logx + 1;             // We add one instead of taking the ceiling
     int number_of_sqrt_terms = logepsilon/logx;

//    int number_of_log_terms = (int)( (log(epsilon) - log(tt))/log(x));
//    int number_of_log_terms_mpfr = (int)ceil(-log(tt)/log(x));
//    int number_of_sqrt_terms = (int)( log(epsilon)/log(x) );

    //number_of_log_terms = 20;
    //number_of_log_terms_mpfr = number_of_log_terms;
    //number_of_sqrt_terms = 20;

    if(verbose::zeta_block_d) {
        cout << "zeta_block_d() called with " << endl;
        cout << "                          v = " << vv << endl;
        cout << "                          K = " << K << endl;
        cout << "                          t = " << tt << endl;
        cout << "                        K/v = " << x << endl;

        cout << "   Number of terms in log taylor expansion is " << number_of_log_terms << endl;
//        cout << "                  Number of terms using mpfr: " << number_of_log_terms_mpfr << endl;
        cout << "   Number of terms in sqrt taylor expansion is " << number_of_sqrt_terms << endl;

        cout << "   Using old code, would have used: " << endl;
        cout << "   Number of terms in log taylor expansion is " << (int)( (log(epsilon) - log(tt))/log(x)) << endl;
//        cout << "                  Number of terms using mpfr: " << (int)ceil(-log(tt)/log(x)) << endl;
        cout << "   Number of terms in sqrt taylor expansion is " <<  (int)( log(epsilon)/log(x) )  << endl;
    }

    Double aa = 53 + ceil((-log(tt)/log(vv))) * log2(K) + 100;      // TODO: think about whether we need to do this.
    if(aa < -log2(epsilon * sqrt(vv)/K)) {
        if(verbose::zeta_block_d) {
            cout << " --Evaluating directly using mpfr" << endl;
            cout << "    aa was " << aa << endl;
        }
        zeta_stats::zeta_block_d_using_mpfr++;
        return zeta_block_mpfr(v, K, t);
    }


    Double a[number_of_log_terms + 1];
    Double b[number_of_sqrt_terms + 1];

    int vsize = mpz_sizeinbase(v, 2);
    int precision = mpfr_get_exp(t) - vsize + 53;       
                                                // We want to accurately compute
                                                // the quantities t/v mod pi. To do so
                                                // it should be enough
                                                // to use log2(t) - log2(v) + 53 bits.
    
    int step_two_precision = mpfr_get_exp(t) + vsize + 53;                  // Later we will want to accurately
                                                                            // compute t log v mod 2pi, so
                                                                            // we will need more precision.
                                                                            // We are overdoing it here a little
                                                                            // bit, and really should just use a
                                                                            // precision equal to log2(t) + log2(log(v)) + 53.
                                                                            //
                                                                            // We will also calculate twopi to this
                                                                            // precision, since we will need it.

    mpfr_t mp_v_power, z, twopi, one_over_v, twopi_l, z1;
    
    mpfr_init2(mp_v_power, precision);
    mpfr_init2(one_over_v, precision);
    mpfr_init2(z, precision);
    mpfr_init2(twopi, step_two_precision);
    mpfr_init2(twopi_l, precision);
    mpfr_init2(z1, 53);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_set_z(one_over_v, v, GMP_RNDN);
    mpfr_ui_div(one_over_v, 1, one_over_v, GMP_RNDN);

    int sign = 1;
    Double one_over_vv = 1.0/vv;
    Double v_power = one_over_vv;
    mpfr_set(mp_v_power, one_over_v, GMP_RNDN);
    mpfr_set(twopi_l, twopi, GMP_RNDN);
    for(int l = 1; l <= number_of_log_terms; l++) {
        //if(l <= number_of_log_terms_mpfr) {
        if(precision >= 53) {
            //
            // The following calls to reduce the precision on each iteration
            // seem to slow things down just a very little bit in the tests I have run, and maybe
            // they shouldn't be there. But I suspect that in cases
            // where v is very large there should be some gain
            // in using them, so I am leaving them for now.
            //
            mpfr_set_prec(z, precision);
            mpfr_prec_round(mp_v_power, precision, GMP_RNDN);
            mpfr_prec_round(twopi_l, precision, GMP_RNDN);
            mpfr_mul(z, t, mp_v_power, GMP_RNDN);
            
            //mpfr_div_si(z, z, l, GMP_RNDN);
            //mpfr_frac(z, z, GMP_RNDN);
            mpfr_fmod(z1, z, twopi_l, GMP_RNDN);
            a[l] = sign * mpfr_get_d(z1, GMP_RNDN)/l;
            mpfr_mul(mp_v_power, mp_v_power, one_over_v, GMP_RNDN);
            v_power = v_power * one_over_vv;
            mpfr_add(twopi_l, twopi_l, twopi, GMP_RNDN);
            precision = precision - vsize;
        }
        else {
            a[l] = tt * sign * v_power/l;
            v_power = v_power * one_over_vv;
        }
        sign = -sign;
    }


    Double s = 1;
    for(int l = 1; l <= number_of_sqrt_terms; l++) {
        s = s * (-.5/vv);
        b[l] = s/l;
    }

    for(int k = 0; k <= K-1; k++) {
        Double k_power = 1;
        Double x = 0;
        Double y = 0;
        for(int l = 1; l <= number_of_log_terms; l++) {
            k_power = k_power * k;
            y = y + a[l] * k_power;
            if(l <= number_of_sqrt_terms)
                x = x + b[l] * k_power;
        }
        S = S + exp(x + I * y);
    }
    S = S / sqrt(vv);

    // Now we need z to have more precision. Enough to store the integer
    // part of t log v + 53 bits

    mpfr_set_prec(z, step_two_precision);

    mpfr_set_z(z, v, GMP_RNDN);
    mpfr_log(z, z, GMP_RNDN);
    mpfr_mul(z, z, t, GMP_RNDN);

    mpfr_fmod(z1, z, twopi, GMP_RNDN); 
    S = S * exp(I * mpfr_get_d(z1, GMP_RNDN));


    if(verbose::zeta_block_d >= 2) {
        cout << "Computed zeta_block_d = " << S << endl;
        Complex z = zeta_block_mpfr(v, K, t);
    }

    if(0) {
        Complex z1 = zeta_block_mpfr(v, K, t);
        Double logerror = log(abs(z1 - S));
        if(logerror > -19) {
            cout << "zeta_block_d() called with " << endl;
            cout << "                          v = " << vv << endl;
            cout << "                          K = " << K << endl;
            cout << "                          t = " << tt << endl;
            cout << "                        K/v = " << x << endl;

            cout << "   Number of terms in log taylor expansion is " << number_of_log_terms << endl;
//            cout << "                  Number of terms using mpfr: " << number_of_log_terms_mpfr << endl;
            cout << "   Number of terms in sqrt taylor expansion is " << number_of_sqrt_terms << endl;
            cout << "Computed zeta_block_d = " << S << endl;
            cout << "      Answer should be: " << z1 << endl;
            cout << "   log of difference is " << logerror << endl;
            cout << "  The logarithm taylor coefficients we computed were: " << endl;
            for(int n = 1; n <= number_of_log_terms; n++) {
                cout << a[n] << "  ";
            }
            cout << endl;
            cout << " Now computing them again, all with mpfr... I get:" << endl;
            mpfr_set_si(mp_v_power, 1, GMP_RNDN);
            int sign = 1;
            for(int l = 1; l <= number_of_log_terms; l++) {
                mpfr_mul_z(mp_v_power, mp_v_power, v, GMP_RNDN);
                mpfr_div(z, t, mp_v_power, GMP_RNDN);
                mpfr_div_si(z, z, l, GMP_RNDN);
                mpfr_fmod(z, z, twopi, GMP_RNDN);
                a[l] = sign * mpfr_get_d(z, GMP_RNDN);
                sign = -sign;
            }
            for(int n = 1; n <= number_of_log_terms; n++) {
                cout << a[n] << "  ";
            }

            cout << endl;
        }
    }

    mpfr_clear(mp_v_power);
    mpfr_clear(z);
    mpfr_clear(twopi);
    mpfr_clear(twopi_l);
    mpfr_clear(one_over_v);
    mpfr_clear(z1);

    return S;
}

Complex zeta_block_d_stupid(mpz_t v, int K, mpfr_t t) {
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);

    Complex S = 0;
    for(int l = 0; l <= K-1; l++) {
        Double n = vv + l;
        S = S + pow(n, -.5 + I * tt);
    }
    return S;
}

Complex initial_zeta_sum_mpfr(mpz_t M, mpfr_t t) {
    mpfr_t x, y;
//    mpfr_init2(x, 53);
    mpfr_init2(y, 53);
//    mpfr_log2(x, t, GMP_RNDN);
    mpfr_set_z(y, M, GMP_RNDN);
    mpfr_log(y, y, GMP_RNDN);
//    mpfr_log2(y, y, GMP_RNDN);
//    mpfr_add(x, x, y, GMP_RNDN);

    int mod_precision = mpfr_get_exp(t) + mpfr_get_exp(y) + 55;
//    int mod_precision = mpfr_get_ui(x, GMP_RNDN) + 55;  // This is the precision that we need when
                                                        // need to calculate the quantity t log n
                                                        // when we mod by 2 pi

//    mpfr_set_z(y, M, GMP_RNDN);
//    mpfr_log2(y, y, GMP_RNDN);
    int n_precision = mpz_sizeinbase(M, 2) + 2;         // This is the precision that we need to exactly
                                                        // represent the largest integer that will occur
                                                        // in the summation.
    
    if(verbose::initial_zeta_sum_mpfr) {
        cout << "In initial_zeta_sum_mpfr using " << mod_precision << " bits of precision for computation." << endl;
    }

    //int precision = mpfr_get_prec(t);
    
    mpfr_t nn, twopi;
    mpfr_init2(nn, n_precision);
    mpfr_init2(twopi, mod_precision);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

 //   mpfr_clear(x);
    mpfr_init2(x, mod_precision);

    mpz_t n;
    mpz_init(n);

    Complex S1 = 0;

    for(mpz_set_si(n, 1); mpz_cmp(n, M) <= 0; mpz_add_ui(n, n, 1)) {
        mpfr_set_z(nn, n, GMP_RNDN);
        mpfr_log(x, nn, GMP_RNDN);           // x = log(n)
        mpfr_mul(x, x, t, GMP_RNDN);         //  x = t log n

        mpfr_fmod(y, x, twopi, GMP_RNDN);

        Double z = mpfr_get_d(y, GMP_RNDN);
        S1 = S1 + exp(I * z)/sqrt(mpz_get_d(n));

//        unsigned int nn = mpfr_get_ui(n, GMP_RNDN);
//        if(nn % 10000 == 0) {
//            cout << nn << endl;
//        }

//        mpfr_sin_cos(b, a, a, GMP_RNDN);    // a + ib = exp(i t log n)
//        mpfr_sqrt(c, n, GMP_RNDN);          // c = sqrt(n)
//        mpfr_div(a, a, c, GMP_RNDN);
//        mpfr_div(b, b, c, GMP_RNDN);
//        mpfr_add(real_part, real_part, a, GMP_RNDN);
//        mpfr_add(imaginary_part, imaginary_part, b, GMP_RNDN);
    }
    
    mpz_clear(n);
    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(twopi);

    return S1;
}

inline void printmp(mpfr_t x) {
    cout << mpfr_get_d(x, GMP_RNDN) << endl;
}

Complex initial_zeta_sum(mpz_t M, mpfr_t t, Double epsilon) {
    Complex S, S1, S2, S3, S4;
    S = 0.0;
    S1 = 0.0;
    S2 = 0.0;
    S3 = 0.0;
    S4 = 0.0;
    
    unsigned int m = 6;
    unsigned int mm = pow(2, m);
    unsigned int K = 0;
    mpz_t M2, M3, r, v, R;
    mpfr_t x;

    mpz_init(M2);
    mpz_init(M3);
    mpz_init(r);
    mpz_init(v);
    mpz_init(R);

    mpfr_init2(x, mpfr_get_prec(t));            //
    mpfr_cbrt(x, t, GMP_RNDN);                  //
    mpfr_div_si(x, x, mm, GMP_RNDN);            //
    mpfr_get_z(M2, x, GMP_RNDD);                // M2 = floor( t^(1/3)/2^m )


    mpz_mul_ui(M3, M2, mm);
    mpz_sub(M3, M, M3);                         // M3 = M - 2^m M2

    mpz_add_ui(R, M3, 1);                       // R = M3 + 1
    mpz_div_ui(R, R, mm);                       // R = floor((M3 + 1)/2^m)

    if(verbose::initial_zeta_sum) {
        cout << "In initial_zeta_sum(): " << endl;
        cout << "        t = " << mpfr_get_d(t, GMP_RNDN) << endl;
        cout << "        M = " << M << endl;
        cout << "       M2 = " << M2 << endl;
        cout << "       M3 = " << M3 << endl;
        cout << "        R = " << R << endl;
    }


    mpz_sub_ui(r, M2, 1);                       //
    S1 = initial_zeta_sum_mpfr(r, t);           // We first compute the sum up to M2 - 1 using mpfr




    // Next we compute the sum in block sizes of length
    // K, where K starts at 1 and goes up by powers of
    // 2 until K = 2^(m-1)
    //
    // (Really we should start at K = 2 instead of K = 1, for efficiency)
    
    S2 = 0;
    K = 1;
    mpz_set(v, M2);
    if(verbose::initial_zeta_sum) {
        cout << "   In initial_zeta_sum(), to start with, v = " << v << endl;
    }
    for(unsigned int l = 0; l <= m-1; l++) {
//        mpz_mul_ui(v, M2, K);        // v = M2 * 2^l
        for(mpz_set_ui(r, 0); mpz_cmp(r, M2) < 0; mpz_add_ui(r, r, 1)) {
            S2 = S2 + zeta_block_d(v, K, t, epsilon);
            //S2 = S2 + zeta_block_mpfr(v, K, t);
            mpz_add_ui(v, v, K);    // v = (M2 + r + 1) * 2^l
        }
        K = K * 2;
    }

    if(verbose::initial_zeta_sum) {
        cout << "                     After first loop,   v = " << v << endl;
        cout << "                                         K = " << K << endl;
    }

    // Now we compute some remainder terms with block size mm = 2^m

    S3 = 0;
    //K = mm;
    //mpz_mul_ui(v, M2, K);            // v = M2 * 2^m to begin with
    for(mpz_set_ui(r, 0); mpz_cmp(r, R) < 0; mpz_add_ui(r, r, 1)) {
        S3 = S3 + zeta_block_d(v, K, t, epsilon);
        //S3 = S3 + zeta_block_mpfr(v, K, t);
        mpz_add_ui(v, v, K);        // now v = (M2 + r + 1) * 2^m, so it will have the correct value on the next iteration of the loop
    }

    if(verbose::initial_zeta_sum) {
        cout << "                    After second loop,   v = " << v << endl;
    }

    // Finally we compute one more block, from v = (M2 + R)2^m
    // (which already has the correct value), to the end, which
    // is M, unless we have already computed the whole sum correctly.

    mpz_sub(M3, M, v);
    mpz_add_ui(M3, M3, 1);      // to include M, the number of terms we have to compute is M - v + 1
    K = mpz_get_ui(M3);

    if(verbose::initial_zeta_sum) {
        cout << "                        for final block, K = " << K << endl;
    }

    S4 = zeta_block_d(v, K, t, epsilon);
    //S4 = zeta_block_mpfr(v, K, t);

    mpz_clear(M2);
    mpz_clear(M3);
    mpz_clear(r);
    mpz_clear(v);
    mpz_clear(R);
    mpfr_clear(x);

    S = S1 + S2 + S3 + S4;

    return S;
}

Complex zeta_sum_basic(mpfr_t t) {
    mpfr_t x;
    mpz_t z;

    mpfr_init2(x, mpfr_get_prec(t));
    mpz_init(z);

    mpfr_const_pi(x, GMP_RNDN);                 // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);             // x = 2 pi
    mpfr_div(x, t, x, GMP_RNDN);                // x = t/2pi
    mpfr_sqrt(x, x, GMP_RNDN);                  // x = sqrt(t/2pi)
    mpfr_floor(x, x);                           // x = floor(sqrt(t/2pi))

    mpfr_get_z(z, x, GMP_RNDN);

    Complex S = initial_zeta_sum(z, t, exp(-20));

    mpfr_clear(x);
    mpz_clear(z);
    return S;

}

Complex zeta_sum_mpfr(mpfr_t t) {
    mpfr_t x;
    mpz_t z;

    mpfr_init2(x, mpfr_get_prec(t));
    mpz_init(z);

    mpfr_const_pi(x, GMP_RNDN);                 // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);             // x = 2 pi
    mpfr_div(x, t, x, GMP_RNDN);                // x = t/2pi
    mpfr_sqrt(x, x, GMP_RNDN);                  // x = sqrt(t/2pi)
    mpfr_floor(x, x);                           // x = floor(sqrt(t/2pi))

    mpfr_get_z(z, x, GMP_RNDN);

    Complex S = initial_zeta_sum_mpfr(z, t);

    mpfr_clear(x);
    mpz_clear(z);
    return S;
}

Complex zeta_sum(mpfr_t t) {
    Complex S = 0.0;

    mpz_t n1, n2, n3, N2, N3;

    mpz_init(n1);
    mpz_init(n2);
    mpz_init(n3);
    mpz_init(N2);
    mpz_init(N3);

    stage_1_bound(n1, t);
    stage_2_bound(n2, t);
    stage_3_bound(n3, t);

    Complex S1, S2, S3;

    if(mpz_cmp(n2, n3) > 0) {      // if n2 > n3, set n2 = n3
        mpz_set(n2, n3);
    }
    if(mpz_cmp(n1, n3) > 0) {      // if n1 > n3, set n1 = n3
        mpz_set(n1, n3);
    }

    cout << "For t = " << mpfr_get_d(t, GMP_RNDN) << ": " << endl;
    cout << "   Using stage1 up to n = " << n1 << endl;
    cout << "   Using stage2 up to n = " << n2 << endl;
    cout << "   Using stage3 up to n = " << n3 << endl;

    mpz_sub(N2, n2, n1);        // N2 and N3 hold the lengths of the sums for stage 2 and stage 3
    mpz_sub(N3, n3, n2);

    S1 = zeta_sum_stage1(n1, t);
    cout << "Done with stage 1. Sum was: " << S1 << endl;

    S2 = zeta_sum_stage2(n1, N2, t);
    cout << "Done with stage 2. Sum was: " << S2 << endl;
    
    S3 = zeta_sum_stage3(n2, N3, t);
    cout << "Done with stage 3. Sum was: " << S3 << endl;

    mpz_clear(n1);
    mpz_clear(n2);
    mpz_clear(n3);
    mpz_clear(N2);
    mpz_clear(N3);

    S = S1 + S2 + S3;

    return S;
}

Double C(int n, Double *powers_z);

Double rs_remainder(mpfr_t t) {
    //computes 3 remainder term in riemann-siegel formula
    //can be extended to more than 3 terms easily, if needed
    //everything can be safely done in doubles unless otherwise is stated
    //we're using the notation of the ``amortized complexity..." paper
    //
    //This is adapted from Michael Rubinstein's code for the riemann-siegel formula,
    //the lcalc library
    //

    Double remainderTerms;
    Double c0, c1, c2;
    Double z;
    Double *powers_z;
      
    //variables and calculations below should be done in MPFR
    mpfr_t mp_a, mp_p, twopi;
    mpz_t n1;
    mpfr_init2(mp_a, mpfr_get_prec(t));     // TODO: intelligently choose precision?
    mpfr_init2(twopi, mpfr_get_prec(t));     // TODO: intelligently choose precision?
    mpfr_init2(mp_p, 53);
    mpz_init(n1);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);

    mpfr_div(mp_a, t, twopi, GMP_RNDN);         // mp_a = t/twopi
    mpfr_sqrt(mp_a, mp_a, GMP_RNDN);            // mp_a = sqrt(t/twopi)

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    mpfr_get_z(n1, mp_a, GMP_RNDD);                 // n1 = floor(a)

    //a = sqrt(t / (2 * Pi));
    //n1 = floor(a);

    mpfr_frac(mp_p, mp_a, GMP_RNDN);            // mp_p = frac(a)
    //p = a - n1;

    Double p = mpfr_get_d(mp_p, GMP_RNDN);

    //z can be calculated as a double
    z=p-.5;

    //precomputed arrays that probably save time

    powers_z = new Double[51];
    powers_z[0]=1;
    for(int l=1; l<=50; l++)
        powers_z[l] = powers_z[l-1] * z;

    //remainderTerms C(0,p), C(1,p), C(2,p)

    c0=C(0,powers_z);
    c1=C(1,powers_z) * pow(a, -1);
    c2=C(2,powers_z) * pow(a, -2);

    //remainderTerms = C(0,p) + C(1,p)*pow(a,-1) + C(2,p)*pow(a,-2) 

    remainderTerms = c0 + c1 + c2;

    //remainderTerms = remainderTerms * (-1)^(n1+1) / sqrt(a)

    mpz_mod_ui(n1, n1, 2);
    if(mpz_cmp_ui(n1, 1) == 0)
        remainderTerms = remainderTerms * pow(a, -.5);
    else
        remainderTerms = -remainderTerms * pow(a, -.5);

    delete [] powers_z;
    mpz_clear(n1);
    mpfr_clear(twopi);
    mpfr_clear(mp_a);
    mpfr_clear(mp_p);

    return remainderTerms;

}

Double C(int n, Double *powers_z){
    if (n==0) return
        .3826834323650897717284599840304*powers_z[0]
        +1.74896187231008179744118586948533*powers_z[2]
        +2.11802520768549637318456427826417*powers_z[4]
        -.87072166705114807391892407738239*powers_z[6]
        -3.4733112243465167073064116693758*powers_z[8]
        -1.66269473089993244964313630119286*powers_z[10]
        +1.21673128891923213447689352804417*powers_z[12]
        +1.30143041610079757730060538099786*powers_z[14]
        +.03051102182736167242108987123981*powers_z[16]
        -.37558030515450952427981932122934*powers_z[18]
        -.10857844165640659743546975901329*powers_z[20]
        +.05183290299954962337576051067322*powers_z[22]
        +.02999948061990227592040084956912*powers_z[24]
        -.00227593967061256422601994851021*powers_z[26]
        -.00438264741658033830594007013585*powers_z[28]
        -.00040642301837298469930723272116*powers_z[30]
        +.00040060977854221139278910314608*powers_z[32]
        +8.97105799138884129783418195378689e-05*powers_z[34]
        -2.30256500272391071161029452573899e-05*powers_z[36]
        -9.38000660190679248471972940127474e-06*powers_z[38]
        +6.32351494760910750424986123959430e-07*powers_z[40]
        +6.55102281923150166621223123133411e-07*powers_z[42]
        +2.21052374555269725866086890381876e-08*powers_z[44]
        -3.32231617644562883503133517017624e-08*powers_z[46]
        -3.73491098993365608176460476015222e-09*powers_z[48]
        +1.24450670607977391951510000249366e-09*powers_z[50];
    if (n==1) return
        -.05365020525675069405998280791133*powers_z[1]
        +.11027818741081482439896362071917*powers_z[3]
        +1.23172001543152263131956529162206*powers_z[5]
        +1.26349648627994578841755482191213*powers_z[7]
        -1.69510899755950301844944739906596*powers_z[9]
        -2.99987119676501008895548735894141*powers_z[11]
        -.10819944959899208642692257787438*powers_z[13]
        +1.94076629462127126879387632539716*powers_z[15]
        +.78384235615006865328843457488694*powers_z[17]
        -.50548296679003659187902141326173*powers_z[19]
        -.3845072349605797405134273885311*powers_z[21]
        +.03747264646531532067594447494023*powers_z[23]
        +.09092026610973176317258142450576*powers_z[25]
        +.01044923755006450921820113972659*powers_z[27]
        -.01258297965158341649747892224592*powers_z[29]
        -.00339950372115127408505894886137*powers_z[31]
        +.00104109505377148912682954240655*powers_z[33]
        +.00050109490511184868603556526727*powers_z[35]
        -3.95635966900318155954711855696337e-05*powers_z[37]
        -4.76245924535718963865409830268035e-05*powers_z[39]
        -1.85393553380851322734349064569117e-06*powers_z[41]
        +3.19369180800689720404663539343268e-06*powers_z[43]
        +4.09078076085060663265089453677018e-07*powers_z[45]
        -1.54466243325766321284375723273104e-07*powers_z[47]
        -3.46630749176913317222559405934073e-08*powers_z[49];
    if (n==2) return
        .00518854283029316849378458151923*powers_z[0]
        +.00123786335522538984133826974438*powers_z[2]
        -.18137505725166997411491896409414*powers_z[4]
        +.14291492748532126541165603376514*powers_z[6]
        +1.33033917666875653250993329998546*powers_z[8]
        +.35224723534037336775327655505836*powers_z[10]
        -2.4210015958919507237815305433405*powers_z[12]
        -1.67607870225381088533346181492372*powers_z[14]
        +1.36894167233283721842349153807076*powers_z[16]
        +1.55390194302229832214563952655935*powers_z[18]
        -.17221642734729980519582586998918*powers_z[20]
        -.63590680550454309889704902355845*powers_z[22]
        -.09911649873041208105423564341370*powers_z[24]
        +.14033480067387008950738254898316*powers_z[26]
        +.04782352019827292236438803506512*powers_z[28]
        -.01735604064147978079795864709223*powers_z[30]
        -.01022501253402859184447660413126*powers_z[32]
        +.00092741491597948878994270014371*powers_z[34]
        +.00135721943723733853452533619958*powers_z[36]
        +6.41369012029388008996238736394533e-05*powers_z[38]
        -.00012300805698196629883342322937*powers_z[40]
        -1.83135074047892025547675543979621e-05*powers_z[42]
        +7.82162860432262730850139938461872e-06*powers_z[44]
        +2.00875424847599455034985293919157e-06*powers_z[46]
        -3.35327653931857137372749727241453e-07*powers_z[48]
        -1.46160209174182309264510097122760e-07*powers_z[50];
    if (n==3) return
        -.00267943218143891380853967145989*powers_z[1]
        +.02995372109103514963731329491570*powers_z[3]
        -.04257017254182869798501935111688*powers_z[5]
        -.28997965779803887506893209478669*powers_z[7]
        +.48888319992354459725374746407169*powers_z[9]
        +1.23085587639574608119312504336294*powers_z[11]
        -.82975607085274087041796910432976*powers_z[13]
        -2.24976353666656686652045012659903*powers_z[15]
        +.07845139961005471379365473620184*powers_z[17]
        +1.74674928008688940039198666645219*powers_z[19]
        +.45968080979749935109237306173169*powers_z[21]
        -.66193534710397749464339040008983*powers_z[23]
        -.31590441036173634578979632973316*powers_z[25]
        +.12844792545207495988511847476209*powers_z[27]
        +.10073382716626152300969450207513*powers_z[29]
        -.00953018384882526775950465984230*powers_z[31]
        -.01926442168751408889840098069714*powers_z[33]
        -.00124646371587692917124790716458*powers_z[35]
        +.00242439696411030857397215245841*powers_z[37]
        +.00043764769774185701827561290396*powers_z[39]
        -.00020714032687001791275913078304*powers_z[41]
        -6.27434450418651556052610958029804e-05*powers_z[43]
        +1.15753438145956693483789208989316e-05*powers_z[45]
        +5.88385492454037978388597885697078e-06*powers_z[47]
        -3.12467740069633622086961449076033e-07*powers_z[49];
    if (n==4) return
        .00046483389361763381853630462560*powers_z[0]
        -.00402264294613618830391153989145*powers_z[2]
        +.00384717705179612688359130685272*powers_z[4]
        +.06581175135809486002088309200741*powers_z[6]
        -.19604124343694449117695528448205*powers_z[8]
        -.20854053686358853244400012794494*powers_z[10]
        +.95077541851417509458477574151058*powers_z[12]
        +.53415353129148739760517592459894*powers_z[14]
        -1.67634944117634007959116448203404*powers_z[16]
        -1.07674715787512899278784663510432*powers_z[18]
        +1.23533930165659698528788361189251*powers_z[20]
        +1.02578253400572757718348949577914*powers_z[22]
        -.40124095793988544378728137523313*powers_z[24]
        -.50366639951083034479585257591604*powers_z[26]
        +.03573487795502744985807080163387*powers_z[28]
        +.14431763086785416624285239495844*powers_z[30]
        +.01509152741790346941712677290432*powers_z[32]
        -.02609887477919436131761773965448*powers_z[34]
        -.00612662837951926174904909908948*powers_z[36]
        +.00307750312987084118476787782167*powers_z[38]
        +.00115624789340887523161201204220*powers_z[40]
        -.00022775966758472127472807733953*powers_z[42]
        -.00014189637118181444432681579894*powers_z[44]
        +7.46486030795591945312240984450313e-06*powers_z[46]
        +1.24797016454091166174449988846871e-05*powers_z[48]
        +4.86394518400209461907998084746180e-07*powers_z[50];
    if (n==5) return
        .00022686811845737363176557957245*powers_z[1]
        +.00110812468537183880897586725284*powers_z[3]
        -.01621857925555009106408484258686*powers_z[5]
        +.05276503405398741662724126665649*powers_z[7]
        +.02570880200903323999290010111095*powers_z[9]
        -.38058660440806397264435991848146*powers_z[11]
        +.22531987892642315322976926989838*powers_z[13]
        +1.03445733164952217211304499657389*powers_z[15]
        -.55282576970508137898888475296735*powers_z[17]
        -1.52877126410780729962736571714169*powers_z[19]
        +.32828366427719583672031669394059*powers_z[21]
        +1.22911021854008706238425001239677*powers_z[23]
        +.04093693938311529830689289790902*powers_z[25]
        -.55860404726420193442735876775644*powers_z[27]
        -.11241976368059115396788439789609*powers_z[29]
        +.15212677711795591829295940144809*powers_z[31]
        +.05173718845528038784023625510664*powers_z[33]
        -.02561227689700728294043343196050*powers_z[35]
        -.01296367251404617794428713962277*powers_z[37]
        +.00254555748186116327806192744188*powers_z[39]
        +.00211933195108777752885073213414*powers_z[41]
        -9.19139194515677754051761292342159e-05*powers_z[43]
        -.00024413466533855272657049552509*powers_z[45]
        -1.36979826922833871237722416376381e-05*powers_z[47]
        +2.06207850332842380642607669485206e-05*powers_z[49];

    else return 0;
}


Complex rs_rotation(mpfr_t t) {
    //returns rotation factor exp(-i theta(t)) in riemann-siegel formula
    //
    //every thing here is done in doubles unless otherwise is stated
    //
    //

    Double theta;
     
    //first two terms of asymptotic expansion for theta(t)
    //should be done in MPFR, temp should be declated as MPFR variable

    mpfr_t temp, twopi;
    mpfr_init2(twopi, mpfr_get_prec(t));        // TODO: intelligently choose precision?
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_init2(temp, mpfr_get_prec(t));         // TODO: intelligently choose precision?
    
    mpfr_div(temp, t, twopi, GMP_RNDN);         // temp = t/2pi
    mpfr_log(temp, temp, GMP_RNDN);             // temp = log(t/2pi)
    mpfr_sub_ui(temp, temp, 1, GMP_RNDN);       // temp = log(t/2pi) - 1
    mpfr_mul(temp, temp, t, GMP_RNDN);          // temp = t log(t/2pi) - t
    mpfr_div_ui(temp, temp, 2, GMP_RNDN);       // temp = t/2 log(t/2pi) - t/2
    mpfr_fmod(temp, temp, twopi, GMP_RNDN);     // temp = (above) mod 2pi

//    temp =  (t / 2.) * log(t / (2. * pi)) - t / 2.
//    temp = fmod(theta, 2 * pi);

    //a few more terms in asymptotic expansion of theta(t)
    //can be safely done in doubles
    
    Double tt = mpfr_get_d(t, GMP_RNDN);

    theta = mpfr_get_d(temp, GMP_RNDN) - PI/8.0 +  1.0/(tt * 48.0) + 7.0/(tt*tt*tt * 5760.0); 

    //calculate exp(-I*theta(t))

    Complex answer = exp(-I * theta);

    mpfr_clear(twopi);
    mpfr_clear(temp);

    return answer;
}

Complex hardy_Z(mpfr_t t, Complex &R) {
    //
    // Compute the hardy_Z function at t and set R to the rotation
    // factor so that zeta(.5 + it) = R * Z(t)
    //
    Complex main_sum;
    Complex rotation_factor;
    Double remainder_terms;

    main_sum = zeta_sum(t);

    rotation_factor = rs_rotation(t);
    R = rotation_factor;

    remainder_terms = rs_remainder(t);

    if(verbose::hardy_Z) {
        cout << "In hardy_Z(): "<< endl;
        cout << "       computed main_sum = " << main_sum;
        cout << "         rotation_factor = " << rotation_factor;
        cout << "         remainder_terms = " << remainder_terms;
    }

    Double Z_value = 2 * real(rotation_factor * main_sum) + remainder_terms;

    if(verbose::hardy_Z) {
        cout << "Computed Z(" << mpfr_get_d(t, GMP_RNDN) << ") = " << Z_value << endl;
    }

    return Z_value;
}

Complex zeta(mpfr_t t) {
    Complex rotation_factor;
    Complex w = hardy_Z(t, rotation_factor);
    return w * rotation_factor;
}
