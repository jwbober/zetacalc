#include <iostream>
#include <cmath>
#include <cstring>
#include "theta_sums.h"
#include "zeta.h"

#include "mpfr.h"
#include "gmp.h"


unsigned int table_precision;

Complex * exp_table;
Complex * exp_table2;

namespace mp_temps {
    mpz_t x, w, z, L;
};

namespace exp_itlogn_stats {
    int bigger_than_one = 0;
    int smaller_than_one = 0;
};

void create_exp_itlogn_table(mpfr_t t) {
    // We want to accurately compute exp(it log n) for
    // large values of n and t. For this, we need
    // to be able to compute t log n mod 2pi to an accuracy
    // of 53 bits, for which we need roughly log_2(t log n)
    // bits. The largest n we will use is sqrt(t/2pi), so
    // we compute everything to a precision of
    // log_2(t) + .5 log_2(log(t/2pi)) + 56 bits
    
    Double tt = mpfr_get_d(t, GMP_RNDN);
    table_precision = (unsigned int)( log2(tt) + .5 * log2(log(tt/(2 * PI)))) + 56;

    // The number of entries in our table will be equal to the precision
    // that we are using. We add 1 to this, and won't use the 0th entry.

    exp_table = new Complex[table_precision + 1];
    exp_table2 = new Complex[table_precision + 1];
    exp_table[0] = -1.0;
    exp_table2[0] = 1.0;

    mpfr_t z, y, x, twopi;
    mpfr_init2(z, table_precision + 2);
    mpfr_init2(y, table_precision + 2);
    mpfr_init2(x, table_precision + 2);
    mpfr_init2(twopi, table_precision + 2);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_set_ui(z, 1, GMP_RNDN);

    for(unsigned int k = 1; k <= table_precision; k++) {
        mpfr_mul_2ui(z, z, 1, GMP_RNDN);    // z = 2^k
        mpfr_sub_ui(y, z, 1, GMP_RNDN);     // y = 2^k - 1
        mpfr_div(y, z, y, GMP_RNDN);        // y = (2^k)/(2^k - 1)
        mpfr_log(y, y, GMP_RNDN);           // y = log( above )
        mpfr_mul(y, y, t, GMP_RNDN);        // y = t log( 2^k / (2^k - 1) )
        mpfr_fmod(y, y, twopi, GMP_RNDN);
        
        exp_table[k] = exp(I * mpfr_get_d(y, GMP_RNDN));
    }

    mpfr_const_log2(z, GMP_RNDN);           // z = log2
    mpfr_mul(z, z, t, GMP_RNDN);            // z = tlog2
    mpfr_set(y, z, GMP_RNDN);               // y = t log 2

    for(unsigned int k = 1; k <= table_precision; k++) {
        mpfr_fmod(x, z, twopi, GMP_RNDN);                   // x = z mod 2pi
        exp_table2[k] = exp(I * mpfr_get_d(x, GMP_RNDN));   // table[k] = exp(I x)
        mpfr_add(z, z, y, GMP_RNDN);                        // z = z + y
    }


    mpz_init2(mp_temps::x, table_precision);
    mpz_init2(mp_temps::z, table_precision);
    mpz_init2(mp_temps::L, table_precision);
    mpz_init2(mp_temps::w, table_precision);

    mpfr_clear(twopi);
    mpfr_clear(z);
    mpfr_clear(y);
    mpfr_clear(x);
}

Complex exp_itlogn3(mpz_t n) {
    // We start by padding n with zeros so that
    // it has exactly table_precision bits.
    unsigned int a = mpz_sizeinbase(n, 2);

    if(a > table_precision) {
        cout << "n was too big." << endl;
        return 0.0/0.0;
    }

    unsigned int N = table_precision;
    int shift = N - a;

    using namespace mp_temps;
    using namespace exp_itlogn_stats;

//    mpz_t x, z, L, w;
//    mpz_init2(x, N);
//    mpz_init2(z, N);
//    mpz_init2(L, N);
//    mpz_init2(w, N);


    mpz_mul_2exp(x, n, shift);
//    mpz_set_ui(L, 1);
//    mpz_mul_2exp(L, L, N - 1);

    Complex y = 1.0;

    unsigned long k = 1u;
    unsigned long current_bit_index = N - (k + 1);
    unsigned long current_limb_index = current_bit_index / GMP_NUMB_BITS;     // even this division shouldn't be necessary.
                                                                              // we can be more clever.
    unsigned long bit_index_in_limb = current_bit_index % GMP_NUMB_BITS;
    mp_srcptr current_limb_ptr = x->_mp_d + current_limb_index;
    mp_limb_t current_limb = *current_limb_ptr;

    unsigned long current_bit_mask = (1u << bit_index_in_limb);

    int index_of_last_limb = current_limb_index;
    int number_of_limbs = index_of_last_limb + 1;
    unsigned long N_minus_1_bit_index = (N - 1) % GMP_NUMB_BITS;
    unsigned long N_minus_1_bit_mask = (1u << N_minus_1_bit_index);

//    cout << "starting" << endl;

    while(1) {
        //while(k <= N - 1 && mpz_tstbit(x, N - (k + 1)) == 0) {
        //    k++;
        //}

        //current_limb_ptr = x->_mp_d + current_limb_index;
        current_limb = *current_limb_ptr;

        //while((k <= N - 1) && (current_limb & current_bit_mask) == 0) {
        while( (current_limb & current_bit_mask) == 0) {
            k++;
            if (__builtin_expect( (current_bit_mask == 1u), 0) ) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_limb_index--;
                current_bit_mask = (1u << (GMP_NUMB_BITS - 1u));
                current_limb = *current_limb_ptr;
            }
            else {
                current_bit_mask = current_bit_mask >> 1u;
            }
//            cout << "Index: " << bit_index_in_limb << endl;
//            cout << "Mask: " << current_bit_mask << endl;
//            cout << "k = " << k << endl;
//            cout << "N = " << N << endl;
//            cout << "current bit should be " << N - (k + 1) << endl;
//            cout << "------" << endl;
        }
        

//        if(k == N) {
//            mpz_clear(x);
//            mpz_clear(z);
//            mpz_clear(L);
//            mpz_clear(w);

//            cout << "here" << endl;
//            return y * exp_table2[a - 1];
//        }

        //mpz_div_2exp(z, x, k);

        //cout << "Want to shift " << x << " by " << k << " to the right." << endl;

        unsigned long number_of_limbs_to_shift = k/GMP_NUMB_BITS;
        unsigned long leftover_shift = k % GMP_NUMB_BITS;
        z->_mp_size = number_of_limbs - number_of_limbs_to_shift;

        //cout << "First shifting " << number_of_limbs_to_shift << " limbs." << endl;
        //cout << "Then shifting an additional " << leftover_shift << " bits."  << endl;

        //cout << "Copying " << sizeof(mp_limb_t) * z->_mp_size << " bytes.";

        //memcpy(z->_mp_d, x->_mp_d + number_of_limbs_to_shift, sizeof(mp_limb_t) * z->_mp_size);                 // note that we are not zeroing any bits of z
                                                                                                                // hopefully this is not necessary, since we
                                                                                                              // reset the size.
        
        //
        // The following seems to be much faster than using memcpy.
        // This is probably because we are only copying a few bytes.
        // It still might be worth investigating ways to possibly do
        // this faster.
        //
        // Also, we don't actually need to copy the full size of z.
        // We know that x is of the form
        //
        //      100..----.001????????...?
        //       <-k zeros-> < N-a bits >
        //
        //  so z will be of the form
        //
        //      00.......00100.......001????..........??
        //      < k zeros > < k zeros > < (N-a-k) bits >
        //
        // We should take advantage of this, but currently we don't.
        //
        mp_limb_t * z_limb_ptr = z->_mp_d;
        mp_limb_t * x_limb_ptr = x->_mp_d + number_of_limbs_to_shift;
        for(int i = 0; i < z->_mp_size; i++) {
            //z->_mp_d[i] = x->_mp_d[number_of_limbs_to_shift + i];
            *z_limb_ptr = *x_limb_ptr;
            z_limb_ptr++;
            x_limb_ptr++;
        }
        //cout << "Now z = " << z << endl;
        if( __builtin_expect(leftover_shift > 0, 1) )
            mpn_rshift(z->_mp_d, z->_mp_d, z->_mp_size, leftover_shift);

        //cout << "Result I get is: " << z << endl;

        mpz_sub(x, x, z);
        //if(mpz_cmp(w, L) >= 1) {
        //if(mpz_tstbit(w, N - 1)) {
        if( __builtin_expect(  __builtin_expect(index_of_last_limb < x->_mp_size, 1) && __builtin_expect( (x->_mp_d[index_of_last_limb] & N_minus_1_bit_mask) != 0, 1), 1) ) {
            //bigger_than_one++; 
            //mpz_swap(x, w);
            
            y = y * exp_table[k];

            k++;
            if (__builtin_expect( (current_bit_mask == 1u), 0) ) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_limb_index--;
                current_bit_mask = (1u << (GMP_NUMB_BITS - 1u));
            }
            else {
                current_bit_mask = current_bit_mask >> 1u;
            }

        }
        else {
            //smaller_than_one++;
            mpz_add(x, x, z);               // We could have done the subtraction x <-- x - z
                                            // before with a temporary variable to avoid this step,
                                            // but this case is very unlikely, and the overhead from
                                            // an mpz_swap call almost every time is larger than
                                            // the overhead from an occasional extra addition.
            //mpz_div_2exp(z, z, 1);
            mpn_rshift(z->_mp_d, z->_mp_d, z->_mp_size, 1);
            k++;
            if( __builtin_expect(current_bit_mask == 1u, 0)) {
                current_limb_ptr--;
                current_bit_mask = (1u << (GMP_NUMB_BITS - 1u));
                current_limb_index--;
            }
            else {
                current_bit_mask = (current_bit_mask >> 1u);
            }
            mpz_sub(x, x, z);
            y = y * exp_table[k];
        }

        //mpz_div_2exp(z, x, k);
    }
}

Complex exp_itlogn2(mpz_t n) {
    // We start by padding n with zeros so that
    // it has exactly table_precision bits.
    unsigned int a = mpz_sizeinbase(n, 2);

    if(a > table_precision) {
        cout << "n was too big." << endl;
        return 0.0/0.0;
    }

    unsigned int N = table_precision;
    int shift = N - a;

    using namespace mp_temps;
    using namespace exp_itlogn_stats;

//    mpz_t x, z, L, w;
//    mpz_init2(x, N);
//    mpz_init2(z, N);
//    mpz_init2(L, N);
//    mpz_init2(w, N);


    mpz_mul_2exp(x, n, shift);
    mpz_set_ui(L, 1);
    mpz_mul_2exp(L, L, N - 1);

    Complex y = 1.0;

    unsigned long k = 1u;
    unsigned long current_bit_index = N - (k + 1);
    unsigned long current_limb_index = current_bit_index / GMP_NUMB_BITS;     // even this division shouldn't be necessary.
                                                                              // we can be more clever.
    unsigned long bit_index_in_limb = current_bit_index % GMP_NUMB_BITS;
    mp_srcptr current_limb_ptr = x->_mp_d + current_limb_index;
    mp_limb_t current_limb = *current_limb_ptr;

    unsigned long current_bit_mask = (1u << bit_index_in_limb);

    int index_of_last_limb = current_limb_index;
    int number_of_limbs = index_of_last_limb + 1;
    unsigned long N_minus_1_bit_index = (N - 1) % GMP_NUMB_BITS;
    unsigned long N_minus_1_bit_mask = (1u << N_minus_1_bit_index);

//    cout << "starting" << endl;

    while(1) {
        //while(k <= N - 1 && mpz_tstbit(x, N - (k + 1)) == 0) {
        //    k++;
        //}

        current_limb_ptr = x->_mp_d + current_limb_index;
        current_limb = *current_limb_ptr;

        while((k <= N - 1) && (current_limb & current_bit_mask) == 0) {
            k++;
            if (__builtin_expect( (current_bit_mask == 1u), 0) ) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_limb_index--;
                current_bit_mask = (1u << (GMP_NUMB_BITS - 1u));
                current_limb = *current_limb_ptr;
            }
            else {
                current_bit_mask = current_bit_mask >> 1u;
            }
//            cout << "Index: " << bit_index_in_limb << endl;
//            cout << "Mask: " << current_bit_mask << endl;
//            cout << "k = " << k << endl;
//            cout << "N = " << N << endl;
//            cout << "current bit should be " << N - (k + 1) << endl;
//            cout << "------" << endl;
        }
        

        if(k == N) {
//            mpz_clear(x);
//            mpz_clear(z);
//            mpz_clear(L);
//            mpz_clear(w);

            return y * exp_table2[a - 1];
        }

        mpz_div_2exp(z, x, k);

        //cout << "Want to shift " << x << " by " << k << " to the right." << endl;

        //unsigned long number_of_limbs_to_shift = k/GMP_NUMB_BITS;
        //unsigned long leftover_shift = k % GMP_NUMB_BITS;
        //z->_mp_size = number_of_limbs - number_of_limbs_to_shift;

        //cout << "First shifting " << number_of_limbs_to_shift << " limbs." << endl;
        //cout << "Then shifting an additional " << leftover_shift << " bits."  << endl;

        //cout << "Copying " << sizeof(mp_limb_t) * z->_mp_size << " bytes.";

        //memcpy(z->_mp_d, x->_mp_d + number_of_limbs_to_shift, sizeof(mp_limb_t) * z->_mp_size);                 // note that we are not zeroing any bits of z
                                                                                                                // hopefully this is not necessary, since we
                                                                                                                    // reset the size.
        //for(int i = 0; i < z->_mp_size; i++) {
        //    z->_mp_d[i] = x->_mp_d[number_of_limbs_to_shift + i];
        //}
        //cout << "Now z = " << z << endl;
        //if(leftover_shift)
        //    mpn_rshift(z->_mp_d, z->_mp_d, z->_mp_size, leftover_shift);

        //cout << "Result I get is: " << z << endl;

        mpz_sub(x, x, z);
        //if(mpz_cmp(w, L) >= 1) {
        //if(mpz_tstbit(w, N - 1)) {
        if( __builtin_expect(  __builtin_expect(index_of_last_limb < x->_mp_size, 1) && __builtin_expect( (x->_mp_d[index_of_last_limb] & N_minus_1_bit_mask) != 0, 1), 1) ) {
            bigger_than_one++; 
            //mpz_swap(x, w);
        }
        else {
            smaller_than_one++;
            mpz_add(x, x, z);               // We could have done the subtraction x <-- x - z
                                            // before with a temporary variable to avoid this step,
                                            // but this case is very unlikely, and the overhead from
                                            // an mpz_swap call almost every time is larger than
                                            // the overhead from an occasional extra addition.
            //mpz_div_2exp(z, z, 1);
            mpn_rshift(z->_mp_d, z->_mp_d, z->_mp_size, 1);
            k++;
            if( __builtin_expect(current_bit_mask == 1u, 0)) {
                current_bit_mask = (1u << (GMP_NUMB_BITS - 1u));
                current_limb_index--;
            }
            else {
                current_bit_mask = (current_bit_mask >> 1u);
            }
            mpz_sub(x, x, z);
        }

        y = y * exp_table[k];
        //mpz_div_2exp(z, x, k);
    }


}

Complex exp_itlogn(mpz_t n) {
    // We start by padding n with zeros so that
    // it has exactly table_precision bits.
    unsigned int a = mpz_sizeinbase(n, 2);
    
    if(a > table_precision) {
        cout << "n was too big." << endl;
        return 0.0/0.0;
    }

    int N = table_precision;
    int shift = N - a;


    mpz_t x, z, L, w;
    mpz_init(x);
    mpz_init(z);
    mpz_init(L);
    mpz_init(w);


    mpz_mul_2exp(x, n, shift);
    mpz_set_ui(L, 1);
    mpz_mul_2exp(L, L, N - 1);

    Complex y = 1.0;

    int k = 1;
    while(1) {
        while(k <= N - 1 && mpz_tstbit(x, N - (k + 1)) == 0) {
            k++;
        }
        if(k == N) {
            mpz_clear(x);
            mpz_clear(z);
            mpz_clear(L);
            mpz_clear(w);

            return y * exp_table2[a - 1];
        }
        mpz_div_2exp(z, x, k);
        mpz_sub(w, x, z);
        //if(mpz_cmp(w, L) >= 1) {
        if(mpz_tstbit(w, N - 1)) {
            mpz_set(x, w);
        }
        else {
            mpz_div_2exp(z, z, 1);
            k++;
            mpz_sub(x, x, z);
        }

        y = y * exp_table[k];
        mpz_div_2exp(z, x, k);
    }
}
