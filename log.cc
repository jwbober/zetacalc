//
// Functions to quickly compute exp(it log n) for fixed t and varying n,
// after a short precomputation. These functions are much faster than using
// mpfr for the same calculation, but not quite as fast as using double precision.
// (Using double precision for these calculations wouldn't give accurate answers, though.)
//
// To use:
//
//      First the function create_exp_itlogn_table(mpfr_t t) should be called with the value
//      of t that we are interested in. After this, the function exp_itlogn(mpz_t n) returns
//      exp(it log n) as a complex<double> to almost 53 bits of precision.
//
//      (More precisely, the function exp_itlogn_table() creates a table of values of
//      exp(i t log(2^k/(2^k - 1))) to full double precision. To compute exp(it log n),
//      the function exp_itlogn() ends up multiplying some of these together, making no
//      attempt to control the small roundoff errors that may occur.)
//
//

// The algorithm used is essentially the one described in exercises 25 and 28 of Section 1.2.2
// of Knuths TAOCP Third Edition.

#include <iostream>
#include <cmath>
#include <cstring>
#include "theta_sums.h"
#include "zeta.h"

#include "mpfr.h"
#include "gmp.h"


unsigned int table_precision;
mp_size_t number_of_limbs_needed;

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


    number_of_limbs_needed = table_precision/GMP_NUMB_BITS + 1;     // There may be one extra limb needed since
                                                                    // the division truncates.

    mpz_init2(mp_temps::x, table_precision + 5);
    mpz_init2(mp_temps::z, table_precision + 5);
    mpz_init2(mp_temps::L, table_precision + 5);
    mpz_init2(mp_temps::w, table_precision + 5);

    mpz_set_si(mp_temps::x, 2);
    mpz_mul_2exp(mp_temps::x, mp_temps::x, table_precision + 4);
    mpz_set(mp_temps::x, mp_temps::x);
    mpz_set(mp_temps::z, mp_temps::x);
    mpz_set(mp_temps::L, mp_temps::x);
    mpz_set(mp_temps::w, mp_temps::x);

    mpfr_clear(twopi);
    mpfr_clear(z);
    mpfr_clear(y);
    mpfr_clear(x);
}

inline void shift_left(mp_ptr output, mp_srcptr input, mp_size_t input_size, unsigned int shift) {
    //
    // leftshift the input by shift bits and put the result in output. We assume that
    // output has enough space to store the result, and that the low order bits of output,
    // which will not be touched, are already 0. (or whatever the caller wishes)
    //

    // We start by shifting limbs to the left, which only
    // gets us a shift which is a multiple of GMP_NUMB_BITS.
    // Afterwards we will shift the rest of the way using mpn_lshift

    mp_size_t number_of_limbs_to_shift = shift/GMP_NUMB_BITS;
    int leftover_shift = shift % GMP_NUMB_BITS;

    //
    // The following seems to be much faster than using memcpy.
    // This is probably because we are only copying a few bytes.
    // It still might be worth investigating ways to possibly do
    // this faster.
    //


    mp_srcptr input_limb_ptr = input;
    mp_ptr output_limb_ptr = output + number_of_limbs_to_shift;
    for(int i = 0; i < input_size; i++) {
        //z->_mp_d[i] = x->_mp_d[number_of_limbs_to_shift + i];
        *output_limb_ptr = *input_limb_ptr;
        output_limb_ptr++;
        input_limb_ptr++;
    }

    mp_limb_t shifted_out_bits = 0;

    // Now we do the leftover shift.
    if( __builtin_expect(leftover_shift > 0, 1) )
        shifted_out_bits = mpn_lshift(output + number_of_limbs_to_shift, output + number_of_limbs_to_shift, input_size, leftover_shift);

    *(output + number_of_limbs_to_shift + input_size) = shifted_out_bits;

    // We don't zero out the low order bits of the output because
    // in the cases that this function will be called, we know that
    // the output is already 0.

    //for(int i = 0; i < number_of_limbs_to_shift; i++) {
    //    *output = 0;
    //    output++;
    //}
}

inline void shift_right(mp_ptr output, mp_srcptr input, mp_size_t input_size, unsigned int shift) {
    //
    // shift the input right by shift bits and put the result in output. We assume that
    // output has enough space to store the result.

    // We start by shifting limbs to the left, which only
    // gets us a shift which is a multiple of GMP_NUMB_BITS.
    // Afterwards we will shift the rest of the way using mpn_lshift

    mp_size_t number_of_limbs_to_shift = shift/GMP_NUMB_BITS;
    int leftover_shift = shift % GMP_NUMB_BITS;

    //
    // The following seems to be much faster than using memcpy.
    // This is probably because we are only copying a few bytes.
    // It still might be worth investigating ways to possibly do
    // this faster.
    //


    mp_srcptr input_limb_ptr = input + number_of_limbs_to_shift;
    mp_ptr output_limb_ptr = output;
    for(int i = 0; i < input_size - number_of_limbs_to_shift; i++) {
        //z->_mp_d[i] = x->_mp_d[number_of_limbs_to_shift + i];
        *output_limb_ptr = *input_limb_ptr;
        output_limb_ptr++;
        input_limb_ptr++;
    }

    // Now we do the leftover shift.
    if( __builtin_expect(leftover_shift > 0, 1) )
        mpn_rshift(output, output, input_size - number_of_limbs_to_shift, leftover_shift);

    // For now we zero out the high order bits of the output, assuming that it has
    // the same size as the input. This shouldn't be necessary if we do things correctly.
 
    output = output + input_size - number_of_limbs_to_shift;
    for(int i = 0; i < number_of_limbs_to_shift; i++) {
        *output = 0;
        output++;
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

    unsigned int N = table_precision;
    int shift = N - a;

    using namespace mp_temps;
    using namespace exp_itlogn_stats;

//    mpz_t x, z, L, w;
//    mpz_init2(x, N);
//    mpz_init2(z, N);
//    mpz_init2(L, N);
//    mpz_init2(w, N);


//    mpz_mul_2exp(x, n, shift);

    mp_limb_t x[number_of_limbs_needed];
    mp_limb_t w[number_of_limbs_needed];
    mp_limb_t z[number_of_limbs_needed];

    for(mp_size_t k = 0; k < number_of_limbs_needed; k++) {
        x[k] = 0;
        w[k] = 0;
        z[k] = 0;
    }

    shift_left(x, n->_mp_d, n->_mp_size, shift);

//    mpz_set_ui(L, 1);
//    mpz_mul_2exp(L, L, N - 1);

    Complex y = 1.0;

    unsigned long k = 1u;
    unsigned long current_bit_index = N - (k + 1);
    mp_size_t current_limb_index = current_bit_index / GMP_NUMB_BITS;     // even this division shouldn't be necessary.
                                                                              // we can be more clever.
    unsigned long bit_index_in_limb = current_bit_index % GMP_NUMB_BITS;
    mp_srcptr current_limb_ptr = x + current_limb_index;
    mp_limb_t current_limb = *current_limb_ptr;

    mp_limb_t current_bit_mask = (1ul << bit_index_in_limb);

    mp_size_t index_of_last_limb = current_limb_index;
    mp_size_t number_of_limbs = index_of_last_limb + 1;
    unsigned long N_minus_1_bit_index = (N - 1) % GMP_NUMB_BITS;
    mp_limb_t N_minus_1_bit_mask = (1ul << N_minus_1_bit_index);

//    cout << "starting" << endl;

    while(1) {
        current_limb = *current_limb_ptr;

        //while((k <= N - 1) && (current_limb & current_bit_mask) == 0) {
        while( (current_limb & current_bit_mask) == 0) {
            k++;
            if (__builtin_expect( (current_bit_mask == 1ul), 0) ) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_limb_index--;
                current_bit_mask = (1ul << (GMP_NUMB_BITS - 1));
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

        shift_right(z, x, number_of_limbs_needed, k);

        //mpz_sub(x, x, z);
        mpn_sub_n(x, x, z, number_of_limbs_needed);

        //if(mpz_cmp(w, L) >= 1) {
        //if(mpz_tstbit(w, N - 1)) {
        if( __builtin_expect(  __builtin_expect(index_of_last_limb < number_of_limbs_needed, 1) && __builtin_expect( (x[index_of_last_limb] & N_minus_1_bit_mask) != 0, 1), 1) ) {
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
                current_bit_mask = (1ul << (GMP_NUMB_BITS - 1u));
            }
            else {
                current_bit_mask = current_bit_mask >> 1u;
            }

        }
        else {
            //smaller_than_one++;
            mpn_add_n(x, x, z, number_of_limbs_needed);             // We could have done the subtraction x <-- x - z
                                                                    // before with a temporary variable to avoid this step,
                                                                    // but this case is very unlikely, and the overhead from
                                                                    // an mpz_swap call almost every time is larger than
                                                                    // the overhead from an occasional extra addition.
            //mpz_div_2exp(z, z, 1);
            mpn_rshift(z, z, number_of_limbs_needed, 1);
            k++;
            if( __builtin_expect(current_bit_mask == 1u, 0)) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_bit_mask = (1ul << (GMP_NUMB_BITS - 1u));
                current_limb_index--;
            }
            else {
                current_bit_mask = (current_bit_mask >> 1u);
            }
            //mpz_sub(x, x, z);
            mpn_sub_n(x, x, z, number_of_limbs_needed);
            y = y * exp_table[k];
        }

        //mpz_div_2exp(z, x, k);
    }
}
