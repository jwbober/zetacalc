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
// NOTE ON MULTITHREADING:
//      It is safe to call exp_itlogn() simultaneously from multiple threads, but
//      the table used is a global variable, so there can only be one t active at any time.
//      In our typical use case, this one thread calls create_exp_itlogn_table() and then
//      many threads call exp_itlogn(). This means that our current implementation
//      supports using many threads to compute a single value of the zeta function,
//      but it does not support computing different values of the zeta function
//      simultaneously.
//

// The algorithm used is essentially the one described in exercises 25 and 28 of Section 1.2.2
// of Knuth's TAOCP Third Edition. It is described further below in the source code.

#include <iostream>
#include <cmath>
#include <cstring>
#include "theta_sums.h"
#include "zeta.h"

#include "mpfr.h"
#include "gmp.h"


unsigned int table_precision;
mp_size_t number_of_limbs_needed;

Complex * exp_table;        // Table to store the values exp(i t log (2^k/(2^k - 1)))
Complex * exp_table2;       // Table to store the values exp(i t log(2^k))

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
    // Compute exp(it log n) to (almost) 53 bits of precision, using whatever
    // t was used as an argument to create_exp_itlogn_table().
    //
    // Conceptually, at least, the algorithm works by first dividing n by a large enough
    // power of 2 to get some x = n/2^(a-1) so that we have 1 <= x < 2. (Thus we will
    // have exp(it log n) = exp(it log 2^(a-1)) * exp(it log x). In fact, a is just
    // the number of binary digits of n.)
    //
    // To compute exp(it log x) we repeatedly divide x by the largest number
    // of the form 2^k/(2^k - 1) which is smaller than x. What we have in the end
    // is
    //
    // exp(it log n) = exp(it log 2^(a-1)) * prod_{some k, possibly with repitition} exp(it log(2^k/(2^k - 1)))
    //
    // Each of these individual quantities is computed to full 53 bit precision in the
    // function create_exp_itlogn_table(), and we ignore any (small) errors that
    // appear when we multiply them together, which is why we say that this
    // is accurate to (almost) 53 bis of precision.
    //
    // One reason that this algorithm work nicely is that dividing x by 2^k/(2^k - 1) is the same as
    // subtracting x/2^k from x, so we can do this division with a shift and a subtraction.

    unsigned int a = mpz_sizeinbase(n, 2);

    if(a > table_precision) {
        cout << "n was too big." << endl;
        return 0.0/0.0;
    }

    unsigned int N = table_precision;
    int shift = N - a;
    
    unsigned long k = 1u;
    mp_limb_t x[number_of_limbs_needed];        // x will be as described above.
    mp_limb_t z[number_of_limbs_needed];        // z is always going to be x/2^k

    for(mp_size_t k = 0; k < number_of_limbs_needed; k++) {
        x[k] = 0;
        z[k] = 0;
    }

    // We start by placing the bits of n into the highest order bits of x,
    // and then we will interpret x as a rational number between 1 and 2.
    //
    // The input, n, has a digits base 2. x will have N digits base 2.
    shift_left(x, n->_mp_d, n->_mp_size, shift);

    // y is an accumulator where we to the product of exp(it log(2^k/(2^k - 1))
    // as we compute it.
    Complex y = 1.0;


    // To check if x > 2^k/(2^k - 1) we will mostly use bit operations. Suppose
    // that we already know that x < 2^n/(2^n - 1) for all n < k. Then, since
    // 2^k/(2^k - 1) = 1 + 2^(-k) + 2^(-2k) + ... we can first examine
    // the kth most significant bit of x to see if it is 1. If not, then
    // x is smaller than 2^k/(2^k - 1). (Otherwise we will need to do a
    // little more work.) This gets a little confusing because the kth
    // most significant bit of x, starting with k = 0, is the (N-(k+1))st bit
    // of x. (All of our bit numberings start with 0.)
    //
    // Thus, we have the following variables to help keep
    // track of these bit operations.

    unsigned long first_bit_index = N - 2;                                  // The index of the 1st most significant bit of x. (Not the 0th).
    mp_size_t first_limb_index = first_bit_index / GMP_NUMB_BITS;           // The index of the limb that this bit is in.
    unsigned long bit_index_in_limb = first_bit_index % GMP_NUMB_BITS;      // The index of the 1st most significant bit inside the most significant limb.

    mp_srcptr current_limb_ptr = x + first_limb_index;                      // Pointer to the limb that contains the kth most significant bit.
    mp_limb_t current_limb = *current_limb_ptr;                             // The actual limb.

    mp_limb_t current_bit_mask = (1ul << bit_index_in_limb);                // The bit mask that we will use to check if the kth most significant bit is 1.

    mp_srcptr most_significant_limb_ptr = x + number_of_limbs_needed - 1;
    unsigned long N_minus_1_bit_index = (N - 1) % GMP_NUMB_BITS;            // To check of x < 2^k/(2^k - 1) we will actually subtract x/2^k from x
    mp_limb_t N_minus_1_bit_mask = (1ul << N_minus_1_bit_index);            // and then simply check if the most significant bit of x is still 1.
                                                                            // This is the bit mask that we will use to check this.


    while(1) {
        // We start by finding the most significant bit of x (after the 0th) that is 1.
        // At this point it should be guaranteed that bits 1 through k-1 (numbered from most significant)
        // are all 0.

        current_limb = *current_limb_ptr;
        while( (current_limb & current_bit_mask) == 0) {
            k++;
            if (__builtin_expect( (current_bit_mask == 1ul), 0) ) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_bit_mask = (1ul << (GMP_NUMB_BITS - 1));
                current_limb = *current_limb_ptr;
            }
            else {
                current_bit_mask = current_bit_mask >> 1u;
            }
        }
        
        // Now we set z to x/2^k and subtract it from x. We don't actually
        // know yet whether the result will be >= 1, but it turns out that
        // it usually is, so it is faster to do the subtraction in place
        // than to use a temporary variable.

        shift_right(z, x, number_of_limbs_needed, k);
        mpn_sub_n(x, x, z, number_of_limbs_needed);

        // Here we actually check if the result is >= 1. Since we are testing >= (instead
        // of just >) we can just check the most significant bit of x and see if it is 1).
        if(__builtin_expect( (*most_significant_limb_ptr & N_minus_1_bit_mask) != 0, 1)) {
            // If we are here, then the result of the subtraction was >= 1
            // so we want to include the factor 2^k/(2^k - 1)
            
            y = y * exp_table[k];

            // In this case, it is impossible for the factor 2^k/(2^k - 1) to be repeated
            // again, so we advance k, and check for the terminating condition.
            k++;
            if (__builtin_expect( (current_bit_mask == 1u), 0) ) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_bit_mask = (1ul << (GMP_NUMB_BITS - 1u));
            }
            else {
                current_bit_mask = current_bit_mask >> 1u;
            }

        }
        else {
            // If we are here, then the result of the subtraction was < 1,
            // so we do not want to include the factor 2^k/(2^k - 1). However
            // if this happens then we know that we want to include the factor
            // 2^(k + 1)/(2^(k + 1) - 1), so we restore x, increase k and shift
            // z one more to the right and subtract it.

            mpn_add_n(x, x, z, number_of_limbs_needed);             // As mentioned above, we could have done the subtraction x <-- x - z
                                                                    // before with a temporary variable to avoid this step,
                                                                    // but this case is very unlikely, and the overhead from
                                                                    // an mpz_swap call almost every time is larger than
                                                                    // the overhead from an occasional extra addition.
            mpn_rshift(z, z, number_of_limbs_needed, 1);
            k++;
            if( __builtin_expect(current_bit_mask == 1u, 0)) {
                if(k == N) {
                    return y * exp_table2[a - 1];
                }
                current_limb_ptr--;
                current_bit_mask = (1ul << (GMP_NUMB_BITS - 1u));
            }
            else {
                current_bit_mask = (current_bit_mask >> 1u);
            }
            mpn_sub_n(x, x, z, number_of_limbs_needed);
            y = y * exp_table[k];
        }
    }
}
