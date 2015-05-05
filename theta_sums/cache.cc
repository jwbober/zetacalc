#include "theta_sums.h"
#include <iomanip>
using namespace std;



theta_cache * build_theta_cache(mpfr_t mp_a, mpfr_t mp_b, int j, int K) {
    //
    // At the beginning of a run of the algorithm we
    // create a cache of a bunch of things that will
    // be needed in many places.
    //
    //
    // This function builds that cache. To use the cache, call
    // the utility functions defined below.
    //
    // The cache itself is of type void*, and it should be thought
    // of as a big binary blob with some structure.
    //
 
    // The beginning of the cache is just a theta_cache struct,
    // but we allocate extra space after the end of it to store
    // extra information.
    //
    // Right now, it should look like
    //                                                      starts at
    //      theta_struct                                    0
    //      (Double) K^-j                                   sizeof(theta_struct)
    //      (Double) K^(-j + 1)
    //      ...
    //      (Double) K^0                                    sizeof(theta_struct) + j * sizeof(Double)
    //      (Double) K^1
    //      ...
    //      (Double) K^j                                    sizeof(theta_struct) + 2j * sizeof(Double)
    //      (Double) (2 PI b)^(-j - 1)/2                    sizeof(theta_struct) + (2j + 1) * sizeof(Double)
    //      (Double) (2 PI b)^-j/2
    //      (Double) (2 PI b)^(-j + 1)/2
    //      ...
    //      (Double) (2 PI b)^0                             sizeof(theta_struct) + (3j + 2) * sizeof(Double)
    //      ...
    //      (Double) (2 PI b)^(j + 1)/2                     sizeof(theta_struct) + (4j + 3) * sizeof(Double)

    // TOTAL SIZE: sizeof(theta_struct) + (2j + 1) sizeof(Double) + (2j + 3) sizeof(Double)

    theta_cache * cache = (theta_cache*)malloc(sizeof(theta_cache) + (2 * j + 1) * sizeof(Double) + (2 * j + 3) * sizeof(Double));
 
    cache->a = mpfr_get_d(mp_a, GMP_RNDN);
    cache->b = mpfr_get_d(mp_b, GMP_RNDN);
    cache->K = K;
    cache->j = j;
    cache->q = int(cache->a + 2 * cache->b * K);
    cache->ExpAK = ExpAK(mp_a, K);
    cache->ExpBK = ExpBK(mp_b, K);
    cache->ExpAK_inverse = conj(cache->ExpAK);
    cache->ExpBK_inverse = conj(cache->ExpBK);
    cache->ExpAB = ExpAB(mp_a, mp_b);
    cache->ExpABK = cache->ExpAK * cache->ExpBK;

    cache->C1 = I * cache->ExpABK;
    cache->C5 = -cache->C1;
    cache->C7 = -cache->C5;
    cache->C8 = -I * cache->ExpBK_inverse;

    Double * K_powers = (Double *)((intptr_t)(cache) + sizeof(theta_cache));
    
    K_powers[0] = pow(K, -j);
    for(int k = 1; k <= 2 * j; k++) {
        K_powers[k] = K_powers[k-1] * K;
    }

    Double * root_2pi_b_powers = (Double *)((intptr_t)(cache) + sizeof(theta_cache) + (2 * j + 1) * sizeof(Double));

    Double root_2pi_b = sqrt(2 * PI * cache->b);
    root_2pi_b_powers[0] = pow(root_2pi_b, -j - 1);
    for(int k = 1; k <= 2 * j + 2; k++) {
        root_2pi_b_powers[k] = root_2pi_b_powers[k - 1] * root_2pi_b;
    }

    // Aside from j, there are only 4 different possible inputs to H_integral_0,
    // and H_integral_0 is likely called with the same inputs over and over again
    //
    // So we should cache these values here.

    return cache;
}

void free_theta_cache(theta_cache * cache) {
    free(cache);
}


