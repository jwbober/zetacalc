typedef double Double;
void create_exp_itlogn_table(mpfr_t t);
Complex exp_itlogn(mpz_t n);

namespace exp_itlogn_stats {
    extern int bigger_than_one;
    extern int smaller_than_one;
};



inline int fastlog2(double x) {
    union {
        double y;
        long long z;
    } a;
    a.y = x;
    (a.z) = ((a.z) >> 52) - 1023;
    
    int d = (int)a.z;

    if(__builtin_expect(d!=0, 1))
        return (int)d;
    else {
        //
        // The above code won't work with denormalized numbers,
        // so we include the following as a fallback in this
        // case. It is still quite fast, I think, and I don't
        // feel like working out this case, since it is somewhat rare.
        //
        int e;
        frexp(x, &e);
        return e - 1;
    }
}



inline int fastlog(double x) {
    //
    // return an approximation to the floor
    // of the log (base e) of x
    //
    // The answer is always exactly equal to the floor of log(x)
    // or else one less than it. In particular it is always <= log(x).
    //
    // It seems to be exact roughly 83.6% of the time, and
    // smaller the other 16.4% of the time.

    int z = fastlog2(x);
    if(z >= 0) {
        z = (z * 726817)/1048576;
    }
    else {
        z = (z * 726817)/1048576 - 1;
    }
    return z;

}
