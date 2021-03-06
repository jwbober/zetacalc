import sys

ranges = {}

range_check = True  # in a quick test, turning off range checking in these functions made no noticable difference
                    # in running time, so it seems best to leave it on. (of course, it could be worth trying to
                    # turn of the range check in the futute to see if it makes a difference.)
ranges['factorial'] = 170
ranges['bernoulli'] = 200
ranges['two_pi_powers'] = 200
ranges['two_pi_power_over_factorial'] = 200
ranges['factorial_over_twopi_power'] = 200
ranges['bernoulli_over_factorial'] = 200
ranges['exp_t_over_N_squared'] = 30
ranges['gamma_s_over_2'] = 330
ranges['binomial'] = 100;
ranges['inverse_binomial'] = 100;
ranges['inverse'] = 500;

def write_all_tables():
    outfile = open("precomputed_tables.h", 'w')

    outfile.write("#include <iostream>\n")
    outfile.write("#include <cstdlib>\n")

    if range_check:
        outfile.write('const bool SKIP_RANGE_CHECK = false;\n')
    else:
        outfile.write('const bool SKIP_RANGE_CHECK = true;\n')

    write_factorial_table(outfile)
    write_bernoulli_table(outfile)
    write_twopi_power_table(outfile)
    write_twopi_over_factorial_table(outfile)
    write_factorial_over_twopi_power_table(outfile)
    write_bernoulli_over_factorial_table(outfile)
    write_exp_table(outfile)
    write_gamma_s_over_2_table(outfile)
    write_binomial_table(outfile)
    write_inverse_binomial_table(outfile)
    write_inverse_table(outfile)

    outfile.close()

def write_factorial_table(outfile):
    N = ranges['factorial']
    R = RealField(100)
    outfile.write("const int factorial_table_range = %s;\n" % (N,))

    outfile.write("const Double factorial_table[factorial_table_range] = {\n")
    for n in srange(N):
        outfile.write( R(factorial(n)).str(truncate=False) )
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double factorial(int n) {
    if(SKIP_RANGE_CHECK || n < factorial_table_range) {
        return factorial_table[n];
    }
    else {
        std::cout << "Warning: the factorial table wasn't large enough. factorial(n) called with n=" << n << ". Exiting. " << std::endl;
        std::exit(1);
        //Double S = factorial_table[factorial_table_range - 1];
        //for(int i = factorial_table_range; i <= n; i++) {
        //    S = S * i;
        //}
        //return S;
    }
}

//inline Double binomial_coefficient(int n, int k) {
//    return factorial(n)/(factorial(k) * factorial(n - k) );
//}
""")



def write_bernoulli_table(outfile):
    N = ranges['bernoulli']
    R = RealField(100)
    outfile.write("const int bernoulli_range = %s;\n" % (N,))

    outfile.write("const Double bernoulli_table[bernoulli_range] = {\n")
    for n in srange(N):
        outfile.write( R(bernoulli(n)).str(truncate=False))
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double bernoulli(int n) {
    if(SKIP_RANGE_CHECK || n < bernoulli_range) {
        return bernoulli_table[n];
    }
    else {
        std::cout << "Warning: the bernoulli table wasn't large enough. bernoulli(n) called with n=" << n << ". Exiting" << std::endl;
        std::exit(1);
    }
}
""")

def write_inverse_table(outfile):
    N = ranges['inverse']
    R = RealField(100)
    outfile.write("const int inverse_range = %s;\n" % (N,))

    outfile.write("const Double inverse_table[inverse_range] = {\n")
    for n in srange(N):
        if(n == 0):
            outfile.write('0.0/0.0')
        else:
            outfile.write( R(1/n).str(truncate=False))
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double inverse(int n) {
    if(n < inverse_range) {
        return inverse_table[n];
    }
    else {
        return 1.0/n;
    }
}
""")





def write_twopi_power_table(outfile):
    N = ranges['two_pi_powers']
    R = RealField(100)
    outfile.write("const int two_pi_power_table_range = %s;\n" % (N,))

    outfile.write("const Double two_pi_powers[two_pi_power_table_range] = {\n")
    for n in srange(N):
        outfile.write( R((2 * pi)^n).str(truncate=False) )
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double two_pi_power(int n) {
    if(SKIP_RANGE_CHECK || n < two_pi_power_table_range) {
        return two_pi_powers[n];
    }
    else {
        std::cout << "Warning: the table of powers of 2pi wasn't large enough." << std::endl;
        std::exit(1);
        //Double S = two_pi_powers[two_pi_power_table_range - 1];
        //for(int i = two_pi_power_table_range; i <= n; i++) {
        //    S = S * (2.0 * PI);
        //}
        //return S;
    }
}
""")


def write_twopi_over_factorial_table(outfile):
    N = ranges['two_pi_power_over_factorial']
    R = RealField(100)
    outfile.write("const int two_pi_over_factorial_table_range = %s;\n" % (N,))

    outfile.write("const Double two_pi_over_factorial_table[two_pi_over_factorial_table_range] = {\n")
    for n in srange(N):
        outfile.write( R( ((2 * pi)^n)/factorial(n)).str(truncate=False)  )
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double two_pi_over_factorial_power(int n) {
    // return (2 pi)^n / n! using table lookup for n < 100
    if(SKIP_RANGE_CHECK || n < two_pi_over_factorial_table_range) {
        return two_pi_over_factorial_table[n];
    }
    else {
        std::cout << "Warning: the table of (2pi)^n / n! wasn't large enough." << std::endl;
        std::exit(1);
        //return two_pi_power(n)/factorial(n);
    }
}
""")



def write_factorial_over_twopi_power_table(outfile):
    N = ranges['factorial_over_twopi_power']
    R = RealField(100)
    outfile.write("const int factorial_over_twopi_power_table_range = %s;\n" % (N,))

    outfile.write("const Double factorial_over_twopi_power_table[factorial_over_twopi_power_table_range] = {\n")
    for n in srange(N):
        outfile.write( R( (factorial(n)/(2 * pi)^n)).str(truncate=False)  )
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double factorial_over_twopi_power(int n) {
    // return (2 pi)^n / n! using table lookup for n < 100
    if(SKIP_RANGE_CHECK || n < factorial_over_twopi_power_table_range) {
        return factorial_over_twopi_power_table[n];
    }
    else {
        std::cout << "Warning: the table of n! / (2 pi )^n wasn't large enough." << std::endl;
        std::exit(1);
        //return two_pi_power(n)/factorial(n);
    }
}
""")





def write_bernoulli_over_factorial_table(outfile):
    N = ranges['bernoulli_over_factorial']
    R = RealField(100)
    outfile.write("const int bernoulli_over_factorial_table_range = %s;\n" % (N,))

    outfile.write("const Double bernoulli_over_factorial_table[bernoulli_over_factorial_table_range] = {\n")
    for n in srange(N):
        outfile.write( R( bernoulli(n)/factorial(n)).str(truncate=False)  )
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double bernoulli_over_factorial(int n) {
    // return B_n / n! using table lookup;
    if(SKIP_RANGE_CHECK || n < bernoulli_over_factorial_table_range) {
        return bernoulli_over_factorial_table[n];
    }
    else {
        std::cout << "bernoulli_over_factorial called out of range." << std::endl;
        std::exit(1);
        return 0.0/0.0;
    }
}
""")

def write_exp_table(outfile):
    M = ranges['exp_t_over_N_squared']
    R = RealField(100)
    outfile.write("const int exp_t_over_N_squared_range = %s;\n" % (M,))

    outfile.write("const Double exp_t_over_N_squared_table[exp_t_over_N_squared_range][exp_t_over_N_squared_range] = {\n")
    for N in srange(M):
        outfile.write("{")
        for t in srange(M):
            if t < N:
                outfile.write( R(exp(-((t/N)^2))).str(truncate=False) )
            else:
                outfile.write("-1.0")
            if(t < M - 1):
                outfile.write(",\n")
        outfile.write("}")
        if(N < M - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double exp_t_over_N_squared(int t, int N) {
    if(SKIP_RANGE_CHECK || N < exp_t_over_N_squared_range) {
        return exp_t_over_N_squared_table[N][t];
    }
    else {
        //std::cout << "exp_t_over_N_squared called out of range." << std::endl;
        //std::exit(1);
        return exp(-(t/(double)N)*(t/(double)N));
    }
}
""")


def write_binomial_table(outfile):
    M = ranges['binomial']
    R = RealField(100)
    outfile.write("const int binomial_range = %s;\n" % (M,))

    outfile.write("const Double binomial_table[binomial_range][binomial_range] = {\n")
    for n in srange(M):
        outfile.write("{")
        for m in srange(M):
            if m <= n:
                outfile.write( R(binomial(n,m)).str(truncate=False) )
            else:
                outfile.write("-1.0")
            if(m < M - 1):
                outfile.write(",\n")
        outfile.write("}")
        if(n < M - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double binomial_coefficient(int n, int m) {
    if(SKIP_RANGE_CHECK || (n < binomial_range && m <= n)) {
        return binomial_table[n][m];
    }
    else {
        std::cout << "binomial_coefficient() called out of range." << std::endl;
        std::cout << "called with n = " << n << std::endl;
        std::cout << "called with m = " << m << std::endl;
        std::exit(1);
    }
}
""")

def write_inverse_binomial_table(outfile):
    M = ranges['inverse_binomial']
    R = RealField(100)
    outfile.write("const int inverse_binomial_range = %s;\n" % (M,))

    outfile.write("const Double inverse_binomial_table[inverse_binomial_range][inverse_binomial_range] = {\n")
    for n in srange(M):
        outfile.write("{")
        for m in srange(M):
            if m <= n:
                outfile.write( R(1/binomial(n,m)).str(truncate=False) )
            else:
                outfile.write("-1.0")
            if(m < M - 1):
                outfile.write(",\n")
        outfile.write("}")
        if(n < M - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double inverse_binomial_coefficient(int n, int m) {
    if(SKIP_RANGE_CHECK || (n < inverse_binomial_range && m <= n)) {
        return binomial_table[n][m];
    }
    else {
        std::cout << "inverse_binomial_coefficient() called out of range." << std::endl;
        std::cout << "called with n = " << n << std::endl;
        std::cout << "called with m = " << m << std::endl;
        std::exit(1);
    }
}
""")






def write_gamma_s_over_2_table(outfile):
    M = ranges['gamma_s_over_2']
    R = RealField(100)
    outfile.write("const int gamma_s_over_2_range = %s;\n" % (M,))

    outfile.write("const Double gamma_s_over_2_table[gamma_s_over_2_range] = {\n")
    for s in srange(M):
        if s == 0:
            outfile.write( "0.0" );
        else:
            outfile.write( R(gamma(s/2)).str(truncate=False) ) 
        if(s < M - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double gamma_s_over_2(int s) {
    if(SKIP_RANGE_CHECK || s < gamma_s_over_2_range) {
        return gamma_s_over_2_table[s];
    }
    else {
        std::cout << "Warning. gamma_s_over_2() called with s out of range." << std::endl;
        std::exit(1);
        //return 0.0/0.0;
    }
}
""")


