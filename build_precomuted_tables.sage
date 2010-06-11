import sys

ranges = {}

ranges['factorial'] = 100
ranges['bernoulli'] = 100
ranges['two_pi_powers'] = 100
ranges['two_pi_power_over_factorial'] = 100
ranges['bernoulli_over_factorial'] = 100

def write_all_tables():
    outfile = open("precomputed_tables.h", 'w')

    outfile.write("#include <iostream>\n")

    write_factorial_table(outfile)
    write_bernoulli_table(outfile)
    write_twopi_power_table(outfile)
    write_twopi_over_factorial_table(outfile)
    write_bernoulli_over_factorial_table(outfile)

    outfile.close()

def write_factorial_table(outfile):
    N = ranges['factorial']
    R = RealField(100)
    outfile.write("const int factorial_table_range = %s;\n" % (N,))

    outfile.write("const Double factorial_table[factorial_table_range] = {\n")
    for n in srange(N):
        outfile.write( str(R(factorial(n))))
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double factorial(int n) {
    if(n < factorial_table_range) {
        return factorial_table[n];
    }
    else {
        std::cout << "Warning: the factorial table wasn't large enough. factorial(n) called with n=" << n << std::endl;
        Double S = factorial_table[factorial_table_range - 1];
        for(int i = factorial_table_range; i <= n; i++) {
            S = S * i;
        }
        return S;
    }
}

inline Double binomial_coefficient(int n, int k) {
    return factorial(n)/(factorial(k) * factorial(n - k) );
}
""")


def write_bernoulli_table(outfile):
    N = ranges['bernoulli']
    R = RealField(100)
    outfile.write("const int bernoulli_range = %s;\n" % (N,))

    outfile.write("const Double bernoulli_table[bernoulli_range] = {\n")
    for n in srange(N):
        outfile.write( str(R(bernoulli(n))))
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double bernoulli(int n) {
    if(n < bernoulli_range) {
        return bernoulli_table[n];
    }
    else {
        std::cout << "Warning: the bernoulli table wasn't large enough. bernoulli(n) called with n=" << n << std::endl;
        return 0.0/0.0;
    }
}
""")


def write_twopi_power_table(outfile):
    N = ranges['two_pi_powers']
    R = RealField(100)
    outfile.write("const int two_pi_power_table_range = %s;\n" % (N,))

    outfile.write("const Double two_pi_powers[two_pi_power_table_range] = {\n")
    for n in srange(N):
        outfile.write( str(R((2 * pi)^n)))
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double two_pi_power(int n) {
    if(n < two_pi_power_table_range) {
        return two_pi_powers[n];
    }
    else {
        std::cout << "Warning: the table of powers of 2pi wasn't large enough." << std::endl;
        Double S = two_pi_powers[two_pi_power_table_range - 1];
        for(int i = two_pi_power_table_range; i <= n; i++) {
            S = S * (2.0 * PI);
        }
        return S;
    }
}
""")


def write_twopi_over_factorial_table(outfile):
    N = ranges['two_pi_power_over_factorial']
    R = RealField(100)
    outfile.write("const int two_pi_over_factorial_table_range = %s;\n" % (N,))

    outfile.write("const Double two_pi_over_factorial_table[two_pi_over_factorial_table_range] = {\n")
    for n in srange(N):
        outfile.write( str(R( ((2 * pi)^n)/factorial(n)  )))
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double two_pi_over_factorial_power(int n) {
    // return (2 pi)^n / n! using table lookup for n < 100
    if(n < two_pi_over_factorial_table_range) {
        return two_pi_over_factorial_table[n];
    }
    else {
        std::cout << "Warning: the table of (2pi)^n / n! wasn't large enough." << std::endl;
        return two_pi_power(n)/factorial(n);
    }
}
""")


def write_bernoulli_over_factorial_table(outfile):
    N = ranges['bernoulli_over_factorial']
    R = RealField(100)
    outfile.write("const int bernoulli_over_factorial_table_range = %s;\n" % (N,))

    outfile.write("const Double bernoulli_over_factorial_table[bernoulli_over_factorial_table_range] = {\n")
    for n in srange(N):
        outfile.write( str(R( bernoulli(n)/factorial(n)  )))
        if(n < N - 1):
            outfile.write(",\n")

    outfile.write("};\n\n")
    
    outfile.write("""
inline Double bernoulli_over_factorial(int n) {
    // return B_n / n! using table lookup for n < 100
    if(n < bernoulli_over_factorial_table_range) {
        return bernoulli_over_factorial_table[n];
    }
    else {
        return 0.0/0.0;
    }
}
""")


