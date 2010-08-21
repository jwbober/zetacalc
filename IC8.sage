def IC8(K, j, a, b):
    z = exp(-2 * pi * i * a^2/(4 * b))
    z = z * 2^(-3 * j/2 - 1) * (b * pi)^( -(j + 1)/2) * K^(-j) * factorial(j) * sqrt(2 * pi) * exp(I * pi/4 + j * 3 * pi * I/4)


    S = 0
    for l in srange(0, j + 1):
        if( (j - l) % 2 == 0 ):
            sign = 0;
            if ((j + l)/2) % 2 == 0:
                sign = 1
            else:
                sign = -1

            S = S + sign * a^l * exp(-3 * pi * I * l/4) * (2 * pi/b)^(l/2)/(factorial(l) * factorial( (j - l)/2))

    S = S * z
    return S
