#include "gmpfrxx/gmpfrxx.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main() {
    mpfr_class x;
    mpfr_class y;

    x.set_prec(400);
    x = 2 * const_pi();

    y.set_prec(x.get_prec());

    y = x;

    cout << setprecision(200);
    cout << y << endl;

    return 0;
}
