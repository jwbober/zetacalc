#include "theta_sums.h"

using namespace std;

int main() {
    for(int j = 0; j < 20; j++) {
        for(double x = 0.0; x < 20; x += .1) {
            for(double y = 0.0; y < 20; y += .03) {
                complex<double> alpha(x,y);
                cout << j << " " << x << " " << y << " " << H(j, alpha, pow(2, -30)) << endl;
            }
        }
    }
    return 0;
}
