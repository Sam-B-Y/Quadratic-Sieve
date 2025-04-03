#include "smoothness_bound.h"
#include <cmath>

unsigned long smoothnessBound(const mpz_class &n) {

    // We can convert mpz_class to double, which is fine as we only need a rough estimate of the log
    double ln_n = log(n.get_d());
    double ln_ln_n = log(ln_n);

    // Calculate the smoothness bound B using chosen formula
    double B_double = exp(0.5 * sqrt(ln_n * ln_ln_n));
    unsigned long B = static_cast<unsigned long>(B_double);

    if (B < 2) {
        B = 2;
    }

    return B;
}
