#include "gcd.h"

mpz_class gcd(const mpz_class &a, const mpz_class &b) {
    if (b == 0) {
        return a;
    }
    return gcd(b, a % b);
}