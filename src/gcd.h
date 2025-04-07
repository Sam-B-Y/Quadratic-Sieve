#ifndef GCD_H
#define GCD_H

#include <gmpxx.h>

// Find the GCD of two numbers using the Euclidean algorithm
mpz_class gcd(const mpz_class &a, const mpz_class &b);

#endif // GCD_H