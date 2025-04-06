#ifndef PROBABLE_PRIME_H
#define PROBABLE_PRIME_H

#include <gmpxx.h>

// Checks if n is a strong probable prime using the Miller-Rabin test
bool isProbablePrime(const mpz_class &n, int reps);

#endif // PROBABLE_PRIME_H