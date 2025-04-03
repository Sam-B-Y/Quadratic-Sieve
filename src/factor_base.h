#ifndef FACTOR_BASE_H
#define FACTOR_BASE_H

#include <vector>
#include <gmpxx.h>

// Generates the primes from 2 to B for which n is a quadratic residue modulo p
std::vector<unsigned long> generateFactorBase(unsigned long B, const mpz_class &n);

#endif // FACTOR_BASE_H
