#ifndef FACTORS_H
#define FACTORS_H

#include <vector>
#include <gmpxx.h>

// Generates the primes from 2 to B for which n is a quadratic residue modulo p
std::pair<std::vector<unsigned long>, std::vector<unsigned long>> generateFactorBase(unsigned long B, const mpz_class &n);

#endif // FACTORS_H
