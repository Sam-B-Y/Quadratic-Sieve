#ifndef FACTORS_H
#define FACTORS_H

#include <vector>
#include <gmpxx.h>

// Generates the primes from 2 to B for which n is a quadratic residue modulo p
std::vector<unsigned long> generatefactors(unsigned long B, const mpz_class &n);

#endif // FACTORS_H
