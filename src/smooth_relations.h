#ifndef SMOOTH_RELATIONS_H
#define SMOOTH_RELATIONS_H

#include <gmpxx.h>
#include <map>

struct Relation
{
    mpz_class x;
    mpz_class Q;                // Q(x) = x^2 - N
    std::vector<int> exponents; // Exponent vector (mod 2)
};

// Function to find B-smooth relations
std::vector<Relation> find_smooth_relations(
    const mpz_class &N,
    const std::vector<unsigned long> &factor_base,
    unsigned long sieve_interval,
    std::vector<Relation> &existing_relations,
    mpz_class &start_x);

// Compute ceil(sqrt(n))
mpz_class isqrt(const mpz_class &n);

#endif // SMOOTH_RELATIONS_H
