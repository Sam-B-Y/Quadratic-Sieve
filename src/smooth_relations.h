#ifndef SMOOTH_RELATIONS_H
#define SMOOTH_RELATIONS_H

#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <map>
#include "smooth_relations.h"

// given a number f, write it as a product of the factor base 
std::map<unsigned long, unsigned int> factorPolynomial(const mpz_class f, const std::vector<unsigned long> factor_base);

struct SmoothRelation {
    mpz_class x;
    mpz_class fx; //x^2 - n
    std::map<unsigned long, unsigned int> factorization;
    bool isNegative;
};

// running for root n, collect b-smooth values
std::vector<SmoothRelation> findSmoothRelations(const mpz_class &n, const std::vector<unsigned long> &factor_base, unsigned long required_relations);
#endif // SMOOTH_RELATIONS_H
