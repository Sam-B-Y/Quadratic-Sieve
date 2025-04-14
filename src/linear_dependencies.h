#ifndef LINEAR_DEPENDENCIES_H
#define LINEAR_DEPENDENCIES_H

#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "smooth_relations.h"

// struct SmoothRelation {
//     mpz_class x;  // Candidate x value.
//     mpz_class fx; // f(x) = x^2 - n (or Q(a) = (a+sqrt(n))^2 - n)
//     std::map<unsigned long, unsigned int> factorization of x^2 - n; Maps prime -> exponent 
//     bool isNegative;
// };

std::vector<bool> findSquareSubset(const std::vector<SmoothRelation>& relations, const std::vector<unsigned long>& factor_base);


#endif // LINEAR_DEPENDENCIES_H