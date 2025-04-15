#ifndef LINEAR_H
#define LINEAR_H

#include <gmpxx.h>
#include <vector>
#include "smooth_relations.h"

struct Relation;

// Gaussian elimination to find all dependencies in a matrix
bool gaussian_elimination_all(const std::vector<std::vector<int>> &M_input, std::vector<std::vector<int>> &dependencies);

// Solve the dependency relation to find a nontrivial factor
mpz_class solve_dependency(const std::vector<Relation> &relations,
                           const std::vector<int> &dep,
                           const mpz_class &N);

#endif // LINEAR_H