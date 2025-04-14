#include "linear_dependencies.h"
#include "smooth_relations.h"



/*
Given input of relations and factor base
want to construct a matrix such that each row corresponds to a relation, and each element in the row is the number of the exponent in the factor base, 
then find a square.

ex:
       2 3 5 7 11 13 17 19
9398 | 0 0 5 0  0  0  0  1
19095| 2 0 1 0  1  1  0  1
1964 | 0 2 0 0  0  3  0  0
17078| 6 2 0 0  1  0  0  0
8077 | 1 0 0 0  0  0  0  1
3397 | 5 0 1 0  0  2  0  0
14262| 0 0 2 2  0  1  0  0

and then find relations such that the sum of the combinations are even for each exponent.
The return should just be a vector of bools, where true means the relation is used in the final factorization. 
*/

std::vector<bool> findSquareSubset(const std::vector<SmoothRelation>& relations, 
    const std::vector<unsigned long>& factor_base) {
size_t numRelations = relations.size();
size_t numPrimes = factor_base.size();

// Build the matrix A where A[i][j] = (exponent of factor_base[j] in relation i) mod 2.
std::vector<std::vector<int>> A(numRelations, std::vector<int>(numPrimes, 0));
for (size_t i = 0; i < numRelations; ++i) {
for (size_t j = 0; j < numPrimes; ++j) {
unsigned long p = factor_base[j];
// Look up the exponent for prime p in the i-th relation.
auto it = relations[i].factorization.find(p);
int expMod2 = (it != relations[i].factorization.end()) ? (it->second % 2) : 0;
A[i][j] = expMod2;
}
}

// Create an identity matrix D (size: numRelations x numRelations) to track row operations.
std::vector<std::vector<int>> D(numRelations, std::vector<int>(numRelations, 0));
for (size_t i = 0; i < numRelations; ++i) {
D[i][i] = 1;
}

// Gaussian elimination modulo 2.
size_t row = 0;
for (size_t col = 0; col < numPrimes && row < numRelations; ++col) {
// Find a pivot row with a 1 in the current column.
size_t pivot = row;
while (pivot < numRelations && A[pivot][col] == 0)
++pivot;
if (pivot == numRelations)
continue;  // No pivot found in this column.

// Swap pivot row into the current row if needed.
if (pivot != row) {
std::swap(A[row], A[pivot]);
std::swap(D[row], D[pivot]);
}

// Eliminate the 1s in the current column from all other rows.
for (size_t i = 0; i < numRelations; ++i) {
if (i != row && A[i][col] == 1) {
for (size_t j = col; j < numPrimes; ++j)
A[i][j] = (A[i][j] + A[row][j]) % 2;
for (size_t j = 0; j < numRelations; ++j)
D[i][j] = (D[i][j] + D[row][j]) % 2;
}
}
row++;
}

// Look for a row in A that is completely zero.
for (size_t i = 0; i < numRelations; ++i) {
bool allZero = true;
for (size_t j = 0; j < numPrimes; ++j) {
if (A[i][j] != 0) {
allZero = false;
break;
}
}
if (allZero) {
// The corresponding row in D is a nontrivial dependency.
std::vector<bool> dependency(numRelations, false);
for (size_t j = 0; j < numRelations; ++j) {
dependency[j] = (D[i][j] != 0);
}
return dependency;
}
}

// If no dependency is found, return an empty vector.
return std::vector<bool>();
}

