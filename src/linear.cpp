#include "linear.h"
#include "smooth_relations.h"
#include <omp.h>

// performs gaussian eliminiation on a mod 2 matrix
bool gaussian_elimination_all(const std::vector<std::vector<int>> &M_input, std::vector<std::vector<int>> &dependencies)
{
    int m = M_input.size();
    if (m == 0)
        return false;
    int n = M_input[0].size();

    // Make a copy
    std::vector<std::vector<int>> M = M_input;

    // T will be an m x m identity matrix that tracks row operations
    // will be used to remake the relationships that have dependencies
    std::vector<std::vector<int>> T(m, std::vector<int>(m, 0));

#pragma omp parallel for
    for (int i = 0; i < m; i++)
        T[i][i] = 1;

    // Keep track of processed rows and columns for faster elimination
    std::vector<bool> processed_rows(m, false);
    std::vector<bool> processed_cols(n, false);

    // Perform Gaussian elimination
    for (int col = 0; col < n; col++)
    {
        // Find a pivot row where M[row][col] == 1
        int pivot_row = -1;
        for (int row = 0; row < m; row++)
        {
            if (!processed_rows[row] && M[row][col] == 1)
            {
                pivot_row = row;
                break;
            }
        }

        if (pivot_row == -1)
        {
            // No pivot found in this column
            continue;
        }

        // Mark this column as processed
        processed_cols[col] = true;

        // Mark this row as processed
        processed_rows[pivot_row] = true;

// Eliminate 1's in this column from all other rows
// We can parallelize this operation since each row operation is independent
#pragma omp parallel for
        for (int row = 0; row < m; row++)
        {
            if (row != pivot_row && M[row][col] == 1)
            {
                // XOR the pivot row with the current row (in binary arithmetic, + is XOR)
                for (int j = col; j < n; j++)
                {
                    M[row][j] = (M[row][j] + M[pivot_row][j]) % 2;
                }

                // Also update the transformation matrix T
                for (int j = 0; j < m; j++)
                {
                    T[row][j] = (T[row][j] + T[pivot_row][j]) % 2;
                }
            }
        }
    }

    // Find dependencies (rows of all zeros in the reduced matrix)
    dependencies.clear();

// We can parallelize the check for zero rows
#pragma omp parallel
    {
        std::vector<std::vector<int>> local_dependencies;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < m; i++)
        {
            if (processed_rows[i])
                continue; // Skip processed rows as they have pivots

            // Check if this row is all zeros
            bool is_zero_row = true;
            for (int j = 0; j < n; j++)
            {
                if (M[i][j] != 0)
                {
                    is_zero_row = false;
                    break;
                }
            }

            if (is_zero_row)
            {
                local_dependencies.push_back(T[i]);
            }
        }

// Also check processed rows - they might have become all zeros during elimination
#pragma omp for schedule(dynamic)
        for (int i = 0; i < m; i++)
        {
            if (!processed_rows[i])
                continue; // Already handled above

            bool is_zero_row = true;
            for (int j = 0; j < n; j++)
            {
                if (M[i][j] != 0)
                {
                    is_zero_row = false;
                    break;
                }
            }

            if (is_zero_row)
            {
                local_dependencies.push_back(T[i]);
            }
        }

// Merge local_dependencies into global dependencies
#pragma omp critical
        {
            dependencies.insert(dependencies.end(), local_dependencies.begin(), local_dependencies.end());
        }
    }

    return !dependencies.empty();
}

// Optimized version for finding non-trivial factors
mpz_class solve_dependency(const std::vector<Relation> &relations,
                           const std::vector<int> &dep,
                           const mpz_class &N)
{
    // Count the number of relations used in this dependency
    int num_relations_used = 0;
    for (size_t i = 0; i < dep.size(); i++)
    {
        if (dep[i] == 1)
        {
            num_relations_used++;
        }
    }

    // Allocate vectors to store the factors efficiently
    std::vector<unsigned long> a_factors;
    a_factors.reserve(num_relations_used * 10); // Estimate for average factors per relation

    mpz_class B = 1;

    // Process each relation in the dependency
    for (size_t i = 0; i < dep.size(); i++)
    {
        if (dep[i] == 1)
        {
            // Multiply B by the relation's x value using efficient squaring/modulo
            mpz_mul(B.get_mpz_t(), B.get_mpz_t(), relations[i].x.get_mpz_t());

            // Take the modulo N after each multiplication to keep numbers small
            B %= N;
        }
    }

    // Compute A by square root of the product of all Q values
    // We'll compute this by first calculating the product
    mpz_class A_sq = 1;

    for (size_t i = 0; i < dep.size(); i++)
    {
        if (dep[i] == 1)
        {
            mpz_mul(A_sq.get_mpz_t(), A_sq.get_mpz_t(), relations[i].Q.get_mpz_t());
        }
    }

    // Calculate square root more efficiently
    mpz_class A = isqrt(A_sq);

    // Try to find a non-trivial factor
    mpz_class diff, factor;

    // First try: gcd(B - A, N)
    diff = B - A;
    mpz_gcd(factor.get_mpz_t(), diff.get_mpz_t(), N.get_mpz_t());

    if (factor == 1 || factor == N)
    {
        // Second try: gcd(B + A, N)
        diff = B + A;
        mpz_gcd(factor.get_mpz_t(), diff.get_mpz_t(), N.get_mpz_t());
    }

    return factor;
}