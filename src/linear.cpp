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

    std::vector<std::vector<int>> M = M_input;

    // T used to remake the relationships that have dependencies
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
        // Try to find pivot row
        int pivot_row = -1;
        for (int row = 0; row < m; row++)
        {
            if (!processed_rows[row] && M[row][col] == 1)
            {
                pivot_row = row; // Found a pivot row
                break;
            }
        }

        if (pivot_row == -1)
        {
            // No pivot found in this column
            continue;
        }

        // Mark the pivot row and column as processed so we don't use them again
        processed_cols[col] = true;
        processed_rows[pivot_row] = true;

// To eliminate all 1's in the column from other rows
#pragma omp parallel for // we can parallelize this operation since each row operation is independent
        for (int row = 0; row < m; row++)
        {
            if (row != pivot_row && M[row][col] == 1)
            {
                for (int j = col; j < n; j++)
                {
                    M[row][j] = (M[row][j] + M[pivot_row][j]) % 2;
                }

                // Keep track of chanes in T
                for (int j = 0; j < m; j++)
                {
                    T[row][j] = (T[row][j] + T[pivot_row][j]) % 2;
                }
            }
        }
    }

    dependencies.clear();

// One again, can parallelize check for zero rows
#pragma omp parallel
    {
        std::vector<std::vector<int>> local_dependencies;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < m; i++)
        {
            if (processed_rows[i])
                continue; // Skip processed rows as they have pivots

            bool is_zero_row = true;
            for (int j = 0; j < n; j++)
            { //checks each column in the row for nonzero value
                if (M[i][j] != 0)
                {
                    is_zero_row = false;
                    break;
                }
            }

            if (is_zero_row)
            { // means we found a dependency
                local_dependencies.push_back(T[i]);
            }
        }

// Need to check rows that were processed since they might have been made zero rows
#pragma omp for schedule(dynamic)
        for (int i = 0; i < m; i++)
        {
            if (!processed_rows[i])
                continue; // just did this

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

// Merge into global dependencies
#pragma omp critical
        {
            dependencies.insert(dependencies.end(), local_dependencies.begin(), local_dependencies.end());
        }
    }

    return !dependencies.empty();
}

mpz_class solve_dependency(const std::vector<Relation> &relations,
                           const std::vector<int> &dep,
                           const mpz_class &N)
{
    int num_relations_used = 0;
    for (size_t i = 0; i < dep.size(); i++)
    {    // Counting each relation in dependency 
        if (dep[i] == 1)
        {
            num_relations_used++;
        }
    }

    std::vector<unsigned long> a_factors;
    a_factors.reserve(num_relations_used * 10); // Estimate for average factors per relation

    mpz_class B = 1;

    // Process each relation in the dependency
    for (size_t i = 0; i < dep.size(); i++)
    {
        if (dep[i] == 1)
        {
            mpz_mul(B.get_mpz_t(), B.get_mpz_t(), relations[i].x.get_mpz_t());
            B %= N;
        }
    }

    // Compute A by square root of the product of all x^2 - n (Q) values
    mpz_class A_sq = 1;

    for (size_t i = 0; i < dep.size(); i++)
    {
        if (dep[i] == 1)
        {
            mpz_mul(A_sq.get_mpz_t(), A_sq.get_mpz_t(), relations[i].Q.get_mpz_t());
        }
    }
    // basic principle
    mpz_class A = isqrt(A_sq);
    mpz_class diff, factor;
    diff = B - A;
    mpz_gcd(factor.get_mpz_t(), diff.get_mpz_t(), N.get_mpz_t());

    if (factor == 1 || factor == N)
    {
        diff = B + A;
        mpz_gcd(factor.get_mpz_t(), diff.get_mpz_t(), N.get_mpz_t());
    }

    return factor;
}