#include <iostream>
#include <vector>
#include <cstdlib>
#include <set>
#include <chrono>
#include <cmath>
#include <gmpxx.h>  // Use GMP library to handle large integers
#include "config.h" // Global configuration file
#include "smoothness_bound.h"
#include "factors.h"
#include "smooth_relations.h"
#include "probable_prime.h"
#include "linear.h"

using namespace std;

int print_factors_set(const set<mpz_class> &factors, const chrono::high_resolution_clock::time_point &start)
{
    if (VERBOSE)
    { // Print out time taken to find factors
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

        cout << "Time taken: " << duration.count() / 1000.0 << " seconds" << endl;

        cout << string(60, '-') << endl;
    }

    cout << "Final Factor(s): ";

    for (const auto &factor : factors)
    {
        cout << factor.get_str() << " ";
    }
    cout << endl;
    return 0;
}

int main()
{
    // Promt user for composite number n
    string nStr;
    cout << "Enter composite number n: ";
    cin >> nStr;

    auto start = chrono::high_resolution_clock::now();

    for (char c : nStr)
    { // Check if the input is a valid number
        if (!isdigit(c))
        {
            cerr << "Error: Invalid input: '" << c << "'. Enter a valid composite number." << endl;
            return EXIT_FAILURE;
        }
    }

    if (nStr.length() > MAX_DIGITS)
    { // Only proceed if the number of digits isn't too long
        cerr << "Error: Number exceeds the maximum allowed digit limit of " << MAX_DIGITS << "." << endl;
        return EXIT_FAILURE;
    }

    if (VERBOSE)
    {
        cout << string(60, '-') << endl;
    }

    mpz_class n(nStr); // n stores the composite number as a GMP integer to handle large numbers

    set<mpz_class> final_factors; // Use a set to store unique factors

    // use sieve of Eratosthenes to find small prime factors up till log(n)

    double n_d = n.get_d(); // we can use double as we don't need precision here as log can be approximate
    unsigned long limit = static_cast<unsigned long>(ceil(log(n_d)));
    if (VERBOSE)
    {
        cout << "Finding factors up to log(n) = " << limit << " using brute force" << endl;
    }

    // actual sieve
    vector<bool> is_prime(limit + 1, true);
    if (limit >= 0)
    {
        is_prime[0] = is_prime[1] = false;
    }

    for (unsigned long i = 2; i * i <= limit; ++i)
    {
        if (is_prime[i])
        {
            for (unsigned long j = i * i; j <= limit; j += i)
            {
                is_prime[j] = false;
            }
        }
    }

    vector<unsigned long> primes;
    for (unsigned long i = 2; i <= limit; ++i)
    {
        if (is_prime[i])
        {
            primes.push_back(i);
        }
    }

    // Remove small prime factors
    for (unsigned long p : primes)
    {
        while (n % p == 0)
        { // keep dividing n by p until it is no longer divisible
            mpz_class prime_factor(p);
            if (final_factors.find(prime_factor) == final_factors.end())
            {
                if (VERBOSE)
                {
                    cout << "Removed factor: " << p << endl;
                }
                final_factors.insert(prime_factor);
            }
            n /= p;
        }
    }

    if (n == 1)
    { // If n is fully factorized
        print_factors_set(final_factors, start);
        return EXIT_SUCCESS;
    }

    if (VERBOSE)
    {
        cout << "Remaining number after removing small factors: " << n << endl;
    }

    // Miller-Rabin strong probable prime test: if n is prime, do not proceed
    if (EXIT_ON_MILLER_RABIN_FAIL && isProbablePrime(n, MAX_ITERATIONS))
    {
        if (final_factors.empty())
        {
            // If n is prime and has no factors, error out
            cerr << "Error: The number is prime. Enter a composite number." << endl;
            return EXIT_FAILURE;
        }

        // If n is prime and has factors, we are done
        if (n != 1)
        {
            final_factors.insert(n); // Add the remaining prime number to the factors
        }

        print_factors_set(final_factors, start);
        return EXIT_SUCCESS;
    }

    if (EXIT_ON_MILLER_RABIN_FAIL && VERBOSE)
    {
        cout << "Miller-Rabin passed: number is likely composite." << endl;
    }

    // TODO add check for n being a power using Newton's method

    unsigned long B = smoothnessBound(n);
    if (VERBOSE)
    {
        cout << "Smoothness bound B: " << B << endl;
    }

    // Generate the factor base
    pair<vector<unsigned long>, vector<unsigned long>> factorBaseAndDividers = generateFactorBase(B, n);
    vector<unsigned long> factorBase = factorBaseAndDividers.first;
    vector<unsigned long> dividers = factorBaseAndDividers.second;

    if (dividers.size() > 0)
    { // In case we found numbers with Legendre symbol 0
        if (VERBOSE)
        {
            cout << "Found small prime factors: ";
            for (unsigned long divider : dividers)
            {
                cout << divider << " ";
            }
            cout << endl;
        }

        // Add the small prime factors to the final factors
        for (unsigned long divider : dividers)
        {
            mpz_class prime_factor(divider);
            final_factors.insert(prime_factor);
            mpz_divexact_ui(n.get_mpz_t(), n.get_mpz_t(), divider);
        }
    }

    if (VERBOSE)
    {
        cout << "Factor Base: ";
        for (unsigned long prime : factorBase)
        {
            cout << prime << " ";
        }
        cout << endl;
    }

    mpz_class sqrt_n;
    mpz_sqrt(sqrt_n.get_mpz_t(), n.get_mpz_t());

    // ceil
    sqrt_n = isqrt(n);

    if (VERBOSE)
    {
        cout << "Starting B-smooth search around x = " << sqrt_n << endl;
    }

    bool found_factor = false;
    unsigned long sieve_interval = SIEVE_INTERVAL;
    int attempt = 0;

    // start searching for smooth relations at sqrt(n)
    vector<Relation> relations;
    mpz_class start_x = sqrt_n;

    // while we haven't found a factor, keep searching by starting at a higher point and increasing the sieve interval
    while (!found_factor)
    {
        attempt++;
        if (VERBOSE)
        {
            cout << "\nAttempt " << attempt << " with sieve interval: " << sieve_interval << endl;
            cout << "Starting search at x = " << start_x << endl;
            cout << "Current relations count: " << relations.size() << endl;
        }

        relations = find_smooth_relations(n, factorBase, sieve_interval, relations, start_x);

        if (VERBOSE)
        {
            cout << "Found " << relations.size() << " smooth relations so far." << endl;
        }

        // Check if we have enough relations to try finding dependencies (need at least pi(B))
        if (relations.size() > factorBase.size())
        {
            vector<vector<int>> matrix;
            for (const auto &rel : relations)
                matrix.push_back(rel.exponents);

            vector<vector<int>> dependencies;
            if (!gaussian_elimination_all(matrix, dependencies))
            {
                cout << "No nontrivial dependency found; need more relations." << endl;
                // We don't need to increase the sieve interval, just continue collecting more relations
                continue;
            }

            cout << "\nFound " << dependencies.size() << " dependency vector(s)." << endl;

            mpz_class factor;
            bool found = false;
            for (size_t k = 0; k < dependencies.size(); k++)
            {
                factor = solve_dependency(relations, dependencies[k], n);
                if (factor != 1 && factor != n)
                {
                    found = true;
                    found_factor = true;
                    cout << "\nDependency vector " << k << " produced a nontrivial factor." << endl;

                    mpz_class other_factor = n / factor;

                    final_factors.insert(factor);
                    final_factors.insert(other_factor);
                    break;
                }
            }

            if (!found)
            {
                cout << "\nNone of the dependency vectors produced a nontrivial factor." << endl;
                // Continue collecting more relations - no need to increase interval
            }
        }
        else
        {
            // If we don't have enough relations yet, continue collecting more
            cout << "Need more relations. Currently have " << relations.size()
                 << " of " << (factorBase.size() + 1) << " required." << endl;
        }

        // Increase the sieve interval after several attempts with no new relations (don't exceed MAX_SIEVE_INTERVAL as program may take too long)
        if (attempt % 5 == 0 && relations.size() < factorBase.size() / 2 && sieve_interval < MAX_SIEVE_INTERVAL)
        {
            sieve_interval *= 10;
            cout << "Increasing sieve interval to " << sieve_interval << endl;
        }
    }

    if (!found_factor)
    {
        cerr << "Failed to find nontrivial factor after " << attempt << " attempts." << endl;
        return EXIT_FAILURE;
    }

    print_factors_set(final_factors, start);
    return EXIT_SUCCESS;
}