#include "smooth_relations.h"
#include <omp.h>

using namespace std;

mpz_class isqrt(const mpz_class &n)
{ // ceil(sqrt(n))
    mpz_class root;
    mpz_sqrt(root.get_mpz_t(), n.get_mpz_t());
    if (root * root < n)
        root += 1; //so that we returnt the ceiling of sqrt n. 
    return root;
}

// Fast modular exponentiation for unsigned long integers
// just uses bitwise shifing and squaring
unsigned long mod_exp(unsigned long base, unsigned long exp, unsigned long mod)
{
    unsigned long result = 1;
    base %= mod;
    while (exp > 0) // iteration over each bit
    {
        if (exp & 1)
            result = (result * base) % mod;
        exp >>= 1;
        base = (base * base) % mod;
    }
    return result;
}

// Tonelli-Shanks algorithm for finding square roots modulo p
// outputs the square roots x where x^2 = a mod p
// as seen in https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
vector<unsigned long> tonelli_shanks(const mpz_class &a_mpz, unsigned long p)
{
    vector<unsigned long> sol;
    if (p == 2)
    {
        sol.push_back(a_mpz.get_ui() % 2);
        return sol;
    }

    // Work with unsigned long: a = a mod p.
    unsigned long a = mpz_class(a_mpz % p).get_ui();

    // Check if a is a quadratic residue mod p. 
    // requirement for the algorithm
    if (mod_exp(a, (p - 1) / 2, p) != 1)
    {
        return sol;
    }

    // 1. Write p-1 as Q * 2^S with Q odd
    unsigned long S = 0;
    unsigned long Q = p - 1;
    while ((Q & 1UL) == 0)
    {
        Q >>= 1;
        S++;
    }
    //2.  Find a quadratic non-residue z
    unsigned long z = 2;
    while (mod_exp(z, (p - 1) / 2, p) == 1)
    {
        z++;
    }

    //3. defintions
    unsigned long c = mod_exp(z, Q, p);
    unsigned long R = mod_exp(a, (Q + 1) / 2, p);
    unsigned long t = mod_exp(a, Q, p);
    unsigned long M = S;


    //4. loop
    while (t != 1)
    {
        // Find the smallest integer i (0 < i < M) such that t^(2^i) ≡ 1 (mod p)
        unsigned long temp = t;
        unsigned long i = 0;
        for (; i < M; i++) 
        {
            if (temp == 1)
            {
                break;
            }
            temp = (temp * temp) % p; //squares have happened 2^i times at this step
        }

        // Compute b = c^(2^(M-i-1)) mod p
        unsigned long exp = 1UL << (M - i - 1); // exp = 2^{M - i - 1} (just faster using bit manipulation)
        unsigned long b = mod_exp(c, exp, p);
        R = (R * b) % p;
        t = (t * b * b) % p;
        c = (b * b) % p;
        M = i;
    }
    sol.push_back(R); //first sol

    // The other solution is p - R
    if (R != 0)
    {

        sol.push_back(p - R);
    }
    return sol;
}

// finds B-smooth values over a given interval
vector<Relation> find_smooth_relations(const mpz_class &N,
                                       const vector<unsigned long> &factor_base,
                                       unsigned long sieve_interval,
                                       vector<Relation> &existing_relations,
                                       mpz_class &start_x)
{
    // Use logarithmic sieving
    vector<double> log_factor_base(factor_base.size());

#pragma omp parallel for // parallelize the initialization of log_factor_base
    for (size_t i = 0; i < factor_base.size(); i++)
    {
        log_factor_base[i] = log(factor_base[i]);
    }

    vector<double> sieve_array(sieve_interval, 0.0); 

// Initialize sieve_array with logarithmic values of x^2 - N
#pragma omp parallel for // parallelize the initialization of sieve_array
    for (unsigned long i = 0; i < sieve_interval; i++)
    {
        mpz_class x = start_x + i;
        mpz_class Q;

        // x^2 - N
        mpz_mul(Q.get_mpz_t(), x.get_mpz_t(), x.get_mpz_t());
        Q -= N;

        // Store the logarithm of the absolute value
        sieve_array[i] = log(abs(Q.get_d()));
    }

    // Store original values for precise checking later
    vector<mpz_class> q_values(sieve_interval);

#pragma omp parallel for
    for (unsigned long i = 0; i < sieve_interval; i++)
    {
        mpz_class x = start_x + i;
        mpz_class Q;

        // x^2 - N
        mpz_mul(Q.get_mpz_t(), x.get_mpz_t(), x.get_mpz_t());
        Q -= N;

        q_values[i] = Q;
    }

// For each prime p in the factor base, subtract log(p) at appropriate positions
#pragma omp parallel for schedule(dynamic)
    for (size_t idx = 0; idx < factor_base.size(); idx++)
    {
        unsigned long p = factor_base[idx];
        double log_p = log_factor_base[idx];

        // handling p = 2 specially
        if (p == 2)
        {
            for (unsigned long i = 0; i < sieve_interval; i++)
            {
                if (mpz_divisible_ui_p(q_values[i].get_mpz_t(), 2))
                {
                    // For each power of 2 that divides q_values[i], subtract log(2)
                    mpz_class temp = q_values[i];
                    while (mpz_divisible_ui_p(temp.get_mpz_t(), 2))
                    {
#pragma omp atomic
                        sieve_array[i] -= log_p;
                        mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), 2);
                    }
                }
            }
            continue;
        }

        // using tonelli shanks to find solutions quickly
        vector<unsigned long> sols = tonelli_shanks(N, p);
        if (sols.empty())
            continue;

        for (unsigned long r : sols)
        {
            unsigned long start_x_mod = mpz_class(start_x % p).get_ui();
            unsigned long offset = (r >= start_x_mod) ? (r - start_x_mod) : (p - (start_x_mod - r));

            // Subtract log(p) for each position divisible by p
            for (unsigned long i = offset; i < sieve_interval; i += p)
            {
                mpz_class temp = q_values[i];
                while (mpz_divisible_ui_p(temp.get_mpz_t(), p))
                {
#pragma omp atomic
                    sieve_array[i] -= log_p;
                    mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), p);
                }
            }
        }
    }

    // Identify potentially smooth candidates
    vector<unsigned long> candidates;
    for (unsigned long i = 0; i < sieve_interval; i++)
    {
        if (sieve_array[i] < 0.1)
        { // include floating point errors
            candidates.push_back(i);
        }
    }

    // Process the candidates to find actual B-smooth relations
    vector<Relation> relations = existing_relations; // Start with existing relations

// Use parallel processing for checking candidates
#pragma omp parallel
    {
        vector<Relation> local_relations;

#pragma omp for schedule(dynamic)
        for (size_t candidate_idx = 0; candidate_idx < candidates.size(); candidate_idx++)
        {
            unsigned long i = candidates[candidate_idx];

// Skip if we already have enough relations
#pragma omp critical
            {
                if (relations.size() >= factor_base.size() + 1)
                    continue;
            }

            // Verify smoothness by trial division
            mpz_class temp = q_values[i];

            for (unsigned long p : factor_base)
            {
                while (mpz_divisible_ui_p(temp.get_mpz_t(), p))
                {
                    mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), p);
                }
            }

            // If temp is 1 or -1, we have a B-smooth number
            if (temp == 1 || temp == -1)
            {
                Relation rel;
                rel.x = start_x + i;
                rel.Q = q_values[i];

                // Build the exponent vector mod 2
                vector<int> vec;

                // For sign: if Q(x) is negative, record a 1 for -1
                vec.push_back(q_values[i] < 0 ? 1 : 0);

                // Work with the absolute value
                temp = q_values[i];
                if (temp < 0)
                    temp = -temp;

                // For each prime in factor_base, count the exponent (mod 2) by trial division
                for (unsigned long p : factor_base)
                {
                    int count = 0;
                    while (mpz_divisible_ui_p(temp.get_mpz_t(), p))
                    {
                        mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), p);
                        count++;
                    }
                    vec.push_back(count % 2); // add the exponent mod 2
                }

                rel.exponents = vec;
                local_relations.push_back(rel);
            }
        }

// Merge current relations into the global relations vector
#pragma omp critical
        {
            relations.insert(relations.end(), local_relations.begin(), local_relations.end());
        }
    }

    // Update the start_x for the next iteration
    start_x = start_x + sieve_interval;

    return relations;
}