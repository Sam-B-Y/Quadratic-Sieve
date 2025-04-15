#include "factors.h"

std::pair<std::vector<unsigned long>, std::vector<unsigned long>> generateFactorBase(unsigned long B, const mpz_class &n)
{
    // Initally mark all numbers from 0 to B as prime
    std::vector<bool> is_prime(B + 1, true);
    is_prime[0] = is_prime[1] = false;

    // Sieve of Eratosthenes to find all primes up to B
    for (unsigned long i = 2; i * i <= B; ++i)
    {
        if (is_prime[i])
        {
            for (unsigned long j = i * i; j <= B; j += i)
            {
                is_prime[j] = false;
            }
        }
    }

    std::vector<unsigned long> factor_base;
    std::vector<unsigned long> dividers;

    if (B >= 2)
    {
        factor_base.push_back(2);
    }

    // Iterate through odd primes
    for (unsigned long i = 3; i <= B; i += 2)
    {
        if (is_prime[i])
        {
            // Check if n is a quadratic residue modulo p
            // Legendre symbol (n/p) = n^(p-1)/2 mod p
            mpz_class legendre_symbol = mpz_legendre(n.get_mpz_t(), mpz_class(i).get_mpz_t());

            if (legendre_symbol == 0)
            {
                if (mpz_divisible_ui_p(n.get_mpz_t(), i))
                {
                    dividers.push_back(i);
                }
            }
            else if (legendre_symbol == 1)
            {
                factor_base.push_back(i);
            }
        }
    }

    return {factor_base, dividers};
}
