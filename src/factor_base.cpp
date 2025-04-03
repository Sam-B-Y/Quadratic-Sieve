#include "factor_base.h"


std::vector<unsigned long> generateFactorBase(unsigned long B, const mpz_class &n) {
    // Initally mark all numbers from 0 to B as prime
    std::vector<bool> is_prime(B + 1, true);
    is_prime[0] = is_prime[1] = false;
    
    // Sieve of Eratosthenes to find all primes up to B
    for (unsigned long i = 2; i * i <= B; ++i) {
        if (is_prime[i]) {
            for (unsigned long j = i * i; j <= B; j += i) {
                is_prime[j] = false;
            }
        }
    }
    
    std::vector<unsigned long> factor_base;
    
    if (B >= 2) {
        factor_base.push_back(2);
    }
    
    // Iterate through odd primes
    for (unsigned long i = 3; i <= B; i += 2) {
        if (is_prime[i]) {

            // Check if n is a quadratic residue modulo p
            mpz_class p(i);
            mpz_class n_mod_p = n % p;
            mpz_class exponent = (p - 1) / 2;
            mpz_class legendre_symbol;

            mpz_powm(legendre_symbol.get_mpz_t(), n_mod_p.get_mpz_t(), exponent.get_mpz_t(), p.get_mpz_t()); // Might want to implement this manually
            
            if (legendre_symbol == 1) {
                factor_base.push_back(i);
            }
        }
    }
    
    return factor_base;
}
