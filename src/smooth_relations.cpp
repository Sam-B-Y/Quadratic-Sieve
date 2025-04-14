#include "smooth_relations.h"
#include "factors.h"

mpz_class computePolynomial(const mpz_class &n, const mpz_class &x) {
    // computing x^2 - n
    return (x*x) - n;
}

std::map<unsigned long, unsigned int> factorPolynomial(const mpz_class f, const std::vector<unsigned long> factor_base) {
    mpz_class rem = f;

    if (rem < 0) {
        rem = -1 * rem;
    }
    
    std::map<unsigned long, unsigned int> factorization;  // exponent vector 

    for (unsigned long p : factor_base) {
        unsigned int exponent = 0;
        while (rem % mpz_class(p) == 0) {
            exponent++;
            rem /= mpz_class(p);
        }
        factorization[p] = exponent;
    }

    // If rem is not 1, then f has a factor outside the factor base.
    if (rem != 1) {
        factorization.clear();  // indicate failure: not B-smooth.
    }
    
    return factorization;
}

//TO DO: Function to find the B-Smooth values from sqrt n onwards.