#include "probable_prime.h"


bool isProbablePrime(const mpz_class &n, int reps) {
    if (n <= 2) return false;

    mpz_class d = n - 1;
    int s = 0;

    // Factor out powers of 2 from n-1
    while (d % 2 == 0) {
        d /= 2;
        ++s;
    }

    for (int i = 0; i < reps; i++) {
        // Pick a random base in the range [2, n - 2]
        mpz_class a = rand() % (n - 4) + 2;
        mpz_class x;
        mpz_powm(x.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());

        if (x == 1 || x == (n - 1)) continue;

        bool foundWitness = false;
        for (int j = 0; j < s - 1; j++) {
            x = (x * x) % n;
            if (x == n - 1) {
                foundWitness = true;
                break;
            }
        }

        if (!foundWitness) return false;
    }

    return true;
}