#include <iostream>
#include <vector>
#include <cstdlib>
#include <unordered_set>
#include <cmath>
#include <gmpxx.h> // Use GMP library to handle large integers
#include "config.h" // Global configuration file
#include "smoothness_bound.h" 
#include "factors.h"
#include "probable_prime.h"
#include "gcd.h"

using namespace std;

int print_factors_set(const unordered_set<unsigned long> &factors) {
    if (VERBOSE){
        cout << string(60, '-') << endl;
    }
    cout << "Final Factor(s): ";

    // sort the factors in ascending order
    vector<unsigned long> sorted_factors(factors.begin(), factors.end());
    sort(sorted_factors.begin(), sorted_factors.end());

    for (const auto &factor : factors) {
        cout << factor << " ";
    }
    cout << endl;
    return 0;
}
 
int main() {
    // Promt user for composite number n
    string nStr;
    cout << "Enter composite number n: ";
    cin >> nStr;

    for (char c : nStr) { // Check if the input is a valid number
        if (!isdigit(c)) {
            cerr << "Error: Invalid input: '" << c << "'. Enter a valid composite number." << endl;
            return EXIT_FAILURE;
        }
    }

    if(nStr.length() > MAX_DIGITS) { // Only proceed if the number of digits isn't too long
        cerr << "Error: Number exceeds the maximum allowed digit limit of " << MAX_DIGITS << "." << endl;
        return EXIT_FAILURE;
    }

    if (VERBOSE) { 
        cout << string(60, '-') << endl;
    }

    mpz_class n(nStr); // n stores the composite number as a GMP integer
    unordered_set<unsigned long> final_factors; // Set to store the final factors (avoiding duplicates)
    
    // use sieve of Eratosthenes to find small prime factors up till log(n)

    double n_d = n.get_d(); // we can use double as we don't need precision here
    unsigned long limit = static_cast<unsigned long>(ceil(log(n_d)));
    if (VERBOSE) { 
        cout << "Finding factors up to log(n) = " << limit << " using brute force" << endl;
    }

    while (n % 2 == 0) {
        final_factors.insert(2);
        n /= 2;
    }

    // sieve
    vector<bool> is_prime(limit + 1, true);
    if (limit >= 0) {
        is_prime[0] = is_prime[1] = false;
    }

    for (unsigned long i = 2; i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (unsigned long j = i * i; j <= limit; j += i) {
                is_prime[j] = false;
            }
        }
    }

    vector<unsigned long> primes;
    for (unsigned long i = 2; i <= limit; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }

    // Remove small prime factors
    for (unsigned long p : primes) {
        while (n % p == 0) { // keep dividing n by p until it is no longer divisible
            if (final_factors.find(p) == final_factors.end()) { 
                if (VERBOSE) { 
                    cout << "Removed factor: " << p << endl;
                }
                final_factors.insert(p);
            }
            n /= p;
        }
    }

    if (n == 1) {
        print_factors_set(final_factors);
        return EXIT_SUCCESS;
    }

    if (VERBOSE) { 
        cout << "Remaining number after removing small factors: " << n << endl;
    }

    // Miller-Rabin strong probable prime test: if n is prime, do not proceed 
    if (isProbablePrime(n, MAX_ITERATIONS)) {
        if (final_factors.empty()) {
            // If n is prime and has no factors, error out
            cerr << "Error: The number is prime. Enter a composite number." << endl;
            return EXIT_FAILURE;
        }

        // If n is prime and has factors, we are done
        if (n != 1) {
            final_factors.insert(n.get_ui()); // Add the remaining prime number to the factors
        }
        
        print_factors_set(final_factors);
        return EXIT_SUCCESS;
    }

    if (VERBOSE) { 
        cout << "Miller-Rabin passed: number is likely composite." << endl;
    }

    unsigned long B = smoothnessBound(n);
    if (VERBOSE) { 
        cout << "Smoothness bound B: " << B << endl;
    }

    vector<unsigned long> potential_factors = generatefactors(B, n);
    if (VERBOSE) { 
        cout << "Potential prime factors: ";
        for (unsigned long prime : potential_factors) {
            cout << prime << " ";
        }
        cout << endl;
    }

    return EXIT_SUCCESS;
}