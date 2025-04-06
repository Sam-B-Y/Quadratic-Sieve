#include <iostream>
#include <vector>
#include <cstdlib>
#include <gmpxx.h> // Use GMP library to handle large integers
#include "config.h" // Global configuration file
#include<unordered_set>
#include "smoothness_bound.h"
#include "factors.h"
#include "probable_prime.h"

using namespace std;

int print_factors_set(const unordered_set<unsigned long> &factors) {
    cout << "Factors: ";
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

    mpz_class n(nStr); // n stores the composite number as a GMP integer
    unordered_set<unsigned long> final_factors; // Set to store the final factors (avoiding duplicates)

    // TODO: Remove small prime factors up to it's logarithm just by trial division (just doing 2 right now)
    while (n % 2 == 0) {
        final_factors.insert(2);
        n /= 2;
    }

    if (n == 1) {
        print_factors_set(final_factors);
        return EXIT_SUCCESS;
    }

    cout << "Remaining composite number after removing small factors: " << n << endl;

    // Miller-Rabin strong probable prime test: if n is prime, do not proceed 
    if (isProbablePrime(n, MAX_ITERATIONS)) {
        
        // If n is prime and has factors, we are done
        if (!final_factors.empty()) {
            if (n != 1) {
                final_factors.insert(n.get_ui()); // Add the remaining composite number to the factors
            }
            print_factors_set(final_factors);
            return EXIT_SUCCESS;
        }

        // If n is prime and has no factors, error out
        cerr << "Error: The number is prime. Enter a composite number." << endl;
        return EXIT_FAILURE;
    }

    cout << "Miller-Rabin passed: number is likely composite." << endl;

    if(nStr.length() > MAX_DIGITS) { // Only proceed if the number of digits isn't too long
        cerr << "Error: Number exceeds the maximum allowed digit limit of " << MAX_DIGITS << "." << endl;
        return EXIT_FAILURE;
    }

    unsigned long B = smoothnessBound(n);
    cout << "Smoothness bound B: " << B << endl;

    std::vector<unsigned long> potential_factors = generatefactors(B, n);
    cout << "Potential prime factors: ";
    for (unsigned long prime : potential_factors) {
        cout << prime << " ";
    }
    cout << endl;

    return EXIT_SUCCESS;
}