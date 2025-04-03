#include <iostream>
#include <vector>
#include <cstdlib>
#include <gmpxx.h> // Use GMP library to handle large integers
#include "config.h" // Global configuration file
#include "smoothness_bound.h"
#include "factor_base.h"

using namespace std;

int main() {
    
    // Promt user for composite number n

    string nStr;
    cout << "Enter composite number n: ";
    cin >> nStr;

    for (char c : nStr) { // Check if the input is a valid number
        if (!isdigit(c)) {
            cerr << "Error: Invalid input: '" << c << "'. Please enter a valid composite number." << endl;
            return EXIT_FAILURE;
        }
    }

    mpz_class n(nStr); // n stores the composite number as a GMP integer

    if(nStr.length() > MAX_DIGITS) { // Only proceed if the number of digits isn't too long
        cerr << "Error: Number exceeds the maximum allowed digit limit of " << MAX_DIGITS << "." << endl;
        return EXIT_FAILURE;
    }

    cout << "Factoring composite number: " << n << endl;

    unsigned long B = smoothnessBound(n);
    cout << "Smoothness bound B: " << B << endl;

    std::vector<unsigned long> factorBase = generateFactorBase(B, n);
    cout << "Factor base: ";
    for (unsigned long prime : factorBase) {
        cout << prime << " ";
    }
    cout << endl;

    return EXIT_SUCCESS;
}