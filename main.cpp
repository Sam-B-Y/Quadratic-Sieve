#include <iostream>
#include <vector>
#include <cstdlib>
#include <gmpxx.h>

using namespace std;

int main() {
    
    // Promt user for composite number n
    // Use GMP library to handle large integers

    string nStr;
    cout << "Enter composite number n: ";
    cin >> nStr;
    mpz_class n(nStr); // n stores the composite number as a GMP integer

    cout << "Factoring composite number: " << n << endl;

}