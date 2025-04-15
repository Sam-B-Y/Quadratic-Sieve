#ifndef CONFIG_H
#define CONFIG_H

// Maximum number of digits allowed for the composite number
#define MAX_DIGITS 100

// Maximum number of iterations for the Miller-Rabin test
#define MAX_ITERATIONS 20

// Print out info messages
#define VERBOSE 1

// Exit if Miller-Rabin test fails
#define EXIT_ON_MILLER_RABIN_FAIL 1

// Minimum smoothness bound
#define MIN_SMOOTHNESS_BOUND 1000

// Sieve Interval
#define SIEVE_INTERVAL 10000
#define MAX_SIEVE_INTERVAL 10000000

#endif // CONFIG_H