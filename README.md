# Quadratic Sieve Implementation

A C++ implementation of the Quadratic Sieve algorithm for integer factorization, optimized with OpenMP parallelization.

## Overview

The Quadratic Sieve is one of the most efficient general-purpose factorization algorithms for integers in the 20-120 digit range. This implementation provides:

- OpenMP parallelization for performance on multi-core systems
- Logarithmic sieving for efficient smooth relation finding
- Optimized Gaussian elimination for large sparse binary matrices
- Dynamic parameter selection based on input size

## Features

- **Customizable Smoothness Bound**: Automatically calculates an optimal smoothness bound based on theoretical results, but allows user customization.
- **Parallel Processing**: Utilizes OpenMP to parallelize both the sieving phase and Gaussian elimination.
- **Small Prime Optimization**: Quickly removes small prime factors before applying the quadratic sieve.
- **Trial Division Avoidance**: Uses the Miller-Rabin primality test to avoid unnecessary work on prime inputs.
- **Memory-Efficient**: Implements optimizations to minimize memory usage during the sieving phase.

## Requirements

- C++11 or higher compiler
- GMP (GNU Multiple Precision Arithmetic Library)
- OpenMP support

## Installation

1. Install the required dependencies:

```bash
# Debian/Ubuntu
sudo apt-get install libgmp-dev
sudo apt-get install libomp-dev

# macOS with Homebrew
brew install gmp
brew install libomp
```

2. Clone this repository:

```bash
git clone https://github.com/yourusername/quadratic-sieve.git
cd quadratic-sieve
```

3. Compile the code:

```bash
make
```

## Usage

Run the program:

```bash
./quadratic_sieve
```

The program will prompt you to enter a composite number and will output the factors.

## Configuration

Edit the `config.h` file to customize algorithm parameters:

- `MAX_DIGITS`: Maximum number of digits allowed for input
- `MAX_ITERATIONS`: Number of Miller-Rabin iterations for primality testing
- `EXIT_ON_MILLER_RABIN_FAIL`: Whether to test for primality before sieving or not
- `MIN_SMOOTHNESS_BOUND`: Minimum value for the smoothness bound
- `SIEVE_INTERVAL`: Initial sieve interval size
- `MAX_SIEVE_INTERVAL`: Maximum sieve interval size
- `VERBOSE`: Set to 1 to enable verbose output

## Technical Details

### Algorithm Overview

The quadratic sieve works by:

1. Searching for values of x where f(x) = x² - n is B-smooth (has all prime factors ≤ B)
2. Collecting enough smooth relations to form a dependency in the exponent vectors mod 2
3. Using this dependency to construct a congruence of squares: a² ≡ b² (mod n)
4. Computing gcd(a-b, n) to find a non-trivial factor

### Implementation Notes

#### Smooth Relation Finding

The `smooth_relations.cpp` module implements:

- Efficient logarithmic sieving
- Tonelli-Shanks algorithm for solving quadratic congruences
- Parallel processing of sieve intervals

#### Linear Algebra

The `linear.cpp` module provides:

- Parallelized Gaussian elimination for solving large binary matrices
- Efficient dependency finding
- Optimization for sparse matrices

#### Polynomial Evaluation

The sieving step evaluates the polynomial f(x) = x² - n over a large interval, using:

- Optimized modular arithmetic with GMP
- Logarithmic approximations for quick smoothness checks

## Performance

Performance varies based on the input size and hardware. Some general guidelines:

- Numbers up to 35 digits: Typically factors in seconds
- Numbers 35-60 digits: May require minutes
- Numbers 60-100 digits: May require hours or days, depending on hardware

The parallel implementation scales well with the number of available CPU cores.

## Acknowledgments

This implementation is based on the work of Carl Pomerance and other researchers in the field of computational number theory. Key references include:

- Pomerance, C. (2008). "Smooth numbers and the quadratic sieve" in Algorithmic Number Theory.
- Crandall, R., & Pomerance, C. (2005). Prime numbers: A computational perspective.
- Davis, J., & Montgomery, P. (1984). "A block Lanczos algorithm for finding dependencies over GF(2)"
