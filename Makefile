# Makefile for Quadratic Sieve using GMP
# Designed for macOS with Homebrew on Apple Silicon

CXX = g++
CXXFLAGS = -std=c++11 -Wall -I/opt/homebrew/include # might need to adjust include path for GMP
LDFLAGS = -L/opt/homebrew/lib -lgmpxx -lgmp # likewise, adjust library path for GMP

SRC = src/main.cpp src/smoothness_bound.cpp src/factors.cpp src/probable_prime.cpp src/gcd.cpp src/smooth_relations.cpp src/linear_dependencies.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = quadratic_sieve

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)