# Makefile for Quadratic Sieve using GMP
# Designed for macOS with Homebrew on Apple Silicon

CXX = g++
CXXFLAGS = -std=c++11 -Wall -I/opt/homebrew/include # might need to adjust include path for GMP
LDFLAGS = -L/opt/homebrew/lib -lgmpxx -lgmp # likewise, adjust library path for GMP

SRC = main.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = quadratic_sieve

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)