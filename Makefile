CXX = g++
CXXFLAGS = -Wall -std=c++17

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=build/%.o)

TARGETS = build/parallel_main build/sequential_main build/tests_main

all: $(TARGETS)

build/sequential_main: src/sequential_main.cpp src/sequential_algorithms.cpp src/matrix.hpp | build/
	$(CXX) $(CXXFLAGS) src/sequential_main.cpp src/sequential_algorithms.cpp -o $@

build/parallel_main: src/parallel_main.cpp src/parallel_algorithms.cpp src/matrix.hpp | build/
	$(CXX) $(CXXFLAGS) -fopenmp src/parallel_main.cpp src/parallel_algorithms.cpp -o $@

build/tests_main: src/tests_main.cpp src/parallel_algorithms.cpp src/sequential_algorithms.cpp src/tests.hpp | build/
	$(CXX) $(CXXFLAGS) -fopenmp src/tests_main.cpp src/parallel_algorithms.cpp src/sequential_algorithms.cpp -o $@

build/:
	mkdir -p build

clean:
	rm -rf build
