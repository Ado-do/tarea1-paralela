CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=build/%.o)

TARGETS = build/parallel_main build/sequential_main

all: $(TARGETS)


build/sequential_main: src/sequential_main.cpp src/matrix.hpp | build/
	$(CXX) $(CXXFLAGS) $< -o $@

build/parallel_main: src/parallel_main.cpp src/matrix.hpp | build/
	$(CXX) $(CXXFLAGS) $< -o $@

build/:
	mkdir -p build

clean:
	rm -rf build
