CXX      := g++
CXXFLAGS := -Wall -std=c++17 -fopenmp
BUILD_DIR := build
SRC_DIR   := src

# Requerimietnos de objects
COMMON_OBJS := $(BUILD_DIR)/sequential_algorithms.o
PARALLEL_OBJS := $(BUILD_DIR)/parallel_algorithms.o $(COMMON_OBJS)

TARGETS := $(BUILD_DIR)/parallel_main \
           $(BUILD_DIR)/sequential_main \
           $(BUILD_DIR)/tests_main

.PHONY: all clean

all: $(TARGETS)

# Linking
$(BUILD_DIR)/parallel_main: $(BUILD_DIR)/parallel_main.o $(PARALLEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/sequential_main: $(BUILD_DIR)/sequential_main.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/tests_main: $(BUILD_DIR)/tests_main.o $(PARALLEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Regla para objects
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
