#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "sequential_algorithms.hpp"

using namespace std;

const string CSV_PATH = "sequential.csv";
const size_t BLOCK_SIZE = 32;
const size_t TEST_SIZE = 2048;
const int REPETITIONS = 4;

void run_sequential_experiment(size_t n, ofstream &csv) {
    int reps = ((n >= 2048)? 1 : REPETITIONS);

    cout << "* Experimento: n = " << n << ", b = " << BLOCK_SIZE << ", reps = " << reps << endl;
    csv << n << ',';

    // Inicializar matrices
    Matrix A(n), B(n), C(n);
    static random_device rd;
    static mt19937_64 gen(rd());
    static uniform_real_distribution<double> dis(0.0, 1.0);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            A(i, j) = dis(gen);
            B(i, j) = dis(gen);
        }
    }

    // 1. Clásica
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < reps; ++i)
        sequential_classic_multiply(A, B, C);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = (end - start) / reps;
    csv << diff.count() << ',';
    cout << "  Clásica:\t\t" << diff.count() << "s" << endl;

    // 2. Cache-friendly
    start = chrono::high_resolution_clock::now();
    for(int i = 0; i < reps; ++i) sequential_cachefriendly_multiply(A, B, C, BLOCK_SIZE);
    end = chrono::high_resolution_clock::now();
    diff = (end - start) / reps;
    csv << diff.count() << ',';
    cout << "  Cache-friendly:\t" << diff.count() << "s" << endl;

    // 3. Strassen
    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < reps; ++i) sequential_strassen_multiply(A, B, C);
    end = chrono::high_resolution_clock::now();
    diff = (end - start) / reps;
    csv << diff.count() << '\n';
    cout << "  Strassen:\t\t" << diff.count() << "s" << endl;
}

int main() {
    vector<size_t> sizes = {128, 256, 512, 1024, 2048, 4096};

    ofstream csv(CSV_PATH);
    csv << "n,classic,cache,strassen\n";

    cout << "** EJECUTANDO EXPERIMENTOS DE ALGORITMOS SECUENCIALES:\n";
    for (size_t n : sizes) run_sequential_experiment(n, csv);

    return 0;
}
