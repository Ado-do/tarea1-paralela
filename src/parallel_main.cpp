#include <iostream>
#include <random>
#include <fstream>
#include <omp.h>

#include "parallel_algorithms.hpp"

using namespace std;

const string CSV_PATH = "parallel.csv";
const size_t BLOCK_SIZE = 128;
const size_t TEST_SIZE = 2048;
const int REPETITIONS = 5;
const int N0 = 32;

void run_parallel_experiment(size_t n, int p, ofstream &csv) {
    int reps = REPETITIONS;
    cout << "* Experimento con n = " << n << " y p = " << p << endl;
    csv << n << ',' << p << ',';
    omp_set_num_threads(p);

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

    // 1. Tiled
    double total_time = 0;
    for (int i = 0; i < reps; i++) {
        double start = omp_get_wtime();
        parallel_tiled_multiply(A, B, C, BLOCK_SIZE);
        double end = omp_get_wtime();
        total_time += (end - start);
    }
    double tiled_time = total_time / reps;
    csv << tiled_time << ',';
    cout << "  Tiled:\t\t" << tiled_time << "s" << endl;

    // 2. Strassen
    total_time = 0;
    for (int i = 0; i < reps; i++) {
        double start = omp_get_wtime();
        parallel_strassen_multiply(A, B, C, N0);
        double end = omp_get_wtime();
        total_time += (end - start);
    }
    double strassen_time = total_time / reps;
    csv << strassen_time << ',';
    cout << "  Strassen:\t\t" << strassen_time << "s" << endl;

    // 3. Hybrid
    total_time = 0;
    for (int i = 0; i < reps; i++) {
        double start = omp_get_wtime();
        parallel_hybrid_multiply(A, B, C, BLOCK_SIZE, N0);
        double end = omp_get_wtime();
        total_time += (end - start);
    }
    double hybrid_time = total_time / reps;
    csv << hybrid_time << ',';
    cout << "  Hybrid:\t\t" << hybrid_time << "s" << endl;
}

int main() {
    vector<size_t> sizes = {256, 512, 1024, 2048, 4096};
    vector<int> threads = {1, 2, 4, 8};

    ofstream csv(CSV_PATH);
    csv << "n,p,tiled,strassen,hybrid\n";

    cout << "** EJECUTANDO EXPERIMENTOS DE ALGORITMOS PARALELOS:\n";
    
    for (size_t n : sizes) {
        for (int p : threads) {
            run_parallel_experiment(n, p, csv);
        }
    }

    return 0;
}
