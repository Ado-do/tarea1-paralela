#include <iostream>
#include <random>
#include <fstream>
#include <omp.h>

#include "parallel_algorithms.hpp"

using namespace std;

// Nombre del CSV
const string csv_path = "parallel.csv";
// Tamaño de los bloques para algoritmo cache-friendly
const size_t b = 32;

void run_parallel_experiment(size_t n, int rep, int p, ofstream &csv) {
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
    for (int i = 0; i < rep; i++) {
        double start = omp_get_wtime();
        parallel_tiled_multiply(A, B, C, b);
        double end = omp_get_wtime();
        total_time += (end - start);
    }
    double tiled_time = total_time / rep;
    csv << tiled_time << ',';
    cout << "  Tiled:\t\t" << tiled_time << "s" << endl;

    // 2. Strassen
    total_time = 0;
    for (int i = 0; i < rep; i++) {
        double start = omp_get_wtime();
        parallel_strassen_multiply(A, B, C);
        double end = omp_get_wtime();
        total_time += (end - start);
    }
    double strassen_time = total_time / rep;
    csv << strassen_time << ',';
    cout << "  Strassen:\t\t" << strassen_time << "s" << endl;

    // 3. Hybrid
    total_time = 0;
    for (int i = 0; i < rep; i++) {
        double start = omp_get_wtime();
        parallel_hybrid_multiply(A, B, C, b);
        double end = omp_get_wtime();
        total_time += (end - start);
    }
    double hybrid_time = total_time / rep;
    csv << hybrid_time << ',';
    cout << "  Hybrid:\t\t" << hybrid_time << "s" << endl;
}

int main() {
    vector<size_t> sizes = {256, 512, 1024, 2048, 4096};
    vector<int> threads = {1, 2, 4, 8};
    int repetitions = 3;

    ofstream csv(csv_path);
    csv << "n,p,tiled,strassen,hybrid\n";

    cout << "** EJECUTANDO EXPERIMENTOS DE ALGORITMOS PARALELOS:\n";
    //for (size_t n : sizes) run_sequential_experiment(n, repetitions, csv);
    run_parallel_experiment(1024, repetitions, 2, csv);

    return 0;
}

// int main() {
//     cout << "Hello World!\n";
//
//     size_t n = 2;
//     Matrix A(n), B(n), C(n);
//     for (size_t i = 0; i < n; i++) {
//         for (size_t j = 0; j < n; j++) {
//             A(i, j) = B(i, j) = (i * n) + j + 1;
//         }
//     }
//
//     cout << "* A:\n";
//     A.print();
//     cout << "* B:\n";
//     B.print();
//
//     // Probar algoritmos
//     //parallel_tiled_multiply(A, B, C);
//
//     cout << "* C:\n";
//     C.print();
//
//     return 0;
// }
