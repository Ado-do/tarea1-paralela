#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "sequential_algorithms.hpp"

using namespace std;

// Nombre del CSV
const string csv_path = "sequential.csv";
// Tamaño de los bloques para algoritmo cache-friendly
const size_t b = 32;


void run_sequential_experiment(size_t n, int repetitions, ofstream &csv) {
    cout << "* Experimento con n = " << n << endl;
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
    for(int i = 0; i < repetitions; ++i) sequential_classic_multiply(A, B, C);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = (end - start) / repetitions;
    csv <<  diff.count() << ',';
    cout << "  Clásica:\t\t" << diff.count() << "s" << endl;

    // 2. Cache-friendly
    start = chrono::high_resolution_clock::now();
    for(int i = 0; i < repetitions; ++i) sequential_cachefriendly_multiply(A, B, C, b);
    end = chrono::high_resolution_clock::now();
    diff = (end - start) / repetitions;
    csv << diff.count() << ',';
    cout << "  Cache-friendly:\t" << diff.count() << "s" << endl;

    // 3. Strassen
    start = chrono::high_resolution_clock::now();
    for(int i = 0; i < repetitions; ++i) sequential_strassen_multiply(A, B, C);
    end = chrono::high_resolution_clock::now();
    diff = (end - start) / repetitions;
    csv << diff.count() << '\n';
    cout << "  Strassen:\t\t" << diff.count() << "s" << endl;
}

int main() {
    vector<size_t> sizes = {256, 512, 1024, 2048, 4096};
    int repetitions = 3;

    ofstream csv(csv_path);
    csv << "n,classic,cache,strassen\n";


    cout << "** EJECUTANDO EXPERIMENTOS DE ALGORITMOS SECUENCIALES:\n";
    //for (size_t n : sizes) run_sequential_experiment(n, repetitions, csv);
    run_sequential_experiment(32, repetitions, csv);

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
//     sequential_classic_multiply(A, B, C);
//
//     cout << "* C:\n";
//     C.print();
//
//     return 0;
// }
