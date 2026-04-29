#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <sys/prctl.h> // para manejar counters de perf

#include "sequential_algorithms.hpp"

using namespace std;

const string CSV_PATH = "sequential.csv";
const size_t BLOCK_SIZE = 128;
const size_t TEST_SIZE = 2048;
const int REPETITIONS = 5;
const int N0 = 32;

void run_sequential_experiment(size_t n, ofstream &csv) {
    // int reps = ((n >= 2048)? 1 : REPETITIONS);
    int reps = REPETITIONS;

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
    //for (int i = 0; i < reps; ++i) sequential_strassen_multiply(A, B, C, N0);
    end = chrono::high_resolution_clock::now();
    diff = (end - start) / reps;
    csv << diff.count() << '\n';
    cout << "  Strassen:\t\t" << diff.count() << "s" << endl;
}

int main(int argc, char* argv[]) {
    // Modo de ejecución único (./sequential_main [n] [classic | cache | strassen])
    if (argc == 3) {
        size_t n = stoul(argv[1]);
        string algo = argv[2];

        Matrix A(n), B(n), C(n);
        random_device rd;
        mt19937_64 gen(rd());
        uniform_real_distribution<double> dis(0.0, 1.0);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                A(i, j) = dis(gen); B(i, j) = dis(gen);
            }
        }

        cout << "** Experimento targetado: n=" << n << " algo=" << algo << endl;
        
        // Permitimos counters de perf
        if (prctl(PR_TASK_PERF_EVENTS_ENABLE) == -1) {
            perror("prctl enable fallado");
        }
        
        if (algo == "classic") {
            for(int i=0; i < REPETITIONS; ++i) sequential_classic_multiply(A, B, C);
        } else if (algo == "cache") {
            for(int i=0; i < REPETITIONS; ++i) sequential_cachefriendly_multiply(A, B, C, BLOCK_SIZE);
        } else if (algo == "strassen") {
            for(int i=0; i < REPETITIONS; ++i) sequential_strassen_multiply(A, B, C, N0);
        }

        // Detenemos counters de perf
        if (prctl(PR_TASK_PERF_EVENTS_DISABLE) == -1) {
            perror("prctl disable fallado");
        }
                
        return 0;
    }



    vector<size_t> sizes = {128, 256, 512, 1024, 2048, 4096};
    //vector<size_t> sizes = {4096};

    ofstream csv(CSV_PATH);
    csv << "n,classic,cache,strassen\n";

    cout << "** EJECUTANDO EXPERIMENTOS DE ALGORITMOS SECUENCIALES:\n";
    for (size_t n : sizes) run_sequential_experiment(n, csv);

    return 0;
}
