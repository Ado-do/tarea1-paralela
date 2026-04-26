#pragma once

#include <cmath>
#include <string>
#include <iostream>
#include <random>
#include <functional>
#include <iomanip>
#include "sequential_algorithms.hpp"

#include "matrix.hpp"
#include "parallel_algorithms.hpp"
#include "parallel_algorithms.hpp"

// Comparación de matrices de floats con epsilon de tolerancia
inline bool compare_matrices(const Matrix& A, const Matrix& B, double epsilon = 1e-9) {
    if (A.size() != B.size()) return false;
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (std::abs(A(i, j) - B(i, j)) > epsilon) {
                return false;
            }
        }
    }
    return true;
}

template <typename Func>
void run_unit_tests(Func to_test, const std::string& name) {
    std::cout << "### Comenzando tests unitarios para: " << name << " ###\n";

    // Función lambda para test case
    auto test_case = [&](size_t n, const std::string& case_name, bool identity_A = false, bool zero_B = false) {
        Matrix A(n), B(n), C_ref(n), C_test(n);
        
        // Setup inputs
        if (identity_A) {
            for (size_t i = 0; i < n; ++i) A(i, i) = 1.0;
        } else {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(-10.0, 10.0);
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    A(i, j) = dis(gen);
                }
            }
        }

        if (!zero_B) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(-10.0, 10.0);
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    B(i, j) = dis(gen);
                }
            }
        }

        // Corremos baseline y target
        sequential_classic_multiply(A, B, C_ref);
        to_test(A, B, C_test);

        bool passed = compare_matrices(C_ref, C_test);
        std::cout << "[Test] n=" << n << " | " << std::left << std::setw(20) << case_name 
                  << " | Estado: " << (passed ? "LOGRADO" : "FALLADO") << "\n";
        return passed;
    };

    // Caso borde 1: 1x1 
    test_case(1, "1x1 Multiplicación escalar");

    // Caso borde 2: Matriz identidad (A * I = A)
    test_case(4, "Multiplicación identidad", true, false);

    // Caso borde 3: Matriz nula (A * 0 = 0)
    test_case(4, "Multiplicación nula", false, true);

    // Algunos casos estándar potencias de 2
    test_case(2, "Pequeña 2x2");
    test_case(8, "Pequeña 8x8");
    test_case(64, "Mediana 64x64");

    std::cout << "### Tests Completados ###\n\n";
}

