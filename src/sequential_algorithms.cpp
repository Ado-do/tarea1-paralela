#include "sequential_algorithms.hpp"

void sequential_classic_multiply(const Matrix& A, const Matrix &B, Matrix &C) {
    size_t n = A.size();

    // iterar filas de A
    for (size_t i = 0; i < n; i++) {
        // iterar columnas de B
        for (size_t j = 0; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += A(i, k) * B(k, j);
            }
            C(i, j) = sum;
        }
    }
}

void sequential_cachefriendly_multiply(const Matrix& A, const Matrix &B, Matrix &C, size_t b) {
    // implementar
}

void sequential_strassen_multiply(const Matrix& A, const Matrix &B, Matrix &C) {
    // implementar
}
