#include "sequential_algorithms.hpp"
#include "parallel_algorithms.hpp"
#include "matrix.hpp"


// Paralelización en base a bloques de la matriz resultado, colapsamos los dos ciclos.
void parallel_tiled_multiply(const Matrix& A, const Matrix& B, Matrix &C, size_t b) {
    size_t n = A.size();
    size_t n_blocks = (n + b - 1) / b;
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t I = 0; I < n_blocks; ++I) {
        for (size_t J = 0; J < n_blocks; ++J) {
            multiply_block_gemm(I, J, b, A, B, C);
        }
    }
}

void parallel_strassen_multiply(const Matrix& A, const Matrix& B, Matrix &C) {
    // implementar
}

void parallel_hybrid_multiply(const Matrix& A, const Matrix& B, Matrix &C, size_t b) {
    // implemetar
}



