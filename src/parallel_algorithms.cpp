#include "parallel_algorithms.hpp"
#include "matrix.hpp"
#include "sequential_algorithms.hpp"


// TODO para matías: Esto probablemente sea lo que tenías que hacer tú.
void multiply_block_gemm(size_t I, size_t J, size_t b, const Matrix& A, const Matrix& B, Matrix& C) {
    size_t n = A.size();
    size_t n_blocks = (n + b - 1) / b;

    // Dimensiones del bloque, considerando que no estamos haciendo padding.
    size_t i0 = I * b;
    size_t j0 = J * b;
    size_t i1 = std::min(i0 + b, n);
    size_t j1 = std::min(j0 + b, n);

    // Para todo lo que viene también se podrían usar views, pero habría que usar padding, creo.
    // MatrixView blockC = MatrixView(C.data() + (I * b * n) + (J * b), b, n);


    for (size_t i = i0; i < i1 ; ++i) {
        for (size_t j = j0; j < j1; ++j) {
            C(i, j) = 0;
        }
    }

    for (size_t K = 0; K < n_blocks; ++K) { // Recorremos por bloques
        size_t k0 = K * b;
        size_t k1 = std::min(k0 + b, n);
                                            //
        //MatrixView blockA = MatrixView(A.data() + (I * b * n) + (K * b), b, n);
        //MatrixView blockB = MatrixView(B.data() + (K * b * n) + (J * b), b, n);

        for (size_t i = i0; i < i1; ++i) {
            for (size_t k = k0; k < k1; ++k) { // Esto reordena los loops; ahora acumulamos sobre C
                double tempA = A(i, k);  
                for (size_t j = j0; j < j1; ++j) {
                    C(i, j) += tempA * B(k, j);
                }
            }
        }
    }


}

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



