#pragma once

#include "matrix.hpp"
#include <omp.h>

template <typename MatA, typename MatB, typename MatC>
void parallel_add_matrices(const MatA& A, const MatB& B, MatC& C) {
    const size_t n = A.size();
    
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C(i, j) = A(i, j) + B(i, j);
        }
    }
}

template <typename MatA, typename MatB, typename MatC>
void parallel_sub_matrices(const MatA& A, const MatB& B, MatC& C) {
    const size_t n = A.size();
    
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C(i, j) = A(i, j) - B(i, j);
        }
    }
}

template <typename MatA, typename MatB, typename MatC>
void parallel_tiled_multiply(const MatA& A, const MatB& B, MatC& C, size_t b);

template <typename MatA, typename MatB, typename MatC>
void parallel_strassen_multiply(const MatA& A, const MatB& B, MatC& C, size_t n0);

template <typename MatA, typename MatB, typename MatC>
void parallel_hybrid_multiply(const MatA& A, const MatB& B, MatC& C, size_t b, size_t n0);

// // Suma de matrices

// template void parallel_add_matrices<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&);
// template void parallel_add_matrices<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&);
// template void parallel_add_matrices<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&);
// template void parallel_add_matrices<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&);
// template void parallel_add_matrices<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&);
// template void parallel_add_matrices<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&);
// template void parallel_add_matrices<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&);
// template void parallel_add_matrices<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&);

// // Resta de matrices

// template void parallel_sub_matrices<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&);
// template void parallel_sub_matrices<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&);
// template void parallel_sub_matrices<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&);
// template void parallel_sub_matrices<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&);
// template void parallel_sub_matrices<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&);
// template void parallel_sub_matrices<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&);
// template void parallel_sub_matrices<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&);
// template void parallel_sub_matrices<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&);
