#pragma once

#include "matrix.hpp"


// Voy a dejar esto acá por conveniencia, ya que es con template.
template <typename MatA, typename MatB, typename MatC>
void sequential_add_matrices(const MatA& A, const MatB& B, MatC& C) {
    const size_t n = A.size();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C(i, j) = A(i, j) + B(i, j);
        }
    }
}

template <typename MatA, typename MatB, typename MatC>
void sequential_sub_matrices(const MatA& A, const MatB& B, MatC& C) {
    const size_t n = A.size();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C(i, j) = A(i, j) - B(i, j);
        }
    }
}


// ----- Algoritmos de multiplicación de matrices -----
template <typename MatA, typename MatB, typename MatC>
void sequential_classic_multiply(const MatA& A, const MatB& B, MatC& C);

void sequential_cachefriendly_multiply(const Matrix& A, const Matrix &B, Matrix &C, size_t b);

template <typename MatA, typename MatB, typename MatC>
void sequential_strassen_multiply(const MatA& A, const MatB &B, MatC &C);
