#pragma once

#include "matrix.hpp"


// Función de multiplicación general de matrices cache-friendly
void multiply_block_gemm(size_t I, size_t J, size_t b, const Matrix& A, const Matrix& B, const Matrix& C);

void parallel_tiled_multiply(const Matrix& A, const Matrix& B, Matrix& C, size_t b);

void parallel_strassen_multiply(const Matrix& A, const Matrix& B, Matrix &C);
void parallel_hybrid_multiply(const Matrix& A, const Matrix& B, Matrix &C, size_t b);
