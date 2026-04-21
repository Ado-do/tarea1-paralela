#pragma once

#include "matrix.hpp"

void parallel_tiled_multiply(const Matrix& A, const Matrix& B, Matrix &C, size_t b);
void parallel_strassen_multiply(const Matrix& A, const Matrix& B, Matrix &C);
void parallel_hybrid_multiply(const Matrix& A, const Matrix& B, Matrix &C, size_t b);
