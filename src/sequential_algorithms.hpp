#pragma once

#include "matrix.hpp"

void sequential_classic_multiply(const Matrix& A, const Matrix &B, Matrix &C);
void sequential_cachefriendly_multiply(const Matrix& A, const Matrix &B, Matrix &C, size_t b);
void sequential_strassen_multiply(const Matrix& A, const Matrix &B, Matrix &C);
