#pragma once

#include "matrix.hpp"


template <typename MatA, typename MatB, typename MatC>
void parallel_tiled_multiply(const MatA& A, const MatB& B, MatC& C, size_t b);

template <typename MatA, typename MatB, typename MatC>
void parallel_strassen_multiply(const MatA& A, const MatB& B, MatC& C);

template <typename MatA, typename MatB, typename MatC>
void parallel_hybrid_multiply(const MatA& A, const MatB& B, MatC& C, size_t b);
