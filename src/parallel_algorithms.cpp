#include "sequential_algorithms.hpp"
#include "parallel_algorithms.hpp"
#include "matrix.hpp"


// Paralelización en base a bloques de la matriz resultado, colapsamos los dos ciclos.
template <typename MatA, typename MatB, typename MatC>
void parallel_tiled_multiply(const MatA& A, const MatB& B, MatC& C, size_t b) {
    size_t n = A.size();
    size_t n_blocks = (n + b - 1) / b;
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t I = 0; I < n_blocks; ++I) {
        for (size_t J = 0; J < n_blocks; ++J) {
            multiply_block_gemm(I, J, b, A, B, C);
        }
    }
}

template <typename MatA, typename MatB, typename MatC>
void _parallel_strassen_multiply(const MatA& A, const MatB& B, MatC& C) {
    size_t n = A.size();
    if (n <= 32) { // n_0 distinto probablemente sea mejor.
        sequential_classic_multiply(A, B, C);
        return;
    }

    size_t m = n / 2;

    // Descomponemos en submatrices 
    MatrixView A11 = A.get_quadrant(0);
    MatrixView A12 = A.get_quadrant(1);
    MatrixView A21 = A.get_quadrant(2);
    MatrixView A22 = A.get_quadrant(3);

    MatrixView B11 = B.get_quadrant(0);
    MatrixView B12 = B.get_quadrant(1);
    MatrixView B21 = B.get_quadrant(2);
    MatrixView B22 = B.get_quadrant(3);

    // Preparamos submatrices (shared) para los siete productos
    Matrix M1(m), M2(m), M3(m), M4(m), M5(m), M6(m), M7(m);

    // M1 = (A11 + A22)(B11 + B22)
    #pragma omp task shared(M1)
    {
        Matrix T1_1(m), T1_2(m);
        sequential_add_matrices(A11, A22, T1_1);
        sequential_add_matrices(B11, B22, T1_2);
        _parallel_strassen_multiply(T1_1, T1_2, M1);
    }

    // M2 = (A21 + A22)B11
    #pragma omp task shared(M2)
    {
        Matrix T2(m);
        sequential_add_matrices(A21, A22, T2);
        _parallel_strassen_multiply(T2, B11, M2);
    }

    // M3 = A11(B12 - B22)
    #pragma omp task shared(M3)
    {
        Matrix T3(m);
        sequential_sub_matrices(B12, B22, T3);
        _parallel_strassen_multiply(A11, T3, M3);
    }

    // M4 = A22(B21 - B11)
    #pragma omp task shared(M4)
    {
        Matrix T4(m);
        sequential_sub_matrices(B21, B11, T4);
        _parallel_strassen_multiply(A22, T4, M4);
    }

    // M5 = (A11 + A12)B22
    #pragma omp task shared(M5)
    {
        Matrix T5(m);
        sequential_add_matrices(A11, A12, T5);
        _parallel_strassen_multiply(T5, B22, M5);
    }

    // M6 = (A21 - A11)(B11 + B12)
    #pragma omp task shared(M6)
    {
        Matrix T6_1(m), T6_2(m);
        sequential_sub_matrices(A21, A11, T6_1);
        sequential_add_matrices(B11, B12, T6_2);
        _parallel_strassen_multiply(T6_1, T6_2, M6);
    }

    // M7 = (A12 - A22)(B21 + B22)
    #pragma omp task shared(M7)
    {
        Matrix T7_1(m), T7_2(m);
        sequential_sub_matrices(A12, A22, T7_1);
        sequential_add_matrices(B21, B22, T7_2);
        _parallel_strassen_multiply(T7_1, T7_2, M7);
    }

    // Sincronización
    #pragma omp taskwait

    // Combinamos los resultados en los cuadrantes de C
    MatrixView C11 = C.get_quadrant(0);
    MatrixView C12 = C.get_quadrant(1);
    MatrixView C21 = C.get_quadrant(2);
    MatrixView C22 = C.get_quadrant(3);

    // C11 = M1 + M4 - M5 + M7
    sequential_add_matrices(M1, M4, C11);
    sequential_sub_matrices(C11, M5, C11);
    sequential_add_matrices(C11, M7, C11);

    // C12 = M3 + M5
    sequential_add_matrices(M3, M5, C12);

    // C21 = M2 + M4
    sequential_add_matrices(M2, M4, C21);

    // C22 = M1 - M2 + M3 + M6
    sequential_sub_matrices(M1, M2, C22);
    sequential_add_matrices(C22, M3, C22);
    sequential_add_matrices(C22, M6, C22);
}

template <typename MatA, typename MatB, typename MatC>
void parallel_strassen_multiply(const MatA& A, const MatB& B, MatC& C) {
    #pragma omp parallel
    {
        #pragma omp single
        {
            _parallel_strassen_multiply(A, B, C);
        }
    }
}

template <typename MatA, typename MatB, typename MatC>
void _parallel_hybrid_multiply(const MatA& A, const MatB& B, MatC& C, size_t b) {
    size_t n = A.size();
    if (n <= 256) { // n_0 distinto probablemente sea mejor.
        parallel_tiled_multiply(A, B, C, b);
        return;
    }

    size_t m = n / 2;

    // Descomponemos en submatrices 
    MatrixView A11 = A.get_quadrant(0);
    MatrixView A12 = A.get_quadrant(1);
    MatrixView A21 = A.get_quadrant(2);
    MatrixView A22 = A.get_quadrant(3);

    MatrixView B11 = B.get_quadrant(0);
    MatrixView B12 = B.get_quadrant(1);
    MatrixView B21 = B.get_quadrant(2);
    MatrixView B22 = B.get_quadrant(3);

    // Preparamos submatrices (shared) para los siete productos
    Matrix M1(m), M2(m), M3(m), M4(m), M5(m), M6(m), M7(m);

    // M1 = (A11 + A22)(B11 + B22)
    #pragma omp task shared(M1)
    {
        Matrix T1_1(m), T1_2(m);
        sequential_add_matrices(A11, A22, T1_1);
        sequential_add_matrices(B11, B22, T1_2);
        _parallel_strassen_multiply(T1_1, T1_2, M1);
    }

    // M2 = (A21 + A22)B11
    #pragma omp task shared(M2)
    {
        Matrix T2(m);
        sequential_add_matrices(A21, A22, T2);
        _parallel_strassen_multiply(T2, B11, M2);
    }

    // M3 = A11(B12 - B22)
    #pragma omp task shared(M3)
    {
        Matrix T3(m);
        sequential_sub_matrices(B12, B22, T3);
        _parallel_strassen_multiply(A11, T3, M3);
    }

    // M4 = A22(B21 - B11)
    #pragma omp task shared(M4)
    {
        Matrix T4(m);
        sequential_sub_matrices(B21, B11, T4);
        _parallel_strassen_multiply(A22, T4, M4);
    }

    // M5 = (A11 + A12)B22
    #pragma omp task shared(M5)
    {
        Matrix T5(m);
        sequential_add_matrices(A11, A12, T5);
        _parallel_strassen_multiply(T5, B22, M5);
    }

    // M6 = (A21 - A11)(B11 + B12)
    #pragma omp task shared(M6)
    {
        Matrix T6_1(m), T6_2(m);
        sequential_sub_matrices(A21, A11, T6_1);
        sequential_add_matrices(B11, B12, T6_2);
        _parallel_strassen_multiply(T6_1, T6_2, M6);
    }

    // M7 = (A12 - A22)(B21 + B22)
    #pragma omp task shared(M7)
    {
        Matrix T7_1(m), T7_2(m);
        sequential_sub_matrices(A12, A22, T7_1);
        sequential_add_matrices(B21, B22, T7_2);
        _parallel_strassen_multiply(T7_1, T7_2, M7);
    }

    // Sincronización
    #pragma omp taskwait

    // Combinamos los resultados en los cuadrantes de C
    MatrixView C11 = C.get_quadrant(0);
    MatrixView C12 = C.get_quadrant(1);
    MatrixView C21 = C.get_quadrant(2);
    MatrixView C22 = C.get_quadrant(3);

    // C11 = M1 + M4 - M5 + M7
    sequential_add_matrices(M1, M4, C11);
    sequential_sub_matrices(C11, M5, C11);
    sequential_add_matrices(C11, M7, C11);

    // C12 = M3 + M5
    sequential_add_matrices(M3, M5, C12);

    // C21 = M2 + M4
    sequential_add_matrices(M2, M4, C21);

    // C22 = M1 - M2 + M3 + M6
    sequential_sub_matrices(M1, M2, C22);
    sequential_add_matrices(C22, M3, C22);
    sequential_add_matrices(C22, M6, C22);
}

template <typename MatA, typename MatB, typename MatC>
void parallel_hybrid_multiply(const MatA& A, const MatB& B, MatC& C, size_t b) {
    #pragma omp parallel
    {
        #pragma omp single
        {
            _parallel_hybrid_multiply(A, B, C, b);
        }
    }
}

// instanciaciones explícitas de los template
// Strassen (8 combinaciones)
template void parallel_strassen_multiply<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&);
template void parallel_strassen_multiply<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&);
template void parallel_strassen_multiply<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&);
template void parallel_strassen_multiply<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&);
template void parallel_strassen_multiply<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&);
template void parallel_strassen_multiply<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&);
template void parallel_strassen_multiply<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&);
template void parallel_strassen_multiply<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&);

template void _parallel_strassen_multiply<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&);
template void _parallel_strassen_multiply<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&);
template void _parallel_strassen_multiply<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&);
template void _parallel_strassen_multiply<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&);
template void _parallel_strassen_multiply<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&);
template void _parallel_strassen_multiply<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&);
template void _parallel_strassen_multiply<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&);
template void _parallel_strassen_multiply<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&);

// Por bloque, partición por bloques
template void parallel_tiled_multiply<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&, size_t);
template void parallel_tiled_multiply<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&, size_t);
template void parallel_tiled_multiply<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&, size_t);
template void parallel_tiled_multiply<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&, size_t);
template void parallel_tiled_multiply<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&, size_t);
template void parallel_tiled_multiply<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&, size_t);
template void parallel_tiled_multiply<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&, size_t);
template void parallel_tiled_multiply<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&, size_t);

// Hibrido
template void parallel_hybrid_multiply<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&, size_t);
template void parallel_hybrid_multiply<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&, size_t);
template void parallel_hybrid_multiply<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&, size_t);
template void parallel_hybrid_multiply<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&, size_t);
template void parallel_hybrid_multiply<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&, size_t);
template void parallel_hybrid_multiply<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&, size_t);
template void parallel_hybrid_multiply<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&, size_t);
template void parallel_hybrid_multiply<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&, size_t);

template void _parallel_hybrid_multiply<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&, size_t);
template void _parallel_hybrid_multiply<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&, size_t);
template void _parallel_hybrid_multiply<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&, size_t);
template void _parallel_hybrid_multiply<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&, size_t);
template void _parallel_hybrid_multiply<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&, size_t);
template void _parallel_hybrid_multiply<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&, size_t);
template void _parallel_hybrid_multiply<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&, size_t);
template void _parallel_hybrid_multiply<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&, size_t);