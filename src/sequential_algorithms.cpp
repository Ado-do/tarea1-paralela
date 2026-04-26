#include "sequential_algorithms.hpp"



template <typename MatA, typename MatB, typename MatC>
void sequential_classic_multiply(const MatA& A, const MatB& B, MatC& C) {
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

template <typename MatA, typename MatB, typename MatC>
void sequential_strassen_multiply(const MatA& A, const MatB &B, MatC &C) {
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

    // Preparamos submatrices para los siete productos
    Matrix T1(m), T2(m);
    
    // M1 = (A11 + A22)(B11 + B22)
    sequential_add_matrices(A11, A22, T1);
    sequential_add_matrices(B11, B22, T2);
    Matrix M1(m); sequential_strassen_multiply(T1, T2, M1);

    // M2 = (A21 + A22)B11
    sequential_add_matrices(A21, A22, T1);
    Matrix M2(m); sequential_strassen_multiply(T1, B11, M2);

    // M3 = A11(B12 - B22)
    sequential_sub_matrices(B12, B22, T2);
    Matrix M3(m); sequential_strassen_multiply(A11, T2, M3);

    // M4 = A22(B21 - B11)
    sequential_sub_matrices(B21, B11, T2);
    Matrix M4(m); sequential_strassen_multiply(A22, T2, M4);

    // M5 = (A11 + A12)B22
    sequential_add_matrices(A11, A12, T1);
    Matrix M5(m); sequential_strassen_multiply(T1, B22, M5);

    // M6 = (A21 - A11)(B11 + B12)
    sequential_sub_matrices(A21, A11, T1);
    sequential_add_matrices(B11, B12, T2);
    Matrix M6(m); sequential_strassen_multiply(T1, T2, M6);

    // M7 = (A12 - A22)(B21 + B22)
    sequential_sub_matrices(A12, A22, T1);
    sequential_add_matrices(B21, B22, T2);
    Matrix M7(m); sequential_strassen_multiply(T1, T2, M7);

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


// instanciaciones explícitas de los template
// Strassen (8 combinaciones)
template void sequential_strassen_multiply<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&);
template void sequential_strassen_multiply<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&);
template void sequential_strassen_multiply<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&);
template void sequential_strassen_multiply<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&);
template void sequential_strassen_multiply<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&);
template void sequential_strassen_multiply<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&);
template void sequential_strassen_multiply<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&);
template void sequential_strassen_multiply<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&);

// Clásica (8 combinaciones)
template void sequential_classic_multiply<Matrix, Matrix, Matrix>(const Matrix&, const Matrix&, Matrix&);
template void sequential_classic_multiply<MatrixView, Matrix, Matrix>(const MatrixView&, const Matrix&, Matrix&);
template void sequential_classic_multiply<Matrix, MatrixView, Matrix>(const Matrix&, const MatrixView&, Matrix&);
template void sequential_classic_multiply<MatrixView, MatrixView, Matrix>(const MatrixView&, const MatrixView&, Matrix&);
template void sequential_classic_multiply<Matrix, Matrix, MatrixView>(const Matrix&, const Matrix&, MatrixView&);
template void sequential_classic_multiply<MatrixView, Matrix, MatrixView>(const MatrixView&, const Matrix&, MatrixView&);
template void sequential_classic_multiply<Matrix, MatrixView, MatrixView>(const Matrix&, const MatrixView&, MatrixView&);
template void sequential_classic_multiply<MatrixView, MatrixView, MatrixView>(const MatrixView&, const MatrixView&, MatrixView&);
