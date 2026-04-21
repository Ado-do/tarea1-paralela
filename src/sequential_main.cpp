#include <iomanip>
#include <iostream>

#include "sequential_algorithms.hpp"

using namespace std;

void printMatrix(const string& name, Matrix& M, size_t n) {
    cout << "* " << name << ":\n";
    for (size_t i = 0; i < n; i++) {
        cout << "[ ";
        for (size_t j = 0; j < n; j++) {
            cout << setw(3) << M(i, j) << " ";
        }
        cout << " ]\n";
    }
    cout << endl;
}

int main() {
    std::cout << "Hello World!\n";

    size_t n = 2;
    Matrix A(n), B(n), C(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            A(i, j) = B(i, j) = (i * n) + j + 1;
        }
    }

    printMatrix("A", A, n);
    printMatrix("B", B, n);

    sequential_classic_multiply(A, B, C);

    printMatrix("C", C, n);

    return 0;
}
