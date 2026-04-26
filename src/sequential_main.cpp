#include <iostream>

#include "sequential_algorithms.hpp"

using namespace std;

int main() {
    std::cout << "Hello World!\n";

    size_t n = 2;
    Matrix A(n), B(n), C(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            A(i, j) = B(i, j) = (i * n) + j + 1;
        }
    }

    cout << "* A:\n";
    A.print();
    cout << "* B:\n";
    B.print();

    sequential_classic_multiply(A, B, C);

    cout << "* C:\n";
    C.print();

    return 0;
}
