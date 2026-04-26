#pragma once

#include <iostream>
#include <iomanip>

class MatrixView {
private:
    double* array;
    size_t n;
    size_t stride; // offset para saltar a siguiente fila

public:
    MatrixView(double* data, size_t size, size_t offset)
        : array(data), n(size), stride(offset) {}

    inline double& operator()(size_t i, size_t j) {
        return array[i * stride + j];
    }

    inline const double& operator()(size_t i, size_t j) const {
        return array[i * stride + j];
    }

    size_t size() const { return n; }

    // Extraer cuadrantes
    // cuadrantes: 0 = top-left, 1 = top-right, 2 = bottom-left, 3 = bottom-right
    MatrixView get_quadrant(int quad) const {
        size_t row_offset = (quad == 2 || quad == 3) ? n/2 : 0;
        size_t col_offset = (quad == 1 || quad == 3) ? n/2 : 0;
        double* data_ptr = array + (row_offset * stride) + col_offset;
        return MatrixView(data_ptr, n/2, stride);
    }
};

class Matrix {
private:
    double *array;
    size_t n;

public:
    Matrix(size_t size) : array(new double[size * size]()), n(size) {};

    inline double& operator()(size_t i, size_t j) {
        return array[i*n + j];
    }
    inline const double& operator()(size_t i, size_t j) const {
        return array[i*n + j];
    }

    size_t size() const { return n; };
    double* data() const { return array; };
    MatrixView view() const { return MatrixView(array, n, n); }
    MatrixView get_quadrant(int quad) const {
        size_t half = n / 2;
        size_t row_offset = (quad == 2 || quad == 3) ? n/2 : 0;
        size_t col_offset = (quad == 1 || quad == 3) ? n/2 : 0;
        return MatrixView(array + (row_offset * n) + col_offset, half, n);
    };

    void print() {
        for (size_t i = 0; i < n; i++) {
            std::cout << "[ ";
            for (size_t j = 0; j < n; j++) {
                std::cout << std::setw(3) << (*this)(i, j) << " ";
            }
            std::cout << " ]\n";
        }
        std::cout << std::endl;
    }
};
