#include <vector>
#include <functional>
#include <tuple> 

#include "matrix.hpp"
#include "parallel_algorithms.hpp"
#include "tests.hpp"

using namespace std;

int main() {
    using gemm_func = function<void(const Matrix&, const Matrix&, Matrix&)>;
    using func_tuple = tuple<gemm_func, string>;

    vector<func_tuple> tasks = {
        {
            [](const Matrix& A, const Matrix& B, Matrix& C) {
                sequential_strassen_multiply(A, B, C);
            },
            "Strassen secuencial"
        },
        {
            [](const Matrix& A, const Matrix& B, Matrix& C) {
                parallel_tiled_multiply(A, B, C, 32);
            },
            "Submatrices paralelo con bloques de 32x32"
        },
        {
            [](const Matrix& A, const Matrix& B, Matrix& C) {
                parallel_tiled_multiply(A, B, C, 64);
            },
            "Submatrices paralelo con bloques de 64x64"
        }
    };
    
    for (auto& task : tasks) { 
        run_unit_tests(get<0>(task), get<1>(task));
    }
    return 0;
}
