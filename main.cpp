#include <iostream>
#include <chrono>
#include "Matrix.h"

using namespace std::chrono;

int main() {
    auto A = Matrix<double>(&std::cin);
    auto T = Matrix<double>(3, 3);

    inverseAsync<double>(A, T);
    T.print();
    (A * T).print();
    return 0;
}
