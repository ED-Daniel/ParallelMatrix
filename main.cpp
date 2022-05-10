#include <iostream>
#include <chrono>
#include "Matrix.h"

using namespace std::chrono;

int main() {
    Matrix<double> A = Matrix<double>(&std::cin);
    Matrix<double> B = Matrix<double>(&std::cin);

    (A - B).print();

    Matrix<int> E = createIdentityMatrix<int>(4, 4);
    E.print();
    std::cout << (E == 1) << std::endl;
    std::cout << (Matrix<int>(3, 3) == 0) << std::endl;

    (!A).print();
    (A * !A).print();
    std::cout << (A * !A == 1) << std::endl;

    return 0;
}
