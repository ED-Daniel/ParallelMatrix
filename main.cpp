#include <iostream>
#include <chrono>
#include "Matrix.h"

using namespace std::chrono;

int main() {
    auto A = Matrix<double>(&std::cin);
    A.print();
    auto B = !A;
    auto C = B;
    auto test = std::move(A);
    C.setElement(-300, 1, 1);
    B.print();
    C.print();
    test.print();

    auto test2 = Matrix<double>({{{0, 0}, 1}, {{0, 1}, 2}, {{1, 0}, 3}, {{7, 8}, 324}});
    test2.print();
}
