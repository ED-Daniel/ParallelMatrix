#include <iostream>
#include <chrono>
#include <fstream>
#include "Matrix.h"

using namespace std::chrono;

int main() {
//    auto A = Matrix<double>(&std::cin);
//    A.print();
//    auto B = !A;
//    auto C = B;
//    auto test = std::move(A);
//    C.setElement(-300, 1, 1);
//    B.print();
//    C.print();
//    test.print();
//
//    auto test2 = Matrix<double>({{{0, 0}, 1}, {{0, 1}, 2}, {{1, 0}, 3}, {{7, 8}, 324}});
//    test2.print();

    auto test10 = Matrix<double>(R"(C:\Users\MSI GS66 Stealth\CLionProjects\ParallelMatrix\test10.txt)");
    std::cout << "READ 10" << std::endl;
    auto test1000 = Matrix<double>(R"(C:\Users\MSI GS66 Stealth\CLionProjects\ParallelMatrix\test1000.txt)");
    std::cout << "READ 1000" << std::endl;
    auto test100 = Matrix<double>(R"(C:\Users\MSI GS66 Stealth\CLionProjects\ParallelMatrix\test100.txt)");
    std::cout << "READ 100" << std::endl;

    std::ofstream out("Output.txt");
    if (out.is_open()) {
        out << "ASYNC VERSION" << std::endl;

        test10.setThreads(2);
        test1000.setThreads(2);
        test100.setThreads(2);

        out << "Test 10: ";
        auto start = std::chrono::high_resolution_clock::now();
        test10.threadMultiplication(test10);
        auto end = std::chrono::high_resolution_clock::now();
        out << std::chrono::duration_cast<std::chrono::microseconds>(end - start) << std::endl;

        out << "Test 1000: ";
        start = std::chrono::high_resolution_clock::now();
        test1000.threadMultiplication(test1000);
        end = std::chrono::high_resolution_clock::now();
        out << std::chrono::duration_cast<std::chrono::milliseconds>(end - start) << std::endl;

        out << "Test 100: ";
        start = std::chrono::high_resolution_clock::now();
        test100 * test100;
        end = std::chrono::high_resolution_clock::now();
        out << std::chrono::duration_cast<std::chrono::milliseconds>(end - start) << std::endl;
        out.close();
    }
}
