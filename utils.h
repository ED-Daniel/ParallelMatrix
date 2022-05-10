//
// Created by MSI GS66 Stealth on 10.05.2022.
//

#ifndef PARALLELMATRIX_UTILS_H
#define PARALLELMATRIX_UTILS_H

#include <iostream>
#include <fstream>

template<typename T>
T readData(std::ifstream *in) {
    T read;
    *in >> read;

    if ((*in).fail()) {
        //std::cout << read << " ";
        (*in).clear();
        throw std::exception("Bad input");
    }

    return read;
}

#endif //PARALLELMATRIX_UTILS_H
