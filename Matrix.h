//
// Created by MSI GS66 Stealth on 08.05.2022.
//

#ifndef PARALLELMATRIX_MATRIX_H
#define PARALLELMATRIX_MATRIX_H

#include<vector>
#include<fstream>
#include<iostream>
#include<thread>
#include<functional>
#include<future>
#include<initializer_list>
#include "utils.h"

template<typename T>
class Matrix;

template<typename T>
static Matrix<T> createIdentityMatrix(size_t rows, size_t cols) {
    Matrix<T> result = Matrix<T>(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (i == j) result.setElement(1, i, j);
            else result.setElement(0, i, j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> operator* (const Matrix<T> & lhs, const Matrix<T> & rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::exception("Columns and rows must be the same");
    }

    Matrix<T> product = Matrix<T>(lhs.rows, rhs.cols);

    for (size_t i = 0; i < product.rows; ++i) {
        for (size_t j = 0; j < product.cols; ++j) {
            T local_sum = 0;
            for (size_t k = 0; k < lhs.cols; ++k) {
                local_sum += lhs.content[i][k] * rhs.content[k][j];
            }
            if (abs(local_sum) < 0.000000001) {
                local_sum = 0;
            }
            product.content[i][j] = local_sum;
        }
    }

    return product;
}

template<typename T>
Matrix<T> operator+ (const Matrix<T> & lhs, const Matrix<T> & rhs) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) throw std::exception("Matrices must be the same dimensions");
    Matrix<T> sum = Matrix<T>(lhs.rows, lhs.cols);
    for (size_t i = 0; i < lhs.rows; ++i) {
        for (size_t j = 0; j < lhs.cols; ++j) {
            sum.content[i][j] = lhs.content[i][j] + rhs.content[i][j];
        }
    }
    return sum;
}

template<typename T>
Matrix<T> operator- (const Matrix<T> & lhs, const Matrix<T> & rhs) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) throw std::exception("Matrices must be the same dimensions");
    Matrix<T> sum = Matrix<T>(lhs.rows, lhs.cols);
    for (size_t i = 0; i < lhs.rows; ++i) {
        for (size_t j = 0; j < lhs.cols; ++j) {
            sum.content[i][j] = lhs.content[i][j] - rhs.content[i][j];
        }
    }
    return sum;
}

template<typename T>
bool operator== (const Matrix<T> & lhs, const Matrix<T> & rhs) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) return false;
    bool trigger = true;
    for (size_t i = 0; i < lhs.rows; ++i) {
        for (size_t j = 0; j < lhs.cols; ++j) {
            if (lhs.content[i][j] != lhs.content[i][j]) {
                trigger = false;
                break;
            }
        }
    }
    return trigger;
}

template<typename T>
bool operator== (const Matrix<T> & lhs, int digit) {
    if (digit == 0 || digit == 1) {
        bool trigger = true;
        if (digit == 0) {
            return lhs == Matrix<T>(lhs.rows, lhs.cols);
        }
        if (digit == 1) {
            return lhs == createIdentityMatrix<T>(lhs.rows, lhs.cols);
        }
    } else throw std::exception("Failed to compare to digit");
}

template<typename T>
void getCofactor(Matrix<T> & matrix, Matrix<T> & temp, size_t exRow, size_t exCol, size_t n) {
    int i = 0, j = 0;

    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            if (row != exRow && col != exCol)
            {
                temp.content[i][j++] = matrix.content[row][col];

                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

template<typename T>
T findDeterminant(Matrix<T> & matrix, size_t n) {
    T determinant{};
    if (n == 1) return matrix.content[0][0];

    Matrix<T> temp = Matrix<T>(matrix.rows, matrix.cols);
    int sign = 1;

    for (size_t f = 0; f < n; ++f) {
        getCofactor<T>(matrix, temp, 0, f, n);
        determinant += sign * matrix.content[0][f] * findDeterminant<T>(temp, n - 1);

        sign = -sign;
    }
    return determinant;
}

template<typename T>
void findAdjoint(Matrix<T> & matrix, Matrix<T> & adj)
{
    if (matrix.cols == 1)
    {
        adj.content[0][0] = 1;
        return;
    }

    int sign = 1;
    Matrix<T> temp = Matrix<T>(matrix.rows, matrix.cols);

    for (size_t i = 0; i < matrix.rows; i++)
    {
        for (size_t j = 0; j < matrix.cols; j++)
        {
            getCofactor(matrix, temp, i, j, matrix.rows);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj.content[j][i] = sign * (findDeterminant<T>(temp, matrix.rows - 1));
        }
    }
}

template<typename T>
T findAdjointAt(Matrix<T> & matrix, size_t i, size_t j) {
    int sign = 1;
    Matrix<T> temp = Matrix<T>(matrix.rows, matrix.cols);
    getCofactor(matrix, temp, i, j, matrix.rows);
    sign = ((i + j) % 2 == 0) ? 1 : -1;
    return sign * (findDeterminant<T>(temp, matrix.rows - 1));
}

template<typename T>
bool inverse(Matrix<T> & matrix, Matrix<T> & inverse)
{
    if (matrix.cols != matrix.rows) {
        std::cout << "Singular matrix, can't find its inverse" << std::endl;
        return false;
    }

    T det = findDeterminant<T>(matrix, matrix.rows);
    if (abs(det) <= 0.0000000001) {
        std::cout << "Singular matrix, can't find its inverse" << std::endl;
        return false;
    }

    Matrix<T> adj = Matrix<T>(matrix.rows, matrix.cols);
    findAdjoint<T>(matrix, adj);

    for (int i = 0; i < matrix.rows; ++i)
        for (int j = 0; j < matrix.cols; ++j)
            inverse.setElement(adj.content[i][j] / T(det), i , j);

    return true;
}

template<typename T>
bool inverseParallel(Matrix<T> & matrix, Matrix<T> & inverse)
{
    if (matrix.cols != matrix.rows) {
        std::cout << "Singular matrix, can't find its inverse" << std::endl;
        return false;
    }

    T det = findDeterminant<T>(matrix, matrix.rows);
    if (abs(det) <= 0.0000000001) {
        std::cout << "Singular matrix, can't find its inverse" << std::endl;
        return false;
    }

    Matrix<T> adj = Matrix<T>(matrix.rows, matrix.cols);

    const size_t threadsCount = 4;
    std::vector<std::thread> threads = std::vector<std::thread>();

    for (size_t i = 0; i < threadsCount; ++i) {
        std::vector<size_t> indexes;

        for (size_t index = i; index < matrix.rows * matrix.cols; index += threadsCount) {
            indexes.push_back(index);
        }

        std::function<void ()> job = [indexes, &matrix, &adj]() {
            if (matrix.cols == 1)
            {
                adj.content[0][0] = 1;
                return;
            }

            int sign = 1;
            Matrix<T> temp = Matrix<T>(matrix.rows, matrix.cols);

            for (auto & index : indexes) {
                size_t i = index / matrix.cols;
                size_t j = index % matrix.cols;

                getCofactor(matrix, temp, i, j, matrix.rows);
                sign = ((i + j) % 2 == 0) ? 1 : -1;
                adj.content[j][i] = sign * (findDeterminant<T>(temp, matrix.rows - 1));
            }
        };

        threads.emplace_back(job);
    }

    for (size_t i = 0; i < threadsCount; ++i) threads[i].join();

    for (int i = 0; i < matrix.rows; ++i)
        for (int j = 0; j < matrix.cols; ++j)
            inverse.setElement(adj.content[i][j] / T(det), i , j);

    return true;
}

template<typename T>
bool inverseAsync(Matrix<T> & matrix, Matrix<T> & inverse)
{
    if (matrix.cols != matrix.rows) {
        std::cout << "Singular matrix, can't find its inverse" << std::endl;
        return false;
    }

    T det = findDeterminant<T>(matrix, matrix.rows);
    if (abs(det) <= 0.0000000001) {
        std::cout << "Singular matrix, can't find its inverse" << std::endl;
        return false;
    }

    Matrix<T> adj = Matrix<T>(matrix.rows, matrix.cols);

    for (size_t i = 0; i < matrix.rows; i++)
    {
        std::vector<std::future<T>> temp_res;
        for (size_t j = 0; j < matrix.cols; j++)
        {
            temp_res.push_back(std::async(std::launch::async, [](Matrix<T> matrix, size_t i, size_t j) {
                int sign = 1;
                Matrix<T> temp = Matrix<T>(matrix.rows, matrix.cols);
                getCofactor(matrix, temp, i, j, matrix.rows);
                sign = ((i + j) % 2 == 0) ? 1 : -1;
                return sign * (findDeterminant<T>(temp, matrix.rows - 1));
                }, matrix, i, j));
        }
        for (size_t j = 0; j < matrix.cols; j++)
        {
            adj.content[j][i] = temp_res[j].get();
        }
        temp_res.clear();
    }

    for (int i = 0; i < matrix.rows; ++i)
        for (int j = 0; j < matrix.cols; ++j)
            inverse.setElement(adj.content[i][j] / T(det), i , j);

    return true;
}

template<typename T>
Matrix<T> operator! (Matrix<T> & matrix) {
    if (matrix.cols != matrix.rows) throw std::exception("Matrix must be a square");
    Matrix<T> temp = Matrix<T>(matrix.rows, matrix.cols);
    if (!inverseParallel<T>(matrix, temp)) {
        throw std::exception("Matrix can't be inverse");
    }
    return temp;
}

template<typename T>
Matrix<T> operator* (Matrix<T> & lhs, int rhs) {
    auto temp = Matrix<T>(lhs.rows, lhs.cols);
    for (int i = 0; i < temp.rows; ++i) {
        for (int j = 0; j < temp.cols; ++j) {
            temp.content[i][j] = lhs.content[i][j] * rhs;
        }
    }
    return temp;
}

template<typename T>
Matrix<T> operator* (int lhs, Matrix<T> & rhs) {
    auto temp = Matrix<T>(rhs.rows, rhs.cols);
    for (int i = 0; i < temp.rows; ++i) {
        for (int j = 0; j < temp.cols; ++j) {
            temp.content[i][j] = rhs.content[i][j] * lhs;
        }
    }
    return temp;
}

template<typename T>
class Matrix {
private:
    size_t cols = 0;
    size_t rows = 0;

    std::vector<std::vector<T>> content;

    friend void getCofactor<T>(Matrix<T> & matrix, Matrix<T> & temp, size_t exRow, size_t exCol, size_t n);
    friend T findDeterminant<T>(Matrix<T> & matrix, size_t n);
    friend void findAdjoint<T>(Matrix<T> & matrix, Matrix<T> & adj);
    friend bool inverse<T>(Matrix<T> & matrix, Matrix<T> & inverse);
    friend bool inverseParallel<T>(Matrix<T> & matrix, Matrix<T> & inverse);
    friend T findAdjointAt<T>(Matrix<T> & matrix, size_t i, size_t j);
    friend bool inverseAsync<T>(Matrix<T> & matrix, Matrix<T> & inverse);

    friend Matrix<T> operator* <T> (const Matrix<T> & lhs, const Matrix<T> & rhs);
    friend Matrix<T> operator* <T> (Matrix<T> & lhs, int rhs);
    friend Matrix<T> operator* <T> (int lhs, Matrix<T> & rhs);
    friend Matrix<T> operator+ <T> (const Matrix<T> & lhs, const Matrix<T> & rhs);
    friend Matrix<T> operator- <T> (const Matrix<T> & lhs, const Matrix<T> & rhs);
    friend Matrix<T> operator! <T> (Matrix<T> & matrix);
    friend bool operator== <T> (const Matrix<T> & lhs, const Matrix<T> & rhs);
    friend bool operator== <T> (const Matrix<T> & lhs, int digit);
public:
    explicit Matrix(const std::string & path) {
        std::ifstream input(path);

        if (input.is_open()) {
            try {
                rows = readData<size_t>(&input);
                cols = readData<size_t>(&input);
            }
            catch (std::exception & e) {
                std::cout << "Failed to read columns and rows: " << e.what() << std::endl;
            }

            content = std::vector<std::vector<T>>(rows);
            for (size_t i = 0; i < rows; ++i) {
                content[i] = std::vector<T>(cols);
                for (size_t j = 0; j < cols; ++j) {
                    try {
                        content[i][j] = readData<T>(&input);
                    }
                    catch(std::exception & e) {
                        std::cout << "Failed to read value at: " << i << " " << j << " - " << e.what() << std::endl;
                    }
                }
            }
        }
    }

    explicit Matrix(std::istream * istream) {
        (*istream) >> rows;

        if ((*istream).fail()) {
            std::cout << "Failed to read rows" << std::endl;
        }

        (*istream) >> cols;

        if ((*istream).fail()) {
            std::cout << "Failed to read columns" << std::endl;
        }

        content = std::vector<std::vector<T>>(rows);
        for (size_t i = 0; i < rows; ++i) {
            content[i] = std::vector<T>(cols);
            for (size_t j = 0; j < cols; ++j) {
                (*istream) >> content[i][j];

                if ((*istream).fail()) {
                    std::cout << "Failed to read value at: " << i << " " << j << std::endl;
                    exit(1);
                }
            }
        }
    }

    Matrix(size_t rows, size_t cols) {
        this->rows = rows;
        this->cols = cols;
        content = std::vector<std::vector<T>>(rows);
        for (size_t i = 0; i < rows; ++i) { content[i] = std::vector<T>(cols); }
    }

    Matrix(const Matrix& other) {
        rows = other.rows;
        cols = other.cols;
        content = std::vector<std::vector<T>>(rows);
        for (size_t i = 0; i < rows; ++i) content[i] = other.content[i];
    }

    Matrix(Matrix&& other) noexcept {
        rows = other.rows;
        cols = other.cols;
        content = std::move(other.content);
        //for (size_t i = 0; i < rows; ++i) content[i] = std::move(other.content[i]);
        other.rows = 0;
        other.cols = 0;
    }

    Matrix(std::initializer_list<std::tuple<std::tuple<size_t, size_t>, T>> elements) {
        content = std::vector<std::vector<T>>();

        for (auto element : elements) {
            std::tuple<size_t, size_t> pos;
            T value;
            std::tie(pos, value) = element;

            size_t row;
            size_t col;
            std::tie(row, col) = pos;

            if (row >= content.size()) {
                content.resize(row + 1);
                for (size_t i = 0; i < content.size(); i++) content[i].resize(content[0].size());
            }
            if (col >= content[0].size()) for (size_t i = 0; i < content.size(); i++) content[i].resize(col + 1);
            content[row][col] = value;
        }

        rows = content.size();
        cols = content[0].size();
    }

    Matrix() {
        content = std::vector<std::vector<T>>();
    }

    void print() {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << content[i][j] << "\t";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    void writeMatrix(const std::string & path) {
        std::ofstream out(path);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                out << content[i][j] << "\t";
            }
            out << "\n";
        }
        out << std::endl;
    }

    size_t getCols() { return cols; }
    size_t getRows() { return rows; }

    void setCols(size_t newCols) {
        if (newCols <= 0) {
            cols = 1;
        }
    }
    void setRows(size_t newRows) {
        if (newRows <= 0) {
            rows = 1;
        }
    }

    T getElement(size_t row, size_t col) { return content[row][cols]; }

    void setElement(T element, size_t row, size_t col) { content[row][col] = element; }

    Matrix<T> threadMultiplication(const Matrix<T> & other) {
        if (cols != other.rows) throw std::exception("Columns and rows must be the same");

        Matrix<T> product = Matrix<T>(rows, other.cols);

        const size_t threadsCount = 4;
        std::vector<std::thread> threads = std::vector<std::thread>();

        for (size_t i = 0; i < threadsCount; ++i) {
            std::vector<size_t> indexes;

            for (size_t index = i; index < rows * cols; index += threadsCount) {
                indexes.push_back(index);
            }

            std::function<void ()> job = [this, indexes, &other, &product]() {
                for (auto & index : indexes) {
                    size_t row = index / this->cols;
                    size_t col = index % this->cols;

                    T local_sum = 0;
                    for (size_t k = 0; k < cols; ++k) {
                        local_sum += this->content[row][k] * other.content[k][col];
                    }
                    if (abs(local_sum) < 0.000000001) {
                        local_sum = 0;
                    }
                    product.content[row][col] = local_sum;
                }
            };

            threads.emplace_back(job);
        }

        for (size_t i = 0; i < threadsCount; ++i) {
            threads[i].join();
        }

        return product;
    }

    Matrix<T> & operator= (const Matrix<T> & other) {
        other.cols = cols;
        other.rows = rows;
        other.content = std::copy(content.begin(), content.end());
    }
};


#endif //PARALLELMATRIX_MATRIX_H
