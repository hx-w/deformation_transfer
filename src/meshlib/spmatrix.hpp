/**
 * implement a sparse matrix class
 * with multiplication
 */

#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include "matrix.hpp"

namespace MeshLib {

template <class T>
class Triplet {
public:
    Triplet(int row, int col, T value) : row(row), col(col), value(value) {}
    int row;
    int col;
    T value;

    // for sorting
    bool operator<(const Triplet &other) const {
        if (row < other.row) {
            return true;
        } else if (row == other.row) {
            return col < other.col;
        } else {
            return false;
        }
    }
};


// csr format
// Compressed Sparse Row matrix
template <class T>
class SparseMatrix {
public:
    SparseMatrix() {}
    SparseMatrix(int rows, int cols) : m_rows(rows), m_cols(cols) {
        m_row_ind.resize(rows + 1);
        m_row_ind[0] = 0;
    }
    SparseMatrix(int rows, int cols, const std::vector<Triplet<T>> &triplets):
        m_rows(rows), m_cols(cols) {
        // sort triplets
        std::vector<Triplet<T>> sorted_triplets(triplets);
        std::sort(sorted_triplets.begin(), sorted_triplets.end());

        // construct csr format: row_ind, col_ind, values
        m_row_ind.emplace_back(0);
        for (auto i = 0; i < sorted_triplets.size(); ++i) {
            int row = sorted_triplets[i].row;
            int col = sorted_triplets[i].col;
            T value = sorted_triplets[i].value;

            m_col_ind.emplace_back(col);
            m_values.emplace_back(value);

            if (i == sorted_triplets.size() - 1 || row != sorted_triplets[i + 1].row) {
                m_row_ind.emplace_back(m_col_ind.size());
            }
        }
    }

    // copy construct
    SparseMatrix(const SparseMatrix &other) {
        m_rows = other.m_rows;
        m_cols = other.m_cols;
        m_row_ind = other.m_row_ind;
        m_col_ind = other.m_col_ind;
        m_values = other.m_values;
    }

    // move construct
    SparseMatrix(SparseMatrix &&other) {
        m_rows = other.m_rows;
        m_cols = other.m_cols;
        m_row_ind = std::move(other.m_row_ind);
        m_col_ind = std::move(other.m_col_ind);
        m_values = std::move(other.m_values);
    }

    void resize(int rows, int cols) {
        m_rows = rows;
        m_cols = cols;
        m_row_ind.resize(rows + 1);
        m_row_ind[0] = 0;
    }

    // copy assignment
    SparseMatrix &operator=(const SparseMatrix &other) {
        m_rows = other.m_rows;
        m_cols = other.m_cols;
        m_row_ind = other.m_row_ind;
        m_col_ind = other.m_col_ind;
        m_values = other.m_values;
        return *this;
    }

    int rows() const { return m_rows; }
    int cols() const { return m_cols; }

    void to_dense(Matrix<T> &mat) const {
        mat = Matrix<T>(m_rows, m_cols);
        for (auto i = 0; i < m_rows; ++i) {
            for (auto j = m_row_ind[i]; j < m_row_ind[i + 1]; ++j) {
                mat(i, m_col_ind[j]) = m_values[j];
            }
        }
    }

    // save to file, in binary format
    void save(const std::string &filename) const {
        std::ofstream fout(filename, std::ios::binary);
        fout.write((char *)&m_rows, sizeof(int));
        fout.write((char *)&m_cols, sizeof(int));
        int nnz = m_values.size();
        fout.write((char *)&nnz, sizeof(int));
        fout.write((char *)m_row_ind.data(), sizeof(int) * (m_rows + 1));
        fout.write((char *)m_col_ind.data(), sizeof(int) * nnz);
        fout.write((char *)m_values.data(), sizeof(T) * nnz);
        fout.close();
    }

    // load from file
    void load(const std::string &filename) {
        std::ifstream fin(filename, std::ios::binary);
        fin.read((char *)&m_rows, sizeof(int));
        fin.read((char *)&m_cols, sizeof(int));
        int nnz;
        fin.read((char *)&nnz, sizeof(int));
        m_row_ind.resize(m_rows + 1);
        m_col_ind.resize(nnz);
        m_values.resize(nnz);
        fin.read((char *)m_row_ind.data(), sizeof(int) * (m_rows + 1));
        fin.read((char *)m_col_ind.data(), sizeof(int) * nnz);

        fin.read((char *)m_values.data(), sizeof(T) * nnz);
        fin.close();
    }

    // save to txt
    void save_txt(const std::string &filename) const {
        std::ofstream fout(filename);
        fout << m_rows << " " << m_cols << std::endl;
        for (auto i = 0; i < m_rows; ++i) {
            for (auto j = m_row_ind[i]; j < m_row_ind[i + 1]; ++j) {
                fout << i << " " << m_col_ind[j] << " " << m_values[j] << std::endl;
            }
        }
        fout.close();
    }

    // load from txt
    void load_txt(const std::string &filename) {
        std::ifstream fin(filename);
        fin >> m_rows >> m_cols;
        m_row_ind.resize(m_rows + 1);
        m_row_ind[0] = 0;
        int row, col;
        T value;
        while (fin >> row >> col >> value) {
            m_col_ind.emplace_back(col);
            m_values.emplace_back(value);
            if (m_col_ind.size() == m_row_ind[m_row_ind.size() - 1] + 1) {
                m_row_ind.emplace_back(m_col_ind.size());
            }
        }
        fin.close();
    }

    // get the value at (row, col)
    T& operator()(int row, int col) {
        assert(row >= 0 && row < m_rows);
        assert(col >= 0 && col < m_cols);

        int start = m_row_ind[row];
        int end = m_row_ind[row + 1];
        for (auto i = start; i < end; ++i) {
            if (m_col_ind[i] == col) {
                return m_values[i];
            }
        }

        // not found,
        // insert a new element of CSR format
        m_col_ind.emplace(m_col_ind.begin() + end, col);
        m_values.emplace(m_values.begin() + end, 0);
        for (auto i = row + 1; i < m_row_ind.size(); ++i) {
            m_row_ind[i] += 1;
        }

        return m_values[end];
    }

    // get the value at (row, col)
    const T at(int row, int col) const {
        assert(row >= 0 && row < m_rows);
        assert(col >= 0 && col < m_cols);

        int start = m_row_ind[row];
        int end = m_row_ind[row + 1];
        for (auto i = start; i < end; ++i) {
            if (m_col_ind[i] == col) {
                return m_values[i];
            }
        }
        // not found
        return 0;
    }

    // basic operator
    SparseMatrix operator+(const SparseMatrix &other) const {
        assert(m_rows == other.m_rows);
        assert(m_cols == other.m_cols);

        SparseMatrix result(m_rows, m_cols);
        // add two csr matrix: m_row_ind, m_col_ind, m_values
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                auto& value = result(i, m_col_ind[j]);
                result(i, m_col_ind[j]) += m_values[j];
            }
        }

        for (auto i = 0; i < other.m_row_ind.size() - 1; ++i) {
            int start = other.m_row_ind[i];
            int end = other.m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                result(i, other.m_col_ind[j]) += other.m_values[j];
            }
        }


        return result;
    }

    SparseMatrix operator-(const SparseMatrix &other) const {
        assert(m_rows == other.m_rows);
        assert(m_cols == other.m_cols);

        std::vector<Triplet<T>> triplets;
        for (auto i = 0; i < m_rows; ++i) {
            for (auto j = 0; j < m_cols; ++j) {
                T value = at(i, j) - other.at(i, j);
                if (value != 0) {
                    triplets.emplace_back(i, j, value);
                }
            }
        }
        return SparseMatrix(m_rows, m_cols, triplets);
    }

    SparseMatrix operator*(const SparseMatrix &other) const {
        assert(m_cols == other.m_rows);
        SparseMatrix result(m_rows, other.m_cols);

        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                int col = m_col_ind[j];
                T value = m_values[j];
                for (auto k = other.m_row_ind[col]; k < other.m_row_ind[col + 1]; ++k) {
                    result(i, other.m_col_ind[k]) += value * other.m_values[k];
                }
            }
        }

        return result;
    }

    // scalar operator
    SparseMatrix operator*(const T &scalar) const {
        SparseMatrix result(m_rows, m_cols);
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                result(i, m_col_ind[j]) = m_values[j] * scalar;
            }
        }
        return result;
    }

    SparseMatrix operator/(const T &scalar) const {
        SparseMatrix result(m_rows, m_cols);
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                result(i, m_col_ind[j]) = m_values[j] / scalar;
            }
        }
        return result;
    }

    // inplace operator
    SparseMatrix& operator+=(const SparseMatrix &other) {
        assert(m_rows == other.m_rows);
        assert(m_cols == other.m_cols);

        for (auto i = 0; i < other.m_row_ind.size() - 1; ++i) {
            int start = other.m_row_ind[i];
            int end = other.m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                (*this)(i, other.m_col_ind[j]) += other.m_values[j];
            }
        }

        return *this;
    }

    SparseMatrix& operator-=(const SparseMatrix &other) {
        assert(m_rows == other.m_rows);
        assert(m_cols == other.m_cols);

        for (auto i = 0; i < other.m_row_ind.size() - 1; ++i) {
            int start = other.m_row_ind[i];
            int end = other.m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                (*this)(i, other.m_col_ind[j]) -= other.m_values[j];
            }
        }

        return *this;
    }

    // transpose
    SparseMatrix transpose() const {
        SparseMatrix result(m_cols, m_rows);
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                result(m_col_ind[j], i) = m_values[j];
            }
        }
        return result;
    }

    // cofactor
    SparseMatrix cofactor() const {
        assert(m_rows == m_cols);
        SparseMatrix result(m_rows, m_cols);
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                result(i, m_col_ind[j]) = std::pow(-1, i + m_col_ind[j]) * at(m_col_ind[j], i);
            }
        }
        return result;
    }

    // submatrix
    SparseMatrix submatrix(int row, int col) const {
        assert(row < m_rows);
        assert(col < m_cols);

        std::vector<Triplet<T>> triplets;
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            if (i == row) {
                continue;
            }
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                if (m_col_ind[j] == col) {
                    continue;
                }
                triplets.emplace_back(i < row ? i : i - 1, m_col_ind[j] < col ? m_col_ind[j] : m_col_ind[j] - 1, m_values[j]);
            }
        }
        return SparseMatrix(m_rows - 1, m_cols - 1, triplets);
    }

    // determinant
    T determinant() const {
        assert(m_rows == m_cols);
        if (m_rows == 1) {
            return at(0, 0);
        }
        if (m_rows == 2) {
            return at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);
        }
        T result = 0;
        for (auto i = 0; i < m_cols; ++i) {
            result += at(0, i) * std::pow(-1, i) * submatrix(0, i).determinant();
        }
        return result;
    }


    // inverse matrix
    SparseMatrix inverse() const {
        assert(m_rows == m_cols);
        SparseMatrix result(m_rows, m_cols);
        
        T det = determinant();
        assert(det != 0);

        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                result(i, m_col_ind[j]) = std::pow(-1, i + m_col_ind[j]) * submatrix(m_col_ind[j], i).determinant() / det;
            }
        }
        return result;
    }

    // LU decompose
    void LU_decompose(SparseMatrix &L, SparseMatrix &U) const {
        assert(m_rows == m_cols);

        // L = SparseMatrix::identity(m_rows);
        L = SparseMatrix::identity(m_rows);
        U = SparseMatrix(m_rows, m_cols);

        // sparse csr matrix LU decompose
        // L is lower triangular matrix
        // U is upper triangular matrix
        // L * U = A
    }

    // Cholesky decompose
    void Cholesky_decompose(SparseMatrix &L) const {
        assert(m_rows == m_cols);
        // sparse csr matrix Cholesky decompose
        // L is lower triangular matrix
        // L * L^T = A
        // A is symmetric positive definite matrix
        L = SparseMatrix::identity(m_rows);

        for (auto i = 0; i < m_rows; ++i) {
            if (i % 3 == 0)
            std::cout << i << " " << m_rows << " " << L.m_values.size() << std::endl;
            for (auto j = 0; j < i; ++j) {
                T sum = 0;
                for (auto k = 0; k < j; ++k) {
                    sum += L.at(i, k) * L.at(j, k);
                }
                auto _v = (at(i, j) - sum) / L.at(j, j);
                if (std::abs(_v) > 1e-16) {
                    L(i, j) = _v;
                }
            }
            T sum = 0;
            for (auto k = 0; k < i; ++k) {
                sum += std::pow(L.at(i, k), 2);
            }
            auto ii = std::sqrt(at(i, i) - sum);
            L(i, i) = ii;
            assert(std::isnan(ii) == false);
        }
    }


    // A shape(m, m)
    // b shape(m, k)
    void solve(const SparseMatrix<T>& b, SparseMatrix<T>& x) {
        assert(m_rows == m_cols);
        assert(m_rows == b.m_rows);

        // solve Ax=b, A is symmetric positive definite matrix
        SparseMatrix L;

        if (std::ifstream(".cache/Cholesky_decompose_L.mat")) {
            L.load(".cache/Cholesky_decompose_L.mat");
            std::cout << "load Cholesky_decompose_L.mat" << std::endl;

            for (int r = 330; r < 340; r++) {
                for (int c = 330; c < 340; c++) {
                    std::cout << L.at(r, c) << " ";
                }
                std::cout << std::endl;
            }
        } else {
            Cholesky_decompose(L);
            L.save(".cache/Cholesky_decompose_L.mat");
        }
        SparseMatrix Lt = L.transpose();

        // check if has nan in L
        int _r = 0, _c = 0;
        for (auto i = 0; i < L.m_row_ind.size() - 1; ++i) {
            int start = L.m_row_ind[i];
            int end = L.m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                if (std::isnan(L.m_values[j])) {
                    _r = i, _c = L.m_col_ind[j];
                    break;
                }
            }
            if (_r != 0) {
                break;
            }
        }
        if (_r != 0) {
            std::cout << "L has nan at " << _r << ", " << _c << std::endl;
        }


        // solve Ly = b
        SparseMatrix y(b.m_rows, b.m_cols);
        for (int i = 0; i < b.m_rows; ++i) {
            if (i % 500 == 0)
            std::cout << "solve: " << i << "/" << b.m_rows << std::endl;
            for (int j = 0; j < b.m_cols; ++j) {
                T sum = 0;
                for (int k = 0; k < i; ++k) {
                    sum += L.at(i, k) * y.at(k, j);
                }
                y(i, j) = (b.at(i, j) - sum) / L.at(i, i);
            }
        }

        // solve L^Tx = y
        x = SparseMatrix(b.m_rows, b.m_cols);
        for (int i = b.m_rows - 1; i >= 0; --i) {
            if (i % 500 == 0)
            std::cout << "solve2: " << i << "/" << b.m_rows << std::endl;
            for (int j = 0; j < b.m_cols; ++j) {
                T sum = 0;
                for (int k = i + 1; k < b.m_rows; ++k) {
                    sum += Lt.at(i, k) * x.at(k, j);
                }
                x(i, j) = (y.at(i, j) - sum) / Lt.at(i, i);
            }
        }

    }

    void solve_(const SparseMatrix<T>& b, SparseMatrix<T>& x) {
        assert(m_rows == m_cols);
        assert(m_rows == b.rows());

        // Ax = b
        // Jacobi iteration
        const auto iter_num = 500;
        x = SparseMatrix<T>(b.rows(), b.cols());

        SparseMatrix D_inv(m_rows, m_cols);
        SparseMatrix R(m_rows, m_cols);
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                if (m_col_ind[j] == i) {
                    D_inv(i, i) = 1.0 / m_values[j];
                }
                else {
                    R(i, m_col_ind[j]) = m_values[j];
                }
            }
        }

        for (auto i = 0; i < iter_num; ++i) {
            x = D_inv * (b - R * x);
        }
        
    }

    bool is_diagonally_dominant() const {
        for (auto i = 0; i < m_rows; ++i) {
            T sum = 0;
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                if (m_col_ind[j] != i) {
                    sum += std::abs(m_values[j]);
                }
            }
            if (std::abs(m_values[start]) < sum) {
                return false;
            }
        }
        return true;
    }

    static SparseMatrix identity(int n) {
        std::vector<Triplet<T>> triplets;
        for (auto i = 0; i < n; ++i) {
            triplets.emplace_back(i, i, 1);
        }
        return SparseMatrix(n, n, triplets);
    }

    // diagonal matrix
    static SparseMatrix diagonal(const std::vector<T>& v) {
        std::vector<Triplet<T>> triplets;
        for (auto i = 0; i < v.size(); ++i) {
            triplets.emplace_back(i, i, v[i]);
        }
        return SparseMatrix(v.size(), v.size(), triplets);
    }

    // slice
    SparseMatrix slice(int row_start, int row_end, int col_start, int col_end) const {
        assert(row_start >= 0);
        assert(row_end <= m_rows);
        assert(col_start >= 0);
        assert(col_end <= m_cols);
        assert(row_start < row_end);
        assert(col_start < col_end);

        std::vector<Triplet<T>> triplets;
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            if (i < row_start || i >= row_end) {
                continue;
            }
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                if (m_col_ind[j] < col_start || m_col_ind[j] >= col_end) {
                    continue;
                }
                triplets.emplace_back(i - row_start, m_col_ind[j] - col_start, m_values[j]);
            }
        }
        return SparseMatrix(row_end - row_start, col_end - col_start, triplets);
    }

    // slice with indices, axis = 1 for column, axis = 0 for row
    SparseMatrix slice(const std::vector<size_t>& indices, size_t axis) const {
        assert(axis == 0 || axis == 1);
        if (axis == 0) {
            assert(indices.size() <= m_rows);
            std::vector<Triplet<T>> triplets;
            for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
                if (std::find(indices.begin(), indices.end(), i) == indices.end()) {
                    continue;
                }
                int start = m_row_ind[i];
                int end = m_row_ind[i + 1];
                for (auto j = start; j < end; ++j) {
                    triplets.emplace_back(i, m_col_ind[j], m_values[j]);
                }
            }
            return SparseMatrix(indices.size(), m_cols, triplets);
        } else {
            assert(indices.size() <= m_cols);
            std::vector<Triplet<T>> triplets;
            for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
                int start = m_row_ind[i];
                int end = m_row_ind[i + 1];
                for (auto j = start; j < end; ++j) {
                    if (std::find(indices.begin(), indices.end(), m_col_ind[j]) == indices.end()) {
                        continue;
                    }
                    triplets.emplace_back(i, m_col_ind[j], m_values[j]);
                }
            }
            return SparseMatrix(m_rows, indices.size(), triplets);
        }
    }

    Vector<T> to_vector() const {
        Vector<T> result(m_rows * m_cols);
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                result(i * m_cols + m_col_ind[j]) = m_values[j];
            }
        }
        return result;
    }

private:
    int m_rows;
    int m_cols;
    std::vector<int> m_row_ind;
    std::vector<int> m_col_ind;
    std::vector<T> m_values;
};

}

#endif
