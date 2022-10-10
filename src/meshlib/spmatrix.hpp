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
        L = SparseMatrix(m_rows, m_cols);
        U = SparseMatrix(m_rows, m_cols);
        
        for (auto i = 0; i < m_row_ind.size() - 1; ++i) {
            int start = m_row_ind[i];
            int end = m_row_ind[i + 1];
            for (auto j = start; j < end; ++j) {
                if (m_col_ind[j] < i) {
                    T sum = 0;
                    for (auto k = m_row_ind[m_col_ind[j]]; k < m_row_ind[m_col_ind[j] + 1]; ++k) {
                        if (m_col_ind[k] < m_col_ind[j]) {
                            sum += L.at(i, m_col_ind[k]) * U.at(m_col_ind[k], m_col_ind[j]);
                        }
                    }
                    L(i, m_col_ind[j]) = (m_values[j] - sum) / U.at(m_col_ind[j], m_col_ind[j]);
                } else {
                    T sum = 0;
                    for (auto k = m_row_ind[i]; k < m_row_ind[i + 1]; ++k) {
                        if (m_col_ind[k] < i) {
                            sum += L.at(i, m_col_ind[k]) * U.at(m_col_ind[k], m_col_ind[j]);
                        }
                    }
                    U(i, m_col_ind[j]) = m_values[j] - sum;
                }
            }
        }
        // diagonal of L is 1
        for (auto i = 0; i < m_rows; ++i) {
            L(i, i) = static_cast<T>(1);
        }
    }

    static SparseMatrix identity(int n) {
        std::vector<Triplet<T>> triplets;
        for (auto i = 0; i < n; ++i) {
            triplets.emplace_back(i, i, 1);
        }
        return SparseMatrix(n, n, triplets);
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
