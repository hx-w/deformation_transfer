// matrix
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <omp.h>

namespace MeshLib {
template <class T> class Vector;
template <class T> class Matrix;

using TKey = std::pair<size_t, size_t>;
using MatrixXd = Matrix<double>;
using MatrixXi = Matrix<size_t>;
using VectorXd = Vector<double>;
using VectorXi = Vector<int>;
// n_dim vector class with math operators
template <class T>
class Vector {
public:
    Vector() : m_data() {}
    Vector(size_t n) : m_data(n) {}
    Vector(size_t n, const T& val) : m_data(n, val) {}
    Vector(const Vector<T>& v) : m_data(v.m_data) {}
    Vector(const std::vector<T>& v) : m_data(v) {}
    Vector(const std::initializer_list<T>& v) : m_data(v) {}
    Vector<T>& operator=(const Vector<T>& v) {
        m_data.assign(v.m_data.begin(), v.m_data.end());
        return *this;
    }
    Vector<T>& operator=(const T* v) {
        m_data.assign(v, v+m_data.size());
        return *this;
    }
    Vector<T>& operator+=(const Vector<T>& v) {
        assert(m_data.size() == v.m_data.size());
        for (int i = 0; i < m_data.size(); ++i)
            m_data[i] += v.m_data[i];
        return *this;
    }
    Vector<T>& operator-=(const Vector<T>& v) {
        assert(m_data.size() == v.m_data.size());
        for (int i = 0; i < m_data.size(); ++i)
            m_data[i] -= v.m_data[i];
        return *this;
    }
    Vector<T>& operator*=(const T& val) {
        for (int i = 0; i < m_data.size(); ++i)
            m_data[i] *= val;
        return *this;
    }
    Vector<T>& operator/=(const T& val) {
        for (int i = 0; i < m_data.size(); ++i)
            m_data[i] /= val;
        return *this;
    }
    Vector<T> operator+(const Vector<T>& v) const {
        Vector<T> r(*this);
        r += v;
        return r;
    }
    Vector<T> operator-(const Vector<T>& v) const {
        Vector<T> r(*this);
        r -= v;
        return r;
    }
    Vector<T> operator*(const T& val) const {
        Vector<T> r(*this);
        r *= val;
        return r;
    }
    Vector<T> operator/(const T& val) const {
        Vector<T> r(*this);
        r /= val;
        return r;
    }

    // get value
    T& operator[](int i) { return m_data.at(i); }
    const T& operator[](int i) const { return m_data.at(i); }
    decltype(auto) data() const { return m_data; }

    // get size
    decltype(auto) size() const { return m_data.size(); }

    // get norm
    double norm() const {
        double n = 0;
        for (auto i = 0; i < m_data.size(); ++i)
            n += m_data[i] * m_data[i];
        return sqrt(n);
    }

    // normalize
    Vector<T> normalize() const {
        Vector<T> r(*this);
        r /= r.norm();
        return r;
    }

    // dot product
    T dot(const Vector<T>& v) const {
        assert(m_data.size() == v.m_data.size());
        T r = 0;
        for (int i = 0; i < m_data.size(); ++i)
            r += m_data[i] * v.m_data[i];
        return r;
    }

    // cross product (only for 3D)
    Vector<T> cross(const Vector<T>& v) const {
        assert(m_data.size() == 3);
        Vector<T> r(3);
        r[0] = m_data[1] * v.m_data[2] - m_data[2] * v.m_data[1];
        r[1] = m_data[2] * v.m_data[0] - m_data[0] * v.m_data[2];
        r[2] = m_data[0] * v.m_data[1] - m_data[1] * v.m_data[0];
        return r;
    }

    // to matrix, axis = 0: column, axis = 1: row
    Matrix<T> to_matrix(int axis = 0) const {
        auto n = m_data.size();
        std::vector<size_t> rc = {n, 1};
        Matrix<T> r(rc[axis], rc[1 - axis]);
        for (auto i = 0; i < n; ++i) {
            if (axis == 0) {
                r(i, 0) = m_data[i];
            }
            else {
                r(0, i) = m_data[i];
            }
        }
        return r;
    }

    // slice
    Vector<T> slice(int start, int end) const {
        assert(start >= 0 && end <= m_data.size());
        Vector<T> r(end - start);
        for (auto i = start; i < end; ++i)
            r[i - start] = m_data[i];
        return r;
    }

private:
    std::vector<T> m_data;
};


template <class T>
class Matrix {
public:
    Matrix() : m_data() {}

    Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols) {}

    Matrix(std::vector<std::vector<T>> data) : m_rows(data.size()), m_cols(data[0].size()) {
        // m_data: map<std::make_pair, T>
        // std::make_pair: pair<size_t, size_t>
        for (auto i = 0; i < m_rows; ++i) {
            for (auto j = 0; j < m_cols; ++j) {
                (*this)(i, j) = data[i][j];
            }
        }
    }

    Matrix(const Matrix& other) : m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.m_data) {}

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            m_rows = other.m_rows;
            m_cols = other.m_cols;
            for (auto [k, v] : other.m_data) {
                m_data[k] = v;
            }
        }
        return *this;
    }

    Matrix(Matrix&& other) : m_rows(other.m_rows), m_cols(other.m_cols), m_data(std::move(other.m_data)) {}

    // fetch values
    const T operator()(size_t row, size_t col) const {
        if (m_data.find(std::make_pair(row, col)) != m_data.end()) {
            return m_data.at(std::make_pair(row, col));
        }
        return static_cast<T>(0);
    }

    T& operator()(size_t row, size_t col) {
        if (m_data.find(std::make_pair(row, col)) != m_data.end()) {
            return m_data[std::make_pair(row, col)];
        }
        set(row, col, static_cast<T>(0));
        return m_data[std::make_pair(row, col)];
    }

    // const fetch
    const T at(size_t row, size_t col) const {
        if (m_data.find(std::make_pair(row, col)) != m_data.end()) {
            return m_data.at(std::make_pair(row, col));
        }
        return static_cast<T>(0);
    }

    void set(size_t row, size_t col, const T& val) {
        m_data[std::make_pair(row, col)] = val;
        m_rows = std::max(m_rows, row + 1);
        m_cols = std::max(m_cols, col + 1);
    }

    // get row and column with index
    Vector<T> row(size_t i) const {
        Vector<T> r(m_cols);
        for (auto j = 0; j < m_cols; ++j) {
            r[j] = (*this)(i, j);
        }
        return r;
    }

    Vector<T> col(size_t j) const {
        Vector<T> r(m_rows);
        for (auto i = 0; i < m_rows; ++i) {
            r[i] = (*this)(i, j);
        }
        return r;
    }

    decltype(auto) data() const { 
        std::vector<T> _values;
        for (auto [_k, _v]: m_data) {
            _values.emplace_back(_v);
        }
        return _values;
    }

    void append_row(const Vector<T>& v) {
        assert(v.size() == m_cols);
        const auto _row = m_rows;
        for (auto j = 0; j < m_cols; ++j) {
            (*this)(_row, j) = v[j];
        }
    }

    void append_col(const Vector<T>& v) {
        assert(v.size() == m_rows);
        const auto _col = m_cols;
        for (auto i = 0; i < m_rows; ++i) {
            (*this)(i, _col) = v[i];
        }
    }

    // axis = 1: column, axis = 0: row
    void append(const Matrix<T>& m, int axis) {
        assert(axis == 0 || axis == 1);
        if (axis == 0) {
            auto old_rows = m_rows;
            assert(m.cols() == m_cols);
            for (auto [_key, _val] : m.m_data) {
                auto i = _key.first, j = _key.second;
                (*this)(i + old_rows, j) = _val;
            }
            m_rows = old_rows + m.rows();
        }
        else {
            auto old_cols = m_cols;
            assert(m.rows() == m_rows);
            for (auto [_key, _val] : m.m_data) {
                auto i = _key.first, j = _key.second;
                (*this)(i, j + old_cols) = _val;
            }
            m_cols = old_cols + m.cols();
        }
    }

    auto rows() const { return m_rows; }
    auto cols() const { return m_cols; }
    auto shape() const { return std::make_pair(m_rows, m_cols); }
    auto shape(size_t dim) const { return dim == 0 ? m_rows : m_cols; }

    // addition
    Matrix operator+(const Matrix& other) const {
        assert(m_rows == other.m_rows && m_cols == other.m_cols);
        Matrix result(m_rows, m_cols);
        for (auto [_key, _val] : m_data) {
            auto i = _key.first, j = _key.second;
            result(i, j) = (*this).at(i, j) + other.at(i, j);
        }

        for (auto [_key, _val] : other.m_data) {
            auto i = _key.first, j = _key.second;
            result(i, j) = (*this).at(i, j) + other.at(i, j);
        }
        return result;
    }

    // subtraction
    Matrix operator-(const Matrix& other) const {
        assert(m_rows == other.m_rows && m_cols == other.m_cols);
        Matrix result(m_rows, m_cols);

        for (auto [_key, _val] : m_data) {
            auto i = _key.first, j = _key.second;
            result(i, j) = (*this).at(i, j) - other.at(i, j);
        }

        for (auto [_key, _val] : other.m_data) {
            auto i = _key.first, j = _key.second;
            result(i, j) = (*this).at(i, j) - other.at(i, j);
        }

        return result;
    }

    // multiplication
    // need change to sparse multiplication
    Matrix operator*(const Matrix& other) const {
        assert(m_cols == other.m_rows);
        Matrix result(m_rows, other.m_cols);
        // m_data: map (i, j) -> value
        // sparse matrix multiplication
        for (auto [_key_f, _val_f] : m_data) {
            auto i = _key_f.first, j = _key_f.second;
// #pragma omp parallel for
//             for (int idx = 0; idx < other.m_data.size(); ++idx) {
//                 auto iter = other.m_data.begin();
//                 iter++;
//                 // auto p = (other.m_data.begin() + idx)->first.first;
//                 // auto q = (other.m_data.begin() + idx)->first.second;
//                 // result(i, q) += (*this).at(i, j) * other.at(p, q);
//             }
            for (auto& [_key_s, _val_s] : other.m_data) {
                auto p = _key_s.first, q = _key_s.second;
            }
        }
        return result;
    }

    // scalar multiplication
    Matrix operator*(const T& scalar) const {
        Matrix result(m_rows, m_cols);
        for (auto [_key, _val] : m_data) {
            auto i = _key.first, j = _key.second;
            result(i, j) = (*this).at(i, j) * scalar;
        }
        return result;
    }

    // slice
    Matrix slice(size_t row_start, size_t row_end, size_t col_start, size_t col_end) const {
        assert(row_start < row_end && col_start < col_end);
        Matrix result(row_end - row_start, col_end - col_start);

        for (auto [_key, _val] : m_data) {
            auto i = _key.first, j = _key.second;
            if (i >= row_start && i < row_end && j >= col_start && j < col_end) {
                result(i - row_start, j - col_start) = (*this).at(i, j);
            }
        }
        return result;
    }

    // slice with indices, axis = 1 for column, axis = 0 for row
    Matrix slice(const std::vector<size_t>& indices, size_t axis) const {
        assert(axis == 0 || axis == 1);
        Matrix result(axis == 0 ? indices.size() : m_rows, axis == 1 ? indices.size() : m_cols);
        map<size_t, size_t> index_map;
        for (auto i = 0; i < indices.size(); ++i) {
            index_map[indices[i]] = i;
        }

        for (auto [_key, _val] : m_data) {
            auto i = _key.first, j = _key.second;
            if (std::find(indices.begin(), indices.end(), axis == 0 ? i : j) != indices.end()) {
                if (axis == 0) {
                    result(index_map[i], j) = (*this).at(i, j);
                }
                else {
                    result(i, index_map[j]) = (*this).at(i, j);
                }
            }
        }

        return result;
    }

    // transpose
    Matrix transpose() const {
        Matrix result(m_cols, m_rows);
        for (auto [_key, _val] : m_data) {
            auto i = _key.first, j = _key.second;
            result(j, i) = (*this).at(i, j);
        }
        return result;
    }

    // determinant
    T determinant() const {
        assert(m_rows == m_cols);
        if (m_rows == 1) {
            return (*this)(0, 0);
        }
        T result = 0;
        for (size_t i = 0; i < m_cols; ++i) {
            result += (*this)(0, i) * cofactor(0, i);
        }
        return result;
    }

    // cofactor
    T cofactor(size_t row, size_t col) const {
        assert(m_rows == m_cols);
        Matrix minor(m_rows - 1, m_cols - 1);
        for (size_t i = 0; i < m_rows; ++i) {
            for (size_t j = 0; j < m_cols; ++j) {
                if (i != row && j != col) {
                    minor(i < row ? i : i - 1, j < col ? j : j - 1) = (*this)(i, j);
                }
            }
        }
        return minor.determinant() * ((row + col) % 2 == 0 ? 1 : -1);
    }

    // inverse
    decltype(auto) inverse() const {
        assert(m_rows == m_cols);
        T det = determinant();
        assert(det != 0);
        Matrix result(m_rows, m_cols);
        for (size_t i = 0; i < m_rows; ++i) {
            for (size_t j = 0; j < m_cols; ++j) {
                result(i, j) = cofactor(i, j) / det;
            }
        }
        return result;
    }

    // to vector only if matrix is 1xN or Nx1
    Vector<T> to_vector() const {
        assert(m_rows == 1 || m_cols == 1);
        // get all values to vector
        std::vector<T> _;
        for (auto& [_k, _v] : m_data) {
            _.emplace_back(_v);
        }
        return Vector<T>(_);
    }

    static decltype(auto) indentity(size_t size) {
        Matrix result(size, size);
        for (size_t i = 0; i < size; ++i) {
            result(i, i) = 1;
        }
        return result;
    }

private:
    size_t m_rows;
    size_t m_cols;
    std::map<TKey, T> m_data;
};



// ===================================================
/**
 * concatenate matrices horizontally if axis = 1, vertically if axis = 0
 */
template <typename T>
decltype(auto) concact_matrices(std::initializer_list<Matrix<T>> matrices, int axis) {
    assert(matrices.size() > 0);
    auto new_rc = matrices.begin()->shape();
    axis == 0 ? new_rc.first = 0 : new_rc.second = 0;
    // check if all matrices have the same row or column count
    for (auto& m : matrices) {
        assert(m.shape(1 - axis) == matrices.begin()->shape(1 - axis));
        (axis == 0 ? new_rc.first : new_rc.second) += m.shape(axis);
    }
    Matrix<T> result(new_rc.first, new_rc.second);

    auto _last = 0;
    for (auto& m: matrices) {
        auto _shape = m.shape();
        for (auto i = 0; i < _shape.first; ++i) {
            for (auto j = 0; j < _shape.second; ++j) {
                if (axis == 0) {
                    result(i + _last, j) = m.at(i, j);
                } else {
                    result(i, j + _last) = m.at(i, j);
                }
            }
        }
        _last += (axis == 0 ? _shape.first : _shape.second);
    }

    return result;
}

// find max element in a vector or matrix
template <class Mx>
decltype(auto) mx_max(const Mx& vm) {
    auto raw_data = vm.data();
    return *std::max_element(raw_data.begin(), raw_data.end());
}

}
#endif