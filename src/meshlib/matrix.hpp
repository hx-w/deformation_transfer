// matrix
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

namespace MeshLib {


template <typename T>
class Matrix {
public:
    Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
        data_.resize(rows_ * cols_);
    }

    Matrix(size_t rows, size_t cols, const T& initial) : rows_(rows_), cols_(cols_) {
        data_.resize(rows_ * cols_, initial);
    }

    Matrix(const Matrix& other) : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows_ = other.rows_;
            cols_ = other.cols_;
            data_.assign(other.data_.begin(), other.data_.end());
        }
        return *this;
    }

    Matrix(Matrix&& other) : rows_(other.rows_), cols_(other.cols_), data_(std::move(other.data_)) {}

    // fetch values
    T& operator()(size_t row, size_t col) {
        assert(row < rows_ && col < cols_);
        return data_[row * cols_ + col];
    }

    const T& operator()(size_t row, size_t col) const {
        assert(row < rows_ && col < cols_);
        return data_[row * cols_ + col];
    }

    auto rows() const { return rows_; }
    auto cols() const { return cols_; }
    auto shape() const { return std::make_pair(rows_, cols_); }
    auto shape(size_t dim) const { return dim == 0 ? rows_ : cols_; }

    // addition
    Matrix operator+(const Matrix& other) const {
        assert(rows_ == other.rows_ && cols_ == other.cols_);
        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }

    // subtraction
    Matrix operator-(const Matrix& other) const {
        assert(rows_ == other.rows_ && cols_ == other.cols_);
        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

    // multiplication
    Matrix operator*(const Matrix& other) const {
        assert(cols_ == other.rows_);
        Matrix result(rows_, other.cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                for (size_t k = 0; k < cols_; ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        return result;
    }

    // scalar multiplication
    Matrix operator*(const T& scalar) const {
        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }
        return result;
    }

    // slice
    Matrix slice(size_t row_start, size_t row_end, size_t col_start, size_t col_end) const {
        assert(row_start < row_end && col_start < col_end);
        Matrix result(row_end - row_start, col_end - col_start);
        for (size_t i = row_start; i < row_end; ++i) {
            for (size_t j = col_start; j < col_end; ++j) {
                result(i - row_start, j - col_start) = (*this)(i, j);
            }
        }
        return result;
    }

    // transpose
    Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    // determinant
    T determinant() const {
        assert(rows_ == cols_);
        if (rows_ == 1) {
            return (*this)(0, 0);
        }
        T result = 0;
        for (size_t i = 0; i < cols_; ++i) {
            result += (*this)(0, i) * cofactor(0, i);
        }
        return result;
    }

    // cofactor
    T cofactor(size_t row, size_t col) const {
        assert(rows_ == cols_);
        Matrix minor(rows_ - 1, cols_ - 1);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                if (i != row && j != col) {
                    minor(i < row ? i : i - 1, j < col ? j : j - 1) = (*this)(i, j);
                }
            }
        }
        return minor.determinant() * ((row + col) % 2 == 0 ? 1 : -1);
    }

    // inverse
    decltype(auto) inverse() const {
        assert(rows_ == cols_);
        T det = determinant();
        assert(det != 0);
        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = cofactor(i, j) / det;
            }
        }
        return result;
    }

private:
    size_t rows_;
    size_t cols_;
    std::vector<T> data_;
};


// ===================================================

decltype(auto) indentity(size_t size) {
    Matrix<double> result(size, size);
    for (size_t i = 0; i < size; ++i) {
        result(i, i) = 1;
    }
    return result;
}

/**
 * concatenate matrices horizontally if axis = 1, vertically if axis = 0
 */
template <typename T>
decltype(auto) concact_matrices(std::initializer_list<Matrix<T>> matrices, size_t axis) {
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
                    result(i + _last, j) = m(i, j);
                } else {
                    result(i, j + _last) = m(i, j);
                }
            }
        }
        _last += (axis == 0 ? _shape.first : _shape.second);
    }

    return result;
}

}

#endif