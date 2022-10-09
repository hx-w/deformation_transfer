// implement solver

#include "solver.h"

#include <iostream>

namespace MeshLib {
using namespace std;

Solver::Solver(const MatrixXd& A) : _A(A) {
    __lu_decompose();
}

/**
 * solve Ax = b
 * A: shape(m, n) LU decomposed
 * x: shape(n, k)
 * b: shape(m, k)
 */
void Solver::solve(const MatrixXd& b, MatrixXd& x) {
    // solve Ly = b
    auto y = MatrixXd(b.rows(), b.cols());
    for (auto i = 0; i < b.rows(); ++i) {
        for (auto j = 0; j < b.cols(); ++j) {
            auto sum = 0.0;
            for (auto k = 0; k < i; ++k) {
                sum += _L(i, k) * y(k, j);
            }
            y(i, j) = (b(i, j) - sum) / _L(i, i);
     
        }
    }
    // solve Ux = y
    x = MatrixXd(b.rows(), b.cols());
    for (int i = b.rows() - 1; i >= 0; --i) { // i must be int
        for (auto j = 0; j < b.cols(); ++j) {
            auto sum = 0.0;
            for (auto k = i + 1; k < b.rows(); ++k) {
                sum += _U(i, k) * x(k, j);
            }
            x(i, j) = (y(i, j) - sum) / _U(i, i);
        }
    }
}

void Solver::__lu_decompose() {
    auto n = _A.rows();
    _L = MatrixXd(n, n);
    _U = MatrixXd(n, n);
    // LU decomposition
    for (auto i = 0; i < n; ++i) {
        for (auto j = 0; j < n; ++j) {
            _U(i, j) = _A(i, j);
            for (auto k = 0; k < i; ++k) {
                _U(i, j) -= _L(i, k) * _U(k, j);
            }
        }
        for (auto j = i; j < n; ++j) {
            _L(j, i) = _A(j, i);
            for (auto k = 0; k < i; ++k) {
                _L(j, i) -= _L(j, k) * _U(k, i);
            }
            _L(j, i) /= _U(i, i);
        }
    }

    // show L and U
    // cout << "L:" << endl;
    // for (auto i = 0; i < n; ++i) {
    //     for (auto j = 0; j < n; ++j) {
    //         cout << _L(i, j) << " ";
    //     }
    //     cout << endl;
    // }

    // cout << "U:" << endl;
    // for (auto i = 0; i < n; ++i) {
    //     for (auto j = 0; j < n; ++j) {
    //         cout << _U(i, j) << " ";
    //     }
    //     cout << endl;
    // }
}

}