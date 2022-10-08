// implement solver

#include "solver.h"

#include <iostream>

namespace MeshLib {
using namespace std;

Solver::Solver(const MatrixXd& A) : _A(A) {
    __lu_decompose();
}

void Solver::solve(const MatrixXd& b, MatrixXd& x) {
    auto n = _A.rows();
    MatrixXd y(n, 1);
    for (auto i = 0; i < n; ++i) {
        y(i, 0) = b(i, 0);
        for (auto j = 0; j < i; ++j) {
            y(i, 0) -= _L(i, j) * y(j, 0);
        }
    }

    for (int i = n - 1; i >= 0; --i) { // i must be signed int
        x(i, 0) = y(i, 0);
        for (auto j = i + 1; j < n; ++j) {
            x(i, 0) -= _U(i, j) * x(j, 0);
        }
        x(i, 0) /= _U(i, i);
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