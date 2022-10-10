// implement solver

#include "solver.h"

#include <iostream>

namespace MeshLib {
using namespace std;

Solver::Solver(const SpMatrixXd& A) : _A(A) {
    __lu_decompose();
}

/**
 * solve Ax = b
 * A: shape(m, n) LU decomposed
 * x: shape(n, k)
 * b: shape(m, k)
 */
void Solver::solve(const SpMatrixXd& b, SpMatrixXd& x) {
    // solve Ly = b
    auto y = SpMatrixXd(b.rows(), b.cols());
    for (auto i = 0; i < b.rows(); ++i) {
        for (auto j = 0; j < b.cols(); ++j) {
            auto sum = 0.0;
            for (auto k = 0; k < i; ++k) {
                sum += _L(i, k) * y(k, j);
            }
            y(i, j) = (b.at(i, j) - sum) / _L(i, i);
     
        }
    }
    // solve Ux = y
    x = SpMatrixXd(b.rows(), b.cols());
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
    _L = SpMatrixXd(n, n);
    _U = SpMatrixXd(n, n);

    _A.LU_decompose(_L, _U);
}

}