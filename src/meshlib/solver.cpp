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
    cout << "solving Ly = b" << endl;
    SpMatrixXd y(b.rows(), b.cols());
    for (int i = 0; i < b.rows(); ++i) {
        for (int j = 0; j < b.cols(); ++j) {
            double sum = 0;
            for (int k = 0; k < i; ++k) {
                sum += _L.at(i, k) * y.at(k, j);
            }
            y(i, j) = b.at(i, j) - sum;
        }
    }

    std::cout << "y: " << std::endl;
    for (auto i = 0; i < 10; ++i) {
        for (auto j = 0; j < y.cols(); ++j) {
            std::cout << y.at(i, j) << " ";
        }
        std::cout << std::endl;
    }

    cout << "solving Ux = y" << endl;
    // solve Ux = y
    x = SpMatrixXd(b.rows(), b.cols());
    for (int i = b.rows() - 1; i >= 0; --i) {
        for (int j = 0; j < b.cols(); ++j) {
            double sum = 0;
            for (int k = i + 1; k < b.rows(); ++k) {
                sum += _U.at(i, k) * x.at(k, j);
            }
            x(i, j) = (y.at(i, j) - sum) / _U.at(i, i);
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