#ifndef SOLVER_H
#define SOLVER_H

/**
 * solve Ax=b equation with LU decompose
 */

// #include "matrix.hpp"
#include "matrix_proxy.hpp"

namespace MeshLib {

class Solver {
public:
    Solver(const SpMatrixXd& A);
    Solver(const Solver&) = delete;

    void solve(const SpMatrixXd& b, SpMatrixXd& x);

private:
    void __lu_decompose();

private:
    SpMatrixXd _A;
    SpMatrixXd _L;
    SpMatrixXd _U;
};
}

#endif
