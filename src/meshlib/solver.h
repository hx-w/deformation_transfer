#ifndef SOLVER_H
#define SOLVER_H

/**
 * solve Ax=b equation with LU decompose
 */

#include "matrix.hpp"

namespace MeshLib {

class Solver {
public:
    Solver(const MatrixXd& A);
    Solver(const Solver&) = delete;

    void solve(const MatrixXd& b, MatrixXd& x);

private:
    void __lu_decompose();

private:
    MatrixXd _A;
    MatrixXd _L;
    MatrixXd _U;
};
}

#endif
