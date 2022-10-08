// test meshlib/solver.h

#include <iostream>

#include "../src/meshlib/solver.h"

using namespace std;
using namespace MeshLib;

int main() {
    MatrixXd A(3, 3);
    A(0, 0) = 1, A(0, 1) = 1, A(0, 2) = 2;
    A(1, 0) = 2, A(1, 1) = 3, A(1, 2) = 6;
    A(2, 0) = 3, A(2, 1) = 8, A(2, 2) = 9;

    MatrixXd b(3, 1);
    b(0, 0) = 1;
    b(1, 0) = 2;
    b(2, 0) = 3;

    MatrixXd x(3, 1);

    cout << "start solver" << endl;
    Solver solver(A);
    solver.solve(b, x);

    cout << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << endl;


    return 0;
}