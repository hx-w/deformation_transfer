#include "../src/meshlib/matrix_proxy.hpp"

#include <iostream>

using namespace std;
using namespace MeshLib;

int main() {
    using TripList = std::vector<Triplet<double>>;
    TripList trip_list;
    trip_list.emplace_back(Triplet{0, 0, 1.0});
    trip_list.emplace_back(Triplet{0, 1, 2.0});
    trip_list.emplace_back(Triplet{1, 0, 3.0});
    trip_list.emplace_back(Triplet{1, 1, 4.0});
    SparseMatrix<double> spm(2, 2, trip_list);

    cout << "SPM: " << endl;
    cout << spm(0, 0) << " " << spm(0, 1) << endl;
    cout << spm(1, 0) << " " << spm(1, 1) << endl;

    auto spm2 = spm + spm;
    cout << "SPM2: " << endl;
    cout << spm2(0, 0) << " " << spm2(0, 1) << endl;
    cout << spm2(1, 0) << " " << spm2(1, 1) << endl;


    auto spm3 = spm2 * spm;
    cout << "SPM3: " << endl;
    cout << spm3(0, 0) << " " << spm3(0, 1) << endl;
    cout << spm3(1, 0) << " " << spm3(1, 1) << endl;

    auto spm4 = spm3 * 2.0 - spm2;
    cout << "SPM4: " << endl;
    cout << spm4(0, 0) << " " << spm4(0, 1) << endl;
    cout << spm4(1, 0) << " " << spm4(1, 1) << endl;

    auto spm_inv = spm.inverse();
    cout << "SPM_INV: " << endl;
    cout << spm_inv(0, 0) << " " << spm_inv(0, 1) << endl;
    cout << spm_inv(1, 0) << " " << spm_inv(1, 1) << endl;

    // LU decompose
    SpMatrixXd A(3, 3);
    A(0, 0) = 4, A(0, 1) = 12, A(0, 2) = -16;
    A(1, 0) = 12, A(1, 1) = 37, A(1, 2) = -43;
    A(2, 0) = -16, A(2, 1) = -43, A(2, 2) = 98;
    // SpMatrixXd L, U;
    // A.LU_decompose(L, U);
    // cout << "L: " << endl;
    // cout << L.at(0, 0) << " " << L.at(0, 1) << " " << L.at(0, 2) << endl;
    // cout << L.at(1, 0) << " " << L.at(1, 1) << " " << L.at(1, 2) << endl;
    // cout << L.at(2, 0) << " " << L.at(2, 1) << " " << L.at(2, 2) << endl;
    // cout << "U: " << endl;
    // cout << U.at(0, 0) << " " << U.at(0, 1) << " " << U.at(0, 2) << endl;
    // cout << U.at(1, 0) << " " << U.at(1, 1) << " " << U.at(1, 2) << endl;
    // cout << U.at(2, 0) << " " << U.at(2, 1) << " " << U.at(2, 2) << endl;


    // auto r = L * U;
    // cout << "R: " << endl;
    // cout << r.at(0, 0) << " " << r.at(0, 1) << " " << r.at(0, 2) << endl;
    // cout << r.at(1, 0) << " " << r.at(1, 1) << " " << r.at(1, 2) << endl;
    // cout << r.at(2, 0) << " " << r.at(2, 1) << " " << r.at(2, 2) << endl;

    // solve
    auto b = SpMatrixXd::identity(3);
    SpMatrixXd x;
    cout << "solve: " << endl;
    A.solve(b, x);
    cout << "X: " << endl;
    cout << x.at(0, 0) << " " << x.at(0, 1) << " " << x.at(0, 2) << endl;
    cout << x.at(1, 0) << " " << x.at(1, 1) << " " << x.at(1, 2) << endl;
    cout << x.at(2, 0) << " " << x.at(2, 1) << " " << x.at(2, 2) << endl;

    cout << "check: " << endl;
    auto A_inv = A.inverse();
    cout << A_inv.at(0, 0) << " " << A_inv.at(0, 1) << " " << A_inv.at(0, 2) << endl;
    cout << A_inv.at(1, 0) << " " << A_inv.at(1, 1) << " " << A_inv.at(1, 2) << endl;
    cout << A_inv.at(2, 0) << " " << A_inv.at(2, 1) << " " << A_inv.at(2, 2) << endl;


    // Jacobi iteration
    SpMatrixXd A2(4, 4);
    A2(0, 0) = 8, A2(0, 1) = -3, A2(0, 2) = 2, A2(0, 3) = -1;
    A2(1, 0) = 4, A2(1, 1) = 11, A2(1, 2) = -1, A2(1, 3) = 3;
    A2(2, 0) = 6, A2(2, 1) = 3, A2(2, 2) = 12, A2(2, 3) = 2;
    A2(3, 0) = 2, A2(3, 1) = -10, A2(3, 2) = 8, A2(3, 3) = 15;
    SpMatrixXd b2(4, 3);
    b2(0, 0) = 20, b2(0, 1) = 33, b2(0, 2) = 36;
    b2(1, 0) = 36, b2(1, 1) = 66, b2(1, 2) = 76;
    b2(2, 0) = 31, b2(2, 1) = 57, b2(2, 2) = 69;
    b2(3, 0) = 29, b2(3, 1) = 32, b2(3, 2) = 38;

    SpMatrixXd x_;
    A2.solve_(b2, x_);
    x_ = x_;
    cout << "X_: " << endl;
    cout << x_.at(0, 0) << " " << x_.at(0, 1) << " " << x_.at(0, 2) << endl;
    cout << x_.at(1, 0) << " " << x_.at(1, 1) << " " << x_.at(1, 2) << endl;
    cout << x_.at(2, 0) << " " << x_.at(2, 1) << " " << x_.at(2, 2) << endl;
    cout << x_.at(3, 0) << " " << x_.at(3, 1) << " " << x_.at(3, 2) << endl;

    auto target = A2 * x_;
    cout << "target: " << endl;
    cout << target.at(0, 0) << " " << target.at(0, 1) << " " << target.at(0, 2) << endl;
    cout << target.at(1, 0) << " " << target.at(1, 1) << " " << target.at(1, 2) << endl;
    cout << target.at(2, 0) << " " << target.at(2, 1) << " " << target.at(2, 2) << endl;
    cout << target.at(3, 0) << " " << target.at(3, 1) << " " << target.at(3, 2) << endl;

    
    // test save/load
    // inv.save("A.obj");
    // SpMatrixXd A2;
    // A2.load("A.obj");
    // cout << "A2: " << endl;
    // cout << A2.at(0, 0) << " " << A2.at(0, 1) << " " << A2.at(0, 2) << endl;
    // cout << A2.at(1, 0) << " " << A2.at(1, 1) << " " << A2.at(1, 2) << endl;
    // cout << A2.at(2, 0) << " " << A2.at(2, 1) << " " << A2.at(2, 2) << endl;

    return 0;
}