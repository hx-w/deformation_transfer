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
    A(0, 0) = 1, A(0, 1) = 2, A(0, 2) = 3;
    A(1, 0) = 2, A(1, 1) = 5, A(1, 2) = 7;
    A(2, 0) = 3, A(2, 1) = 5, A(2, 2) = 3;
    SpMatrixXd L, U;
    A.LU_decompose(L, U);
    cout << "L: " << endl;
    cout << L.at(0, 0) << " " << L.at(0, 1) << " " << L.at(0, 2) << endl;
    cout << L.at(1, 0) << " " << L.at(1, 1) << " " << L.at(1, 2) << endl;
    cout << L.at(2, 0) << " " << L.at(2, 1) << " " << L.at(2, 2) << endl;
    cout << "U: " << endl;
    cout << U.at(0, 0) << " " << U.at(0, 1) << " " << U.at(0, 2) << endl;
    cout << U.at(1, 0) << " " << U.at(1, 1) << " " << U.at(1, 2) << endl;
    cout << U.at(2, 0) << " " << U.at(2, 1) << " " << U.at(2, 2) << endl;


    auto r = L * U;
    cout << "R: " << endl;
    cout << r.at(0, 0) << " " << r.at(0, 1) << " " << r.at(0, 2) << endl;
    cout << r.at(1, 0) << " " << r.at(1, 1) << " " << r.at(1, 2) << endl;
    cout << r.at(2, 0) << " " << r.at(2, 1) << " " << r.at(2, 2) << endl;

    // solve
    auto b = SpMatrixXd::identity(3);
    SpMatrixXd x;
    cout << "solve: " << endl;
    A.solve(b, x);
    cout << "X: " << endl;
    cout << x.at(0, 0) << " " << x.at(0, 1) << " " << x.at(0, 2) << endl;
    cout << x.at(1, 0) << " " << x.at(1, 1) << " " << x.at(1, 2) << endl;
    cout << x.at(2, 0) << " " << x.at(2, 1) << " " << x.at(2, 2) << endl;

    auto inv = A.inverse();
    cout << "INV: " << endl;
    cout << inv.at(0, 0) << " " << inv.at(0, 1) << " " << inv.at(0, 2) << endl;
    cout << inv.at(1, 0) << " " << inv.at(1, 1) << " " << inv.at(1, 2) << endl;
    cout << inv.at(2, 0) << " " << inv.at(2, 1) << " " << inv.at(2, 2) << endl;

    // test save/load
    inv.save("A.obj");
    SpMatrixXd A2;
    A2.load("A.obj");
    cout << "A2: " << endl;
    cout << A2.at(0, 0) << " " << A2.at(0, 1) << " " << A2.at(0, 2) << endl;
    cout << A2.at(1, 0) << " " << A2.at(1, 1) << " " << A2.at(1, 2) << endl;
    cout << A2.at(2, 0) << " " << A2.at(2, 1) << " " << A2.at(2, 2) << endl;

    return 0;
}