// test meshlib/matrix.hpp

#include <iostream>

#include "../src/meshlib/matrix.hpp"

using namespace std;
using namespace MeshLib;

int main() {
    Matrix<double> m(3, 3);

    m(0, 0) = 1;
    m(0, 1) = 1.1;
    m(0, 2) = 3.3;
    m(1, 0) = 4;
    m(1, 1) = 5;
    m(1, 2) = 6;
    m(2, 0) = 7;
    m(2, 1) = 8;
    m(2, 2) = 9;

    cout << m(0, 0) << " " << m(0, 1) << " " << m(0, 2) << endl;
    cout << m(1, 0) << " " << m(1, 1) << " " << m(1, 2) << endl;
    cout << m(2, 0) << " " << m(2, 1) << " " << m(2, 2) << endl;

    auto im = m.inverse();

    cout << im(0, 0) << " " << im(0, 1) << " " << im(0, 2) << endl;
    cout << im(1, 0) << " " << im(1, 1) << " " << im(1, 2) << endl;
    cout << im(2, 0) << " " << im(2, 1) << " " << im(2, 2) << endl;

    // test concatenate
    auto I = Matrix<double>::indentity(3).slice(0, 3, 0, 2);
    cout << I.rows() << " " << I.cols() << endl;

    auto m2 = concact_matrices({m, I}, 1);

    cout << "input concat:" << endl;
    // 3 x 5
    cout << m2(0, 0) << " " << m2(0, 1) << " " << m2(0, 2) << " " << m2(0, 3) << " " << m2(0, 4) << endl;
    cout << m2(1, 0) << " " << m2(1, 1) << " " << m2(1, 2) << " " << m2(1, 3) << " " << m2(1, 4) << endl;
    cout << m2(2, 0) << " " << m2(2, 1) << " " << m2(2, 2) << " " << m2(2, 3) << " " << m2(2, 4) << endl;

    cout << "----- append ------" << endl;
    m.append(I, 1);
    cout << m(0, 0) << " " << m(0, 1) << " " << m(0, 2) << " " << m(0, 3) << " " << m(0, 4) << endl;
    cout << m(1, 0) << " " << m(1, 1) << " " << m(1, 2) << " " << m(1, 3) << " " << m(1, 4) << endl;
    cout << m(2, 0) << " " << m(2, 1) << " " << m(2, 2) << " " << m(2, 3) << " " << m(2, 4) << endl;


    // ===========================
    Vector<double> v({1, 2, 3});
    cout << v[0] << " " << v[1] << " " << v[2] << endl;

    auto v2 = v.normalize();
    cout << v2[0] << " " << v2[1] << " " << v2[2] << endl;

    auto v3 = v2.cross(v);
    cout << v3[0] << " " << v3[1] << " " << v3[2] << endl << endl;;


    auto vm = v.to_matrix(0);
    cout << vm(0, 0) << " " << vm(1, 0) << " " << vm(2, 0) << endl;

    auto mv = vm.to_vector();
    cout << mv[0] << " " << mv[1] << " " << mv[2] << endl;

    cout << mx_max(vm) << endl;

    return 0;
}