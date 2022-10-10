#include "../src/meshlib/spmatrix.hpp"

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



    return 0;
}