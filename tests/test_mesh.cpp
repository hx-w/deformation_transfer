// test meshlib/mesh.h

#include <iostream>

#include "../src/meshlib/mesh.h"

using namespace std;
using namespace MeshLib;

int main() {
    Mesh mesh;
    mesh.load("../static/source_mesh/reference.obj");

    mesh.to_4d();
    // mesh.to_3d();
    vector<MatrixXd> inv_hat_list;
    mesh.get_inv_hat(inv_hat_list);

    return 0;
}