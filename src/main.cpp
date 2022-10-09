/**
 * 1. Correspondence 
 * 2. Deformation Transfer
 */

#include <iostream>
#include "meshlib/mesh.h"
#include "correspondace.h"

using namespace std;
using namespace MeshLib;

int main() {
    Mesh src_mesh("../static/source_mesh/reference.obj");
    Mesh tgt_mesh("../static/target_mesh/reference.obj");

    MatrixXi markers;
    read_markers("../static/markers.txt", markers);

    MatrixXi mappings;
    compute_correspondence(src_mesh, tgt_mesh, markers, mappings);

    return 0;
}