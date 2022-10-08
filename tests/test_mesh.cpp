// test meshlib/mesh.h

#include <iostream>

#include "../src/meshlib/mesh.h"

using namespace std;
using namespace MeshLib;

int main() {
    Mesh mesh;
    mesh.load("../static/source_mesh/reference.obj");

    return 0;
}