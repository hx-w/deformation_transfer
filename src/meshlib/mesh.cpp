// implement of mesh.h

#include "mesh.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace MeshLib {
using namespace std;

Mesh::Mesh(const string& filename) {
    load(filename);
}

// only .obj file is supported
bool Mesh::load(const string& filename) {
    if (filename.find(".obj") != string::npos) {
        return __load_obj(filename);
    }
    return false;
}

bool Mesh::__load_obj(const string& filename) {
    using DimXd = vector<double>;
    using DimXi = vector<int>;
    vector<DimXd> _verts;
    vector<DimXi> _faces;

    ifstream ifs(filename);
    if (!ifs.is_open()) {
        return false;
    }

    string line = "";
    while (getline(ifs, line)) {
        if (line[0] == 'v' && line[1] == ' ') {
            double x = 0., y = 0., z = 0.;
            sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
            _verts.emplace_back(DimXd{x, y, z});
        } else if (line[0] == 'f') {
            vector<string> words;
            __split_words(line, words, ' ');
            DimXi _f{};
            for (auto i = 1; i < words.size(); ++i) {
                int idx = 0;
                sscanf(words[i].c_str(), "%d", &idx);
                _f.emplace_back(idx - 1);
            }
            _faces.emplace_back(_f);
        }
    }

    // vertices: (n, 3)
    // faces: (m, 3)
    m_vertices = MatrixXd(_verts);
    m_faces = MatrixXi(_faces);

    return true;
}

void Mesh::__split_words(const string& line, vector<string>& words, const char delim) {
    stringstream ss(line);
    string word = "";
    while (getline(ss, word, delim)) {
        words.emplace_back(word);
    }
}


}