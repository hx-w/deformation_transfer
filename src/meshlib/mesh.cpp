// implement of mesh.h

#include "mesh.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace MeshLib {
using namespace std;

using TEdge = pair<int, int>;

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

bool Mesh::save(const string& filename) {
    if (filename.find(".obj") == string::npos) {
        return false;
    }

    ofstream fout(filename);
    if (!fout.is_open()) {
        return false;
    }

    // only save v and f
    for (auto i = 0; i < m_vertices.rows(); ++i) {
        fout << "v " << m_vertices(i, 0) << " " << m_vertices(i, 1) << " " << m_vertices(i, 2) << endl;
    }

    for (auto i = 0; i < m_faces.rows(); ++i) {
        fout << "f " << m_faces(i, 0) + 1 << " " << m_faces(i, 1) + 1 << " " << m_faces(i, 2) + 1 << endl;
    }

    fout.close();
    return true;
}

bool Mesh::__load_obj(const string& filename) {
    using DimXd = vector<double>;
    using DimXi = vector<size_t>;
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
            split_words(line, words, ' ');
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

    ifs.close();
    return true;
}

void Mesh::split_words(const string& line, vector<string>& words, const char delim) {
    stringstream ss(line);
    string word = "";
    while (getline(ss, word, delim)) {
        words.emplace_back(word);
    }
}


void Mesh::to_4d() {
    if (m_faces.cols() == 4) {
        return;
    }
    // append new vertices(v4), [0, old_verts_num) is old vertices
    auto old_verts_num = m_vertices.rows();
    // for all triangles, compute v4 = (v2 - v1) x (v3 - v1)
    for (auto i = 0; i < m_faces.rows(); ++i) {
        auto v1 = m_vertices.row(m_faces.at(i, 0));
        auto v2 = m_vertices.row(m_faces.at(i, 1));
        auto v3 = m_vertices.row(m_faces.at(i, 2));
        auto v4 = (v2 - v1).cross(v3 - v1).normalize();
        m_vertices.append_row(v4);
    }

    // change faces to 4d by concat [old_verts_num, m_vertices.rows()) to m_faces by column
    vector<size_t> new_verts_idx;
    for (auto i = old_verts_num; i < m_vertices.rows(); ++i) {
        new_verts_idx.emplace_back(i);
    }
    m_faces.append_col(new_verts_idx);
}

void Mesh::to_3d() {
    if (m_faces.cols() == 3) {
        return;
    }
    const auto last_verts_num = m_faces.rows();
    // remove last_verts_num vertices
    m_vertices = m_vertices.slice(0, m_vertices.rows() - last_verts_num, 0, m_vertices.cols());
    m_faces = m_faces.slice(0, m_faces.rows(), 0, m_faces.cols() - 1);
}

void Mesh::get_inv_hat(vector<MatrixXd>& inv_hat_list) const {
    if (m_faces.cols() != 4) {
        return;
    }
    inv_hat_list.clear();

    // for each triangle, V = [v2 - v1, v3 - v1, v4 - v1]
    // V^hat-1 = C * V^{-1}
    for (auto i = 0; i < m_faces.rows(); ++i) {
        auto v1 = m_vertices.row(m_faces.at(i, 0));
        auto v2 = m_vertices.row(m_faces.at(i, 1));
        auto v3 = m_vertices.row(m_faces.at(i, 2));
        auto v4 = m_vertices.row(m_faces.at(i, 3));
        auto V = concact_matrices({
            (v2 - v1).to_matrix(),
            (v3 - v1).to_matrix(),
            (v4 - v1).to_matrix(),
        }, 1);
        auto inv_V = V.inverse();
        // compute sum of each row of inv_V
        vector<double> _sum_cols{};
        for (auto j = 0; j < 3; ++j) {
            _sum_cols.emplace_back(inv_V.at(0, j) + inv_V.at(1, j) + inv_V.at(2, j));
        }
        auto inv_hat_V = MatrixXd(m_vertices.rows(), 3);
        for (auto j = 0; j < 3; ++j) {
            inv_hat_V(m_faces.at(i, 0), j) = -_sum_cols[j];
            inv_hat_V(m_faces.at(i, 1), j) = inv_V.at(0, j);
            inv_hat_V(m_faces.at(i, 2), j) = inv_V.at(1, j);
            inv_hat_V(m_faces.at(i, 3), j) = inv_V.at(2, j);
        }

        inv_hat_list.emplace_back(inv_hat_V);
    }
}

void Mesh::get_triangle_adj(vector<vector<size_t>>& adj_list) const {
    if (m_faces.cols() != 3) {
        return;
    }
    adj_list.clear();
    adj_list.resize(m_faces.rows());

    // edge to triangle idx
    map<TEdge, vector<size_t>> edge_to_tri;

    for (auto i = 0; i < m_faces.rows(); ++i) {
        // sort triangle vertices idx
        vector<size_t> _tri{m_faces.at(i, 0), m_faces.at(i, 1), m_faces.at(i, 2)};
        sort(_tri.begin(), _tri.end());

        auto e1 = make_pair(_tri[0], _tri[1]);
        auto e2 = make_pair(_tri[0], _tri[2]);
        auto e3 = make_pair(_tri[1], _tri[2]);

        edge_to_tri[e1].emplace_back(i);
        edge_to_tri[e2].emplace_back(i);
        edge_to_tri[e3].emplace_back(i);
    }

    for (auto& [edge, tri_idx_list] : edge_to_tri) {
        // if size > 2, geometry is not a manifold,
        // if size == 1, it is a boundary edge, ignore it
        if (tri_idx_list.size() == 2) {
            adj_list[tri_idx_list[0]].emplace_back(tri_idx_list[1]);
            adj_list[tri_idx_list[1]].emplace_back(tri_idx_list[0]);
        }
    }
}


}