// implement correspondences

#include <iostream>
#include <fstream>
#include <string>

#include "correspondace.h"

using namespace std;
using namespace MeshLib;

void compute_correspondence(
    Mesh& src_mesh, Mesh& tgt_mesh,
    const MatrixXi& markers, MatrixXi& mappings
) {
    // to 4d
    src_mesh.to_4d();
    tgt_mesh.to_4d();

    auto face_adj_list = vector<vector<int>>();
    src_mesh.get_triangle_adj(face_adj_list);

    auto inv_hat_list = vector<MatrixXd>();
    src_mesh.get_inv_hat(inv_hat_list);

    MatrixXd AEi, Bi, AEs, Bs, AEc, Bc;

    construct_indentity_cost(inv_hat_list, AEi, Bi);
    construct_smoothness_cost(inv_hat_list, face_adj_list, AEs, Bs);
}


void construct_indentity_cost(
    const vector<MatrixXd>& inv_hat_V,
    MatrixXd& AEi, MatrixXd& Bi
) {
    // transpose each element in inv_hat_V, and concat along column
    AEi = inv_hat_V[0].transpose();
    Bi = MatrixXd::indentity(3);
    for (auto i = 1; i < inv_hat_V.size(); ++i) {
        AEi.append(inv_hat_V[i].transpose(), 0);
        Bi.append(MatrixXd::indentity(3), 0);
    }
    cout << "AEi: " << AEi.rows() << " " << AEi.cols() << endl;
    cout << "Bi: " << Bi.rows() << " " << Bi.cols() << endl;
    AEi.transpose() * AEi;
}


void construct_smoothness_cost(
    const vector<MatrixXd>& inv_hat_V,
    const vector<vector<int>>& face_adj_list,
    MatrixXd& AEs, MatrixXd& Bs
) {

}


void read_markers(const string& filename, MatrixXi& markers) {
    // read markers format:
    // source_markers:target_markers
    vector<vector<int>> markers_list;
    ifstream fin(filename);
    if (!fin.is_open()) {
        cout << "Error: cannot open file " << filename << endl;
        return;
    }
    string line = "";

    while (getline(fin, line)) {
        auto words = vector<string>();
        Mesh::split_words(line, words, ':');

        markers_list.emplace_back(vector<int>{stoi(words[0]), stoi(words[1])});
    }

    markers = MatrixXi(markers_list);
    cout << "Makers loaded: " << markers.rows() << ", " << markers.cols() << endl; 
}