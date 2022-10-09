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
    vector<vector<size_t>> face_adj_list;
    src_mesh.get_triangle_adj(face_adj_list);

    // to 4d
    src_mesh.to_4d();
    tgt_mesh.to_4d();

    vector<MatrixXd> inv_hat_list;
    src_mesh.get_inv_hat(inv_hat_list);

    MatrixXd AEi, Bi, AEs, Bs, AEc, Bc;

    construct_indentity_cost(inv_hat_list, AEi, Bi);
    construct_smoothness_cost(inv_hat_list, face_adj_list, AEs, Bs);

    // substract markers from AE, B
    apply_markers(AEi, Bi, tgt_mesh, markers);
    apply_markers(AEs, Bs, tgt_mesh, markers);
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
}


void construct_smoothness_cost(
    const vector<MatrixXd>& inv_hat_V,
    const vector<vector<size_t>>& face_adj_list,
    MatrixXd& AEs, MatrixXd& Bs
) {
    AEs = MatrixXd(0, inv_hat_V[0].rows());

    for (auto i = 0; i < inv_hat_V.size(); ++i) {
        auto adj_faces = face_adj_list[i];
        for (const auto& j : adj_faces) {
            if (i == j) continue;
            auto diff = inv_hat_V[i] - inv_hat_V[j];
            AEs.append(diff.transpose(), 0);
        }
    }
    Bs = MatrixXd(AEs.rows(), 3);
    cout << "AEs: " << AEs.rows() << " " << AEs.cols() << endl;
    cout << "Bs: " << Bs.rows() << " " << Bs.cols() << endl;
}


void read_markers(const string& filename, MatrixXi& markers) {
    // read markers format:
    // source_markers:target_markers
    vector<vector<size_t>> markers_list;
    ifstream fin(filename);
    if (!fin.is_open()) {
        cout << "Error: cannot open file " << filename << endl;
        return;
    }
    string line = "";

    while (getline(fin, line)) {
        auto words = vector<string>();
        Mesh::split_words(line, words, ':');

        markers_list.emplace_back(vector<size_t>{stoul(words[0]), stoul(words[1])});
    }

    markers = MatrixXi(markers_list);
    cout << "Makers loaded: " << markers.rows() << ", " << markers.cols() << endl; 
}


void apply_markers(MatrixXd& AE, MatrixXd& B, Mesh& tgt_mesh, const MatrixXi& markers) {
    auto src_marked = markers.col(0).data();
    auto tgt_marked = markers.col(1).data();

    auto vdim = AE.cols();
    // compute src_unmarked
    auto src_unmarked = vector<size_t>();
    for (auto vi = 0; vi < vdim; ++vi) {
        if (find(src_marked.begin(), src_marked.end(), vi) == src_marked.end()) {
            src_unmarked.emplace_back(vi);
        }
    }

    cout << "source unmarked: " << src_unmarked.size() << endl;

    auto AEu = AE.slice(src_unmarked, 1);
    auto AEm = AE.slice(src_marked, 1);

    auto Xm = tgt_mesh.vertices().slice(tgt_marked, 0);

    cout << "AEm: " << AEm.rows() << " " << AEm.cols() << endl;
    cout << "Xm: " << Xm.rows() << " " << Xm.cols() << endl;

    B = B - AEm * Xm;
}
