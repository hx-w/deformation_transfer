// implement correspondences

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>

#include "correspondace.h"
#include "meshlib/solver.h"
#include "meshlib/matrix_proxy.hpp"

using namespace std;
using namespace MeshLib;

void compute_correspondence(
    Mesh& src_mesh, Mesh& tgt_mesh,
    const MatrixXi& markers, MatrixXi& mappings
) {
    // coffecient
    auto Wi = 0.001;
    auto Ws = 1.0;

    vector<vector<size_t>> face_adj_list;
    src_mesh.get_triangle_adj(face_adj_list);

    // to 4d
    src_mesh.to_4d();
    tgt_mesh.to_4d();

    // if file find, load it
    SpMatrixXd AtA, AtB;
    if (ifstream(".cache/AtA.mat")) {
        AtA.load(".cache/AtA.mat");
        AtB.load(".cache/AtB.mat");

        cout << "load from cache: AtA" << endl;
    }
    else {
        SpMatrixXd sp_AE, sp_B, sp_At;
        if (ifstream(".cache/AE.mat")) {
            sp_AE.load(".cache/AE.mat");
            sp_B.load(".cache/B.mat");
            sp_At.load(".cache/At.mat");
            cout << "load from cache: AE" << endl;
        }
        else {
            MatrixXd AEi, Bi, AEs, Bs, AEc, Bc;

            vector<MatrixXd> inv_hat_list;
            src_mesh.get_inv_hat(inv_hat_list);
        

            construct_identity_cost(inv_hat_list, AEi, Bi);
            construct_smoothness_cost(inv_hat_list, face_adj_list, AEs, Bs);

            // substract markers from AE, B
            apply_markers(AEi, Bi, tgt_mesh, markers);
            apply_markers(AEs, Bs, tgt_mesh, markers);

            // ignore AEc
            auto AE = concact_matrices({AEi * Wi, AEs * Ws}, 0);
            auto B = concact_matrices({Bi * Wi, Bs * Ws}, 0);
            to_sparse(AE, sp_AE);
            to_sparse(B, sp_B);
            sp_AE.save(".cache/AE.mat");
            sp_B.save(".cache/B.mat");

            sp_At = sp_AE.transpose();
            sp_At.save(".cache/At.mat");
        }

        auto print = [](auto& mat, auto& name) {
            cout << name << ": " << mat.rows() << "x" << mat.cols() << endl;
            for (auto r = 0; r < 10; ++r) {
                for (auto c = 0; c < 10; ++c) {
                    cout << mat.at(r, c) << " ";
                }
                cout << endl;
            }
        };
        print(sp_AE, "AE");
        cout << sp_AE.at(186844, 0) << endl;

        // print time now: %H:%M:%S
        auto now = chrono::system_clock::now();
        auto in_time_t = chrono::system_clock::to_time_t(now);
        cout << "start multi: " << ctime(&in_time_t);

        AtA = sp_At * sp_AE;

        AtB = sp_At * sp_B;

        AtA.save(".cache/AtA.mat");
        AtB.save(".cache/AtB.mat");
    }

    if (AtA.is_diagonally_dominant()) {
        cout << "diagonally dominant" << endl;
    }
    else {
        cout << "not diagonally dominant" << endl;
    }

    for (int i  = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            cout << AtA(i, j) << " ";
        }
        cout << endl;
    }

    SpMatrixXd sp_X;
    AtA.solve(AtB, sp_X);

    // 10 rows print
    for (auto i = 0; i < 10; i++) {
        cout << sp_X.at(i, 0) << " " << sp_X.at(i, 1) << " " << sp_X.at(i, 2)  << endl;
    }

    // // print time now: %H:%M:%S
    auto now = chrono::system_clock::now();
    auto in_time_t = chrono::system_clock::to_time_t(now);
    cout << "end solving: " << ctime(&in_time_t);
    sp_X.save(".cache/correspondence.mat");
    // sp_X.save_txt(".cache/correspondence.txt");

    // revert markers
    MatrixXd ds_X;
    to_dense(sp_X, ds_X);
    MatrixXd X_revert;
    revert_markers(ds_X, tgt_mesh, markers, X_revert);

    // print time now: %H:%M:%S
    now = chrono::system_clock::now();
    in_time_t = chrono::system_clock::to_time_t(now);
    cout << "end revert: " << ctime(&in_time_t);

    // to 3d
    src_mesh.vertices() = X_revert;
    src_mesh.to_3d();
    src_mesh.save("test.obj");
}


void construct_identity_cost(
    const vector<MatrixXd>& inv_hat_V,
    MatrixXd& AEi, MatrixXd& Bi
) {
    // transpose each element in inv_hat_V, and concat along column
    AEi = inv_hat_V[0].transpose();
    Bi = MatrixXd::identity(3);
    for (auto i = 1; i < inv_hat_V.size(); ++i) {
        AEi.append(inv_hat_V[i].transpose(), 0);
        Bi.append(MatrixXd::identity(3), 0);
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

    auto AEm = AE.slice(src_marked, 1);
    AE = AE.slice(src_unmarked, 1);

    auto Xm = tgt_mesh.vertices().slice(tgt_marked, 0);

    cout << "AEm: " << AEm.rows() << " " << AEm.cols() << endl;
    cout << "Xm: " << Xm.rows() << " " << Xm.cols() << endl;

    B = B - AEm * Xm;
}


void revert_markers(
    const MatrixXd& X,
    Mesh& tgt_mesh,
    const MatrixXi& markers,
    MatrixXd& X_reverted
) {
    auto src_marked = markers.col(0).data();
    auto tgt_marked = markers.col(1).data();
    auto vdim = X.rows() + src_marked.size();

    auto src_unmarked = vector<size_t>();
    for (auto vi = 0; vi < vdim; ++vi) {
        if (find(src_marked.begin(), src_marked.end(), vi) == src_marked.end()) {
            src_unmarked.emplace_back(vi);
        }
    }

    cout << "unmarked: " << src_unmarked.size() << " X.rows " << X.rows() << endl;
    assert(src_unmarked.size() == X.rows());
    X_reverted = MatrixXd(vdim, 3);
    for (int i = 0; i < src_unmarked.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            X_reverted(src_unmarked[i], j) = X.at(i, j);
        }
    }

    for (int i = 0; i < src_marked.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            X_reverted(src_marked[i], j) = tgt_mesh.vertices().at(tgt_marked[i], j);
        }
    }
}