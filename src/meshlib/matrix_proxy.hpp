#ifndef MATRIX_PROXY_HPP
#define MATRIX_PROXY_HPP

#include "matrix.hpp"
#include "spmatrix.hpp"

namespace MeshLib {

using MatrixXd = Matrix<double>;
using MatrixXi = Matrix<size_t>;
using SpMatrixXd = SparseMatrix<double>;
using SpMatrixXi = SparseMatrix<size_t>;


template <class T>
void to_sparse(const Matrix<T>& mat, SparseMatrix<T>& spmat) {
    spmat.resize(mat.rows(), mat.cols());
    const auto& _raw = mat.raw();
    std::vector<Triplet<T>> trip_list;
    for (const auto& [_key, _val] : _raw) {
        trip_list.emplace_back(Triplet<T>{
            static_cast<int>(_key.first),
            static_cast<int>(_key.second),
            _val
        });
    }
    spmat = SparseMatrix<T>(mat.rows(), mat.cols(), trip_list);
    std::cout << "to_sparse finished" << std::endl;
}

template <class T>
void to_dense(const SparseMatrix<T>& spmat, Matrix<T>& mat) {
    spmat.to_dense(mat);
    std::cout << "to_dense finished" << std::endl;
}

}

#endif
