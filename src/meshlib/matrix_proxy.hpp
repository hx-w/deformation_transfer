#ifndef MATRIX_PROXY_HPP
#define MATRIX_PROXY_HPP

#include "matrix.hpp"
#include "spmatrix.hpp"

namespace MeshLib {

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
    std::cout << "trip_list size: " << trip_list.size() << std::endl;
    spmat = SparseMatrix<T>(mat.rows(), mat.cols(), trip_list);
}


}

#endif
