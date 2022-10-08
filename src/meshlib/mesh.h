// mesh
#ifndef MESH_H
#define MESH_H

#include <string>
#include "matrix.hpp"

namespace MeshLib {

class Mesh {
public:
    Mesh() = default;
    Mesh(const std::string& filename);

    bool load(const std::string& filename);

    void to_4d();
    void to_3d();

    // get V^{\hat}^{-1} list, ordered by face index
    void get_inv_hat(std::vector<MatrixXd>& inv_hat_list) const;

private:
    bool __load_obj(const std::string& filename);
    void __split_words(
        const std::string& line,
        std::vector<std::string>& words,
        const char delim
    );

private:
    MatrixXd m_vertices;
    MatrixXi m_faces;    
};

}

#endif
