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

    bool save(const std::string& filename);
    bool load(const std::string& filename);


    void to_4d();
    void to_3d();

    // get V^{\hat}^{-1} list, ordered by face index
    void get_inv_hat(std::vector<MatrixXd>& inv_hat_list) const;

    void get_triangle_adj(std::vector<std::vector<size_t>>& adj_list) const;

    static void split_words(
        const std::string& line,
        std::vector<std::string>& words,
        const char delim
    );

    MatrixXd& vertices() { return m_vertices; }
    MatrixXi& faces() { return m_faces; }

private:
    bool __load_obj(const std::string& filename);

private:
    MatrixXd m_vertices;
    MatrixXi m_faces;    
};

}

#endif
