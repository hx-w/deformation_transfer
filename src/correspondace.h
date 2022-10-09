// method of correspondence

#include "meshlib/mesh.h"


/**
 * compute triangle mappings between two meshes with markders
 * @param src_mesh source mesh
 * @param tgt_mesh target mesh
 * @param markers shape (k, 2), k is the number of markers
 * @return mappings shape (m, 2), m is the number of triangles in target mesh
 */
void compute_correspondence(
    MeshLib::Mesh& src_mesh,
    MeshLib::Mesh& tgt_mesh,
    const MeshLib::MatrixXi& markers,
    MeshLib::MatrixXi& mappings
);

void construct_indentity_cost(
    // MeshLib::Mesh& mesh,
    const std::vector<MeshLib::MatrixXd>& inv_hat_V,
    MeshLib::MatrixXd& AEi,
    MeshLib::MatrixXd& Bi
);

void construct_smoothness_cost(
    const std::vector<MeshLib::MatrixXd>& inv_hat_V,
    const std::vector<std::vector<int>>& face_adj_list,
    MeshLib::MatrixXd& AEs,
    MeshLib::MatrixXd& Bs
);

void read_markers(
    const std::string& filename,
    MeshLib::MatrixXi& markers
);
