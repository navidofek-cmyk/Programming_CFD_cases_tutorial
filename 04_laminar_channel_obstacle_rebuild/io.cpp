#include "io.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>

// Write a cell-centred scalar field to CSV.
// Rows are j (top to bottom in file → y decreasing), columns are i (x increasing).
void writeCSV(const std::string& filename,
              const Mesh& mesh,
              const Field2D& field)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot open " + filename);
    f << std::scientific << std::setprecision(6);

    f << "y\\x";
    for (int i = 0; i < mesh.nx; ++i) f << ',' << mesh.xc[i];
    f << '\n';

    for (int j = mesh.ny-1; j >= 0; --j) {
        f << mesh.yc[j];
        for (int i = 0; i < mesh.nx; ++i) f << ',' << field(i,j);
        f << '\n';
    }
}

// Write u, v averaged to cell centres from staggered faces.
// u_face: (nx+1)×ny,  v_face: nx×(ny+1)
void writeUV(const std::string& u_file, const std::string& v_file,
             const Mesh& mesh,
             const Field2D& uf, const Field2D& vf)
{
    const int nx = mesh.nx, ny = mesh.ny;

    // Build cell-centred u and v
    Field2D uc(nx, ny), vc(nx, ny);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            uc(i,j) = 0.5*(uf(i,j) + uf(i+1,j));  // average west+east u-faces
            vc(i,j) = 0.5*(vf(i,j) + vf(i,j+1));  // average south+north v-faces
        }

    writeCSV(u_file, mesh, uc);
    writeCSV(v_file, mesh, vc);
}

void writeMask(const std::string& filename, const Mesh& mesh)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot open " + filename);
    f << "y\\x";
    for (int i = 0; i < mesh.nx; ++i) f << ',' << mesh.xc[i];
    f << '\n';
    for (int j = mesh.ny-1; j >= 0; --j) {
        f << mesh.yc[j];
        for (int i = 0; i < mesh.nx; ++i)
            f << ',' << (mesh.isSolid(i,j) ? 1 : 0);
        f << '\n';
    }
}
