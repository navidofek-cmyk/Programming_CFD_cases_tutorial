#include "io.hpp"

#include <fstream>
#include <iomanip>

void write_vtk(const std::string& filename,
               const Mesh2D& mesh,
               const PrimitiveFields& primitive) {
    std::ofstream out(filename);
    out << std::setprecision(16);

    out << "# vtk DataFile Version 3.0\n";
    out << "2D Euler quadrant Riemann problem\n";
    out << "ASCII\n";
    out << "DATASET RECTILINEAR_GRID\n";
    out << "DIMENSIONS " << mesh.nx + 1 << ' ' << mesh.ny + 1 << " 1\n";

    out << "X_COORDINATES " << mesh.nx + 1 << " double\n";
    for (double x : mesh.xf) {
        out << x << '\n';
    }

    out << "Y_COORDINATES " << mesh.ny + 1 << " double\n";
    for (double y : mesh.yf) {
        out << y << '\n';
    }

    out << "Z_COORDINATES 1 double\n0\n";

    out << "CELL_DATA " << mesh.nx * mesh.ny << '\n';

    out << "SCALARS density double 1\nLOOKUP_TABLE default\n";
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            out << primitive.rho(i, j) << '\n';
        }
    }

    out << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            out << primitive.p(i, j) << '\n';
        }
    }

    out << "SCALARS mach double 1\nLOOKUP_TABLE default\n";
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            out << primitive.M(i, j) << '\n';
        }
    }

    out << "VECTORS velocity double\n";
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            out << primitive.u(i, j) << ' ' << primitive.v(i, j) << " 0\n";
        }
    }
}
