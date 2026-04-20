#include "mesh.hpp"

Mesh1D build_mesh(int nx, double x_min, double x_max) {
    Mesh1D mesh;
    mesh.nx = nx;
    mesh.x_min = x_min;
    mesh.x_max = x_max;
    mesh.dx = (x_max - x_min) / static_cast<double>(nx - 1);
    mesh.x.resize(static_cast<std::size_t>(nx), 0.0);

    for (int i = 0; i < nx; ++i) {
        mesh.x[static_cast<std::size_t>(i)] = x_min + static_cast<double>(i) * mesh.dx;
    }

    return mesh;
}
