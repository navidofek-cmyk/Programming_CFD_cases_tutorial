#include "mesh.hpp"

Mesh2D build_mesh(int nx, int ny, double x_min, double x_max, double y_min, double y_max) {
    Mesh2D mesh;
    mesh.nx = nx;
    mesh.ny = ny;
    mesh.x_min = x_min;
    mesh.x_max = x_max;
    mesh.y_min = y_min;
    mesh.y_max = y_max;
    mesh.dx = (x_max - x_min) / static_cast<double>(nx);
    mesh.dy = (y_max - y_min) / static_cast<double>(ny);

    mesh.xc.resize(static_cast<std::size_t>(nx), 0.0);
    mesh.yc.resize(static_cast<std::size_t>(ny), 0.0);
    mesh.xf.resize(static_cast<std::size_t>(nx + 1), 0.0);
    mesh.yf.resize(static_cast<std::size_t>(ny + 1), 0.0);

    for (int i = 0; i <= nx; ++i) {
        mesh.xf[static_cast<std::size_t>(i)] = x_min + static_cast<double>(i) * mesh.dx;
    }
    for (int j = 0; j <= ny; ++j) {
        mesh.yf[static_cast<std::size_t>(j)] = y_min + static_cast<double>(j) * mesh.dy;
    }
    for (int i = 0; i < nx; ++i) {
        mesh.xc[static_cast<std::size_t>(i)] = x_min + (static_cast<double>(i) + 0.5) * mesh.dx;
    }
    for (int j = 0; j < ny; ++j) {
        mesh.yc[static_cast<std::size_t>(j)] = y_min + (static_cast<double>(j) + 0.5) * mesh.dy;
    }

    return mesh;
}
