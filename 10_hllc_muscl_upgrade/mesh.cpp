#include "mesh.hpp"

#include <cmath>

Mesh2D build_mesh(int nx,
                  int ny,
                  double x_min,
                  double x_max,
                  double y_min,
                  double y_max,
                  double ramp_start_x,
                  double ramp_angle_deg) {
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
    mesh.solid_mask.resize(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0);

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

    const double ramp_angle_rad = ramp_angle_deg * 3.14159265358979323846 / 180.0;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const double x = mesh.xc[static_cast<std::size_t>(i)];
            const double y = mesh.yc[static_cast<std::size_t>(j)];
            const double wall_height =
                (x >= ramp_start_x) ? std::tan(ramp_angle_rad) * (x - ramp_start_x) : 0.0;
            mesh.solid_mask[static_cast<std::size_t>(mesh.idx(i, j))] = (y < wall_height) ? 1 : 0;
        }
    }

    return mesh;
}
