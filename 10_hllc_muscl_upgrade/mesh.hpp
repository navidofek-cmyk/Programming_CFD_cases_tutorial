#pragma once

#include <vector>

struct Mesh2D {
    int nx = 0;
    int ny = 0;
    double x_min = 0.0;
    double x_max = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    std::vector<double> xc;
    std::vector<double> yc;
    std::vector<double> xf;
    std::vector<double> yf;
    std::vector<int> solid_mask;

    [[nodiscard]] int idx(int i, int j) const { return j * nx + i; }
    [[nodiscard]] bool is_inside(int i, int j) const { return i >= 0 && i < nx && j >= 0 && j < ny; }
    [[nodiscard]] bool is_solid(int i, int j) const { return solid_mask[static_cast<std::size_t>(idx(i, j))] != 0; }
    [[nodiscard]] bool is_fluid(int i, int j) const { return !is_solid(i, j); }
};

Mesh2D build_mesh(int nx,
                  int ny,
                  double x_min,
                  double x_max,
                  double y_min,
                  double y_max,
                  double ramp_start_x,
                  double ramp_angle_deg);
