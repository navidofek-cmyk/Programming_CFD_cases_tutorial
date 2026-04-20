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

    [[nodiscard]] int idx(int i, int j) const { return j * nx + i; }
};

Mesh2D build_mesh(int nx, int ny, double x_min, double x_max, double y_min, double y_max);
