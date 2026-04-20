#pragma once

#include <vector>

struct Mesh1D {
    int nx = 0;
    double x_min = 0.0;
    double x_max = 0.0;
    double dx = 0.0;
    std::vector<double> x;
};

Mesh1D build_mesh(int nx, double x_min, double x_max);
