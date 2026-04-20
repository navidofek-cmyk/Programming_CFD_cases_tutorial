#pragma once

#include <vector>

struct Field2D {
    int nx = 0;
    int ny = 0;
    std::vector<double> values;

    Field2D() = default;

    Field2D(int nx_in, int ny_in, double value = 0.0)
        : nx(nx_in), ny(ny_in), values(static_cast<std::size_t>(nx_in) * static_cast<std::size_t>(ny_in), value) {}

    [[nodiscard]] double& operator()(int i, int j) {
        return values[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i)];
    }

    [[nodiscard]] const double& operator()(int i, int j) const {
        return values[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i)];
    }
};
