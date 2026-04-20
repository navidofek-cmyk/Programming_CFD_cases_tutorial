#pragma once

#include <vector>

struct ScalarField1D {
    std::vector<double> values;

    explicit ScalarField1D(int n = 0, double value = 0.0)
        : values(static_cast<std::size_t>(n), value) {}

    [[nodiscard]] double& operator()(int i) { return values[static_cast<std::size_t>(i)]; }
    [[nodiscard]] const double& operator()(int i) const { return values[static_cast<std::size_t>(i)]; }
    [[nodiscard]] int size() const { return static_cast<int>(values.size()); }
};
