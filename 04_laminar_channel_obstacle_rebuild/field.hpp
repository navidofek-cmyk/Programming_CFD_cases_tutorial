#pragma once
#include <vector>
#include <stdexcept>

// 2D cell-centered scalar field stored in row-major order (i = x-index, j = y-index).
// Access: f(i, j) where i in [0, nx), j in [0, ny).
class Field2D {
public:
    int nx, ny;

    Field2D() : nx(0), ny(0) {}
    Field2D(int nx_, int ny_, double init = 0.0)
        : nx(nx_), ny(ny_), data_(static_cast<std::size_t>(nx_ * ny_), init) {}

    double& operator()(int i, int j)       { return data_[idx(i,j)]; }
    double  operator()(int i, int j) const { return data_[idx(i,j)]; }

    void fill(double v) { std::fill(data_.begin(), data_.end(), v); }

    const std::vector<double>& raw() const { return data_; }
          std::vector<double>& raw()       { return data_; }

private:
    std::vector<double> data_;
    int idx(int i, int j) const { return j * nx + i; }
};
