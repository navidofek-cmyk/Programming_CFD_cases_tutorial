#include "geometry.hpp"

ScalarField1D build_area_field(const Mesh1D& mesh, double throat_position) {
    ScalarField1D area(mesh.nx, 1.0);

    for (int i = 0; i < mesh.nx; ++i) {
        const double x = mesh.x[static_cast<std::size_t>(i)];
        area(i) = 1.0 + 2.2 * (x - throat_position) * (x - throat_position);
    }

    return area;
}
