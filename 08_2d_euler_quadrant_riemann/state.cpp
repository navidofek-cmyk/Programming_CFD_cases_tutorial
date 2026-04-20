#include "state.hpp"

#include <algorithm>
#include <cmath>

PrimitiveState conservative_to_primitive(double U1, double U2, double U3, double U4, double gamma_value) {
    const double rho = std::max(1.0e-8, U1);
    const double u = U2 / rho;
    const double v = U3 / rho;
    const double kinetic = 0.5 * rho * (u * u + v * v);
    const double p = std::max(1.0e-8, (gamma_value - 1.0) * (U4 - kinetic));
    return {rho, u, v, p};
}

void primitive_to_conservative(const PrimitiveState& state,
                               double gamma_value,
                               double& U1,
                               double& U2,
                               double& U3,
                               double& U4) {
    U1 = state.rho;
    U2 = state.rho * state.u;
    U3 = state.rho * state.v;
    U4 = state.p / (gamma_value - 1.0)
         + 0.5 * state.rho * (state.u * state.u + state.v * state.v);
}

double sound_speed(const PrimitiveState& state, double gamma_value) {
    return std::sqrt(gamma_value * state.p / state.rho);
}
