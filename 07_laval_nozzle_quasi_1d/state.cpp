#include "state.hpp"

#include <cmath>

double pressure_from_state(double rho, double u, double energy, double gamma_value) {
    return (gamma_value - 1.0) * (energy - 0.5 * rho * u * u);
}

double sound_speed(double p, double rho, double gamma_value) {
    return std::sqrt(gamma_value * p / rho);
}
