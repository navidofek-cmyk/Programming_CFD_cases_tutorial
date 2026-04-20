#pragma once

namespace config {

constexpr double gamma_value = 1.4;
constexpr double gas_constant = 1.0;

constexpr double x_min = 0.0;
constexpr double x_max = 3.0;
constexpr int nx = 121;

constexpr int n_steps = 2500;
constexpr double cfl = 0.15;
constexpr int report_interval = 100;
constexpr double artificial_viscosity = 0.02;

// Simple nozzle shape: throat at x = 1.5
// A(x) = 1 + 2.2 * (x - 1.5)^2
constexpr double throat_position = 1.5;

// Default initial state: mild subsonic field that relaxes toward nozzle flow.
constexpr double rho_initial = 1.0;
constexpr double u_initial = 0.1;
constexpr double p_initial = 1.0;

// Boundary conditions:
// inlet: prescribed stagnation state, with velocity extrapolated from interior
// outlet: prescribed static pressure if subsonic, extrapolated if supersonic
constexpr double T0_inlet = 1.0;
constexpr double p0_inlet = 1.0;

constexpr double p_exit = 0.6784;

}  // namespace config
