#pragma once

#include <cmath>

namespace config {

constexpr double gamma_value = 1.4;

constexpr double x_min = 0.0;
constexpr double x_max = 1.2;
constexpr double y_min = 0.0;
constexpr double y_max = 0.8;

constexpr int nx = 240;
constexpr int ny = 160;

constexpr double cfl = 0.35;
constexpr int n_steps = 500;
constexpr int report_interval = 20;
constexpr int write_interval = 125;

constexpr double ramp_start_x = 0.20;
constexpr double ramp_angle_deg = 15.0;

constexpr double rho_inlet = 1.0;
constexpr double p_inlet = 1.0;
constexpr double mach_inlet = 2.0;
constexpr double u_inlet = mach_inlet * std::sqrt(gamma_value * p_inlet / rho_inlet);
constexpr double v_inlet = 0.0;

}  // namespace config
