#pragma once

namespace config {

constexpr double gamma_value = 1.4;

constexpr double x_min = 0.0;
constexpr double x_max = 1.0;
constexpr double y_min = 0.0;
constexpr double y_max = 1.0;

constexpr int nx = 160;
constexpr int ny = 160;

constexpr double cfl = 0.35;
constexpr int n_steps = 300;
constexpr int report_interval = 20;
constexpr int write_interval = 100;

constexpr double x_split = 0.5;
constexpr double y_split = 0.5;

// Lax-Liu style four-state 2D Riemann problem.
constexpr double rho_ul = 0.5323;
constexpr double u_ul = 1.206;
constexpr double v_ul = 0.0;
constexpr double p_ul = 0.3;

constexpr double rho_ur = 1.5;
constexpr double u_ur = 0.0;
constexpr double v_ur = 0.0;
constexpr double p_ur = 1.5;

constexpr double rho_ll = 0.138;
constexpr double u_ll = 1.206;
constexpr double v_ll = 1.206;
constexpr double p_ll = 0.029;

constexpr double rho_lr = 0.5323;
constexpr double u_lr = 0.0;
constexpr double v_lr = 1.206;
constexpr double p_lr = 0.3;

}  // namespace config
