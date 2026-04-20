#pragma once

// ─── Physical parameters ───────────────────────────────────────────────────
inline constexpr double H      = 1.0;           // obstacle height [m]
inline constexpr double U_IN   = 1.0;           // inlet x-velocity [m/s]
inline constexpr double RE     = 60.0;          // Reynolds number
inline constexpr double NU     = U_IN * H / RE; // kinematic viscosity [m²/s]
inline constexpr double RHO    = 1.0;           // density (incompressible, constant)

// ─── Domain ────────────────────────────────────────────────────────────────
inline constexpr double LX = 32.0 * H;
inline constexpr double LY =  8.0 * H;

// Obstacle extents (cell masking — stair-step on Cartesian mesh)
inline constexpr double OBS_X0 =  8.0 * H;
inline constexpr double OBS_X1 =  9.0 * H;
inline constexpr double OBS_Y0 =  3.5 * H;
inline constexpr double OBS_Y1 =  4.5 * H;

// ─── Mesh ──────────────────────────────────────────────────────────────────
inline constexpr int NX = 240;
inline constexpr int NY = 120;

// ─── Solver ────────────────────────────────────────────────────────────────
inline constexpr int    MAX_ITER = 5000;
inline constexpr double TOL      = 1e-4;   // residual convergence threshold

// Under-relaxation factors (SIMPLE stability)
inline constexpr double ALPHA_U = 0.5;
inline constexpr double ALPHA_V = 0.5;
inline constexpr double ALPHA_P = 0.3;

// ─── Output ────────────────────────────────────────────────────────────────
inline constexpr int PRINT_EVERY = 100;
