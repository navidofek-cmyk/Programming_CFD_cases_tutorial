# CFD Solver Implementations — Incompressible & Compressible Flow

Twelve standalone C++ solvers covering the core numerical methods used in industrial CFD codes. Each case targets a specific physical regime and numerical scheme, written for clarity and correctness rather than production-scale performance.

No external libraries. Standard C++17 only. Each case compiles and runs independently.

---

## Why These Methods Matter in Practice

The numerical kernels implemented here are the same ones that underpin commercial solvers — OpenFOAM, ANSYS Fluent, Star-CCM+:

- **SIMPLE / PISO** — pressure–velocity coupling in every steady and transient incompressible solver used in process equipment, heat exchangers, and duct systems
- **Riemann-based flux schemes (Rusanov, HLLC)** — foundation of compressible solvers for gas turbines, relief valve sizing, safety blow-down, and supersonic nozzle design
- **MUSCL + slope limiters** — the standard route to second-order accuracy in shock-capturing without spurious oscillations; present in any production compressible solver
- **RANS closure** — industrial standard for turbulent flows in pipelines, separators, and heat transfer equipment

Knowing how these work at the implementation level improves solver setup, boundary condition choices, convergence troubleshooting, and mesh sensitivity interpretation in production work.

---

## Cases

### Incompressible Flow

| Case | Physical Setup | Numerical Method |
|---|---|---|
| `01_cavity_case` | Lid-driven cavity | Finite differences, projection/pressure-correction |
| `02_laminar_channel_obstacle` | Steady channel with bluff-body obstacle | Cell-centred FVM, SIMPLE pseudo-time |
| `03_channel_obstacle_unsteady` | Transient wake behind obstacle | Cell-centred FVM, PISO corrections |
| `04_laminar_channel_obstacle_rebuild` | Steady obstacle flow (rebuilt) | Staggered MAC layout, SIMPLE/SOR |
| `05_laminar_channel_simple` | Poiseuille channel reference | Staggered channel solver |

**Physical observables:** pressure–velocity coupling behaviour, bluff-body wake development, convergence rates of segregated solvers.

### Compressible Flow

| Case | Physical Setup | Numerical Method |
|---|---|---|
| `06_exact_riemann_solver_1d_euler` | 1D shock-tube / Riemann problem | Exact star-region solve |
| `07_laval_nozzle_quasi_1d` | Quasi-1D Laval nozzle, choked flow | MacCormack predictor-corrector, area source |
| `08_2d_euler_quadrant_riemann` | 2D quadrant Riemann initial data | Cell-centred FV, Rusanov flux |
| `09_2d_euler_oblique_shock` | Supersonic ramp, oblique shock | Cell-centred FV, Rusanov flux, CFL stepping |
| `10_hllc_muscl_upgrade` | Oblique shock — improved resolution | HLLC Riemann solver, MUSCL, minmod limiter |
| `11_rae2822_trans` | RAE2822 airfoil — transonic regime | 2D FV, transonic flow conditions |
| `12_rae2822_rans` | RAE2822 airfoil — turbulent flow | RANS closure |

**Physical observables:** shock standoff and angle, Mach number distribution, pressure coefficient along airfoil, RANS vs. inviscid comparison.

---

## Numerical Methods Reference

| Category | Methods implemented |
|---|---|
| Pressure–velocity coupling | Explicit projection, SIMPLE, PISO |
| Flux evaluation | Upwind, Rusanov (LLF), HLLC |
| Reconstruction | First-order, MUSCL with minmod limiter |
| Grid types | Structured Cartesian, staggered MAC |
| Time integration | Explicit Euler, MacCormack predictor-corrector |
| Obstacle handling | Cartesian mask (immersed solid) |

---

## Output and Post-processing

All cases export results to CSV and/or VTK format. VTK files open directly in ParaView for field visualisation (velocity components, pressure, density, Mach number, temperature — depending on case).

```bash
# Example: steady obstacle case
cd 02_laminar_channel_obstacle
make
./channel_obstacle
# outputs: fields.csv, pressure.csv

# Example: 2D oblique shock
cd 09_2d_euler_oblique_shock
make
./oblique_shock
paraview output_*.vtk
```

Python scripts for CSV plotting (matplotlib) are included where applicable.

---

## Build

C++17, no external dependencies.

```bash
cd <case_directory>
make
./solver_name
```

Or compile directly:

```bash
g++ -std=c++17 -O2 main.cpp -o solver
./solver
```

Each case directory contains a short note on physical setup, expected flow behaviour, and result interpretation.

**Requirements:**
- GCC 9+ or Clang 10+
- Make
- ParaView (optional, for VTK output)
- Python 3.8+ with matplotlib (optional, for CSV plots)
