# 09_2d_euler_oblique_shock

Educational 2D compressible Euler case with a supersonic inlet and a wedge / ramp.

## Purpose

This case models:

- supersonic inflow from the left boundary
- a wedge / ramp represented as a stair-step solid mask on the Cartesian grid
- oblique shock formation above the ramp

## Numerics

- 2D compressible Euler equations
- cell-centered finite volume method
- uniform Cartesian mesh
- explicit Euler time stepping
- local Rusanov flux in `x` and `y`
- supersonic inlet
- transmissive top / outlet
- reflective wall treatment at the ramp

## Output

The solver writes VTK snapshots into:

- `output/oblique_shock_00125.vtk`
- `output/oblique_shock_00250.vtk`
- ...

Fields include:

- density
- pressure
- Mach number
- velocity vector
- solid mask

## Build

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/09_2d_euler_oblique_shock
make
```

## Run

```bash
./oblique_shock
```

or:

```bash
make run
```

## Notes

- This is an educational oblique-shock case, not a high-resolution production solver.
- The wedge is represented by a stair-step solid mask rather than a body-fitted mesh.
- Rusanov flux is robust and simple, but diffusive.
