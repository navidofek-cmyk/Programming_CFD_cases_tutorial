# 08_2d_euler_quadrant_riemann

Educational first 2D compressible Euler case.

## Purpose

This case solves a classic 2D quadrant Riemann problem:

- four constant primitive states
- discontinuities aligned with `x = x_split` and `y = y_split`
- resulting 2D interacting compressible wave pattern

It is a good first 2D compressible example because:

- the geometry is simple
- there are no curved walls or masks
- it still produces genuinely 2D wave interactions

## Numerics

- 2D compressible Euler equations
- cell-centered finite volume method
- uniform Cartesian mesh
- explicit Euler time stepping
- local Rusanov flux in `x` and `y`
- transmissive boundaries

## Output

The solver writes VTK snapshots into:

- `output/euler2d_00100.vtk`
- `output/euler2d_00200.vtk`
- ...

Fields include:

- density
- pressure
- Mach number
- velocity vector

## Build

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/08_2d_euler_quadrant_riemann
make
```

## Run

```bash
./euler2d
```

or:

```bash
make run
```

## Notes

- This is an educational first 2D compressible Euler solver, not a high-resolution shock-capturing code.
- Rusanov flux is robust and simple, but diffusive.
- It is a good stepping stone before moving to:
  - HLL / HLLC
  - MUSCL reconstruction
  - wedge or oblique-shock geometries
