# 10_hllc_muscl_upgrade

Numerical upgrade of the oblique-shock ramp case using HLLC and MUSCL reconstruction.

## Purpose

This case keeps the same physical setup as `09_2d_euler_oblique_shock`, but improves the numerics.

## Numerics

- 2D compressible Euler equations
- cell-centered finite volume method
- uniform Cartesian mesh
- explicit Euler time stepping
- HLLC approximate Riemann solver
- MUSCL reconstruction
- minmod slope limiter
- supersonic inlet
- transmissive top / outlet
- reflective wall treatment at the ramp

## Output

The solver writes VTK snapshots into:

- `output/oblique_shock_hllc_00125.vtk`
- `output/oblique_shock_hllc_00250.vtk`
- ...

Fields include:

- density
- pressure
- Mach number
- velocity vector
- solid mask

## Build

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/10_hllc_muscl_upgrade
make
```

## Run

```bash
./oblique_shock_hllc
```

or:

```bash
make run
```

## Notes

- This is still an educational compressible solver, not a production shock-capturing code.
- The wedge is represented by a stair-step solid mask rather than a body-fitted mesh.
- HLLC + MUSCL should preserve contact and shock structure better than the first-order Rusanov version in case `09`.
