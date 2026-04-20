# 05 – Laminar Channel Flow (SIMPLE, Staggered Grid)

Steady incompressible laminar flow through a 2D channel solved with the
SIMPLE algorithm on a staggered (MAC) grid.  At full development the
velocity profile is parabolic (Hagen–Poiseuille); the analytic maximum
is `u_max = 1.5 · U_in`.

![CFD domain and staggered cell](figures/domain.png)

---

## Numerical method

| Item | Choice |
|---|---|
| Grid | Staggered (MAC): p at cell centres, u/v at face centres |
| Algorithm | SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) |
| Convection | First-order upwind (UDS) |
| Diffusion | Central differences |
| Pressure Poisson | Conjugate Gradient (exact for SPD system) |
| Under-relaxation | α_u = α_v = 0.7, α_p = 0.3 |

The staggered grid eliminates the checkerboard pressure instability
without Rhie–Chow interpolation.  The pressure-correction equation is
derived directly from the discrete continuity requirement and is SPD,
so Conjugate Gradient converges in O(√κ) iterations.

### Boundary conditions

| Boundary | u | v | p |
|---|---|---|---|
| Inlet (x = 0) | U_in (uniform) | 0 at x = 0 (ghost cell) | zero-gradient |
| Outlet (x = Lx) | zero-gradient | zero-gradient | p = 0 (Dirichlet) |
| Walls (y = 0, Ly) | no-slip (ghost) | 0 | zero-gradient |

---

## Parameters (`config.hpp`)

```
Domain        Lx = 6,  Ly = 1
Grid          Nx = 120,  Ny = 40
Re            100   (nu = U_in · Ly / Re = 0.01)
Inlet         U_in = 1.0
pseudo_dt     0.05
alpha_u/v     0.7
alpha_p       0.3
max_outer     4000
tol_u         1e-7  (RMS change in u between outer iterations)
```

---

## Build & run

```bash
make
./channel_flow
```

Outputs:
- `result.csv` — cell-centred x, y, u, v, p
- `result.vtk` — VTK legacy ASCII for ParaView

---

## Results (Re = 100, 120×40 grid)

```
Converged at iter ~784
u_max(outlet) = 1.491   analytic = 1.5   (error < 1 %)
div_interior  ~ 1e-8    (machine precision)
```

---

## Visualisation in ParaView

1. `File → Open → result.vtk → Apply`
2. Colour by **p** to see the linear pressure drop along the channel
3. Colour by **speed** to see the Poiseuille parabolic profile
4. `Filters → Stream Tracer` on field **velocity** for streamlines
