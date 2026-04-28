*[English](#english) · [Čeština](#čeština)*

---

## English

# CFD Solver Implementations — Incompressible & Compressible Flow

Twelve standalone C++ solvers covering the core numerical methods used in industrial CFD codes. Each case targets a specific physical regime and numerical scheme, written for clarity and correctness rather than production-scale performance.

No external libraries. Standard C++17 only. Each case compiles and runs independently.

### Why These Methods Matter in Practice

The numerical kernels implemented here are the same ones that underpin commercial solvers — OpenFOAM, ANSYS Fluent, Star-CCM+:

- **SIMPLE / PISO** — pressure–velocity coupling in every steady and transient incompressible solver used in process equipment, heat exchangers, and duct systems
- **Riemann-based flux schemes (Rusanov, HLLC)** — foundation of compressible solvers for gas turbines, relief valve sizing, safety blow-down, and supersonic nozzle design
- **MUSCL + slope limiters** — the standard route to second-order accuracy in shock-capturing without spurious oscillations
- **RANS closure** — industrial standard for turbulent flows in pipelines, separators, and heat transfer equipment

Knowing how these work at the implementation level improves solver setup, boundary condition choices, convergence troubleshooting, and mesh sensitivity interpretation in production simulations.

### Cases

#### Incompressible Flow

| Case | Physical Setup | Numerical Method |
|---|---|---|
| `01_cavity_case` | Lid-driven cavity | Finite differences, projection/pressure-correction |
| `02_laminar_channel_obstacle` | Steady channel with bluff-body obstacle | Cell-centred FVM, SIMPLE pseudo-time |
| `03_channel_obstacle_unsteady` | Transient wake behind obstacle | Cell-centred FVM, PISO corrections |
| `04_laminar_channel_obstacle_rebuild` | Steady obstacle flow (rebuilt) | Staggered MAC layout, SIMPLE/SOR |
| `05_laminar_channel_simple` | Poiseuille channel reference | Staggered channel solver |

#### Compressible Flow

| Case | Physical Setup | Numerical Method |
|---|---|---|
| `06_exact_riemann_solver_1d_euler` | 1D shock-tube / Riemann problem | Exact star-region solve |
| `07_laval_nozzle_quasi_1d` | Quasi-1D Laval nozzle, choked flow | MacCormack predictor-corrector, area source |
| `08_2d_euler_quadrant_riemann` | 2D quadrant Riemann initial data | Cell-centred FV, Rusanov flux |
| `09_2d_euler_oblique_shock` | Supersonic ramp, oblique shock | Cell-centred FV, Rusanov flux, CFL stepping |
| `10_hllc_muscl_upgrade` | Oblique shock — improved resolution | HLLC Riemann solver, MUSCL, minmod limiter |
| `11_rae2822_trans` | RAE2822 airfoil — transonic regime | 2D FV, transonic flow conditions |
| `12_rae2822_rans` | RAE2822 airfoil — turbulent flow | RANS closure |

### Numerical Methods Reference

| Category | Methods implemented |
|---|---|
| Pressure–velocity coupling | Explicit projection, SIMPLE, PISO |
| Flux evaluation | Upwind, Rusanov (LLF), HLLC |
| Reconstruction | First-order, MUSCL with minmod limiter |
| Grid types | Structured Cartesian, staggered MAC |
| Time integration | Explicit Euler, MacCormack predictor-corrector |
| Obstacle handling | Cartesian mask (immersed solid) |

### Output and Post-processing

All cases export results to CSV and/or VTK format. VTK files open directly in ParaView.

```bash
cd 09_2d_euler_oblique_shock
make
./oblique_shock
paraview output_*.vtk
```

### Build

```bash
cd <case_directory>
make
./solver_name
```

Or directly:

```bash
g++ -std=c++17 -O2 main.cpp -o solver
./solver
```

**Requirements:** GCC 9+ or Clang 10+, Make, ParaView (optional), Python 3.8+ with matplotlib (optional)

---

## Čeština

# Implementace CFD solverů — nestlačitelné a stlačitelné proudění

Dvanáct samostatných C++ solverů pokrývajících základní numerické metody používané v průmyslových CFD kódech. Každý případ cílí na konkrétní fyzikální režim a numerické schéma — psáno pro čitelnost a správnost, ne pro produkční výkon.

Žádné externí knihovny. Pouze standardní C++17. Každý případ se kompiluje a spouští samostatně.

### Proč jsou tyto metody důležité v praxi

Numerické jádra implementovaná zde jsou totožná s těmi, na kterých stojí komerční solvery — OpenFOAM, ANSYS Fluent, Star-CCM+:

- **SIMPLE / PISO** — vazba tlak–rychlost v každém ustáleném i přechodovém nestlačitelném solveru pro procesní zařízení, výměníky tepla a potrubní systémy
- **Riemannova schémata (Rusanov, HLLC)** — základ stlačitelných solverů pro plynové turbíny, dimenzování pojistných ventilů, odtlakování a nadzvukový design trysek
- **MUSCL + limitery sklonu** — standardní cesta k přesnosti druhého řádu při zachycení rázových vln bez parazitních oscilací
- **RANS uzávěr** — průmyslový standard pro turbulentní proudění v potrubí, separátorech a teplosměnných zařízeních

### Případy

#### Nestlačitelné proudění

| Případ | Fyzikální nastavení | Numerická metoda |
|---|---|---|
| `01_cavity_case` | Lid-driven cavity | Konečné diference, projekční metoda |
| `02_laminar_channel_obstacle` | Ustálené proudění kanálem s překážkou | MKO cell-centred, SIMPLE pseudo-čas |
| `03_channel_obstacle_unsteady` | Přechodová stopa za překážkou | MKO cell-centred, PISO korekce |
| `04_laminar_channel_obstacle_rebuild` | Ustálené proudění s překážkou (přepsáno) | Staggered MAC, SIMPLE/SOR |
| `05_laminar_channel_simple` | Referenční Poiseuilleovo proudění | Staggered channel solver |

#### Stlačitelné proudění

| Případ | Fyzikální nastavení | Numerická metoda |
|---|---|---|
| `06_exact_riemann_solver_1d_euler` | 1D trubice s rázovou vlnou | Přesné řešení hvězdové oblasti |
| `07_laval_nozzle_quasi_1d` | Kvazi-1D Lavalova tryska, ucpané proudění | MacCormack prediktor-korektor, zdrojový člen plochy |
| `08_2d_euler_quadrant_riemann` | 2D kvadrantová Riemannova počáteční data | MKO cell-centred, Rusanovův tok |
| `09_2d_euler_oblique_shock` | Nadzvuková rampa, šikmá rázová vlna | MKO cell-centred, Rusanov, CFL krokování |
| `10_hllc_muscl_upgrade` | Šikmá rázová vlna — vyšší rozlišení | HLLC Riemannův solver, MUSCL, minmod limiter |
| `11_rae2822_trans` | Profil RAE2822 — transonický režim | 2D MKO, transonické podmínky |
| `12_rae2822_rans` | Profil RAE2822 — turbulentní proudění | RANS uzávěr |

### Přehled numerických metod

| Kategorie | Implementované metody |
|---|---|
| Vazba tlak–rychlost | Explicitní projekce, SIMPLE, PISO |
| Vyhodnocení toků | Upwind, Rusanov (LLF), HLLC |
| Rekonstrukce | První řád, MUSCL s minmod limiterem |
| Typy sítí | Strukturovaná kartézská, staggered MAC |
| Časová integrace | Explicitní Euler, MacCormack prediktor-korektor |
| Zpracování překážek | Kartézská maska (immersed solid) |

### Výstup a post-processing

Všechny případy exportují výsledky do CSV a/nebo VTK formátu. VTK soubory se otevírají přímo v ParaView.

### Kompilace

```bash
cd <adresář_případu>
make
./solver
```

**Požadavky:** GCC 9+ nebo Clang 10+, Make, ParaView (volitelné), Python 3.8+ s matplotlib (volitelné)
