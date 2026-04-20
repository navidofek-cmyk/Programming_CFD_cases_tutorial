#include "solver.hpp"
#include "boundary.hpp"
#include "config.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// Staggered MAC grid SIMPLE algorithm (educational)
//
// Grid layout (mesh.hpp):
//   u(i,j): x-velocity at east face of cell (i-1,j) / west face of cell (i,j)
//           → x = i*dx,       y = (j+0.5)*dy   size (nx+1)×ny
//   v(i,j): y-velocity at north face of cell (i,j-1) / south face of cell (i,j)
//           → x = (i+0.5)*dx, y = j*dy          size nx×(ny+1)
//   p(i,j): pressure at cell centre              size nx×ny
//
// SIMPLE outer loop:
//   1. N_MOM_INNER GS sweeps for x-momentum → u*, store aP_u per face
//   2. N_MOM_INNER GS sweeps for y-momentum → v*, store aP_v per face
//   3. Compute b(i,j) = (u*(i+1,j)-u*(i,j))*dy + (v*(i,j+1)-v*(i,j))*dx
//   4. Compute aE,aW,aN,aS for pressure-correction Poisson
//   5. N_PP_SWEEPS alternating TDMA sweeps → pp
//   6. p  += ALPHA_P * pp
//   7. u* -= (pp_right-pp_left)/dx * dy/aP_u  (full SIMPLE correction)
//      v* -= (pp_top-pp_bot)/dy   * dx/aP_v
//   8. Apply BCs; check convergence (max |u-u_old|).
// ─────────────────────────────────────────────────────────────────────────────

static constexpr int N_PP_SWEEPS  = 100;
static constexpr int N_MOM_INNER  = 5;    // inner GS sweeps per outer iteration
// Optimal SOR omega for 240×120 Poisson:  ω ≈ 2/(1+sin(π/max(Nx,Ny))) ≈ 1.96
// 100 SOR sweeps → ~98% Poisson residual reduction vs ~8% with unaccelerated ADI.
static constexpr double OMEGA_SOR = 1.93; // slightly below optimal for robustness


int runSolver(const Mesh& mesh, Field2D& u, Field2D& v, Field2D& p)
{
    const int    nx = mesh.nx;
    const int    ny = mesh.ny;
    const double dx = mesh.dx;
    const double dy = mesh.dy;

    const double aP_diff = 2.0*(NU*dy/dx + NU*dx/dy);
    Field2D aPu(nx+1, ny,   aP_diff);
    Field2D aPv(nx,   ny+1, aP_diff);

    Field2D pp(nx, ny, 0.0);
    Field2D bfield(nx, ny, 0.0);
    Field2D ppAE(nx, ny, 0.0), ppAW(nx, ny, 0.0);
    Field2D ppAN(nx, ny, 0.0), ppAS(nx, ny, 0.0);

    // Save fields from previous outer iteration for convergence check & under-relaxation
    Field2D u_old(nx+1, ny,   0.0);
    Field2D v_old(nx,   ny+1, 0.0);

    // ── Initial conditions ────────────────────────────────────────────────────
    // Velocity: plug flow (avoids zero-velocity transient)
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i)
            u(i, j) = mesh.uBlocked(i, j) ? 0.0 : U_IN;
    v.fill(0.0);
    // Pressure: linear Poiseuille gradient — eliminates long-wave filling transient.
    // dp/dx = -12*NU*U_IN / LY² for laminar channel (H_channel = LY)
    {
        const double dpdx = -12.0*NU*U_IN / (LY*LY);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                p(i, j) = dpdx * (mesh.xc[i] - LX);  // p=0 at outlet x=LX
    }
    applyBC(mesh, u, v, p);

    int    iter  = 0;
    double res_u = 1.0, res_v = 1.0, res_p = 1.0;

    for (iter = 1; iter <= MAX_ITER; ++iter) {

        // Snapshot for convergence check (max change per outer iteration)
        u_old = u;
        v_old = v;

        // ── Steps 1 & 2: momentum GS sweeps ──────────────────────────────────
        for (int inner = 0; inner < N_MOM_INNER; ++inner) {

            // x-momentum
            for (int j = 0; j < ny; ++j) {
                for (int i = 1; i < nx; ++i) {
                    if (mesh.uBlocked(i,j)) { u(i,j)=0.0; aPu(i,j)=aP_diff; continue; }

                    double u_e = 0.5*(u(i,j)+u(i+1,j));
                    double u_w = 0.5*(u(i-1,j)+u(i,j));
                    double v_n = 0.5*(v(i-1,j+1)+v(i,j+1));
                    double v_s = 0.5*(v(i-1,j  )+v(i,j  ));

                    double Fe = u_e*dy, Fw = u_w*dy;
                    double Fn = v_n*dx, Fs = v_s*dx;

                    double aE = NU*dy/dx + std::max(-Fe, 0.0);
                    double aW = NU*dy/dx + std::max( Fw, 0.0);
                    double aN = NU*dx/dy + std::max(-Fn, 0.0);
                    double aS = NU*dx/dy + std::max( Fs, 0.0);

                    if (j == 0)    aS += NU*dx/(0.5*dy);
                    if (j == ny-1) aN += NU*dx/(0.5*dy);

                    double aP0 = aE + aW + aN + aS;
                    double sp  = -(p(i,j) - p(i-1,j)) * dy;

                    double uN_nb = (j < ny-1) ? u(i,j+1) : 0.0;
                    double uS_nb = (j > 0)    ? u(i,j-1) : 0.0;

                    // Use u_old for under-relaxation reference (consistent across inner iters)
                    u(i,j) = (aE*u(i+1,j)+aW*u(i-1,j)+aN*uN_nb+aS*uS_nb+sp
                              + (1.0-ALPHA_U)*aP0*u_old(i,j)) / (aP0/ALPHA_U);
                    aPu(i,j) = aP0;
                }
            }

            // y-momentum
            for (int j = 1; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (mesh.vBlocked(i,j)) { v(i,j)=0.0; aPv(i,j)=aP_diff; continue; }

                    double v_n = 0.5*(v(i,j)+v(i,j+1));
                    double v_s = 0.5*(v(i,j-1)+v(i,j));
                    double u_e = 0.5*(u(i+1,j-1)+u(i+1,j));
                    double u_w = 0.5*(u(i,  j-1)+u(i,  j));

                    double Fe = u_e*dy, Fw = u_w*dy;
                    double Fn = v_n*dx, Fs = v_s*dx;

                    double aE = NU*dy/dx + std::max(-Fe, 0.0);
                    double aW = NU*dy/dx + std::max( Fw, 0.0);
                    double aN = NU*dx/dy + std::max(-Fn, 0.0);
                    double aS = NU*dx/dy + std::max( Fs, 0.0);

                    double aP0 = aE + aW + aN + aS;
                    double sp  = -(p(i,j) - p(i,j-1)) * dx;

                    double vE = (i < nx-1) ? v(i+1,j) : v(i,j);
                    double vW = (i > 0)    ? v(i-1,j) : v(i,j);

                    v(i,j) = (aE*vE+aW*vW+aN*v(i,j+1)+aS*v(i,j-1)+sp
                              + (1.0-ALPHA_V)*aP0*v_old(i,j)) / (aP0/ALPHA_V);
                    aPv(i,j) = aP0;
                }
            }
        } // end inner momentum loop

        // ── Step 3a: mass imbalance and Poisson coefficients ─────────────────
        res_p = 0.0;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (mesh.isSolid(i,j)) {
                    bfield(i,j) = 0.0;
                    ppAE(i,j) = ppAW(i,j) = ppAN(i,j) = ppAS(i,j) = 0.0;
                    continue;
                }
                double b = (u(i+1,j)-u(i,j))*dy + (v(i,j+1)-v(i,j))*dx;
                bfield(i,j) = b;
                res_p += std::abs(b);

                ppAE(i,j) = (i<nx-1 && !mesh.uBlocked(i+1,j)) ? dy*dy/(dx*aPu(i+1,j)) : 0.0;
                ppAW(i,j) = (i>0    && !mesh.uBlocked(i,  j)) ? dy*dy/(dx*aPu(i,  j)) : 0.0;
                ppAN(i,j) = (!mesh.vBlocked(i,j+1))            ? dx*dx/(dy*aPv(i,j+1)) : 0.0;
                ppAS(i,j) = (!mesh.vBlocked(i,j  ))            ? dx*dx/(dy*aPv(i,j  )) : 0.0;
            }
        }

        // ── Step 3b: SOR for pressure-correction Poisson equation ────────────
        // Unaccelerated ADI-TDMA: ~8% residual reduction per 100 sweep-pairs.
        // SOR with ω≈1.93: ~97% residual reduction per 100 sweeps (far superior).
        // BCs: pp(nx-1,j)=0 (Dirichlet outlet); Neumann elsewhere (aW=0 at inlet).
        pp.fill(0.0);

        for (int sw = 0; sw < N_PP_SWEEPS; ++sw) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx-1; ++i) {   // outlet column pp=0 fixed
                    if (mesh.isSolid(i,j)) { pp(i,j) = 0.0; continue; }
                    double aE=ppAE(i,j), aW=ppAW(i,j), aN=ppAN(i,j), aS=ppAS(i,j);
                    double aP_ = aE+aW+aN+aS;
                    if (aP_ < 1e-30) continue;
                    double ppE = (i<nx-2 && !mesh.isSolid(i+1,j)) ? pp(i+1,j) : 0.0;
                    double ppW = (i>0    && !mesh.isSolid(i-1,j)) ? pp(i-1,j) : 0.0;
                    double ppN = (j<ny-1 && !mesh.isSolid(i,j+1)) ? pp(i,j+1) : 0.0;
                    double ppS = (j>0    && !mesh.isSolid(i,j-1)) ? pp(i,j-1) : 0.0;
                    double pp_gs = (aE*ppE + aW*ppW + aN*ppN + aS*ppS - bfield(i,j)) / aP_;
                    pp(i,j) += OMEGA_SOR * (pp_gs - pp(i,j));
                }
            }
        }

        // ── Steps 4 & 5: pressure update and velocity correction ──────────────
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx-1; ++i)
                if (!mesh.isSolid(i,j))
                    p(i,j) += ALPHA_P * pp(i,j);

        for (int j = 0; j < ny; ++j)
            for (int i = 1; i < nx; ++i)
                if (!mesh.uBlocked(i,j)) {
                    double ppR = pp(i,   j);
                    double ppL = (i>0) ? pp(i-1,j) : pp(0,j);
                    u(i,j) -= (ppR-ppL)/dx * dy / aPu(i,j);
                }

        for (int j = 1; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                if (!mesh.vBlocked(i,j)) {
                    v(i,j) -= (pp(i,j)-pp(i,j-1))/dy * dx / aPv(i,j);
                }

        applyBC(mesh, u, v, p);

        // ── Convergence: max change in velocity between outer iterations ───────
        res_u = 0.0;
        res_v = 0.0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i <= nx; ++i)
                if (!mesh.uBlocked(i,j))
                    res_u = std::max(res_u, std::abs(u(i,j) - u_old(i,j)));
        for (int j = 1; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                if (!mesh.vBlocked(i,j))
                    res_v = std::max(res_v, std::abs(v(i,j) - v_old(i,j)));

        // ── Diagnostics ───────────────────────────────────────────────────────
        if (iter % PRINT_EVERY == 0 || iter == 1) {
            double umax = 0.0;
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    if (mesh.isSolid(i,j)) continue;
                    double uc = 0.5*(u(i,j)+u(i+1,j));
                    double vc = 0.5*(v(i,j)+v(i,j+1));
                    umax = std::max(umax, std::hypot(uc,vc));
                }
            std::cout << "iter=" << iter
                      << "  du=" << res_u
                      << "  dv=" << res_v
                      << "  res_p=" << res_p
                      << "  Umax=" << umax << std::endl;
        }

        if (res_u < TOL && res_v < TOL) {
            std::cout << "Converged at iter=" << iter << std::endl;
            break;
        }
    }

    return iter;
}
