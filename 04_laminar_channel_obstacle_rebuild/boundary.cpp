#include "boundary.hpp"
#include "config.hpp"

void applyBC(const Mesh& mesh, Field2D& u, Field2D& v, Field2D& p)
{
    const int nx = mesh.nx;
    const int ny = mesh.ny;

    // ── Inlet (u-face i=0): Dirichlet u=U_IN, p zero-gradient ─────────────
    for (int j = 0; j < ny; ++j) {
        u(0, j) = U_IN;
        p(0, j) = p(1, j);  // zero-gradient (Neumann) pressure at inlet
    }

    // ── Outlet (u-face i=nx): zero-gradient velocity, p=0 ──────────────────
    for (int j = 0; j < ny; ++j) {
        u(nx, j) = u(nx-1, j);   // extrapolate (zero-gradient)
        p(nx-1, j) = 0.0;        // reference pressure at last cell
    }

    // ── Top/bottom wall: v-faces are already zero via vBlocked; also fix
    //    the end u-faces at the physical boundary (j=-½, j=ny-½ ghost).
    //    For no-slip of the tangential (u) component: wall face u = 0.
    //    In this staggered scheme u(i,j) lives at y=(j+0.5)*dy, NOT at y=0.
    //    No-slip is enforced via the diffusion term coefficient (half-cell);
    //    here we just zero all blocked v-faces.
    for (int i = 0; i < nx; ++i) {
        v(i, 0)  = 0.0;  // bottom wall (j=0 face at y=0)
        v(i, ny) = 0.0;  // top wall    (j=ny face at y=Ly)
    }

    // ── Obstacle faces: zero all blocked u/v faces ──────────────────────────
    for (int j = 0; j < ny; ++j)
        for (int i = 1; i < nx; ++i)   // interior u-faces (skip inlet i=0 and outlet i=nx)
            if (mesh.uBlocked(i, j)) u(i, j) = 0.0;

    for (int j = 1; j < ny; ++j)       // interior v-faces (top/bottom already done)
        for (int i = 0; i < nx; ++i)
            if (mesh.vBlocked(i, j)) v(i, j) = 0.0;

    // ── Pressure at solid cells: average from fluid neighbours ─────────────
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (mesh.isSolid(i, j)) {
                double s = 0.0; int cnt = 0;
                auto nb = [&](int ii, int jj){
                    if (!mesh.isSolid(ii, jj)) { s += p(ii,jj); ++cnt; }
                };
                if (i > 0)    nb(i-1, j);
                if (i < nx-1) nb(i+1, j);
                if (j > 0)    nb(i, j-1);
                if (j < ny-1) nb(i, j+1);
                p(i,j) = cnt>0 ? s/cnt : 0.0;
            }
}
