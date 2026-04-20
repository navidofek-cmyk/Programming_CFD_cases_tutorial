#pragma once
#include <vector>
#include "config.hpp"

// ─── Staggered MAC grid ────────────────────────────────────────────────────
//
//  p(i,j)  : cell centre at (xc[i], yc[j])
//  u(i,j)  : x-velocity at west/east face of cell i   → x = i*dx, y = (j+0.5)*dy
//             size (nx+1) × ny   (i = 0..nx, j = 0..ny-1)
//  v(i,j)  : y-velocity at south/north face of cell j → x = (i+0.5)*dx, y = j*dy
//             size nx × (ny+1)   (i = 0..nx-1, j = 0..ny)
//
//  Obstacle cells: solid → all faces TOUCHING a solid cell are blocked (u=v=0).
//  Blocked u-face (i,j): cell (i-1,j) or cell (i,j) is solid
//  Blocked v-face (i,j): cell (i,j-1) or cell (i,j) is solid

struct Mesh {
    int nx, ny;
    double dx, dy;

    std::vector<double> xc; // cell-centre x coords  (size nx)
    std::vector<double> yc; // cell-centre y coords  (size ny)

    std::vector<bool> solidCell;   // solid[j*nx+i]  — solid cell mask

    Mesh(int nx_, int ny_);

    bool isSolid(int i, int j) const {
        if (i<0||i>=nx||j<0||j>=ny) return true; // out-of-domain = wall
        return solidCell[j*nx+i];
    }

    // True if u-face (i,j) is blocked (touches a solid or domain boundary)
    bool uBlocked(int i, int j) const {
        // i=0: inlet (free flow, NOT blocked here — BCs handled elsewhere)
        // i=nx: outlet (free flow)
        if (j<0||j>=ny) return true;
        bool left  = (i==0) ? false : isSolid(i-1,j);
        bool right = (i==nx)? false : isSolid(i,  j);
        return left || right;
    }

    // True if v-face (i,j) is blocked (touches a solid or domain boundary)
    bool vBlocked(int i, int j) const {
        if (i<0||i>=nx) return true;
        bool below = (j==0)  ? true  : isSolid(i,j-1); // bottom wall always blocked
        bool above = (j==ny) ? true  : isSolid(i,j);   // top wall always blocked
        return below || above;
    }
};
