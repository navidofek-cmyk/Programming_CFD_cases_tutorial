#include "mesh.hpp"

Mesh::Mesh(int nx_, int ny_)
    : nx(nx_), ny(ny_),
      dx(LX/nx_), dy(LY/ny_),
      xc(nx_), yc(ny_),
      solidCell(static_cast<std::size_t>(nx_*ny_), false)
{
    for (int i = 0; i < nx; ++i) xc[i] = (i+0.5)*dx;
    for (int j = 0; j < ny; ++j) yc[j] = (j+0.5)*dy;

    // Stair-step obstacle: mark cells whose centre falls inside the obstacle box
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (xc[i] >= OBS_X0 && xc[i] <= OBS_X1 &&
                yc[j] >= OBS_Y0 && yc[j] <= OBS_Y1)
                solidCell[j*nx+i] = true;
}
