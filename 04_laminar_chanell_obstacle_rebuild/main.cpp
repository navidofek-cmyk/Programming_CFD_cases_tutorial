#include <iostream>
#include "config.hpp"
#include "mesh.hpp"
#include "field.hpp"
#include "solver.hpp"
#include "io.hpp"

int main()
{
    std::cout << "=== 2D Laminar Channel + Obstacle (staggered SIMPLE) ===\n";
    std::cout << "  Re=" << RE << "  nu=" << NU
              << "  Nx=" << NX << "  Ny=" << NY << "\n\n";

    Mesh mesh(NX, NY);

    // Staggered MAC fields:
    //   u: (NX+1)×NY — x-velocity at west/east faces
    //   v: NX×(NY+1) — y-velocity at south/north faces
    //   p: NX×NY     — pressure at cell centres
    Field2D u(NX+1, NY), v(NX, NY+1), p(NX, NY);

    int iters = runSolver(mesh, u, v, p);

    std::cout << "\nWriting output files...\n";
    writeUV("u.csv", "v.csv", mesh, u, v);
    writeCSV("p.csv", mesh, p);
    writeMask("mask.csv", mesh);

    std::cout << "Done.  Completed " << iters << " iterations.\n";
    return 0;
}
