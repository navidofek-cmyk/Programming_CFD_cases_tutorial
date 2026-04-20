#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "config.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "solver.hpp"

int main() {
    const Mesh2D mesh =
        build_mesh(config::nx, config::ny, config::x_min, config::x_max, config::y_min, config::y_max);

    PrimitiveFields primitive{
        Field2D(mesh.nx, mesh.ny, 0.0),
        Field2D(mesh.nx, mesh.ny, 0.0),
        Field2D(mesh.nx, mesh.ny, 0.0),
        Field2D(mesh.nx, mesh.ny, 0.0),
        Field2D(mesh.nx, mesh.ny, 0.0)
    };

    ConservativeFields conservative{
        Field2D(mesh.nx, mesh.ny, 0.0),
        Field2D(mesh.nx, mesh.ny, 0.0),
        Field2D(mesh.nx, mesh.ny, 0.0),
        Field2D(mesh.nx, mesh.ny, 0.0)
    };

    initialize_quadrant_problem(mesh,
                                config::x_split,
                                config::y_split,
                                {config::rho_ul, config::u_ul, config::v_ul, config::p_ul},
                                {config::rho_ur, config::u_ur, config::v_ur, config::p_ur},
                                {config::rho_ll, config::u_ll, config::v_ll, config::p_ll},
                                {config::rho_lr, config::u_lr, config::v_lr, config::p_lr},
                                config::gamma_value,
                                primitive,
                                conservative);

    std::filesystem::create_directories("output");

    std::cout << "2D Euler quadrant Riemann problem\n";
    std::cout << "nx=" << config::nx
              << " ny=" << config::ny
              << " steps=" << config::n_steps
              << " CFL=" << config::cfl << '\n';

    for (int step = 1; step <= config::n_steps; ++step) {
        const StepReport report =
            advance_euler_step(mesh, config::gamma_value, config::cfl, primitive, conservative);

        if (step % config::report_interval == 0 || step == 1 || step == config::n_steps) {
            std::cout << "step " << step
                      << "  dt=" << report.dt
                      << "  drho_l2=" << report.rho_change
                      << "  dp_l2=" << report.p_change
                      << "  maxM=" << report.max_mach << '\n';
        }

        if (step % config::write_interval == 0 || step == config::n_steps) {
            std::ostringstream name;
            name << "output/euler2d_" << std::setw(5) << std::setfill('0') << step << ".vtk";
            write_vtk(name.str(), mesh, primitive);
        }
    }

    std::cout << "Wrote VTK snapshots in output/\n";
    return 0;
}
