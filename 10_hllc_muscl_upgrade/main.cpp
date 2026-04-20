#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "config.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "solver.hpp"

int main() {
    const Mesh2D mesh = build_mesh(config::nx,
                                   config::ny,
                                   config::x_min,
                                   config::x_max,
                                   config::y_min,
                                   config::y_max,
                                   config::ramp_start_x,
                                   config::ramp_angle_deg);

    PrimitiveFields primitive{
        Field2D(mesh.nx, mesh.ny, 0.0),
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

    initialize_oblique_shock_case(mesh,
                                  {config::rho_inlet, config::u_inlet, config::v_inlet, config::p_inlet},
                                  config::gamma_value,
                                  primitive,
                                  conservative);

    std::filesystem::create_directories("output");

    std::cout << "2D Euler oblique shock over a ramp (HLLC + MUSCL)\n";
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
            name << "output/oblique_shock_hllc_" << std::setw(5) << std::setfill('0') << step << ".vtk";
            write_vtk(name.str(), mesh, primitive);
        }
    }

    std::cout << "Wrote VTK snapshots in output/\n";
    return 0;
}
