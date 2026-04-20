#include <iostream>

#include "config.hpp"
#include "geometry.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "solver.hpp"

int main() {
    const Mesh1D mesh = build_mesh(config::nx, config::x_min, config::x_max);
    const ScalarField1D area = build_area_field(mesh, config::throat_position);

    PrimitiveFields primitive{
        ScalarField1D(mesh.nx, 0.0),
        ScalarField1D(mesh.nx, 0.0),
        ScalarField1D(mesh.nx, 0.0),
        ScalarField1D(mesh.nx, 0.0),
        ScalarField1D(mesh.nx, 0.0)
    };

    ConservativeFields conservative{
        ScalarField1D(mesh.nx, 0.0),
        ScalarField1D(mesh.nx, 0.0),
        ScalarField1D(mesh.nx, 0.0)
    };

    initialize_solution(mesh,
                        area,
                        config::gamma_value,
                        config::gas_constant,
                        config::rho_initial,
                        config::u_initial,
                        config::p_initial,
                        primitive,
                        conservative);


    std::cout << "Quasi-1D Laval nozzle solver\n";
    std::cout << "nx=" << config::nx
              << " steps=" << config::n_steps
              << " CFL=" << config::cfl << '\n';

    for (int step = 1; step <= config::n_steps; ++step) {
        const StepDiagnostics diagnostics =
            advance_maccormack(mesh,
                               area,
                               config::gamma_value,
                               config::gas_constant,
                               config::cfl,
                               config::artificial_viscosity,
                               config::rho_inlet,
                               config::u_inlet,
                               config::p_inlet,
                               config::p_exit,
                               primitive,
                               conservative);

        if (step % config::report_interval == 0 || step == 1 || step == config::n_steps) {
            std::cout << "step " << step
                      << "  dt=" << diagnostics.dt
                      << "  drho_l2=" << diagnostics.rho_change
                      << "  du_l2=" << diagnostics.u_change
                      << "  dp_l2=" << diagnostics.p_change
                      << "  maxM=" << diagnostics.max_mach << '\n';
        }
    }

    write_solution_csv("solution.csv", mesh, area, primitive);
    std::cout << "Wrote solution.csv\n";
    return 0;
}
