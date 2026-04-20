#include "solver.hpp"

#include <algorithm>
#include <cmath>

namespace {

constexpr double rho_floor = 1.0e-8;
constexpr double p_floor = 1.0e-8;

void update_conservative_from_primitive(const ScalarField1D& area,
                                        double gamma_value,
                                        const PrimitiveFields& primitive,
                                        ConservativeFields& conservative) {
    for (int i = 0; i < primitive.rho.size(); ++i) {
        const double rho = primitive.rho(i);
        const double u = primitive.u(i);
        const double p = primitive.p(i);
        const double total_energy = p / (gamma_value - 1.0) + 0.5 * rho * u * u;

        conservative.U1(i) = rho * area(i);
        conservative.U2(i) = rho * u * area(i);
        conservative.U3(i) = total_energy * area(i);
    }
}

void update_primitive_from_conservative(const ScalarField1D& area,
                                        double gamma_value,
                                        double gas_constant,
                                        const ConservativeFields& conservative,
                                        PrimitiveFields& primitive) {
    for (int i = 0; i < primitive.rho.size(); ++i) {
        const double rho = std::max(rho_floor, conservative.U1(i) / area(i));
        const double u = conservative.U2(i) / conservative.U1(i);
        const double energy_density = conservative.U3(i) / area(i);
        const double p = std::max(p_floor, pressure_from_state(rho, u, energy_density, gamma_value));
        const double T = p / (rho * gas_constant);
        const double a = sound_speed(p, rho, gamma_value);

        primitive.rho(i) = rho;
        primitive.u(i) = u;
        primitive.p(i) = p;
        primitive.T(i) = T;
        primitive.M(i) = u / a;
    }
}

void apply_boundary_conditions(const ScalarField1D& area,
                               double gamma_value,
                               double gas_constant,
                               double rho_inlet,
                               double u_inlet,
                               double p_inlet,
                               double p_exit,
                               PrimitiveFields& primitive,
                               ConservativeFields& conservative) {
    const int n = primitive.rho.size();

    // Simple subsonic inlet: fix density and pressure, extrapolate velocity
    // from the interior. This is less aggressive than fixing all primitives.
    primitive.rho(0) = rho_inlet;
    primitive.p(0) = p_inlet;
    primitive.u(0) = std::max(0.0, 2.0 * primitive.u(1) - primitive.u(2));

    // Simple subsonic outlet: extrapolate rho and u, fix static pressure.
    primitive.rho(n - 1) = std::max(rho_floor, 2.0 * primitive.rho(n - 2) - primitive.rho(n - 3));
    primitive.u(n - 1) = 2.0 * primitive.u(n - 2) - primitive.u(n - 3);
    primitive.p(n - 1) = p_exit;

    (void)u_inlet;

    for (int i : {0, n - 1}) {
        primitive.T(i) = primitive.p(i) / (primitive.rho(i) * gas_constant);
        primitive.M(i) = primitive.u(i) / sound_speed(primitive.p(i), primitive.rho(i), gamma_value);
    }

    update_conservative_from_primitive(area, gamma_value, primitive, conservative);
}

double compute_time_step(const Mesh1D& mesh,
                         double gamma_value,
                         double cfl,
                         const PrimitiveFields& primitive) {
    double dt = 1.0e30;

    for (int i = 0; i < mesh.nx; ++i) {
        const double a = sound_speed(primitive.p(i), primitive.rho(i), gamma_value);
        const double spectral_radius = std::abs(primitive.u(i)) + a;
        dt = std::min(dt, cfl * mesh.dx / spectral_radius);
    }

    return dt;
}

void stabilize_conservative_fields(const Mesh1D& mesh,
                                   double gamma_value,
                                   double artificial_viscosity,
                                   ConservativeFields& conservative) {
    if (artificial_viscosity <= 0.0) {
        return;
    }

    ConservativeFields original = conservative;

    for (int i = 1; i < mesh.nx - 1; ++i) {
        conservative.U1(i) += artificial_viscosity
                              * (original.U1(i + 1) - 2.0 * original.U1(i) + original.U1(i - 1));
        conservative.U2(i) += artificial_viscosity
                              * (original.U2(i + 1) - 2.0 * original.U2(i) + original.U2(i - 1));
        conservative.U3(i) += artificial_viscosity
                              * (original.U3(i + 1) - 2.0 * original.U3(i) + original.U3(i - 1));

        conservative.U1(i) = std::max(rho_floor * 1.0, conservative.U1(i));
        const double kinetic = 0.5 * conservative.U2(i) * conservative.U2(i) / conservative.U1(i);
        const double min_total_energy = kinetic + p_floor / (gamma_value - 1.0);
        conservative.U3(i) = std::max(min_total_energy, conservative.U3(i));
    }
}

}  // namespace

void initialize_solution(const Mesh1D& mesh,
                         const ScalarField1D& area,
                         double gamma_value,
                         double gas_constant,
                         double rho_initial,
                         double u_initial,
                         double p_initial,
                         PrimitiveFields& primitive,
                         ConservativeFields& conservative) {
    primitive.rho = ScalarField1D(mesh.nx, rho_initial);
    primitive.u = ScalarField1D(mesh.nx, u_initial);
    primitive.p = ScalarField1D(mesh.nx, p_initial);
    primitive.T = ScalarField1D(mesh.nx, p_initial / (rho_initial * gas_constant));
    primitive.M = ScalarField1D(mesh.nx, 0.0);

    conservative.U1 = ScalarField1D(mesh.nx, 0.0);
    conservative.U2 = ScalarField1D(mesh.nx, 0.0);
    conservative.U3 = ScalarField1D(mesh.nx, 0.0);

    for (int i = 0; i < mesh.nx; ++i) {
        // Gentle smooth perturbation helps the explicit solver evolve
        // without starting from an overly aggressive profile.
        const double x = mesh.x[static_cast<std::size_t>(i)];
        const double shape = 1.0 - 0.05 * (x - 1.5) * (x - 1.5);
        primitive.rho(i) = rho_initial * shape;
        primitive.u(i) = u_initial;
        primitive.p(i) = p_initial * shape;
        primitive.T(i) = primitive.p(i) / (primitive.rho(i) * gas_constant);
        primitive.M(i) = primitive.u(i) / sound_speed(primitive.p(i), primitive.rho(i), gamma_value);
    }

    update_conservative_from_primitive(area, gamma_value, primitive, conservative);
}

StepDiagnostics advance_maccormack(const Mesh1D& mesh,
                                   const ScalarField1D& area,
                                   double gamma_value,
                                   double gas_constant,
                                   double cfl,
                                   double artificial_viscosity,
                                   double rho_inlet,
                                   double u_inlet,
                                   double p_inlet,
                                   double p_exit,
                                   PrimitiveFields& primitive,
                                   ConservativeFields& conservative) {
    StepDiagnostics diagnostics{};
    const int n = mesh.nx;

    PrimitiveFields primitive_old = primitive;
    ConservativeFields conservative_old = conservative;

    diagnostics.dt = compute_time_step(mesh, gamma_value, cfl, primitive);

    ScalarField1D F1(n, 0.0);
    ScalarField1D F2(n, 0.0);
    ScalarField1D F3(n, 0.0);
    ScalarField1D J2(n, 0.0);

    for (int i = 0; i < n; ++i) {
        const double rho = primitive.rho(i);
        const double u = primitive.u(i);
        const double p = primitive.p(i);
        const double energy = conservative.U3(i) / area(i);

        F1(i) = conservative.U2(i);
        F2(i) = (rho * u * u + p) * area(i);
        F3(i) = u * (energy + p) * area(i);
        const double area_right = area(i + (i < n - 1 ? 1 : 0));
        const double area_left = area(i - (i > 0 ? 1 : 0));
        J2(i) = p * (area_right - area_left) / (2.0 * mesh.dx);
    }

    ConservativeFields predictor = conservative;

    for (int i = 1; i < n - 1; ++i) {
        predictor.U1(i) = conservative.U1(i) - diagnostics.dt / mesh.dx * (F1(i + 1) - F1(i));
        predictor.U2(i) = conservative.U2(i) - diagnostics.dt / mesh.dx * (F2(i + 1) - F2(i)) + diagnostics.dt * J2(i);
        predictor.U3(i) = conservative.U3(i) - diagnostics.dt / mesh.dx * (F3(i + 1) - F3(i));
    }

    PrimitiveFields primitive_predictor = primitive;
    update_primitive_from_conservative(area, gamma_value, gas_constant, predictor, primitive_predictor);
    apply_boundary_conditions(area, gamma_value, gas_constant, rho_inlet, u_inlet, p_inlet, p_exit, primitive_predictor, predictor);

    ScalarField1D F1p(n, 0.0);
    ScalarField1D F2p(n, 0.0);
    ScalarField1D F3p(n, 0.0);
    ScalarField1D J2p(n, 0.0);

    for (int i = 0; i < n; ++i) {
        const double rho = primitive_predictor.rho(i);
        const double u = primitive_predictor.u(i);
        const double p = primitive_predictor.p(i);
        const double energy = predictor.U3(i) / area(i);

        F1p(i) = predictor.U2(i);
        F2p(i) = (rho * u * u + p) * area(i);
        F3p(i) = u * (energy + p) * area(i);
        const double area_right = area(i + (i < n - 1 ? 1 : 0));
        const double area_left = area(i - (i > 0 ? 1 : 0));
        J2p(i) = p * (area_right - area_left) / (2.0 * mesh.dx);
    }

    for (int i = 1; i < n - 1; ++i) {
        conservative.U1(i) =
            0.5 * (conservative_old.U1(i) + predictor.U1(i)
                   - diagnostics.dt / mesh.dx * (F1p(i) - F1p(i - 1)));

        conservative.U2(i) =
            0.5 * (conservative_old.U2(i) + predictor.U2(i)
                   - diagnostics.dt / mesh.dx * (F2p(i) - F2p(i - 1))
                   + diagnostics.dt * J2p(i));

        conservative.U3(i) =
            0.5 * (conservative_old.U3(i) + predictor.U3(i)
                   - diagnostics.dt / mesh.dx * (F3p(i) - F3p(i - 1)));
    }

    stabilize_conservative_fields(mesh, gamma_value, artificial_viscosity, conservative);
    update_primitive_from_conservative(area, gamma_value, gas_constant, conservative, primitive);
    apply_boundary_conditions(area, gamma_value, gas_constant, rho_inlet, u_inlet, p_inlet, p_exit, primitive, conservative);

    double rho_sum = 0.0;
    double u_sum = 0.0;
    double p_sum = 0.0;
    diagnostics.max_mach = 0.0;

    for (int i = 0; i < n; ++i) {
        rho_sum += (primitive.rho(i) - primitive_old.rho(i)) * (primitive.rho(i) - primitive_old.rho(i));
        u_sum += (primitive.u(i) - primitive_old.u(i)) * (primitive.u(i) - primitive_old.u(i));
        p_sum += (primitive.p(i) - primitive_old.p(i)) * (primitive.p(i) - primitive_old.p(i));
        diagnostics.max_mach = std::max(diagnostics.max_mach, std::abs(primitive.M(i)));
    }

    diagnostics.rho_change = std::sqrt(rho_sum / static_cast<double>(n));
    diagnostics.u_change = std::sqrt(u_sum / static_cast<double>(n));
    diagnostics.p_change = std::sqrt(p_sum / static_cast<double>(n));

    return diagnostics;
}
