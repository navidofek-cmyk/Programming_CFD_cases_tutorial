#include "solver.hpp"

#include <algorithm>
#include <cmath>

#include "flux.hpp"

namespace {

PrimitiveState primitive_at(const PrimitiveFields& primitive, int i, int j) {
    return {
        primitive.rho(i, j),
        primitive.u(i, j),
        primitive.v(i, j),
        primitive.p(i, j)
    };
}

PrimitiveState reflect_vertical(const PrimitiveState& state) {
    return {state.rho, -state.u, state.v, state.p};
}

PrimitiveState reflect_horizontal(const PrimitiveState& state) {
    return {state.rho, state.u, -state.v, state.p};
}

void update_primitive_from_conservative(const Mesh2D& mesh,
                                        double gamma_value,
                                        const ConservativeFields& conservative,
                                        PrimitiveFields& primitive) {
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_solid(i, j)) {
                primitive.rho(i, j) = 0.0;
                primitive.u(i, j) = 0.0;
                primitive.v(i, j) = 0.0;
                primitive.p(i, j) = 0.0;
                primitive.M(i, j) = 0.0;
                primitive.mask(i, j) = 1.0;
                continue;
            }

            const PrimitiveState state =
                conservative_to_primitive(conservative.U1(i, j),
                                          conservative.U2(i, j),
                                          conservative.U3(i, j),
                                          conservative.U4(i, j),
                                          gamma_value);

            primitive.rho(i, j) = state.rho;
            primitive.u(i, j) = state.u;
            primitive.v(i, j) = state.v;
            primitive.p(i, j) = state.p;
            primitive.M(i, j) = std::sqrt(state.u * state.u + state.v * state.v) / sound_speed(state, gamma_value);
            primitive.mask(i, j) = 0.0;
        }
    }
}

double compute_time_step(const Mesh2D& mesh,
                         double gamma_value,
                         double cfl,
                         const PrimitiveFields& primitive) {
    double dt = 1.0e30;

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_solid(i, j)) {
                continue;
            }
            const PrimitiveState s = primitive_at(primitive, i, j);
            const double a = sound_speed(s, gamma_value);
            const double spectral =
                (std::abs(s.u) + a) / mesh.dx
                + (std::abs(s.v) + a) / mesh.dy;
            dt = std::min(dt, cfl / spectral);
        }
    }

    return dt;
}

PrimitiveState face_left_state(const Mesh2D& mesh,
                               const PrimitiveFields& primitive,
                               const PrimitiveState& inlet_state,
                               int i,
                               int j) {
    if (i < 0) {
        return inlet_state;
    }
    if (mesh.is_solid(i, j)) {
        return reflect_vertical(primitive_at(primitive, i + 1, j));
    }
    return primitive_at(primitive, i, j);
}

PrimitiveState face_right_state(const Mesh2D& mesh,
                                const PrimitiveFields& primitive,
                                int i,
                                int j) {
    if (i >= mesh.nx) {
        return primitive_at(primitive, mesh.nx - 1, j);
    }
    if (mesh.is_solid(i, j)) {
        return reflect_vertical(primitive_at(primitive, i - 1, j));
    }
    return primitive_at(primitive, i, j);
}

PrimitiveState face_bottom_state(const Mesh2D& mesh,
                                 const PrimitiveFields& primitive,
                                 int i,
                                 int j) {
    if (j < 0) {
        return reflect_horizontal(primitive_at(primitive, i, 0));
    }
    if (mesh.is_solid(i, j)) {
        return reflect_horizontal(primitive_at(primitive, i, j + 1));
    }
    return primitive_at(primitive, i, j);
}

PrimitiveState face_top_state(const Mesh2D& mesh,
                              const PrimitiveFields& primitive,
                              int i,
                              int j) {
    if (j >= mesh.ny) {
        return primitive_at(primitive, i, mesh.ny - 1);
    }
    if (mesh.is_solid(i, j)) {
        return reflect_horizontal(primitive_at(primitive, i, j - 1));
    }
    return primitive_at(primitive, i, j);
}

}  // namespace

void initialize_oblique_shock_case(const Mesh2D& mesh,
                                   const PrimitiveState& inlet_state,
                                   double gamma_value,
                                   PrimitiveFields& primitive,
                                   ConservativeFields& conservative) {
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_solid(i, j)) {
                primitive.mask(i, j) = 1.0;
                primitive.rho(i, j) = 0.0;
                primitive.u(i, j) = 0.0;
                primitive.v(i, j) = 0.0;
                primitive.p(i, j) = 0.0;
                primitive.M(i, j) = 0.0;
                conservative.U1(i, j) = 1.0e-8;
                conservative.U2(i, j) = 0.0;
                conservative.U3(i, j) = 0.0;
                conservative.U4(i, j) = 1.0e-8 / (gamma_value - 1.0);
                continue;
            }

            const PrimitiveState state = inlet_state;
            primitive.rho(i, j) = state.rho;
            primitive.u(i, j) = state.u;
            primitive.v(i, j) = state.v;
            primitive.p(i, j) = state.p;
            primitive.M(i, j) = std::sqrt(state.u * state.u + state.v * state.v) / sound_speed(state, gamma_value);
            primitive.mask(i, j) = 0.0;

            primitive_to_conservative(state,
                                      gamma_value,
                                      conservative.U1(i, j),
                                      conservative.U2(i, j),
                                      conservative.U3(i, j),
                                      conservative.U4(i, j));
        }
    }
}

StepReport advance_euler_step(const Mesh2D& mesh,
                              double gamma_value,
                              double cfl,
                              PrimitiveFields& primitive,
                              ConservativeFields& conservative) {
    StepReport report{};
    report.dt = compute_time_step(mesh, gamma_value, cfl, primitive);
    const PrimitiveState inlet_state = primitive_at(primitive, 0, mesh.ny / 2);

    PrimitiveFields primitive_old = primitive;
    ConservativeFields old = conservative;

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_solid(i, j)) {
                conservative.U1(i, j) = old.U1(i, j);
                conservative.U2(i, j) = old.U2(i, j);
                conservative.U3(i, j) = old.U3(i, j);
                conservative.U4(i, j) = old.U4(i, j);
                continue;
            }

            const PrimitiveState left_face =
                face_left_state(mesh, primitive, inlet_state, i - 1, j);
            const PrimitiveState center = primitive_at(primitive, i, j);
            const PrimitiveState right_face =
                face_right_state(mesh, primitive, i + 1, j);

            const PrimitiveState bottom_face =
                face_bottom_state(mesh, primitive, i, j - 1);
            const PrimitiveState top_face =
                face_top_state(mesh, primitive, i, j + 1);

            double Fe1, Fe2, Fe3, Fe4;
            double Fw1, Fw2, Fw3, Fw4;
            double Gn1, Gn2, Gn3, Gn4;
            double Gs1, Gs2, Gs3, Gs4;

            rusanov_flux_x(center, right_face, gamma_value, Fe1, Fe2, Fe3, Fe4);
            rusanov_flux_x(left_face, center, gamma_value, Fw1, Fw2, Fw3, Fw4);
            rusanov_flux_y(center, top_face, gamma_value, Gn1, Gn2, Gn3, Gn4);
            rusanov_flux_y(bottom_face, center, gamma_value, Gs1, Gs2, Gs3, Gs4);

            conservative.U1(i, j) = old.U1(i, j)
                                    - report.dt / mesh.dx * (Fe1 - Fw1)
                                    - report.dt / mesh.dy * (Gn1 - Gs1);
            conservative.U2(i, j) = old.U2(i, j)
                                    - report.dt / mesh.dx * (Fe2 - Fw2)
                                    - report.dt / mesh.dy * (Gn2 - Gs2);
            conservative.U3(i, j) = old.U3(i, j)
                                    - report.dt / mesh.dx * (Fe3 - Fw3)
                                    - report.dt / mesh.dy * (Gn3 - Gs3);
            conservative.U4(i, j) = old.U4(i, j)
                                    - report.dt / mesh.dx * (Fe4 - Fw4)
                                    - report.dt / mesh.dy * (Gn4 - Gs4);
        }
    }

    update_primitive_from_conservative(mesh, gamma_value, conservative, primitive);

    double rho_sum = 0.0;
    double p_sum = 0.0;
    report.max_mach = 0.0;
    const double count = static_cast<double>(mesh.nx * mesh.ny);

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_solid(i, j)) {
                continue;
            }
            const double drho = primitive.rho(i, j) - primitive_old.rho(i, j);
            const double dp = primitive.p(i, j) - primitive_old.p(i, j);
            rho_sum += drho * drho;
            p_sum += dp * dp;
            report.max_mach = std::max(report.max_mach, primitive.M(i, j));
        }
    }

    report.rho_change = std::sqrt(rho_sum / count);
    report.p_change = std::sqrt(p_sum / count);
    return report;
}
