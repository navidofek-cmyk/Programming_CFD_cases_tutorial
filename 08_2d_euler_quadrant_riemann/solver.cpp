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

void update_primitive_from_conservative(const Mesh2D& mesh,
                                        double gamma_value,
                                        const ConservativeFields& conservative,
                                        PrimitiveFields& primitive) {
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
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

PrimitiveState transmissive_x(const PrimitiveFields& primitive, int i, int j, int nx) {
    const int ii = std::clamp(i, 0, nx - 1);
    return primitive_at(primitive, ii, j);
}

PrimitiveState transmissive_y(const PrimitiveFields& primitive, int i, int j, int ny) {
    const int jj = std::clamp(j, 0, ny - 1);
    return primitive_at(primitive, i, jj);
}

}  // namespace

void initialize_quadrant_problem(const Mesh2D& mesh,
                                 double x_split,
                                 double y_split,
                                 const PrimitiveState& upper_left,
                                 const PrimitiveState& upper_right,
                                 const PrimitiveState& lower_left,
                                 const PrimitiveState& lower_right,
                                 double gamma_value,
                                 PrimitiveFields& primitive,
                                 ConservativeFields& conservative) {
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            const double x = mesh.xc[static_cast<std::size_t>(i)];
            const double y = mesh.yc[static_cast<std::size_t>(j)];

            PrimitiveState state{};
            if (x <= x_split && y >= y_split) {
                state = upper_left;
            } else if (x > x_split && y >= y_split) {
                state = upper_right;
            } else if (x <= x_split && y < y_split) {
                state = lower_left;
            } else {
                state = lower_right;
            }

            primitive.rho(i, j) = state.rho;
            primitive.u(i, j) = state.u;
            primitive.v(i, j) = state.v;
            primitive.p(i, j) = state.p;
            primitive.M(i, j) = std::sqrt(state.u * state.u + state.v * state.v) / sound_speed(state, gamma_value);

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

    PrimitiveFields primitive_old = primitive;
    ConservativeFields old = conservative;

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            const PrimitiveState left_face =
                transmissive_x(primitive, i - 1, j, mesh.nx);
            const PrimitiveState center = primitive_at(primitive, i, j);
            const PrimitiveState right_face =
                transmissive_x(primitive, i + 1, j, mesh.nx);

            const PrimitiveState bottom_face =
                transmissive_y(primitive, i, j - 1, mesh.ny);
            const PrimitiveState top_face =
                transmissive_y(primitive, i, j + 1, mesh.ny);

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
