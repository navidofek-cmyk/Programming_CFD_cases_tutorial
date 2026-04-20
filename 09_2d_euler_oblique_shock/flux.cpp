#include "flux.hpp"

#include <algorithm>

namespace {

void physical_flux_x(const PrimitiveState& s,
                     double gamma_value,
                     double& F1,
                     double& F2,
                     double& F3,
                     double& F4) {
    double U1, U2, U3, U4;
    primitive_to_conservative(s, gamma_value, U1, U2, U3, U4);

    F1 = U2;
    F2 = U2 * s.u + s.p;
    F3 = U2 * s.v;
    F4 = (U4 + s.p) * s.u;
}

void physical_flux_y(const PrimitiveState& s,
                     double gamma_value,
                     double& G1,
                     double& G2,
                     double& G3,
                     double& G4) {
    double U1, U2, U3, U4;
    primitive_to_conservative(s, gamma_value, U1, U2, U3, U4);

    G1 = U3;
    G2 = U3 * s.u;
    G3 = U3 * s.v + s.p;
    G4 = (U4 + s.p) * s.v;
}

}  // namespace

void rusanov_flux_x(const PrimitiveState& left,
                    const PrimitiveState& right,
                    double gamma_value,
                    double& F1,
                    double& F2,
                    double& F3,
                    double& F4) {
    double FL1, FL2, FL3, FL4;
    double FR1, FR2, FR3, FR4;
    physical_flux_x(left, gamma_value, FL1, FL2, FL3, FL4);
    physical_flux_x(right, gamma_value, FR1, FR2, FR3, FR4);

    double UL1, UL2, UL3, UL4;
    double UR1, UR2, UR3, UR4;
    primitive_to_conservative(left, gamma_value, UL1, UL2, UL3, UL4);
    primitive_to_conservative(right, gamma_value, UR1, UR2, UR3, UR4);

    const double lambda =
        std::max(std::abs(left.u) + sound_speed(left, gamma_value),
                 std::abs(right.u) + sound_speed(right, gamma_value));

    F1 = 0.5 * (FL1 + FR1) - 0.5 * lambda * (UR1 - UL1);
    F2 = 0.5 * (FL2 + FR2) - 0.5 * lambda * (UR2 - UL2);
    F3 = 0.5 * (FL3 + FR3) - 0.5 * lambda * (UR3 - UL3);
    F4 = 0.5 * (FL4 + FR4) - 0.5 * lambda * (UR4 - UL4);
}

void rusanov_flux_y(const PrimitiveState& bottom,
                    const PrimitiveState& top,
                    double gamma_value,
                    double& G1,
                    double& G2,
                    double& G3,
                    double& G4) {
    double GB1, GB2, GB3, GB4;
    double GT1, GT2, GT3, GT4;
    physical_flux_y(bottom, gamma_value, GB1, GB2, GB3, GB4);
    physical_flux_y(top, gamma_value, GT1, GT2, GT3, GT4);

    double UB1, UB2, UB3, UB4;
    double UT1, UT2, UT3, UT4;
    primitive_to_conservative(bottom, gamma_value, UB1, UB2, UB3, UB4);
    primitive_to_conservative(top, gamma_value, UT1, UT2, UT3, UT4);

    const double lambda =
        std::max(std::abs(bottom.v) + sound_speed(bottom, gamma_value),
                 std::abs(top.v) + sound_speed(top, gamma_value));

    G1 = 0.5 * (GB1 + GT1) - 0.5 * lambda * (UT1 - UB1);
    G2 = 0.5 * (GB2 + GT2) - 0.5 * lambda * (UT2 - UB2);
    G3 = 0.5 * (GB3 + GT3) - 0.5 * lambda * (UT3 - UB3);
    G4 = 0.5 * (GB4 + GT4) - 0.5 * lambda * (UT4 - UB4);
}
