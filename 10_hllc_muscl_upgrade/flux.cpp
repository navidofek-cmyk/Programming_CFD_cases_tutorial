#include "flux.hpp"

#include <algorithm>
#include <cmath>

namespace {

struct EulerState1D {
    double rho = 1.0;
    double un = 0.0;
    double ut = 0.0;
    double p = 1.0;
    double E = 0.0;
};

EulerState1D rotate_to_x(const PrimitiveState& s, double gamma_value) {
    double U1, U2, U3, U4;
    primitive_to_conservative(s, gamma_value, U1, U2, U3, U4);
    return {s.rho, s.u, s.v, s.p, U4};
}

EulerState1D rotate_to_y(const PrimitiveState& s, double gamma_value) {
    double U1, U2, U3, U4;
    primitive_to_conservative(s, gamma_value, U1, U2, U3, U4);
    return {s.rho, s.v, s.u, s.p, U4};
}

void physical_flux_1d(const EulerState1D& s,
                      double& F1,
                      double& F2,
                      double& F3,
                      double& F4) {
    F1 = s.rho * s.un;
    F2 = s.rho * s.un * s.un + s.p;
    F3 = s.rho * s.un * s.ut;
    F4 = (s.E + s.p) * s.un;
}

void hllc_flux_1d(const EulerState1D& left,
                  const EulerState1D& right,
                  double gamma_value,
                  double& F1,
                  double& F2,
                  double& F3,
                  double& F4) {
    double FL1, FL2, FL3, FL4;
    double FR1, FR2, FR3, FR4;
    physical_flux_1d(left, FL1, FL2, FL3, FL4);
    physical_flux_1d(right, FR1, FR2, FR3, FR4);

    const double aL = std::sqrt(gamma_value * left.p / left.rho);
    const double aR = std::sqrt(gamma_value * right.p / right.rho);

    const double SL = std::min(left.un - aL, right.un - aR);
    const double SR = std::max(left.un + aL, right.un + aR);

    const double numerator =
        right.p - left.p
        + left.rho * left.un * (SL - left.un)
        - right.rho * right.un * (SR - right.un);
    const double denominator =
        left.rho * (SL - left.un)
        - right.rho * (SR - right.un);
    const double SM = numerator / denominator;

    const double rho_star_L = left.rho * (SL - left.un) / (SL - SM);
    const double rho_star_R = right.rho * (SR - right.un) / (SR - SM);

    const double E_star_L =
        ((SL - left.un) * left.E - left.p * left.un + left.p * SM) / (SL - SM);
    const double E_star_R =
        ((SR - right.un) * right.E - right.p * right.un + right.p * SM) / (SR - SM);

    const double UstarL1 = rho_star_L;
    const double UstarL2 = rho_star_L * SM;
    const double UstarL3 = rho_star_L * left.ut;
    const double UstarL4 = E_star_L;

    const double UstarR1 = rho_star_R;
    const double UstarR2 = rho_star_R * SM;
    const double UstarR3 = rho_star_R * right.ut;
    const double UstarR4 = E_star_R;

    const double UL1 = left.rho;
    const double UL2 = left.rho * left.un;
    const double UL3 = left.rho * left.ut;
    const double UL4 = left.E;

    const double UR1 = right.rho;
    const double UR2 = right.rho * right.un;
    const double UR3 = right.rho * right.ut;
    const double UR4 = right.E;

    if (0.0 <= SL) {
        F1 = FL1; F2 = FL2; F3 = FL3; F4 = FL4;
    } else if (SL <= 0.0 && 0.0 <= SM) {
        F1 = FL1 + SL * (UstarL1 - UL1);
        F2 = FL2 + SL * (UstarL2 - UL2);
        F3 = FL3 + SL * (UstarL3 - UL3);
        F4 = FL4 + SL * (UstarL4 - UL4);
    } else if (SM <= 0.0 && 0.0 <= SR) {
        F1 = FR1 + SR * (UstarR1 - UR1);
        F2 = FR2 + SR * (UstarR2 - UR2);
        F3 = FR3 + SR * (UstarR3 - UR3);
        F4 = FR4 + SR * (UstarR4 - UR4);
    } else {
        F1 = FR1; F2 = FR2; F3 = FR3; F4 = FR4;
    }
}

}  // namespace

void hllc_flux_x(const PrimitiveState& left,
                 const PrimitiveState& right,
                 double gamma_value,
                 double& F1,
                 double& F2,
                 double& F3,
                 double& F4) {
    hllc_flux_1d(rotate_to_x(left, gamma_value),
                 rotate_to_x(right, gamma_value),
                 gamma_value,
                 F1, F2, F3, F4);
}

void hllc_flux_y(const PrimitiveState& bottom,
                 const PrimitiveState& top,
                 double gamma_value,
                 double& G1,
                 double& G2,
                 double& G3,
                 double& G4) {
    hllc_flux_1d(rotate_to_y(bottom, gamma_value),
                 rotate_to_y(top, gamma_value),
                 gamma_value,
                 G1, G3, G2, G4);
}
