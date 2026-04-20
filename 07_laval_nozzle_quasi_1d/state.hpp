#pragma once

#include "field.hpp"

struct PrimitiveFields {
    ScalarField1D rho;
    ScalarField1D u;
    ScalarField1D p;
    ScalarField1D T;
    ScalarField1D M;
    ScalarField1D mdot;
};

struct ConservativeFields {
    ScalarField1D U1;
    ScalarField1D U2;
    ScalarField1D U3;
};

struct StepDiagnostics {
    double dt = 0.0;
    double rho_change = 0.0;
    double u_change = 0.0;
    double p_change = 0.0;
    double max_mach = 0.0;
    double throat_mach = 0.0;
    double mass_flow_variation = 0.0;
};

double pressure_from_state(double rho, double u, double energy, double gamma_value);
double sound_speed(double p, double rho, double gamma_value);
