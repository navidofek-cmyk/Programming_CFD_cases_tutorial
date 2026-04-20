#pragma once

#include "field.hpp"

struct PrimitiveState {
    double rho = 1.0;
    double u = 0.0;
    double v = 0.0;
    double p = 1.0;
};

struct ConservativeFields {
    Field2D U1;
    Field2D U2;
    Field2D U3;
    Field2D U4;
};

struct PrimitiveFields {
    Field2D rho;
    Field2D u;
    Field2D v;
    Field2D p;
    Field2D M;
    Field2D mask;
};

struct StepReport {
    double dt = 0.0;
    double rho_change = 0.0;
    double p_change = 0.0;
    double max_mach = 0.0;
};

PrimitiveState conservative_to_primitive(double U1, double U2, double U3, double U4, double gamma_value);
void primitive_to_conservative(const PrimitiveState& state,
                               double gamma_value,
                               double& U1,
                               double& U2,
                               double& U3,
                               double& U4);
double sound_speed(const PrimitiveState& state, double gamma_value);
