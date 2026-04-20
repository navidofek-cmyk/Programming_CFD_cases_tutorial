#pragma once

#include "state.hpp"

void hllc_flux_x(const PrimitiveState& left,
                 const PrimitiveState& right,
                 double gamma_value,
                 double& F1,
                 double& F2,
                 double& F3,
                 double& F4);

void hllc_flux_y(const PrimitiveState& bottom,
                 const PrimitiveState& top,
                 double gamma_value,
                 double& G1,
                 double& G2,
                 double& G3,
                 double& G4);
