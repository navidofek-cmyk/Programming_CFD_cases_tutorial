#pragma once

#include "mesh.hpp"
#include "state.hpp"

void initialize_oblique_shock_case(const Mesh2D& mesh,
                                   const PrimitiveState& inlet_state,
                                   double gamma_value,
                                   PrimitiveFields& primitive,
                                   ConservativeFields& conservative);

StepReport advance_euler_step(const Mesh2D& mesh,
                              double gamma_value,
                              double cfl,
                              PrimitiveFields& primitive,
                              ConservativeFields& conservative);
