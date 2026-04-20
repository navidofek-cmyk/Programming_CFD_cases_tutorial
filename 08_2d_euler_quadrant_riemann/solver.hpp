#pragma once

#include "mesh.hpp"
#include "state.hpp"

void initialize_quadrant_problem(const Mesh2D& mesh,
                                 double x_split,
                                 double y_split,
                                 const PrimitiveState& upper_left,
                                 const PrimitiveState& upper_right,
                                 const PrimitiveState& lower_left,
                                 const PrimitiveState& lower_right,
                                 double gamma_value,
                                 PrimitiveFields& primitive,
                                 ConservativeFields& conservative);

StepReport advance_euler_step(const Mesh2D& mesh,
                              double gamma_value,
                              double cfl,
                              PrimitiveFields& primitive,
                              ConservativeFields& conservative);
