#pragma once

#include "field.hpp"
#include "mesh.hpp"
#include "state.hpp"

void initialize_solution(const Mesh1D& mesh,
                         const ScalarField1D& area,
                         double gamma_value,
                         double gas_constant,
                         double rho_initial,
                         double u_initial,
                         double p_initial,
                         PrimitiveFields& primitive,
                         ConservativeFields& conservative);

StepDiagnostics advance_maccormack(const Mesh1D& mesh,
                                   const ScalarField1D& area,
                                   double gamma_value,
                                   double gas_constant,
                                   double cfl,
                                   double artificial_viscosity,
                                   double T0_inlet,
                                   double p0_inlet,
                                   double p_exit,
                                   PrimitiveFields& primitive,
                                   ConservativeFields& conservative);
