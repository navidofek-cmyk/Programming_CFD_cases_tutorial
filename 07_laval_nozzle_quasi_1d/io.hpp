#pragma once

#include <string>

#include "field.hpp"
#include "mesh.hpp"
#include "state.hpp"

void write_solution_csv(const std::string& filename,
                        const Mesh1D& mesh,
                        const ScalarField1D& area,
                        const PrimitiveFields& primitive);
