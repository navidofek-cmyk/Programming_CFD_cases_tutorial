#pragma once

#include "field.hpp"
#include "mesh.hpp"

ScalarField1D build_area_field(const Mesh1D& mesh, double throat_position);
