#pragma once

#include <string>

#include "mesh.hpp"
#include "state.hpp"

void write_vtk(const std::string& filename,
               const Mesh2D& mesh,
               const PrimitiveFields& primitive);
