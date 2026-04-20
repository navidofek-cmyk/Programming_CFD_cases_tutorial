#pragma once
#include "mesh.hpp"
#include "field.hpp"

// Apply boundary conditions on the staggered fields.
// u: (nx+1)×ny,  v: nx×(ny+1),  p: nx×ny
void applyBC(const Mesh& mesh, Field2D& u, Field2D& v, Field2D& p);
