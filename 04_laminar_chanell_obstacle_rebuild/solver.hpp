#pragma once
#include "mesh.hpp"
#include "field.hpp"

// Run the staggered SIMPLE iteration to steady state.
// u: (nx+1)×ny,  v: nx×(ny+1),  p: nx×ny
// All fields are initialised inside; returns number of iterations performed.
int runSolver(const Mesh& mesh, Field2D& u, Field2D& v, Field2D& p);
