#pragma once
#include <string>
#include "mesh.hpp"
#include "field.hpp"

// Write a cell-centred scalar field to CSV.
void writeCSV(const std::string& filename,
              const Mesh& mesh,
              const Field2D& field);

// Write u and v (staggered face fields) averaged to cell centres.
// uf: (nx+1)×ny,  vf: nx×(ny+1)
void writeUV(const std::string& u_file, const std::string& v_file,
             const Mesh& mesh,
             const Field2D& uf, const Field2D& vf);

// Write obstacle mask (1=solid).
void writeMask(const std::string& filename, const Mesh& mesh);
