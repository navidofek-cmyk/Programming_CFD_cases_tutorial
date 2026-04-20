#include "io.hpp"

#include <fstream>
#include <iomanip>

void write_solution_csv(const std::string& filename,
                        const Mesh1D& mesh,
                        const ScalarField1D& area,
                        const PrimitiveFields& primitive) {
    std::ofstream out(filename);
    out << std::setprecision(16);
    out << "x,A,rho,u,p,T,M\n";

    for (int i = 0; i < mesh.nx; ++i) {
        out << mesh.x[static_cast<std::size_t>(i)] << ','
            << area(i) << ','
            << primitive.rho(i) << ','
            << primitive.u(i) << ','
            << primitive.p(i) << ','
            << primitive.T(i) << ','
            << primitive.M(i) << '\n';
    }
}
