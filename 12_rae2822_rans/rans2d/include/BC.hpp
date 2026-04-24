#pragma once
#include "Solver.hpp"

void bc_wall        (Solver& s);  // j=-1 ghost: inviscid slip (Stage 2)
void bc_wall_noslip (Solver& s);  // j=-1 ghost: no-slip adiabatic (Stage 3+)
void bc_farfield    (Solver& s);  // j=ncj ghost: Riemann characteristic
void bc_wake_cut    (Solver& s);  // j=-1 wake ghosts: mirror across cut
void bc_i_exit      (Solver& s);  // i=-1, i=nci ghosts: zeroth-order extrapolation
