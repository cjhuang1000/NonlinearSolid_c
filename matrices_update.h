#ifndef MATRICES_UPDATE_H
#define MATRICES_UPDATE_H 1

#include "Field_s.h"
#include "boun_func.h"
#include "nrutil.h"
# include <sdtlib.h>

void compute_matricesNonlinearStructure_update(Matrices_S* ptr_ms, Index_S* ptr_i, Grid_S* ptr_g, Solid* ptr_s);

#endif  /* MATRICES_UPDATE_H */
