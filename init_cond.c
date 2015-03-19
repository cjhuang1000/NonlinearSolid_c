#include "init_cond.h"
#include <string.h>

void set_initial(Grid_S* ptr_g, Solid* ptr_s, Field_S* ptr_fs)
{

    int     i,j,k;
    int     nx = ptr_g->Nx, ny = ptr_g->Ny;


    /* set initial condition*/
    if (strcmp(ptr_s->initial,"rest")==0)
    {
        for(i=0; i<nx; i++)
            for(j=0; j<ny; j++){

                k = i+ nx*j;
                ptr_fs->xi_x[k]   = 0.0;
                ptr_fs->xi_y[k]   = 0.0;
                ptr_fs->dxi_x[k]  = 0.0;
                ptr_fs->dxi_y[k]  = 0.0;
                ptr_fs->ddxi_x[k] = 0.0;
                ptr_fs->ddxi_y[k] = 0.0;

        }
    }

}
