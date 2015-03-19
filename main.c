#include "stdio.h"
#include "stdlib.h"

#include "Field_s.h"
#include "user_param.h"
#include "init_setting.h"
#include "boun_func.h"
#include "init_cond.h"

//#include "boun_func.c"



int main(){

    int     i;

	Grid_S	        grid;
	Solid	        solid;
	TimeMarching	timem;
	Index_S	        ind;
	Constraint_S    const_s;
    Field_S         field_s;
    Matrices_S      mat_s;

    /* ------- parameter setting-------------- */
	user_param(&grid,&solid,&timem,&const_s);
	set_index(&ind,&grid,solid.boundary_sign,const_s.fsineumanndirichlet);
    set_boundfunc(grid.dx);

    /*allocation @@ to be changed in the future for PETSC*/
    mat_s.MS_xx   = dmatrix(0,ind.xix_N-1,0,ind.xix_N-1);
    mat_s.MS_yy   = dmatrix(0,ind.xiy_N-1,0,ind.xiy_N-1);
    mat_s.KlS_xx  = dmatrix(0,ind.xix_N-1,0,ind.xix_N-1);
    mat_s.KlS_xy  = dmatrix(0,ind.xix_N-1,0,ind.xiy_N-1);
    mat_s.KlS_yy  = dmatrix(0,ind.xiy_N-1,0,ind.xiy_N-1);
    mat_s.KnlS_xx = dmatrix(0,ind.xix_N-1,0,ind.xix_N-1);
    mat_s.KnlS_yy = dmatrix(0,ind.xiy_N-1,0,ind.xiy_N-1);

    mat_s.NS_xc = dmatrix(0,ind.xix_N-1,0,ind.xix_Ncell_Neumann-1);
    mat_s.NS_yc = dmatrix(0,ind.xiy_N-1,0,ind.xiy_Ncell_Neumann-1);
    mat_s.DS_xc = dmatrix(0,ind.xix_N-1,0,ind.xix_Ncell_Dirichlet-1);
    mat_s.DS_yc = dmatrix(0,ind.xiy_N-1,0,ind.xiy_Ncell_Dirichlet-1);

    mat_s.FS_x  = dvector(0,ind.xix_N-1);
    mat_s.FS_y  = dvector(0,ind.xiy_N-1);

    /* ------- constructing governing matrices-------------- */
    compute_matricesNonlinearStructure(&mat_s, &ind, &grid, &solid, const_s.fsineumanndirichlet);
    reduce_system(&grid,&mat_s,&ind);

    /* ------- constructing governing matrices-------------- */
    field_s.xi_x   = dvector(0,grid.N-1);
    field_s.xi_y   = dvector(0,grid.N-1);
    field_s.dxi_x  = dvector(0,grid.N-1);
    field_s.dxi_y  = dvector(0,grid.N-1);
    field_s.ddxi_x = dvector(0,grid.N-1);
    field_s.ddxi_y = dvector(0,grid.N-1);

    set_initial(&grid, &solid, &field_s);


    /*
    for(i=0;i<grid.Nx;i++)
        for(j=0;j<grid.Ny;j++)
            printf("%d %d %f\n",i,j,solid.boundary_value[i][j]);
*/

   for(i=0;i<grid.N;i++)
        printf("%d %d\n",i,ind.xix.C2c_neumann[i]);

	return 0;
}


