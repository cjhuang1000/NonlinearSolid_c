# include "matrices_update.h"


boundfunc_cell   compute_matricesSubcell(double bound_value[4], char bd_type, int* ptr_type);
void compute_stress_update();


void reduce_system(Grid_S* ptr_g,Matrices_S* ptr_m,Index_S* ptr_i)
{
    /*% REDUCE_SYSTEM reduces the initial system of equation  A * t[u, v, p] = F by removing:
    - irrelevant unknowns in u, v or p which have no-influence on the domain
    - by getting rid of the arbitrary solid-body-motion*/

    int         n_xix =0,n_xiy =0, n_cellx_neu =0,n_celly_neu =0,n_cellx_dir =0,n_celly_dir=0;
    int         i,j;
    int         *index_keep_xix,*index_keep_xiy;
    int         *index_keep_stagx_Neumann,   *index_keep_stagy_Neumann;
    int         *index_keep_stagx_Dirichlet, *index_keep_stagy_Dirichlet;
    Index_S     ind_keep;
    Matrices_S  mat_keep;

    index_keep_xix              = (int*) calloc(sizeof(int),ptr_i->xix_N);
    index_keep_xiy              = (int*) calloc(sizeof(int),ptr_i->xiy_N);
    index_keep_stagx_Neumann    = (int*) calloc(sizeof(int),ptr_i->xix_Ncell_Neumann);
    index_keep_stagy_Neumann    = (int*) calloc(sizeof(int),ptr_i->xiy_Ncell_Neumann);
    index_keep_stagx_Dirichlet  = (int*) calloc(sizeof(int),ptr_i->xix_Ncell_Dirichlet );
    index_keep_stagy_Dirichlet  = (int*) calloc(sizeof(int),ptr_i->xiy_Ncell_Dirichlet );

    /*List of the indices of variables which are not outside the domain*/
    for(i=0; i<ptr_i->xix_N;i++)
    {
        if( ( ptr_m->MS_xx[i][i] != 0 ) || ( ptr_m->KlS_xx[i][i] != 0 ))
        {
            index_keep_xix[n_xix] = i;
            n_xix +=1;
        }
    }

    for(i=0; i<ptr_i->xiy_N;i++)
    {
        if( ( ptr_m->MS_yy[i][i] != 0 ) || ( ptr_m->KlS_yy[i][i] != 0 ))
        {
            index_keep_xiy[n_xiy] = i;
            n_xiy +=1;
        }
    }

   /*List of the indices of cells acting on boundary*/
    for(i=0; i<ptr_i->xix_N;i++)
        for(j=0;j< ptr_i->xix_Ncell_Neumann; j++)
            if (ptr_m->NS_xc[i][j] != 0)
            {
                index_keep_stagx_Neumann[n_cellx_neu] = i;
                n_cellx_neu +=1;
                break;
            }

    for(i=0; i<ptr_i->xiy_N;i++)
        for(j=0;j< ptr_i->xiy_Ncell_Neumann; j++)
            if (ptr_m->NS_yc[i][j] != 0)
            {
                index_keep_stagy_Neumann[n_celly_neu] = i;
                n_celly_neu +=1;
                break;
            }
    for(i=0; i<ptr_i->xix_N;i++)
        for(j=0;j< ptr_i->xix_Ncell_Dirichlet; j++)
            if (ptr_m->DS_xc[i][j] != 0)
            {
                index_keep_stagx_Dirichlet[n_cellx_dir] = i;
                n_cellx_dir +=1;
                break;
            }

    for(i=0; i<ptr_i->xiy_N;i++)
        for(j=0;j< ptr_i->xiy_Ncell_Dirichlet; j++)
            if (ptr_m->DS_yc[i][j] != 0)
            {
                index_keep_stagy_Dirichlet[n_celly_dir] = i;
                n_celly_dir +=1;
                break;
            }

    /*----------------------Update matrices-------------------------*/

    mat_keep.MS_xx   = dmatrix(0,n_xix-1,0,n_xix-1);
    mat_keep.MS_yy   = dmatrix(0,n_xiy-1,0,n_xiy-1);
    mat_keep.KlS_xx  = dmatrix(0,n_xix-1,0,n_xix-1);
    mat_keep.KlS_xy  = dmatrix(0,n_xix-1,0,n_xiy-1);
    mat_keep.KlS_yy  = dmatrix(0,n_xiy-1,0,n_xiy-1);
    mat_keep.KnlS_xx = dmatrix(0,n_xix-1,0,n_xix-1);
    mat_keep.KnlS_yy = dmatrix(0,n_xiy-1,0,n_xiy-1);

    mat_keep.NS_xc = dmatrix(0,n_xix-1,0,n_cellx_neu-1);
    mat_keep.NS_yc = dmatrix(0,n_xiy-1,0,n_celly_neu-1);
    mat_keep.DS_xc = dmatrix(0,n_xix-1,0,n_cellx_dir-1);
    mat_keep.DS_yc = dmatrix(0,n_xiy-1,0,n_cellx_dir-1);

    mat_keep.FS_x  = dvector(0,n_xix-1);
    mat_keep.FS_y  = dvector(0,n_xiy-1);

    for(i=0;i<n_xix;i++)
        for(j=0;j<n_xix;j++)
        {
            mat_keep.MS_xx[i][j]  = ptr_ms->MS_xx[index_keep_xix[i]][index_keep_xix[j]];
            mat_keep.KlS_xx[i][j] = ptr_ms->KlS_xx[index_keep_xix[i]][index_keep_xix[j]];
        }

    for(i=0;i<n_xiy;i++)
        for(j=0;j<n_xiy;j++)
        {
            mat_keep.MS_yy[i][j]  = ptr_ms->MS_yy[index_keep_xiy[i]][index_keep_xiy[j]];
            mat_keep.KlS_yy[i][j] = ptr_ms->KlS_yy[index_keep_xiy[i]][index_keep_xiy[j]];
        }

    for(i=0;i<n_xix;i++)
        for(j=0;j<n_xiy;j++)
            mat_keep.KlS_xy[i][j] = ptr_ms->KlS_xy[index_keep_xix[i]][index_keep_xiy[j]];

    for(i=0;i<n_xix;i++)
    {
        for(j=0;j<n_cellx_neu;j++)
            mat_keep.NS_xc[i][j]  = ptr_ms->NS_xc[index_keep_xix[i]][index_keep_stagx_Neumann[j]];
        for(j=0;j<n_cellx_dir;j++)
            mat_keep.DS_xc[i][j]  = ptr_ms->DS_xc[index_keep_xix[i]][index_keep_stagx_Dirichlet[j]];
    }
    for(i=0;i<n_xiy;i++)
    {
        for(j=0;j<n_celly_neu;j++)
            mat_keep.NS_yc[i][j]  = ptr_ms->NS_yc[index_keep_xiy[i]][index_keep_stagy_Neumann[j]];
        for(j=0;j<n_celly_dir;j++)
            mat_keep.DS_yc[i][j]  = ptr_ms->DS_yc[index_keep_xiy[i]][index_keep_stagy_Dirichlet[j]];
    }

    /* Update the G2g and C2c arrays vvvv*/

    free memory
    have new index and matrices

}

