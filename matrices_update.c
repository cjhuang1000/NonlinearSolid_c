# include "matrices_update.h"


boundfunc_cell   compute_matricesSubcell(double bound_value[4], char bd_type, int* ptr_type);
void compute_stress_update();

void reduce_vector(short int** a, int* v, int n_b, int n_a);
void reduce_matrix(double*** a, int* vx, int* vy, int nx_b, int nx_a ,int ny_b, int ny_a);
void reduce_FSvector(double** a, int* v, int n_b, int n_a);

void reduce_system(Grid_S* ptr_g,Matrices_S* ptr_m,Index_S* ptr_i)
// reduce_system.m
{
    /*% REDUCE_SYSTEM reduces the initial system of equation  A * t[u, v, p] = F by removing:
    - irrelevant unknowns in u, v or p which have no-influence on the domain
    - by getting rid of the arbitrary solid-body-motion*/

    int         n_xix =0,n_xiy =0, n_cellx_neu =0,n_celly_neu =0,n_cellx_dir =0,n_celly_dir=0;
    int         i,j;
    int         *index_keep_xix,*index_keep_xiy;
    int         *index_keep_stagx_Neumann,   *index_keep_stagy_Neumann;
    int         *index_keep_stagx_Dirichlet, *index_keep_stagy_Dirichlet;

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

    /*update matrices*/
    reduce_matrix(&(ptr_m->MS_xx),index_keep_xix,index_keep_xix,ptr_i->xix_N,n_xix,ptr_i->xix_N,n_xix);
    reduce_matrix(&(ptr_m->MS_yy),index_keep_xiy,index_keep_xiy,ptr_i->xiy_N,n_xiy,ptr_i->xiy_N,n_xiy);

    reduce_matrix(&(ptr_m->KlS_xx),index_keep_xix,index_keep_xix,ptr_i->xix_N,n_xix,ptr_i->xix_N,n_xix);
    reduce_matrix(&(ptr_m->KlS_yy),index_keep_xiy,index_keep_xiy,ptr_i->xiy_N,n_xiy,ptr_i->xiy_N,n_xiy);
    reduce_matrix(&(ptr_m->KlS_xy),index_keep_xix,index_keep_xix,ptr_i->xix_N,n_xix,ptr_i->xiy_N,n_xiy);

    reduce_matrix(&(ptr_m->KnlS_xx),index_keep_xix,index_keep_xix,ptr_i->xix_N,n_xix,ptr_i->xix_N,n_xix);
    reduce_matrix(&(ptr_m->KnlS_yy),index_keep_xiy,index_keep_xiy,ptr_i->xiy_N,n_xiy,ptr_i->xiy_N,n_xiy);

    reduce_matrix(&(ptr_m->NS_xc),index_keep_xix,index_keep_stagx_Neumann,ptr_i->xix_N,n_xix,ptr_i->xix_Ncell_Neumann,n_cellx_neu);
    reduce_matrix(&(ptr_m->NS_xc),index_keep_xiy,index_keep_stagy_Neumann,ptr_i->xiy_N,n_xiy,ptr_i->xiy_Ncell_Neumann,n_celly_neu);

    reduce_matrix(&(ptr_m->DS_xc),index_keep_xix,index_keep_stagx_Dirichlet,ptr_i->xix_N,n_xix,ptr_i->xix_Ncell_Dirichlet,n_cellx_dir);
    reduce_matrix(&(ptr_m->DS_xc),index_keep_xiy,index_keep_stagy_Dirichlet,ptr_i->xiy_N,n_xiy,ptr_i->xiy_Ncell_Dirichlet,n_celly_dir);

    reduce_FSvector(&(ptr_m->FS_x),index_keep_xix,ptr_i->xix_N,n_xix);
    reduce_FSvector(&(ptr_m->FS_y),index_keep_xiy,ptr_i->xiy_N,n_xiy);

    /* --------Update the G2g and C2c arrays ------------*/

    //---------------preallocate ind_new and copy value into it

    /*add value to g2G,G2g_before and keep before freeing the memory in the original index*/
    ptr_i->xix.g2G_before = (short int *) malloc(sizeof(short int)*ptr_i->xix_N);
	ptr_i->xix.G2g_before = (short int *) malloc(sizeof(short int)*ptr_g->N);

	ptr_i->xiy.g2G_before = (short int *) malloc(sizeof(short int)*ptr_i->xiy_N);
	ptr_i->xiy.G2g_before = (short int *) malloc(sizeof(short int)*ptr_g->N);

	ptr_i->xix.keep       = (short int *) malloc(sizeof(short int)*n_xix);
	ptr_i->xiy.keep       = (short int *) malloc(sizeof(short int)*n_xiy);

    for(i=0;i<ptr_g->N;i++)
    {
        ptr_i->xix.G2g_before[i] = ptr_i->xix.G2g[i];
        ptr_i->xiy.G2g_before[i] = ptr_i->xiy.G2g[i];
    }

	for(i=0;i<n_xix;i++)
    {
        ptr_i->xix.keep[i] = index_keep_xix[i];
    }

	for(i=0;i<n_xiy;i++)
    {
        ptr_i->xiy.keep[i] = index_keep_xiy[i];
    }

    reduce_vector(&(ptr_i->xix.g2G),index_keep_xix,ptr_i->xix_N,n_xix);
    reduce_vector(&(ptr_i->xiy.g2G),index_keep_xiy,ptr_i->xiy_N,n_xiy);

    reduce_vector(&(ptr_i->xix.c2C_neumann),index_keep_stagx_Neumann,ptr_i->xix_Ncell_Neumann,n_cellx_neu);
    reduce_vector(&(ptr_i->xiy.c2C_neumann),index_keep_stagy_Neumann,ptr_i->xiy_Ncell_Neumann,n_celly_neu);

    reduce_vector(&(ptr_i->xix.c2C_dirichlet),index_keep_stagx_Dirichlet,ptr_i->xix_Ncell_Dirichlet,n_cellx_dir);
    reduce_vector(&(ptr_i->xiy.c2C_dirichlet),index_keep_stagy_Dirichlet,ptr_i->xiy_Ncell_Dirichlet,n_celly_dir);

    for(i=0;i<ptr_g->N;i++)
    {
        ptr_i->xix.G2g[i] = -1;
        ptr_i->xiy.G2g[i] = -1;
        ptr_i->xix.C2c_dirichlet[i] = -1;
        ptr_i->xiy.C2c_dirichlet[i] = -1;
        ptr_i->xix.C2c_neumann[i] = -1;
        ptr_i->xiy.C2c_neumann[i] = -1;
    }

    for(i=0;i<n_xix;i++)
        ptr_i->xix.G2g[ptr_i->xix.g2G[i]]=i;
    for(i=0;i<n_xiy;i++)
        ptr_i->xiy.G2g[ptr_i->xiy.g2G[i]]=i;
    for(i=0;i<n_cellx_neu;i++)
        ptr_i->xix.C2c_neumann[ptr_i->xix.c2C_neumann[i]]=i;
    for(i=0;i<n_celly_neu;i++)
        ptr_i->xiy.C2c_neumann[ptr_i->xiy.c2C_neumann[i]]=i;
    for(i=0;i<n_cellx_dir;i++)
        ptr_i->xix.C2c_dirichlet[ptr_i->xix.c2C_dirichlet[i]]=i;
    for(i=0;i<n_celly_dir;i++)
        ptr_i->xiy.C2c_dirichlet[ptr_i->xiy.c2C_dirichlet[i]]=i;

    ptr_i->xix_N = n_xix;
    ptr_i->xiy_N = n_xiy;
    ptr_i->xix_Ncell_Neumann = n_cellx_neu;
    ptr_i->xiy_Ncell_Neumann = n_celly_neu;
    ptr_i->xix_Ncell_Dirichlet = n_cellx_dir;
    ptr_i->xiy_Ncell_Dirichlet = n_celly_dir;

}

// reduce the size of the array
void reduce_vector(short int** a, int* v, int n_b, int n_a)
{
    int i;
    short int* temp;

    temp = (short int *) malloc(sizeof(short int)*n_b);
    for(i=0;i<n_b;i++)
        temp[i] = (*a)[i];

    free(*a);
    *a = (short int *) malloc(sizeof(short int)*n_a);
    for(i=0;i<n_a;i++)
        (*a)[i] = temp[v[i]];

}

// reduce the size of the governing matrices
void reduce_matrix(double*** a, int* vx, int* vy, int nx_b, int nx_a ,int ny_b, int ny_a)
{
    int i,j;
    double** temp;
    temp = dmatrix(0,nx_b-1,0,ny_b-1);

    for(i=0;i<nx_b;i++)
        for(j=0;j<ny_b;j++)
            temp[i][j] = (*a)[i][j];

    free_dmatrix(*a,0,nx_b-1,0,ny_b-1);
    *a = dmatrix(0,nx_a-1,0,ny_a-1);

     for(i=0;i<nx_a;i++)
        for(j=0;j<ny_a;j++)
            (*a)[i][j] = temp[vx[i]][vy[j]];

}

void reduce_FSvector(double** a, int* v, int n_b, int n_a)
{
    int i;
    double* temp;
    temp = dvector(0,n_b-1);

    for(i=0;i<n_b;i++)
            temp[i] = (*a)[i];

    free_dvector(*a,0,n_b-1);
    *a = dvector(0,n_a-1);

     for(i=0;i<n_a;i++)
            (*a)[i] = temp[v[i]];

}

