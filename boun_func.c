#include "boun_func.h"
#include "init_setting.h"
#include <math.h>
#include "nrutil.h"

double boundFunct_MS_alpha(int k, int l,double x,double y,double alpha);
double boundFunct_MS_gamma(int k, int l,double x,double y,double gamma);

double boundFunct_KS_xx_alpha(int k, int l,double y1,double y2,double alpha, double beta);
double boundFunct_KS_xx_gamma(int k, int l,double x1,double x2,double gamma, double delta);
double boundFunct_KS_yy_alpha(int k, int l,double y1,double y2,double alpha, double beta);
double boundFunct_KS_yy_gamma(int k, int l,double x1,double x2,double gamma, double delta);
double boundFunct_KS_xy_alpha(int k, int l,double y1,double y2,double alpha, double beta);
double boundFunct_KS_xy_gamma(int k, int l,double x1,double x2,double gamma, double delta);

double boundFunct_FS_x_alpha(int k,double y1,double y2,double alpha, double beta);
double boundFunct_FS_x_gamma(int k,double x1,double x2,double gamma, double delta);
double boundFunct_FS_y_alpha(int k,double y1,double y2,double alpha, double beta);
double boundFunct_FS_y_gamma(int k,double x1,double x2,double gamma, double delta);

double boundFunct_NS_alpha(int l, double x, double y, double alpha);
double boundFunct_NS_gamma(int l, double x, double y, double gamma);

static double dx;
static double dx2;

struct boundfunc_edge{

	double MS_xx[4][4];
	double MS_yy[4][4];
	double KS_hxxhxx[4][4];
	double KS_hyyhyy[4][4];
	double KS_hxyhxy[4][4];
	double KS_hyxhyx[4][4];
	double KS_hxxhyx[4][4];
	double KS_hxxhyy[4][4];
	double KS_hxyhyx[4][4];
	double KS_hxyhyy[4][4];
	double KS_hxxhxy[4][4];
	double KS_hyyhyx[4][4];
	double FS_hxx[4];
	double FS_hxy[4];
	double FS_hyy[4];
	double FS_hyx[4];
};

struct boundfunc_edgeNS{

	double NS_xc[4];
	double NS_yc[4];
};

struct boundfunc_edge   compute_matricesEdge   (double start_coord[2], double end_coord[2]);
struct boundfunc_edgeNS compute_matricesEdgeNS (double start_coord[2], double end_coord[2]);
int              type_determination     (double x[4]);
boundfunc_cell   compute_matricesCell_bound(int cell_i, double sdf_value[4], char bd_type, int* ptr_bd);
boundfunc_cell   compute_matricesSubcell(double bound_value[4], char bd_type, int* ptr_type);



void set_boundfunc(double ddx){

    int     i,j,subcell_k;
    double temp1[6][6],temp2[6][6];

    dx = ddx;
    dx2 = ddx*ddx;

    /* Set the Cell Elemental Matrices with. Adding dx^2 and denominators  */
    for(i=0;i<6;i++){
        for(j=0;j<6;j++)
        {
            boundfunc_Mcell0.MS_xx[i][j]     *= (dx2/144);
            boundfunc_Mcell0.MS_yy[i][j]     *= (dx2/144);
            boundfunc_Mcell0.KS_hxxhxx[i][j] /= 24;
            boundfunc_Mcell0.KS_hyyhyy[i][j] /= 24;
            boundfunc_Mcell0.KS_hxyhxy[i][j] /= 12;
            boundfunc_Mcell0.KS_hyxhyx[i][j] /= 12;
            boundfunc_Mcell0.KS_hxxhyx[i][j] /= 96;
            boundfunc_Mcell0.KS_hxxhyy[i][j] /= 64;
            boundfunc_Mcell0.KS_hxyhyx[i][j] /= 64;
            boundfunc_Mcell0.KS_hxyhyy[i][j] /= 96;
            boundfunc_Mcell0.KS_hxxhxy[i][j] /= 16;
            boundfunc_Mcell0.KS_hyyhyx[i][j] /= 16;

            boundfunc_subcell0.MS_xx[i][j]     *= (dx2/576);
            boundfunc_subcell0.MS_yy[i][j]     *= (dx2/576);
            boundfunc_subcell0.KS_hxxhxx[i][j] /= 96;
            boundfunc_subcell0.KS_hyyhyy[i][j] /= 96;
            boundfunc_subcell0.KS_hxyhxy[i][j] /= 96;
            boundfunc_subcell0.KS_hyxhyx[i][j] /= 96;
            boundfunc_subcell0.KS_hxxhyx[i][j] /= 96;
            boundfunc_subcell0.KS_hxxhyy[i][j] /= 64;
            boundfunc_subcell0.KS_hxyhyx[i][j] /= 64;
            boundfunc_subcell0.KS_hxyhyy[i][j] /= 96;
            boundfunc_subcell0.KS_hxxhxy[i][j] /= 64;
            boundfunc_subcell0.KS_hyyhyx[i][j] /= 64;

        }

        boundfunc_Mcell0.FS_hxx[i] *= dx/8;
        boundfunc_Mcell0.FS_hxy[i] *= dx/4;
        boundfunc_Mcell0.FS_hyy[i] *= dx/8;
        boundfunc_Mcell0.FS_hyx[i] *= dx/4;

        boundfunc_subcell0.FS_hxx[i] *= dx/16;
        boundfunc_subcell0.FS_hxy[i] *= dx/16;
        boundfunc_subcell0.FS_hyy[i] *= dx/16;
        boundfunc_subcell0.FS_hyx[i] *= dx/16;
    }

    /*Constructing the subcells for subcell 0~3*/

    for(subcell_k=0; subcell_k<4;subcell_k++)
    {
        for(i=0;i<6;i++){
            for(j=0;j<6;j++)
            {

                boudfunc_subcellint[subcell_k].KS_hxxhxx[i][j]
                = boundfunc_subcell0.KS_hxxhxx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                boudfunc_subcellint[subcell_k].KS_hyyhyy[i][j]
                = boundfunc_subcell0.KS_hyyhyy[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];
                boudfunc_subcellint[subcell_k].KS_hxyhxy[i][j]
                = boundfunc_subcell0.KS_hxyhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                boudfunc_subcellint[subcell_k].KS_hyxhyx[i][j]
                = boundfunc_subcell0.KS_hyxhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];

                boudfunc_subcellint[subcell_k].KS_hxxhyx[i][j]
                = boundfunc_subcell0.KS_hxxhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];
                boudfunc_subcellint[subcell_k].KS_hxxhyy[i][j]
                = boundfunc_subcell0.KS_hxxhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                boudfunc_subcellint[subcell_k].KS_hxyhyx[i][j]
                = boundfunc_subcell0.KS_hxyhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;;
                boudfunc_subcellint[subcell_k].KS_hxyhyy[i][j]
                = boundfunc_subcell0.KS_hxyhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];

                boudfunc_subcellint[subcell_k].KS_hxxhxy[i][j]
                = boundfunc_subcell0.KS_hxxhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                boudfunc_subcellint[subcell_k].KS_hyyhyx[i][j]
                = boundfunc_subcell0.KS_hyyhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
            }
            boudfunc_subcellint[subcell_k].FS_hxx[i] = boundfunc_subcell0.FS_hxx[boundfunc.sym[subcell_k].x[i]]* boundfunc.sym[subcell_k].signx ;
            boudfunc_subcellint[subcell_k].FS_hxy[i] = boundfunc_subcell0.FS_hxy[boundfunc.sym[subcell_k].x[i]]* boundfunc.sym[subcell_k].signy ;
            boudfunc_subcellint[subcell_k].FS_hyy[i] = boundfunc_subcell0.FS_hyy[boundfunc.sym[subcell_k].y[i]]* boundfunc.sym[subcell_k].signy ;
            boudfunc_subcellint[subcell_k].FS_hyx[i] = boundfunc_subcell0.FS_hyx[boundfunc.sym[subcell_k].y[i]]* boundfunc.sym[subcell_k].signx ;
        }

        for(i=0;i<6;i++)
            for(j=0;j<6;j++){
                temp1[i][j] = boudfunc_subcellint[subcell_k].KS_hxxhxy[j][i];
                temp2[i][j] = boudfunc_subcellint[subcell_k].KS_hyyhyx[j][i];
            }

         for(i=0;i<6;i++)
            for(j=0;j<6;j++){
                boudfunc_subcellint[subcell_k].KS_hxxhxy[i][j] += temp1[i][j];
                boudfunc_subcellint[subcell_k].KS_hyyhyx[i][j] += temp2[i][j];
            }
        }
}


/*first time to calculate governing matrices for reducing system: assuming undeformed state*/
/* Kl, M, DS, NS are constructed */
void compute_matricesNonlinearStructure(Matrices_S* ptr_ms, Index_S* ptr_i, Grid_S* ptr_g, Solid* ptr_s, char* fnd){

    double          m2pl = 2*ptr_s->mu + ptr_s->lambda;
    double          sdf[4];
    int             nx   = ptr_g->Nx;
    int             i,j,k,cell_i,nsubcell,ii,jj,bdcell = 0;
    int             index_xix[6], index_xiy[6], index_stagx[2],index_stagy[2];
    struct invol    index_cell;
    boundfunc_cell  Mcell;

    /*initialize the matrices @@ */
    for(i=0;i<ptr_i->xix_N;i++)
        for(j=0;j<ptr_i->xix_N;j++)
        {
            ptr_ms->MS_xx [i][j] = 0;
            ptr_ms->KlS_xx[i][j] = 0;
        }

    for(i=0;i<ptr_i->xiy_N;i++)
        for(j=0;j<ptr_i->xiy_N;j++)
        {
            ptr_ms->KlS_yy[i][j] = 0;
            ptr_ms->MS_yy[i][j] = 0;
        }

    for(i=0;i<ptr_i->xix_N;i++)
        for(j=0;j<ptr_i->xiy_N;j++)
            ptr_ms->KlS_xy[i][j] = 0;

    for(i=0;i<ptr_i->xix_N;i++)
    {
        for(j=0;j<ptr_i->xix_Ncell_Neumann;j++)
            ptr_ms->NS_xc[i][j]=0;
        for(j=0;j<ptr_i->xix_Ncell_Dirichlet;j++)
            ptr_ms->DS_xc[i][j]=0;
    }
    for(i=0;i<ptr_i->xiy_N;i++)
    {
        for(j=0;j<ptr_i->xiy_Ncell_Neumann;j++)
            ptr_ms->NS_yc[i][j]=0;
        for(j=0;j<ptr_i->xiy_Ncell_Dirichlet;j++)
            ptr_ms->DS_yc[i][j]=0;
    }

    /* =========  for interior cells =========== @@ */
    for(i=0;i<ptr_i->cell_N_interior;i++)
    {
       cell_i = (int) ptr_i->cell_interior[i];
       index_cell = involvedIndices_grid(cell_i, nx);

       for(j=0;j<6;j++)
       {
            index_xix[j] = ptr_i->xix.G2g[index_cell.xix[j]];
            index_xiy[j] = ptr_i->xiy.G2g[index_cell.xiy[j]];
       }

       for(j=0;j<6;j++)
        for(k=0;k<6;k++)
        {
            ptr_ms->MS_xx[index_xix[j]][index_xix[k]] += boundfunc_Mcell0.MS_xx[j][k];
            ptr_ms->MS_xx[index_xix[j]][index_xix[k]] += boundfunc_Mcell0.MS_xx[j][k];

            ptr_ms->KlS_xx[index_xix[j]][index_xix[k]]+= m2pl*boundfunc_Mcell0.KS_hxxhxx[j][k]
                                                         + ptr_s->mu*boundfunc_Mcell0.KS_hxyhxy[j][k];
            ptr_ms->KlS_xy[index_xix[j]][index_xiy[k]]+= ptr_s->lambda*boundfunc_Mcell0.KS_hxxhyy[j][k]
                                                         + ptr_s->mu*boundfunc_Mcell0.KS_hxyhyx[j][k];
            ptr_ms->KlS_yy[index_xiy[j]][index_xiy[k]]+= m2pl*boundfunc_Mcell0.KS_hyyhyy[j][k]
                                                         + ptr_s->mu*boundfunc_Mcell0.KS_hyxhyx[j][k];
        }
    }

    /*Determine how many subcells are intersecting with the boundaries @@*/
    nsubcell = 0;
    shapefunc.info = imatrix(0,ptr_g->N-1,0,3);

    for(i=0;i<ptr_i->cell_N_boundary ;i++)
    {
        cell_i = (int) ptr_i->cell_boundary[i];
        ii = cell_i%nx;
        jj = (int) floor((double)cell_i/nx);

        sdf[0] = ptr_s->boundary_value[ii][jj];
        sdf[1] = ptr_s->boundary_value[ii+1][jj];
        sdf[2] = ptr_s->boundary_value[ii][jj+1];
        sdf[3] = ptr_s->boundary_value[ii+1][jj+1];

        nsubcell+= type_determination(sdf);
    }

    nsubcell -= 1;
    shapefunc.cutcell=(boundfunc_cutcell*) malloc(nsubcell*sizeof(boundfunc_cutcell));

    /* =========  for boundary cells =========== @@ */
    for(i=0;i<ptr_i->cell_N_boundary; i++)
    {
       cell_i = (int) ptr_i->cell_boundary[i];
       ii = cell_i%nx;
       jj = (int) floor((double)cell_i/nx);

       index_cell = involvedIndices_grid(cell_i, nx);

       for(j=0;j<6;j++)
       {
            index_xix[j] = ptr_i->xix.G2g[index_cell.xix[j]];
            index_xiy[j] = ptr_i->xiy.G2g[index_cell.xiy[j]];
        }

       sdf[0] = ptr_s->boundary_value[ii][jj];
       sdf[1] = ptr_s->boundary_value[ii+1][jj];
       sdf[2] = ptr_s->boundary_value[ii][jj+1];
       sdf[3] = ptr_s->boundary_value[ii+1][jj+1];

       Mcell = compute_matricesCell_bound(cell_i,sdf, fnd[cell_i] ,&bdcell);

       for(j=0;j<6;j++)
        for(k=0;k<6;k++)
        {
            ptr_ms->MS_xx[index_xix[j]][index_xix[k]] += boundfunc_Mcell0.MS_xx[j][k];
            ptr_ms->MS_xx[index_xix[j]][index_xix[k]] += boundfunc_Mcell0.MS_xx[j][k];

            ptr_ms->KlS_xx[index_xix[j]][index_xix[k]]+= m2pl*boundfunc_Mcell0.KS_hxxhxx[j][k]
                                                         + ptr_s->mu*boundfunc_Mcell0.KS_hxyhxy[j][k];
            ptr_ms->KlS_xy[index_xix[j]][index_xiy[k]]+= ptr_s->lambda*boundfunc_Mcell0.KS_hxxhyy[j][k]
                                                         + ptr_s->mu*boundfunc_Mcell0.KS_hxyhyx[j][k];
            ptr_ms->KlS_yy[index_xiy[j]][index_xiy[k]]+= m2pl*boundfunc_Mcell0.KS_hyyhyy[j][k]
                                                         + ptr_s->mu*boundfunc_Mcell0.KS_hyxhyx[j][k];

        }

       if (fnd[cell_i] == 'D')
       {
           for(j=0;j<2;j++)
            {
                index_stagx[j] = ptr_i->xix.C2c_dirichlet[index_cell.stagx_cell[j]];
                index_stagy[j] = ptr_i->xiy.C2c_dirichlet[index_cell.stagy_cell[j]];
            }
            for(j=0;j<6;j++)
                for(k=0;k<2;k++)
                {
                    ptr_ms->DS_xc[index_xix[j]][index_stagx[k]] +=  Mcell.DS_xc[j][k];
                    ptr_ms->DS_yc[index_xiy[j]][index_stagy[k]] +=  Mcell.DS_yc[j][k];
                }
       }
       else if (fnd[cell_i] == 'N')
       {
            for(j=0;j<2;j++)
            {
                index_stagx[j] = ptr_i->xix.C2c_neumann[index_cell.stagx_cell[j]];
                index_stagy[j] = ptr_i->xiy.C2c_neumann[index_cell.stagy_cell[j]];
            }
           for(j=0;j<6;j++)
                for(k=0;k<2;k++)
                {
                    ptr_ms->NS_xc[index_xix[j]][index_stagx[k]] +=  Mcell.NS_xc[j][k];
                    ptr_ms->NS_yc[index_xiy[j]][index_stagy[k]] +=  Mcell.NS_yc[j][k];
                }
       }
    }
}

// input the sdf of four corners and return how many subcells is intersecting with boundaries
int type_determination(double x[4])
{
    int         number=0,i;
    double      sdf[9] = {x[0],0,x[1],0,0,0,x[2],0,x[3]};
    int         y[4][4] = {{1,2,4,5},{2,3,5,6},{4,5,7,8},{5,6,8,9}},ssdf[9];
    _Bool       isext,isint;

    sdf[2] = (sdf[1]+sdf[3])/2;
    sdf[4] = (sdf[1]+sdf[7])/2;
    sdf[6] = (sdf[3]+sdf[9])/2;
    sdf[8] = (sdf[7]+sdf[9])/2;
    sdf[5] = (sdf[1]+sdf[3]+sdf[7]+sdf[9])/4;

    for(i=0;i<9;i++)
        ssdf[i] = (int) SIGN(1,sdf[i]);

    for(i=0;i<4;i++)
    {
        isext = ((ssdf[y[i][0]]!= -1) && (ssdf[y[i][1]]!= -1) && (ssdf[y[i][2]]!= -1) && (ssdf[y[i][3]]!= -1));
        isint = ((ssdf[y[i][0]]+ssdf[y[i][1]]+ssdf[y[i][2]]+ssdf[y[i][3]])<=-3);
        if (!isext && !isint)
            number +=1;
    }
    return number;
}

boundfunc_cell compute_matricesCell_bound(int cell_i, double sdf_value[4], char bd_type, int* ptr_bd)
{
    boundfunc_cell  Mcell={0};
    int             subcell_k,i,j;
    int             sym[4][4] = {{1,2,5,4},{2,3,6,5},{4,5,8,7},{5,6,9,8}};
    double          sdf[9] = {sdf_value[0],0,sdf_value[1],0,0,0,sdf_value[2],0,sdf_value[3]};
    double          bound_value[4];
    boundfunc_cell  subcell;
    int             type,n;

    sdf[2] = (sdf[1]+sdf[3])/2;
    sdf[4] = (sdf[1]+sdf[7])/2;
    sdf[6] = (sdf[3]+sdf[9])/2;
    sdf[8] = (sdf[7]+sdf[9])/2;
    sdf[5] = (sdf[1]+sdf[3]+sdf[7]+sdf[9])/4;

    for(subcell_k=0;subcell_k<4; subcell_k++ )
    {
        for(i=0; i<4; i++)
            bound_value[i]=sdf[sym[subcell_k][i]];

        subcell=compute_matricesSubcell(bound_value,bd_type,&type);

        for(i=0; i<6; i++)
            for(j=0; j<6; j++)
            {
                Mcell.MS_xx[i][j]+=subcell.MS_xx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                Mcell.MS_yy[i][j]+=subcell.MS_yy[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];

                Mcell.KS_hxxhxx[i][j] +=subcell.KS_hxxhxx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                Mcell.KS_hxyhxy[i][j] +=subcell.KS_hxyhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];

                Mcell.KS_hxxhyy[i][j] +=subcell.KS_hxxhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                Mcell.KS_hxyhyx[i][j] +=subcell.KS_hxyhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;;

                Mcell.KS_hyyhyy[i][j] +=subcell.KS_hyyhyy[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];
                Mcell.KS_hyxhyx[i][j] +=subcell.KS_hyxhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];
            }

        if(bd_type == 'N')
        {
            for(i=0; i<6; i++)
            {
                Mcell.NS_xc[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].cellx] += subcell.NS_xc[i][0];
                Mcell.NS_yc[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].celly] += subcell.NS_yc[i][0];
            }

        }
        else if (bd_type == 'D')
        {
             for(i=0; i<6; i++)
            {
                Mcell.DS_xc[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].cellx] += subcell.DS_xc[i][0];
                Mcell.DS_yc[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].celly] += subcell.DS_yc[i][0];
            }
        }

        if (type == 1)
        {
            n = *ptr_bd;
            for(i=0; i<6; i++)
                for(j=0; j<6; j++)
                {
                    shapefunc.cutcell[n].KS_hxxhxx[i][j]=subcell.KS_hxxhxx[i][j];
                    shapefunc.cutcell[n].KS_hyyhyy[i][j]=subcell.KS_hyyhyy[i][j];
                    shapefunc.cutcell[n].KS_hxyhxy[i][j]=subcell.KS_hxyhxy[i][j];
                    shapefunc.cutcell[n].KS_hyxhyx[i][j]=subcell.KS_hyxhyx[i][j];
                    shapefunc.cutcell[n].KS_hxxhyx[i][j]=subcell.KS_hxxhyx[i][j];
                    shapefunc.cutcell[n].KS_hxxhyy[i][j]=subcell.KS_hxxhyy[i][j];
                    shapefunc.cutcell[n].KS_hxyhyx[i][j]=subcell.KS_hxyhyx[i][j];
                    shapefunc.cutcell[n].KS_hxyhyy[i][j]=subcell.KS_hxyhyy[i][j];
                    shapefunc.cutcell[n].KS_hxxhxy[i][j]=subcell.KS_hxxhxy[i][j];
                    shapefunc.cutcell[n].KS_hyyhyx[i][j]=subcell.KS_hyyhyx[i][j];
                }

            for(i=0; i<6; i++)
            {
                shapefunc.cutcell[n].FS_hxx[i] = subcell.FS_hxx[i];
                shapefunc.cutcell[n].FS_hxx[i] = subcell.FS_hxx[i];
            }
            shapefunc.info[cell_i][subcell_k]=type;
            *ptr_bd += 1;
        }
        else
            shapefunc.info[cell_i][subcell_k]=type;
    }
    return Mcell;
}


boundfunc_cell compute_matricesSubcell(double bound_value[4], char bd_type, int* ptr_type)
{
    /*% We define the variable CORNERS containing the x-coordinates,
    % y-coordinates and signed-function value for each of the corners.
    % The coordinates are with respect to the left-lower corner origin, and
    %normalized by the cell size dim.h*/
    int     i,j;
    int     bound_sign[4],side,side_next,ncorner1 =4,ncorner2;
    int     ind1[4] = {1,2,3,4}, ind2[4]={1,2,4,5};
    double  corners[4][8] = {{0,0.5,0.5,0},{0,0,0.5,0.5},{bound_value[0],bound_value[1],bound_value[2],bound_value[3]},{0,1,2,3}};
    double  new_corner[2],start_coord[2],end_coord[2];
    struct boundfunc_edge      edge;
    struct boundfunc_edgeNS    edgeNS;
    boundfunc_cell             subcell={0};

    for(i=0;i<4;i++)
        bound_sign[i]=SIGN(1,bound_value[i]);

    // Compute the boundary corners
    for(side=3;side>-1;side--)
    {
        side_next = (side+1)%4;
        if( corners[2][side]*corners[2][side_next]<0)
        {
            ncorner1 +=1;

            for(i=0;i<2;i++)
                new_corner[i] = corners[i][side]+(corners[i][side_next]-corners[i][side])
                            *corners[2][side]/(corners[2][side]-corners[2][side_next]);

            for(i=0; i<4;i++)
                for(j=ncorner1-1; j>side+1;j--)
                    corners[i][j] = corners[i][j-1];

            corners[0][side+1] = new_corner[0];
            corners[1][side+1] = new_corner[1];
            corners[2][side+1] = 0;
            corners[3][side+1] = side;
        }
    }
    ncorner2 = ncorner1;

    // Delete the external corners
    for(side=ncorner1-1;side>-1;side--)
    {
        if (corners[2][side] >0)
        {
            ncorner2 -= 1;
            for(i=side; i<ncorner2;i++)
                for(j=0; j<4; j++)
                    corners[i][j] = corners[i][j+1];
        }
    }

    /*Compute for boundary subcells contribution, edge-by-edge*/
    // Add the contribution of each edge into the line integral
    for(side=0; side<ncorner2;side++)
    {
        side_next = (side+1)%4;
        start_coord[0] = corners[0][side];
        start_coord[1] = corners[1][side];
        end_coord[0]   = corners[0][side_next];
        end_coord[1]   = corners[1][side_next];

        if ( corners[1][side] != corners[1][side_next])
        {
            edge = compute_matricesEdge( start_coord, end_coord);

            for(i=0; i<4; i++)
            {
                for(j=0; j<4; j++)
                {
                    subcell.MS_xx[ind1[i]][ind1[j]]     += edge.MS_xx[i][j];
                    subcell.MS_yy[ind2[i]][ind2[j]]     += edge.MS_yy[i][j];

                    subcell.KS_hxxhxx[ind1[i]][ind1[j]] += edge.KS_hxxhxx[i][j];
                    subcell.KS_hyyhyy[ind2[i]][ind2[j]] += edge.KS_hyyhyy[i][j];
                    subcell.KS_hxyhxy[ind1[i]][ind1[j]] += edge.KS_hxyhxy[i][j];
                    subcell.KS_hyxhyx[ind2[i]][ind2[j]] += edge.KS_hyxhyx[i][j];

                    subcell.KS_hxxhyx[ind1[i]][ind2[j]] += edge.KS_hxxhyx[i][j];
                    subcell.KS_hxxhyy[ind1[i]][ind2[j]] += edge.KS_hxxhyy[i][j];
                    subcell.KS_hxyhyx[ind1[i]][ind2[j]] += edge.KS_hxyhyx[i][j];
                    subcell.KS_hxyhyy[ind1[i]][ind2[j]] += edge.KS_hxyhyy[i][j];

                    subcell.KS_hxxhxy[ind1[i]][ind1[j]] += edge.KS_hxxhxy[i][j];
                    subcell.KS_hyyhyx[ind2[i]][ind2[j]] += edge.KS_hyyhyx[i][j];
                }
                subcell.FS_hxx[ind1[i]] += edge.FS_hxx[i];
                subcell.FS_hxy[ind1[i]] += edge.FS_hxy[i];
                subcell.FS_hyy[ind2[i]] += edge.FS_hyy[i];
                subcell.FS_hyx[ind2[i]] += edge.FS_hyx[i];
            }
        }

        if ((corners[2][side]==0) && (corners[2][side_next]==0))
        {
            edgeNS = compute_matricesEdgeNS( start_coord, end_coord);

            /*if (bd_type == 'F')
            {
                subcell.coord = [ corners(1:2,side) , corners(1:2,side_next) ];
                subcell.TS0.xc( : , end+1 ) = - edge.NS.xc;
                subcell.TS0.yc( : , end+1 ) = - edge.NS.yc;
            }*/
            if (bd_type == 'N')
            {
                for(i=0; i<4; i++)
                {
                    subcell.NS_xc[ind1[i]][0] += edgeNS.NS_xc[i];
                    subcell.NS_yc[ind2[i]][0] += edgeNS.NS_yc[i];
                }
            }
            else if (bd_type == 'D')
            {
                for(i=0; i<4; i++)
                {
                    subcell.DS_xc[ind1[i]][0] -= edgeNS.NS_xc[i];
                    subcell.DS_yc[ind2[i]][0] -= edgeNS.NS_yc[i];
                }
            }
        }
    }
    return subcell;
}

struct boundfunc_edge compute_matricesEdge(double start_coord[2], double end_coord[2])
{
    struct boundfunc_edge  edge={0};
    double          alpha, beta,delta;
    int             i,j;

    alpha    = ( end_coord[0]-start_coord[0] ) / ( end_coord[1]-start_coord[1]);

    if (abs(alpha) <= 1)
        beta  = start_coord[0] - alpha*start_coord[1];
    else
        delta = start_coord[1] - start_coord[0]/alpha;

    //Calculation using alpha-formula
    if ( abs(alpha) <= 1)
    {
        for(i=0;i<4;i++)
        {
           for (j=0;j<4;j++)
          {
            //MS.xx coefficients
            edge.MS_xx[i][j] = boundFunct_MS_alpha(  i,  j,  end_coord[0],  end_coord[1], alpha)
                             - boundFunct_MS_alpha(  i,  j,start_coord[0],start_coord[1], alpha);
            //MS.yy coefficients
            edge.MS_yy[i][j] = boundFunct_MS_alpha(4+i,4+j,  end_coord[0],  end_coord[1], alpha)
                             - boundFunct_MS_alpha(4+i,4+j,start_coord[0],start_coord[1], alpha);

            //KS coefficients
            edge.KS_hxxhxx[i][j] = boundFunct_KS_xx_alpha(  i,  j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hyyhyy[i][j] = boundFunct_KS_yy_alpha(4+i,4+j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxyhxy[i][j] = boundFunct_KS_yy_alpha(  i,  j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hyxhyx[i][j] = boundFunct_KS_xx_alpha(4+i,4+j,  start_coord[1],  end_coord[1], alpha, beta);

            edge.KS_hxxhyx[i][j] = boundFunct_KS_xx_alpha(  i,4+j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxxhyy[i][j] = boundFunct_KS_xy_alpha(j+4,  i,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxyhyx[i][j] = boundFunct_KS_xy_alpha(  i,j+4,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxyhyy[i][j] = boundFunct_KS_yy_alpha( j+4, i,  start_coord[1],  end_coord[1], alpha, beta);

            edge.KS_hxxhxy[i][j] = boundFunct_KS_xy_alpha(   i,  j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hyyhyx[i][j] = boundFunct_KS_xy_alpha( j+4,i+4,  start_coord[1],  end_coord[1], alpha, beta);
          }

         edge.FS_hxx[i] = boundFunct_FS_x_alpha(  i,   start_coord[1],  end_coord[1], alpha, beta);
         edge.FS_hxy[i] = boundFunct_FS_y_alpha(  i,   start_coord[1],  end_coord[1], alpha, beta);
         edge.FS_hyy[i] = boundFunct_FS_y_alpha(4+i,   start_coord[1],  end_coord[1], alpha ,beta);
         edge.FS_hyx[i] = boundFunct_FS_x_alpha(4+i,   start_coord[1],  end_coord[1], alpha ,beta);
        }
    }

    //Calculation using gamma-formula
    else
    {
        for(i=0;i<4;i++)
        {
           for (j=0;j<4;j++)
          {
            //MS.xx coefficients
            edge.MS_xx[i][j] = boundFunct_MS_gamma(  i,  j,  end_coord[0],  end_coord[1], 1/alpha)
                             - boundFunct_MS_gamma(  i,  j,start_coord[0],start_coord[1], 1/alpha);
            //MS.yy coefficients
            edge.MS_yy[i][j] = boundFunct_MS_gamma(4+i,4+j,  end_coord[0],  end_coord[1], 1/alpha)
                             - boundFunct_MS_gamma(4+i,4+j,start_coord[0],start_coord[1], 1/alpha);
            //KS coefficients
            edge.KS_hxxhxx[i][j] = boundFunct_KS_xx_gamma(  i,  j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hyyhyy[i][j] = boundFunct_KS_yy_gamma(4+i,4+j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxyhxy[i][j] = boundFunct_KS_yy_gamma(  i,  j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hyxhyx[i][j] = boundFunct_KS_xx_gamma(4+i,4+j,  start_coord[0],  end_coord[0], 1/alpha, delta);

            edge.KS_hxxhyx[i][j] = boundFunct_KS_xx_gamma(  i,4+j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxxhyy[i][j] = boundFunct_KS_xy_gamma(j+4,  i,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxyhyx[i][j] = boundFunct_KS_xy_gamma(  i,j+4,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxyhyy[i][j] = boundFunct_KS_yy_gamma( j+4, i,  start_coord[0],  end_coord[0], 1/alpha, delta);

            edge.KS_hxxhxy[i][j] = boundFunct_KS_xy_gamma(   i,  j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hyyhyx[i][j] = boundFunct_KS_xy_gamma( j+4,i+4,  start_coord[0],  end_coord[0], 1/alpha, delta);

          }
         edge.FS_hxx[i] = boundFunct_FS_x_gamma(  i,   start_coord[0],  end_coord[0], 1/alpha, delta);
         edge.FS_hxy[i] = boundFunct_FS_y_gamma(  i,   start_coord[0],  end_coord[0], 1/alpha, delta);
         edge.FS_hyy[i] = boundFunct_FS_y_gamma(4+i,   start_coord[0],  end_coord[0], 1/alpha ,delta);
         edge.FS_hyx[i] = boundFunct_FS_x_gamma(4+i,   start_coord[0],  end_coord[0], 1/alpha ,delta);
        }
    }
    return edge;
}

struct boundfunc_edgeNS compute_matricesEdgeNS(double start_coord[2], double end_coord[2])
{

    struct boundfunc_edgeNS  edge={0};
    double            alpha, beta,delta;
    int               i;

    alpha    = ( end_coord[0]-start_coord[0] ) / ( end_coord[1]-start_coord[1]);

    if (abs(alpha) <= 1)
        beta  = start_coord[0] - alpha*start_coord[1];
    else
        delta = start_coord[1] - start_coord[0]/alpha;

    //Calculation using alpha-formula
    if ( abs(alpha) <= 1)
    {
        for(i=0;i<4;i++)
        {
         //Calculation of the NS coefficients
         edge.NS_xc[i] = SIGN( 1, end_coord[1]-start_coord[1] ) *
                           ( boundFunct_NS_alpha(   i,  end_coord[0],  end_coord[1], alpha)
                            -boundFunct_NS_alpha(   i,start_coord[0],start_coord[1], alpha) );
         edge.NS_yc[i] = SIGN( 1, end_coord[1]-start_coord[1] ) *
                           ( boundFunct_NS_alpha( 4+i,  end_coord[0],  end_coord[1], alpha)
                            -boundFunct_NS_alpha( 4+i,start_coord[0],start_coord[1], alpha) );
        }
    }
    else
    {
        for(i=0;i<4;i++)
        {
         //Calculation of the NS coefficients
         edge.NS_xc[i] = SIGN( 1,end_coord[0]-start_coord[0] ) *
                           ( boundFunct_NS_gamma(   i,  end_coord[0],  end_coord[1], 1/alpha)
                            -boundFunct_NS_gamma(   i,start_coord[0],start_coord[1], 1/alpha) );

         edge.NS_yc[i] = SIGN( 1,end_coord[0]-start_coord[0] ) *
                           ( boundFunct_NS_gamma( 4+i,  end_coord[0],  end_coord[1], 1/alpha)
                            -boundFunct_NS_gamma( 4+i,start_coord[0],start_coord[1], 1/alpha) );
        }
    }

    return edge;
}

double boundFunct_MS_alpha(int k, int l,double x,double y,double alpha){

     double z;
     double y3 = y*y*y, y4 = y3*y, y5= y4*y, y6 = y5*y;

     z = dx2* boundfunc.coeff[0][k]*boundfunc.coeff[0][l]*(
         (2*(x*x*x)/6       + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x*x)/2  + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*x)
         * (2*(y3)/6        + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y*y)/2  + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*y)
         + (-alpha)*(2*(x*x)/2 + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])* x    + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]  )
         * (2*(y4)/24 + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y3)/6         + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*(y*y)/2 )
         + alpha*alpha*(2* x                                                         + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])                                                        )
         * (2*(y5)/120+ (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y4)/24        + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*(y3)/6 )
         + alpha*alpha*alpha*2
         * (2*(y6)/720+ (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y5)/120       + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*(y4)/24) );

    return z;
}

double boundFunct_MS_gamma(int k, int l,double x,double y,double gamma){

     double z;
     double x3 = x*x*x, x4 = x3*x, x5= x4*x, x6 = x5*x;

     z = dx2 * boundfunc.coeff[0][k] * boundfunc.coeff[0][l]* gamma
         *(   (2*(x4)/24 + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x3)/6  + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*(x*x)/2 )
            * (2*(y*y)/2 + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*y + boundfunc.coeff[2][k]*boundfunc.coeff[2][l])
            + (-gamma)  *(2*(x5)/120+ (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x4)/24 + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*(x3)/6 )
            * (2* y      + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])                                                )
            + gamma*gamma*(2*(x6)/720+ (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x5)/120+ boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*(x4)/24)
            * 2);
     return z;
}

double boundFunct_KS_xx_alpha(int k, int l,double y1,double y2,double alpha, double beta){

    double z;
    double y13 = y1*y1*y1, y14 = y13*y1, y23= y2*y2*y2, y24 = y23*y2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]
        *( alpha *(y24-y14)/4
         + (beta + alpha* (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])) *(y23-y13)/3
         + ( boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*alpha
         + beta*( boundfunc.coeff[2][k] + boundfunc.coeff[2][l]))        *(y2*y2-y1*y1)/2
         + beta*boundfunc.coeff[2][k] * boundfunc.coeff[2][l]            *(y2-y1) );

    return z;
}

double boundFunct_KS_xx_gamma(int k, int l,double x1,double x2,double gamma, double delta){

    double z;
    double x13 = x1*x1*x1, x14 = x13*x1, x23= x2*x2*x2, x24 = x23*x2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]
        *( gamma*gamma*gamma                                                         *(x24-x14)/4
           + gamma*gamma*( boundfunc.coeff[2][k] + boundfunc.coeff[2][l] + 2*delta ) *(x23-x13)/3
           + gamma*(delta + boundfunc.coeff[2][k] )* (delta + boundfunc.coeff[2][l]) *(x2*x2-x1*x1)/2 );

    return z;
}

double boundFunct_KS_yy_alpha(int k, int l,double y1,double y2,double alpha, double beta){

    double z;
    double y13 = y1*y1*y1, y14 = y13*y1, y23= y2*y2*y2, y24 = y23*y2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]
        *(   alpha*alpha*alpha                                                         *(y24-y14)/12
           + alpha*alpha*  ( beta + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])/2 )  *(y23-y13)/3
           + alpha*(beta + boundfunc.coeff[1][k])*(beta + boundfunc.coeff[1][l])       *(y2*y2-y1*y1)/2
           + (beta*beta*beta/3 + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])/2*beta*beta + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*beta ) *(y2-y1) );

    return z;
}

double boundFunct_KS_yy_gamma(int k, int l,double x1,double x2,double gamma, double delta){

    double z;
    double x13 = x1*x1*x1, x14 = x13*x1, x23= x2*x2*x2, x24 = x23*x2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]* gamma
        *( (x24-x14)/12
            + (boundfunc.coeff[1][k] + boundfunc.coeff[1][l]) *(x23-x13)/6
            +  boundfunc.coeff[1][k]*boundfunc.coeff[1][l]    *(x2*x2-x1*x1)/2 );

    return z;
}

double boundFunct_KS_xy_alpha(int k, int l,double y1,double y2,double alpha, double beta){

    double z;
    double y13 = y1*y1*y1, y14 = y13*y1, y23= y2*y2*y2, y24 = y23*y2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]/2
        *( alpha*alpha                                                                 *(y24-y14)/4
           +  alpha* (alpha*boundfunc.coeff[2][l] +2*(beta + boundfunc.coeff[1][k]) )  *(y23-y13)/3
           + (2*boundfunc.coeff[2][l]*alpha+ beta + boundfunc.coeff[1][k])*(beta + boundfunc.coeff[1][k])      *(y2*y2-y1*y1)/2
           + boundfunc.coeff[2][l]*(beta + boundfunc.coeff[1][k])*(beta + boundfunc.coeff[1][k])               *(y2-y1)  );

    return z;
}

double boundFunct_KS_xy_gamma(int k, int l,double x1,double x2,double gamma, double delta){

    double z;
    double x13 = x1*x1*x1, x14 = x13*x1, x23= x2*x2*x2, x24 = x23*x2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]/2
        *( gamma*gamma                                                                *(x24-x14)/4
           + gamma* (2*gamma*boundfunc.coeff[1][k] + delta + boundfunc.coeff[2][l] )  *(x23-x13)/3
           + (2*boundfunc.coeff[1][k]*gamma*( delta +boundfunc.coeff[2][l]) + boundfunc.coeff[1][k]*boundfunc.coeff[1][k]*gamma*gamma)    *(x2*x2-x1*x1)/2
           + boundfunc.coeff[1][k]*boundfunc.coeff[1][k]*gamma*(delta + boundfunc.coeff[2][l])                                            *(x2-x1) );

    return z;
}

double boundFunct_FS_x_alpha(int k,double y1,double y2,double alpha, double beta){

    double z;
    double y12 = y1*y1, y13 = y12*y1, y22 = y2*y2, y23= y22*y2;

    z = dx * boundfunc.coeff[0][k]
        *( alpha                                    *(y23-y13)/3
           + (alpha*boundfunc.coeff[2][k] + beta )  *(y22-y12)/2
           + beta*boundfunc.coeff[2][k]             *(y2-y1) );

    return z;
}

double boundFunct_FS_x_gamma(int k,double x1,double x2,double gamma, double delta){

    double z;
    double x12 = x1*x1, x13 = x12*x1, x22= x2*x2, x23 = x22*x2;

    z = dx * boundfunc.coeff[0][k]
        *( gamma*gamma                          *(x23-x13)/3
        + gamma*(delta + boundfunc.coeff[2][k]) *(x22-x12)/2 );

    return z;
}

double boundFunct_FS_y_alpha(int k,double y1,double y2,double alpha, double beta){

    double z;
    double y12 = y1*y1, y13 = y12*y1, y22 = y2*y2, y23= y22*y2;

    z = dx * boundfunc.coeff[0][k]/2
        *( alpha*alpha                              *(y23-y13)/3
           + alpha*(boundfunc.coeff[1][k] + beta )  *(y22-y12)
           + (boundfunc.coeff[1][k] + beta )*(boundfunc.coeff[1][k] + beta )      *(y2-y1) );

    return z;
}

double boundFunct_FS_y_gamma(int k,double x1,double x2,double gamma, double delta){

    double z;

    z = dx * boundfunc.coeff[0][k]* gamma
        *( ( (x2+boundfunc.coeff[1][k])*(x2+boundfunc.coeff[1][k])*(x2+boundfunc.coeff[1][k])
            -(x1+boundfunc.coeff[1][k])*(x1+boundfunc.coeff[1][k])*(x1+boundfunc.coeff[1][k]) )/6 );

    return z;
}

double boundFunct_NS_alpha(int l, double x, double y, double alpha){

    double z;
    z = dx * boundfunc.coeff[0][l] * sqrt(alpha*alpha+1)
        *(   (x + boundfunc.coeff[1][l] )  *  ( (y*y)/2 + boundfunc.coeff[2][l]* y)
             + (-alpha)  *( (y*y*y)/6 + boundfunc.coeff[2][l]*(y*y)/2 ) );

    return z;
}

double boundFunct_NS_gamma(int l, double x, double y, double gamma){

    double z;
    z = dx * boundfunc.coeff[0][l] * sqrt(gamma*gamma+1)
        *(  ( y + boundfunc.coeff[2][l] )  *  ( (x*x)/2 + boundfunc.coeff[1][l]* x )
              + (-gamma)  *( (x*x*x)/6 + boundfunc.coeff[1][l]*(x*x)/2 ) );
    return z;
}
