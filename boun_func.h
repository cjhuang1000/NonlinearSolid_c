#ifndef BOUN_FUNC_H
#define BOUN_FUNC_H 1

#include "Field_s.h"

typedef struct {

	double KS_hxxhxx[6][6];
	double KS_hyyhyy[6][6];
	double KS_hxyhxy[6][6];
	double KS_hyxhyx[6][6];
	double KS_hxxhyx[6][6];
	double KS_hxxhyy[6][6];
	double KS_hxyhyx[6][6];
	double KS_hxyhyy[6][6];
	double KS_hxxhxy[6][6];
	double KS_hyyhyx[6][6];
	double FS_hxx[6];
	double FS_hxy[6];
	double FS_hyy[6];
	double FS_hyx[6];

} boundfunc_cutcell;

typedef struct {

	/* % Internal Cell Elemental Matrices*/
	// integration within the internal grid
	double MS_xx[6][6];
	double MS_yy[6][6];
	double KS_hxxhxx[6][6];
	double KS_hyyhyy[6][6];
	double KS_hxyhxy[6][6];
	double KS_hyxhyx[6][6];
	double KS_hxxhyx[6][6];
	double KS_hxxhyy[6][6];
	double KS_hxyhyx[6][6];
	double KS_hxyhyy[6][6];
	double KS_hxxhxy[6][6];
	double KS_hyyhyx[6][6];
	double FS_hxx[6];
	double FS_hxy[6];
	double FS_hyy[6];
	double FS_hyx[6];
	double NS_xc[6][2];
	double NS_yc[6][2];
	double DS_xc[6][2];
	double DS_yc[6][2];

} boundfunc_cell;

typedef struct {

	/*Symetry consideration*/
	short int cellx;
	short int celly;
	short int signx;
	short int signy;
	short int x[6];
	short int y[6];
	short int subcell[4];

} bound_func_sym;

typedef struct {

	/* Elementary functions coefficients */
	/* % Coefficients which characterize the bilinear elementary functions N^xi
	 involved on the subcell omega1 (N_x1, N_x2, N_x3, N_x4, N_y7, N_y8,
	 N_y10, N_y11), such that, for (x,y) in omega1:
                       N(x,y) = eps*(x+a)*(y+b)
	The array below contains the [eps; a; b] coefficients for each elementary
	 functions N*/
	double coeff[3][8];
    bound_func_sym sym[4];

}boundary_function;

struct shapefunction{

    int               **info;
    boundfunc_cutcell *cutcell;

};
/*Internal Cell Elemental Matrices*/
static boundfunc_cell boundfunc_Mcell0
= {
    {{2,1,4,2,0,0},{0,2,2,4,0,0},{0,0,28,14,4,2},{0,0,0,28,2,4},{0,0,0,0,2,1},{0,0,0,0,0,2}},
    {{2,4,0,1,2,0},{0,28,4,2,14,2},{0,0,2,0,2,1},{0,0,0,2,4,0},{0,0,0,0,28,4},{0,0,0,0,0,2}},
    {{1,-1,2,-2,0,0},{-1,1,-2,2,0,0},{2,-2,14,-14,2,-2},{-2,2,-14,14,-2,2},{0,0,2,-2,1,-1},{0,0,-2,2,-1,1}},
    {{1,2,0,-1,-2,0},{2,14,2,-2,-14,-2},{0,2,1,0,-2,-1},{-1,-2,0,1,2,0},{-2,-14,-2,2,14,2},{0,-2,-1,0,2,1}},
    {{2,1,-2,-1,0,0},{1,2,-1,-2,0,0},{-2,-1,4,2,-2,-1},{-1,-2,2,4,-1,-2},{0,0,-2,-1,2,1},{0,0,-1,-2,1,2}},
    {{2,-2,0,1,-1,0},{-2,4,-2,-1,2,-1},{0,-2,2,0,-1,1},{1,-1,0,2,-2,0},{-1,2,-1,-2,4,-2},{0,-1,1,0,-2,2}},
    {{5,-5,18,-18,1,-1},{0,0,0,0,0,0},{-5,5,-18,18,-1,1},{1,-1,18,-18,5,-5},{0,0,0,0,0,0},{-1,1,-18,18,-5,5}},
    {{1,-1,6,-6,1,-1},{6,-6,36,-36,6,-6},{1,-1,6,-6,1,-1},{-1,1,-6,6,-1,1},{-6,6,-36,36,-6,6},{-1,1,-6,6,-1,1}},
    {{9,3,-6,-2,-3,-1},{-6,6,4,-4,2,-2},{-3,-9,2,6,1,3},{3,1,6,2,-9,-3},{-2,2,-4,4,6,-6},{-1,-3,-2,-6,3,9}},
    {{5,1,0,0,-5,-1},{18,18,0,0,-18,-18},{1,5,0,0,-1,-5},{-5,-1,0,0,5,1},{-18,-18,0,0,18,18},{-1,-5,0,0,1,5}},
    {{1,-1,3,-3,0,0},{1,-1,3,-3,0,0},{-1,1,0,0,1,-1},{-1,1,0,0,1,-1},{0,0,-3,3,-1,1},{0,0,-3,3,-1,1}},
    {{1,3,0,-1,-3,0},{-1,0,1,1,0,-1},{0,-3,-1,0,3,1},{1,3,0,-1,-3,0},{-1,0,1,1,0,-1},{0,-3,-1,0,3,1}},
    {-1, 1, -6, 6,-1, 1},
    {-1, -1, 0, 0, 1, 1},
    {-1, -6, -1, 1, 6, 1},
    {-1, 0,  1, -1, 0, 1}

};

/*Internal Subcell Elemental Matrices*/
static     boundfunc_cell  boundfunc_subcell0
    ={
        {{7,2,14,4,0,0},{2,1,4,2,0,0},{14,4,49,14,0,0},{4,2,14,7,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{7,14,0,2,4,0},{14,49,0,4,14,0},{0,0,0,0,0,0},{2,4,0,1,2,0},{4,14,0,2,7,0},{0,0,0,0,0,0}},
        {{2,-2,4,-4,0,0},{-2,2,-4,4,0,0},{4,-4,14,-14,0,0},{-4,4,-14,14,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{2,4,0,-2,-4,0},{4,14,0,-4,-14,0},{0,0,0,0,0,0},{-2,-4,0,2,4,0},{-4,-14,0,4,14,0},{0,0,0,0,0,0}},
        {{14,4,-14,-4,0,0},{4,2,-4,-2,0,0},{-14,-4,14,4,0,0},{-4,-2,4,2,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{14,-14,0,4,-4,0},{-14,14,0,-4,4,0},{0,0,0,0,0,0},{4,-4,0,2,-2,0},{-4,4,0,-2,2,0},{0,0,0,0,0,0}},
        {{5,-5,13,-13,0,0},{-5,5,-13,13,0,0},{0,0,0,0,0,0},{1,-1,5,-5,0,0},{-1,1,-5,5,0,0},{0,0,0,0,0,0}},
        {{1,-1,3,-3,0,0},{3,-3,9,-9,0,0},{0,0,0,0,0,0},{-1,1,-3,3,0,0},{-3,3,-9,9,0,0},{0,0,0,0,0,0}},
        {{9,3,-9,-3,0,0},{-9,-3,9,3,0,0},{0,0,0,0,0,0},{3,1,-3,-1,0,0},{-3,-1,3,1,0,0},{0,0,0,0,0,0}},
        {{5,1,-5,-1,0,0},{13,5,-13,-5,0,0},{0,0,0,0,0,0},{-5,-1,5,1,0,0},{-13,-5,13,5,0,0},{0,0,0,0,0,0}},
        {{3,-3,9,-9,0,0},{1,-1,3,-3,0,0},{-3,3,-9,9,0,0},{-1,1,-3,3,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{3,9,0,-3,-9,0},{-3,-9,0,3,9,0},{0,0,0,0,0,0},{1,3,0,-1,-3,0},{-1,-3,0,1,3,0},{0,0,0,0,0,0}},
        {-1, 1, -3, 3, 0, 0},
        {-3,-1,  3, 1, 0, 0},
        {-3, 3,  0,-1, 1, 0},
        {-1,-3,  0, 1, 3, 0}
    };

static boundary_function boundfunc
= {
    .coeff={ {  1, -1, -1,  1,  1, -1, -1,  1},
             { -1,  0, -1,  0,-.5, .5,-.5, .5},
             {-.5,-.5, .5, .5, -1, -1,  0,  0}},
    .sym[0]= {0,0,1,1,{0,1,2,3,4,5},{0,1,2,3,4,5},{0,1,2,3}},
    .sym[1]= {0,1,-1,1,{1,0,3,2,5,4},{2,1,0,5,4,3},{1,0,3,2}},
    .sym[2]= {1,0,1,-1,{4,5,2,3,0,1},{3,4,5,0,1,2},{2,3,0,1}},
    .sym[3]= {1,1,-1,-1,{5,4,3,2,1,0},{5,4,3,2,1,0},{3,2,1,0}}
    // sym has modified since c is starting from 0
    };

static boundfunc_cutcell boudfunc_subcellint[4];
struct shapefunction shapefunc;

void set_boundfunc(double dx);
void compute_matricesNonlinearStructure(Matrices_S* ptr_ms, Index_S* ptr_i, Grid_S* ptr_g, Solid* ptr_s, char* fnd);

#endif  /* BOUN_FUNC_H */
