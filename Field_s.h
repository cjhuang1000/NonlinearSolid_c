#ifndef FIELD_S_H
#define FIELD_S_H 1

typedef struct
{

	short int *G2g,			*g2G;
	short int *G2g_before,	*g2G_before;
	short int *C2c_fsi,		*c2C_fsi;
	short int *C2c_neumann,	*c2C_neumann;
	short int *C2c_dirichlet,*c2C_dirichlet;
    short int *keep;

} Index_xixy;

typedef struct {

  /* list of the internal/boundary cells indices*/
  unsigned int *cell_interior;
  unsigned int *cell_boundary;
  unsigned int *cell_fsi;
  unsigned int *cell_neumann;
  unsigned int *cell_dirichlet;

  /* list of indices of involved variables*/
  Index_xixy xix;

  /* Position of the involved variables in the list*/
  Index_xixy xiy;

  /* Number of involved cells/variables*/
  int cell_N_interior,cell_N_FSI,cell_N_boundary,cell_N_Neumann,cell_N_Dirichlet;
  int xix_N         ,xiy_N;
  int xix_N_before  ,xiy_N_before;
  int xix_Ncell_FSI ,xiy_Ncell_FSI;
  int xix_Ncell_Neumann     ,xiy_Ncell_Neumann;
  int xix_Ncell_Dirichlet   ,xiy_Ncell_Dirichlet;

} Index_S;



typedef struct {

	/* displacement,velocity */
	double *xi_x,       *xi_y;
	double *dxi_x,      *dxi_y;
    double *ddxi_x,     *ddxi_y;
    double *inc_dxi_x,  *inc_dxi_y;
    double *tau_x,      *tau_y;

} Field_S;


typedef struct {

	/* mass, stiffness matrices, load vectors */
	double **MS_xx,   **MS_yy;
	double **KlS_xx,  **KlS_xy,  **KlS_yy;
    double **KnlS_xx, **KnlS_yy;

    double *FS_x,     *FS_y;

    double **NS_xc,   **NS_yc;
    double **DS_xc,   **DS_yc;

} Matrices_S;


/* boundary type at the solid boundary */
typedef struct {

	// body force
	double *body_funct_x,  *body_funct_y;

	// FSI Surface Forces on the staggerred grid
	double *fsi_func_x, *fsi_func_y;

	// Neumann Suface Forces Definition
    double *neumann_funct_x, *neumann_funct_y;

	// Dirichlet Suface Constraints Definition
	double *dirichlet_dxi_x,*dirichlet_dxi_y;

	// Projection of the surface specification on the main grid
	char *fsineumanndirichlet;

} Constraint_S;

typedef struct {

  /* global domain dimension */
  double	char_lengh;  // characteristic length of the problem
  double	domain_x[2],domain_y[2]; // domain interval
  double	domain_size[2];

  double	res,dx;   // resolution and grid spacing

  int		 Nx,Ny,N;     // number of grid in x- and y- direction, total grid number
  double	*x_grid, *y_grid;

} Grid_S;


typedef struct {

	double  length;
	double	rho;
	double	mu;
	double	lambda;

	char	type[80];
	char	initial[80];
	char	boundcond;

	int		**boundary_sign;
	double	**boundary_value;

	double  threshold;

}Solid;


typedef struct {

	double	initial;
	double	dt;
	int		nstep;
	int		nstep_output;
	double	delta;
	double	theta;

}TimeMarching;

#endif  /* FIELD_S_H */
