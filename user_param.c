/* parameters that is set by the users*/
#include "user_param.h"


void set_grid (Grid_S* g);
void set_solid(Solid* s, Grid_S* g);
void set_constraint(Constraint_S *c,Grid_S* g,char bc);

void user_param (Grid_S* ptr_g, Solid* ptr_s, TimeMarching* ptr_t,Constraint_S* ptr_c)
{

	int		solid_resolution	= 5;
	double	domain[4]			={-1,1,-1,1};
	double  char_length			=1.0;

	double  density				=10.0;
	double	nu					=0.4;
	double	E					=6;
	char	shape[]				="Beam"; //"Beam","Disk","Elliptic"

	char	boundarycondition	='N'; // N: all Neumann, F:all FSI, D: all Dirichlet
									  // M: all mixed, P: Neumann+Dirichlet
    double  solid_threshold     =1e-4;

	double	timestep			=1e-2;
	double	finaltime			=6;
	double	Newmark_theta		=1/4;
	int		plot_per_steps		=1000;

	char	initial[]			="rest";

	/* setting grid information*/
	ptr_g->char_lengh	= char_length;
	ptr_g->domain_x[0]	= domain[0];
	ptr_g->domain_x[1]	= domain[1];
	ptr_g->domain_y[0]	= domain[2];
	ptr_g->domain_y[1]	= domain[3];

	ptr_g->domain_size[0] = domain[1]-domain[0];
	ptr_g->domain_size[1] = domain[3]-domain[2];

	ptr_g->res = solid_resolution;

	if (ptr_g->domain_size[0] > ptr_g->domain_size[1])
		ptr_g->dx = ptr_g->domain_size[1]/solid_resolution;
	else
		ptr_g->dx = ptr_g->domain_size[0]/solid_resolution;
	set_grid(ptr_g);

	/* setting structure information*/
	ptr_s->length	= char_length;
	ptr_s->rho		= density;
	ptr_s->lambda	= E*nu/(1+nu)/(1-2*nu);
	ptr_s->mu		= E/2/(1+nu);
	ptr_s->threshold= solid_threshold;
	strcpy( ptr_s->type, shape);
	strcpy( ptr_s->initial, initial);
	ptr_s->boundcond = boundarycondition;
	set_solid(ptr_s,ptr_g);

	/* setting time marching information*/
	ptr_t->initial	= 0;
	ptr_t->dt		= timestep;
	ptr_t->nstep	= (int) ceil(finaltime/timestep);
	ptr_t->nstep_output	= plot_per_steps;
	ptr_t->delta	= 0.5;
	ptr_t->theta	= Newmark_theta;

	/* setting straint of*/
	set_constraint(ptr_c,ptr_g,boundarycondition);
}


void set_grid(Grid_S* g){

	int i;
	int n;
	double dx = g->dx;

	g->Nx = (int) ceil(g->domain_size[0]/dx);
	g->Ny = (int) ceil(g->domain_size[1]/dx);
	g->N  = g->Nx*g->Ny;

	g->x_grid = dvector(0,g->Nx-1);
	g->y_grid = dvector(0,g->Ny-1);

	/* setting the x- y- coordinates of the grid*/
	for(i=0; i<g->Nx; i++)
		g->x_grid[i] = g->domain_x[0] + i*dx;

	for(i=0; i<g->Ny; i++)
		g->y_grid[i] = g->domain_y[0] + i*dx;

}

/* conpute the sdf of the solid*/
void set_solid(Solid* s, Grid_S* g){

	int		m = g->Nx;
	int		n = g->Ny;
	int		i,j;
	double	center_x,center_y,r,hl_ratio,Lx,Ly;
	double  dx2 = g->dx/2;

	s->boundary_value = dmatrix(0,m-1,0,n-1);
	s->boundary_sign  = imatrix(0,m-1,0,n-1);

	// Calculate the siged distance function on each grid point
	// Need to specify the shape parameter
	if(strcmp(s->type, "Beam")==0){

		center_x = 0;
		center_y = 0;
		hl_ratio = 1;
		Lx = s->length;
		Ly = Lx*hl_ratio;

		for(i=0;i<m;i++)
			for(j=0;j<n;j++)
				s->boundary_value[i][j] = - DMIN(
				DMIN( (g->x_grid[i]+dx2-center_x)+ Lx/2,  Lx/2 - (g->x_grid[i]+dx2-center_x) ) ,
                DMIN( (g->y_grid[j]+dx2-center_y)+ Ly/2,  Ly/2 - (g->y_grid[j]+dx2-center_y) ) );

	}
	else if (strcmp(s->type, "Disk")==0){

		center_x = 2.0;
		center_y = 2.0;
		r = 1.0;

		for(i=0;i<m;i++)
			for(j=0;j<n;j++)
				s->boundary_value[i][j] = -r + sqrt((g->x_grid[i]+dx2-center_x)*(g->x_grid[i]+dx2-center_x)
				+(g->y_grid[j]+dx2 - center_y)*(g->y_grid[j]+dx2 - center_y));

	}
	else
		printf("Please enter a new shape of solid");

	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			s->boundary_sign[i][j] = (int)((s->boundary_value[i][j] > 0) - (s->boundary_value[i][j] < 0));

}

void set_constraint(Constraint_S *c,Grid_S* g,char bc){

	int		m = g->Nx,n = g->Ny;
	int		i,j;
	int		toln = m*n;

	// body force
	c->body_funct_x = dvector(0,toln-1);
	c->body_funct_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			c->body_funct_x[j*m+i]=0;
			c->body_funct_y[j*m+i]=0;
		}

	// fsi force
	c->fsi_func_x = dvector(0,toln-1);
	c->fsi_func_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			c->fsi_func_x[j*m+i]=0;
			c->fsi_func_y[j*m+i]=0;
		}

	// Neumann
	c->neumann_funct_x = dvector(0,toln-1);
	c->neumann_funct_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			// uni-direction
			if (g->y_grid[j] > 0.4)
				c->fsi_func_y[j*m+i]=1;
			else if (g->y_grid[j] < 0.4)
				c->fsi_func_y[j*m+i]=-1;
			else
				c->fsi_func_y[j*m+i]=0;

			// zero force
			c->fsi_func_x[j*m+i]=0;
			//c->fsi_func_y[j*m+i]=0;
		}

	//Dirichlet
	c->dirichlet_dxi_x = dvector(0,toln-1);
	c->dirichlet_dxi_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			c->dirichlet_dxi_x[j*m+i]=0;
			c->dirichlet_dxi_y[j*m+i]=0;
		}

	//Surface specification
	c->fsineumanndirichlet=(char *)malloc((size_t) ((m*n)*sizeof(char)));

	switch (bc)
	{
		case 'N':
			for(i=0; i<m; i++)
				for(j=0;j<n;j++){
					c->fsineumanndirichlet[j*m+i]='N';
			}
			break;

		case 'P':
			for(i=0; i<m; i++)
				for(j=0;j<n;j++){
					if ((g->y_grid[j] > 0.4) || (g->y_grid[j] < -0.4))
						c->fsineumanndirichlet[j*m+i]='N';
					else
						c->fsineumanndirichlet[j*m+i]='D';
			}
			break;

		case 'D':
			for(i=0; i<m; i++)
				for(j=0;j<n;j++){
					c->fsineumanndirichlet[j*m+i]='D';
			}
			break;

	}
}
