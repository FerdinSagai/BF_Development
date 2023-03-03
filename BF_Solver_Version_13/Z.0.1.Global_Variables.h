/*
The purpose of this file is to declare all the variables that will be used in the entire code.
These variables and their purpose is defined comprehensively, elsewhere. Only a brief definition is provided here. 

NOTE: Single variables are defined in a line where necessary to provide clarity

APPENDIX A
Raceway_Voidage is the measure of the voidage inside the raceway
	: 0.85 for Correlation based approximations, This value matches our experimental findings well.
	Similar matching methods have been employed by Apte et. al.(1998?) and Kamble et.al. (2020)
Domain_Voidage is the measure of the voidage inside the Packed bed
	: 0.40 for A loosely packed bed of spheres. Value determined by Yu et. al (1991)
*/

/*
//=================================================================================
//				FUTURE SCOPE
//=================================================================================
	=> 3D implementation
	=> Uniform Cartesian is implemented at present. Other grids must be understood and coded.
	=> Multi-Tuyere systems have not been implemented
	=> Axi-symmetric implementation
*/
#define NO_NEIGHBOUR (-1)		// A macro for indicating a non-populated index
#define MAX_NEIGHBOUR 8			// Maximum number of neighbours that a particle can have
#define VOID_SURR_NUM 15		// Maximum number of particle that can form a void
#define MAXNUM 5000				// Maximum number of DEM particles
#define NVOIDS 2*MAXNUM			// Maximum number of Voids in the packed bed
#define NWALLS 4				// Maximum number of Walls in the domain for DEM particle interaction
#define MaxBinPack 3

// SOME CONSTANTS
#define PI acos(-1.0)
#define gravity 9.81			// 	Acceleration due to gravity
#define CD 0.47
#define inv_12 1/12.0

#define ni 40
#define nj 56
#define ncx 500
#define ncy 500



double Surf_Xi[NWALLS], Surf_Yi[NWALLS];
double Surf_Xf[NWALLS], Surf_Yf[NWALLS];
double Surf_dx[NWALLS], Surf_dy[NWALLS];
double Surf_Angle[NWALLS];
double Surf_Normal_Angle[NWALLS];
double Surf_n[NWALLS][3];
double Surf_len[NWALLS];

double acc_x[MAXNUM], acc_y[MAXNUM];		// Linear Acceleration of Bed Particle
double acc_ang[MAXNUM];						// Angular Acceleration of Bed Particle

bool restart_flag, restart_read;
int Rivulet_Update;
int Droplet_Update;
int Rivulet_Break_Update;
int Previous_ITER_LIQUID;
int N_Particle_type;
double insertion_min_X, insertion_min_Y, insertion_max_X, insertion_max_Y;
double deletion_min_X, deletion_min_Y, deletion_max_X, deletion_max_Y;

int insertion_rate;
int discharge_rate;			// Mass of particles exiting the bed

double Area_RC	= 0.00;
double Area_RC_display	= 0.00;
double Growth_RC_Area	= 0.00;

int Gas_Transient_output_count;
int Fines_Transient_output_count;


double Residue_P_old;

double mass_max, d_max;

double Growth_Fines;
int NN_Model;
double max_deformation;
double max_gasforce;
double Port[20];
int ITER_LIQUID;
//==============================================================================================================================
//			BEGINNING OF GENERIC DEFINITIONS
//==============================================================================================================================
bool status_gas_phase;			// A flag to indicate if the GAS phase is active
bool status_solid_phase;		// A flag to indicate if the GAS phase is active
bool status_liquid_phase;		// A flag to indicate if the GAS phase is active
bool status_fines_phase;		// A flag to indicate if the GAS phase is active

bool init_gas_phase;			// A flag to indicate if the GAS phase has been initialized
bool init_solid_phase;			// A flag to indicate if the SOLID phase has been initialized
bool init_liquid_phase;			// A flag to indicate if the LIQUID phase has been initialized
bool init_fines_phase;			// A flag to indicate if the FINES phase has been initialized

bool steady_state_gas;			// A flag to indicate if the GAS phase has attained steady-state
bool steady_state_solid;		// A flag to indicate if the SOLID phase has attained steady-state
bool steady_state_liquid;		// A flag to indicate if the LIQUID phase has attained steady-state
bool steady_state_fines;		// A flag to indicate if the FINES phase has attained steady-state


int Gas_Convergence_count;
int Fines_Convergence_count;

int Transient_output_count;
int Particle_position_count;
int Convergence_output_count;
int Performance_scheme_count;
int Performance_method_count;
//==============================================================================================================================
//==============================================================================================================================
//==============================================================================================================================
int GFS_iter;				// Iteration counter for the gas computation loop
double solution_time;			// Real time of the simulation
//==============================================================================================================================
// ----- GAS PHASE
clock_t time_Ugb,time_Uga;		// CPU Time measured before and after Gas X velocity computation
clock_t time_Vgb,time_Vga;		// CPU Time measured before and after Gas Y velocity computation
clock_t time_Pb,time_Pa;		// CPU Time measured before and after Pressure computation
clock_t time_KEgb,time_KEga;		// CPU Time measured before and after Gas Turbulence Energy computation
clock_t time_degb,time_dega;		// CPU Time measured before and after Gas Turbulence Dissipation computation
double cpu_time_U_gas;			// CPU Time taken for Gas X velocity computation
double cpu_time_V_gas;			// CPU Time taken for Gas Y velocity computation
double cpu_time_Pressure;		// CPU Time taken for Pressure computation
double cpu_time_KE_gas;			// CPU Time taken for Gas Turbulence Energy computation
double cpu_time_de_gas;			// CPU Time taken for Gas Turbulence Dissipation computation
int Conv_U;				// No. of internal iterations needed for convergence of Gas X velocity
int Conv_V;				// No. of internal iterations needed for convergence of Gas Y velocity
int Conv_P;				// No. of internal iterations needed for convergence of Pressure
int Conv_KE;				// No. of internal iterations needed for convergence of Gas Turbulence Energy
int Conv_de;				// No. of internal iterations needed for convergence of Gas Turbulence Dissipation
double max_Ug;				// Maximum value of X Gas velocity after the computation converges
double max_Vg;				// Maximum value of Y Gas velocity after the computation converges
double max_P;				// Maximum value of Pressure after the computation converges
double max_KE;				// Maximum value of Gas Turbulence Energy after the computation converges
double max_de;				// Maximum value of Gas Turbulence Dissipation after the computation converges
double min_Stress;			// Minimum value of effective stress
// ----- FINES PHASE
clock_t time_Ufb,time_Ufa;		// CPU Time measured before and after Fines X velocity computation
clock_t time_Vfb,time_Vfa;		// CPU Time measured before and after Fines Y velocity computation
clock_t time_hdb,time_hda;		// CPU Time measured before and after Fines dynamic holdup computation
clock_t time_hsb,time_hsa;		// CPU Time measured before and after Fines static holdup computation
double cpu_time_U_fines;		// CPU Time taken for Fines X velocity computation
double cpu_time_V_fines;		// CPU Time taken for Fines Y velocity computation
double cpu_time_dynamic_holdup;		// CPU Time taken for Fines dynamic holdup computation
double cpu_time_static_holdup;		// CPU Time taken for Fines static holdup computation
int Conv_Uf;				// No. of internal iterations needed for convergence of Fines X velocity
int Conv_Vf;				// No. of internal iterations needed for convergence of Fines Y velocity
int Conv_epfd;				// No. of internal iterations needed for convergence of Fines dynamic holdup
double max_Uf;				// Maximum value of X Fines velocity after the computation converges
double max_Vf;				// Maximum value of Y Fines velocity after the computation converges
double max_epfs;			// Maximum value of Fines Static holdup after the computation converges
double max_epfd;			// Maximum value of Fines Dynamic holdup after the computation converges
// ----- SOLID PHASE
clock_t time_NNb, time_NNa;		// CPU Time measured before and after Nearest neighbour determination
clock_t time_FCb, time_FCa;		// CPU Time measured before and after Force Computation
clock_t time_DEMb, time_DEMa;		// CPU Time measured before and after DEM Operations (NN, FC, PU)
double cpu_time_NN;			// CPU Time taken for Nearest neighbour determination
double cpu_time_FC;			// CPU Time taken for Force Computation
double cpu_time_DEM;			// CPU Time taken for DEM Operations (NN, FC, PU)
// ----- COUPLING FORCES
double max_Fgp;				// Maximum value of interaction force between Gas and Solid particles
double max_Fgl;				// Maximum value of interaction force between Gas and Liquid phase
double max_Fgf;				// Maximum value of interaction force between Gas and Fines particles
double max_Fpl;				// Maximum value of interaction force between Solid particles and Liquid phase
double max_Fpf;				// Maximum value of interaction force between Solid particles and Fines phase
double max_Flf;				// Maximum value of interaction force between Liquid and Fines particles
//==============================================================================================================================
// PARALLELIZATION-RELATED VARIABLES
int rank;				// Rank of the MPI processor
int nprocs;				// Total number of MPI processor
int ierr;				// Error flag of the MPI process
int tid;				// Id of the OpenMP thread
int nthreads;				// Total number of OpenMP threads
//==============================================================================================================================
int liquid_update_count;

double LIMIT_FOR_INPUT;

double X_INLET;
double Y_INLET;
double YP_MAX;

double CA;
double surface_tension;
int NUM_OF_ROTAMETERS;
double initial_liquid_volume;
double contact_angle;
double SOURCE_OF_LIQUID[10][2];

int head_count;
double INITIAL_LIQ_VOL_VAL;
double CONTACT_ANGLE;
double SURFACE_TENSION_LIQUID;
double ROTAMETER_X;
double ROTAMETER_Y;
double ROTAMETER_OPENING;
double ROTAMETER_FLOWRATE;  
double CUT_OFF_LIQUID_VELOCITY;
double CUT_OFF_LIQUID_VOLUME;
double total_time = 0.0;



int N_Pressure_ports;
double Port_X[20], Port_Y[20];

//==============================================================================================================================
// DETECTION OF VOIDS IN PACKED BED
int NUM_OF_VOIDS;											// Number of voids detected in the particle bed
int VOID_ENCIRCLED_BY[NVOIDS][VOID_SURR_NUM];				// Provides the particle ID surrounding void (1st index) at position (2nd index)
double VOID_POS[NVOIDS][6];									// Stores void data. 
/*
VOID_POS[NVOIDS][0] :: Centroid X
VOID_POS[NVOIDS][1] :: Centroid Y
VOID_POS[NVOIDS][2] :: Void Volume
VOID_POS[NVOIDS][3] :: Volume ratio (Solid fraction)
VOID_POS[NVOIDS][4] :: Effective diameter
*/
	
int NUM_PARTICLE_VOIDS[NVOIDS];						// Number of particle associated with void(1st index)
int NUM_VOIDS_PARTICLE[MAXNUM];						// Number of voids associated with particle(1st index)

int Num_suspended_particles;
int suspended_particles[MAXNUM];					// Stores id of suspended particles
double suspended_solid_volume[NVOIDS];					// Volume of particle suspended in void


//  UNSURE
int VOIDP[NVOIDS][VOID_SURR_NUM][2];					// VOIDP[i][j][k] :: Numbe of sources to void i from particle 2nd index
int LEFT_NEIGH_VOIDP[NVOIDS][VOID_SURR_NUM][2];			// LEFT_NEIGH_VOIDP[i][j][k] :: 
int RIGHT_NEIGH_VOIDP[NVOIDS][VOID_SURR_NUM][2];		// RIGHT_NEIGH_VOIDP[i][j][k] :: 

//==============================================================================================================================
// DISCRETE LIQUID FLOW
double LIQ_IN_VOID[6*MAXNUM];
double LIQ_IN_VOID_OLD[6*MAXNUM];

double LIQUID_ON_PARTICLE_OLD[MAXNUM][NVOIDS];
double LIQUID_ON_PARTICLE[MAXNUM][NVOIDS];
double LIQ_POS_X[MAXNUM][NVOIDS];
double LIQ_POS_Y[MAXNUM][NVOIDS];


double GOING_RIGHT[MAXNUM][NVOIDS];
double LIQVELX[MAXNUM][NVOIDS];
double LIQVELY[MAXNUM][NVOIDS];
double LIQVELX_OLD[MAXNUM][NVOIDS];
double LIQVELY_OLD[MAXNUM][NVOIDS];
double GOING_RIGHT_OLD[MAXNUM][NVOIDS];

int PAR_VOID[NVOIDS][6];

//---------------------------------------------------------------------------------
// PHYSICAL PROPERTIES
double Lx, Ly, Lz;			// Length, Breadth and Height of the packed bed
int IMAX, JMAX, KMAX;			// Number of cells in the X, Y, and Z direction
double max_x, max_y, min_x, min_y;	// Max and Min size of cells along the X and Y directions
double CFL_x, CFL_y;			// CFL number along the X and Y directions
double rho_g;				// Density of gas phase
double mu_g;				// Viscosity of the gas phase
double dp;				// Size of the packed bed particles
double Gf;				// Mass flux of the fines 
double Gg;				// Mass flux of the gas
double rho_l;				// Density of liquid phase
double mu_l;				// Viscosity of the liquid phase
double rho_f;				// Density of gas phase
double mu_f;				// Viscosity of the gas phase
double dpf;				// Size of the packed bed particles
double phi_f;				// Shape factor of the fines particles
double loading;
double phi_s;				// Shape factor of the packed bed particles
double rho_s;				// Density of Packed bed particles and fines
double gas_flowrate;			// Gas flow rate through tuyere
double vin_x, vin_y, vin;		// Gas velocity at inlet and its components in X and Y (Blast velocity)
double v_super;				// Gas superficial velocity at inlet
double Re;				// Reynolds number of the gas
double dt1;				// Appropritate time step to advance the simulation
//---------------------------------------------------------------------------------
// GRID SPECIFICATION
int grid_type;				// Type of Grid (Uniform Cartesian, Nonuniform Cartesian, Unstructured)
int State[ncx+1][ncy+1];		// Status of each grid point. Solid or void space
double x[ncx+1];			// X co-ordinate of the grid point
double y[ncy+1];			// Y co-ordinate of the grid point
double dx[ncx+1];			// Cell size in the X direction
double dy[ncy+1];			// Cell size in the Y direction
double vfrac[ncx+1][ncy+1];		// Value of void fraction at each grid point
double sfrac[ncx+1][ncy+1];		// Value of solid fraction at each grid point (1 - void fraction))
//---------------------------------------------------------------------------------
// GAS FLOW COMPUTATION
double u[ncx+1][ncy+1];			// Corrected Gas velocity in the X direction at each grid point
double v[ncx+1][ncy+1];			// Corrected Gas velocity in the Y direction at each grid point
double f[ncx+1][ncy+1];			// Projected Gas velocity in the X direction at each grid point
double g[ncx+1][ncy+1];			// Projected Gas velocity in the Y direction at each grid point
double p[ncx+1][ncy+1];			// Gas pressure computed by a Poisson-type equation
double vorticity_g[ncx+1][ncy+1];	// Gas vorticity at each grid point
double ke[ncx+1][ncy+1];		// Turbulent Kinetic Energy of the gas (at previous time step)
double ke1[ncx+1][ncy+1];		// Turbulent Kinetic Energy of the gas (at present time step)
double de[ncx+1][ncy+1];		// Dissipation of the Turbulent Kinetic Energy of the gas (at previous time step)
double de1[ncx+1][ncy+1];		// Dissipation of the Turbulent Kinetic Energy of the gas (at present time step)
double un[ncx+1][ncy+1];		// Corrected Gas velocity in the X direction at each grid point (at previous time step)
double vn[ncx+1][ncy+1];		// Corrected Gas velocity in the Y direction at each grid point (at previous time step)
double pn[ncx+1][ncy+1];		// Corrected Pressure at each grid point (at previous time step)
//---------------------------------------------------------------------------------
// ------------------PARTICLE PHASE COMPUTATION
int Condition;				// Flag to indicate increasing or decreasing condition of the gas flow (raceway hysteresis
int Flag_Raceway_Model;			// Flag to indicate the raceway model used
double Raceway_Voidage;			// Voidage inside the Raceway (Refer Comments A)
double Domain_Voidage;			// Voidage inside the Packed bed (Refer Comment A)
double us[ncx+1][ncy+1];		// Solid velocity in the X direction at each grid point
double vs[ncx+1][ncy+1];		// Solid velocity in the Y direction at each grid point
// ------------------STRESS-BASED METHOD VARIABLES
int m = 9;				// ?? 
int n = 0;				// ?? 
double sigma[ncx][ncy];			// ?? 
double sai[ncx][ncy];			// ?? 
double sigma_init[ncx][ncy];		// ?? 
double sigma_eff[ncx][ncy];		// ?? 
double sigma_xx[ncx][ncy];		// ?? 
double sigma_yy[ncx][ncy];		// ?? 
double sigma_alpha[ncx][ncy];		// ?? 
double sigma_beta[ncx][ncy];		// ?? 
double tau_alpha[ncx][ncy];		// ?? 
double tau_beta[ncx][ncy];		// ?? 
double sigma_xx1[ncx][ncy];		// ?? 
double sigma_yy1[ncx][ncy];		// ?? 
double tau_xy[ncx][ncy];		// ?? 
double xa[ni][nj];			// ?? 
double ya[ni][nj];			// ?? 
double iso_stress;			// Value of stress at the Raceway boundary
// ------------------DEM-RELATED VARIABLES
int Nx, Ny, Nz;				// Grid resolution for the DEM computations
int NUM;				// Number of DEM particles used in packing

double xp[MAXNUM],yp[MAXNUM];		// Position for the Bed Particles
double up[MAXNUM],vp[MAXNUM];		// Velocity of Bed Particle
double mass[MAXNUM];			// Mass of the particle in kg
double MMOI[MAXNUM];			// Moment of inertia 

double d[MAXNUM];			// Diameter of the particle
double omega[MAXNUM];			// Angular Velocity of the bed particles
double angle[MAXNUM][MAX_NEIGHBOUR];	// Angle between colliding pairs
int CNUM[MAXNUM][MAX_NEIGHBOUR];	// Neighbour list for a single particle
			
double Kn,Kt;				// Normal and Tangential spring constants of a single particle
double CnVAL, CtVAL;
//double Ct;				// Damping coefficient in the Normal and Tangential direction for a single particle (Magnitude only)
double KnWall,KtWall;			// Normal and Tangential spring constants of the wall
double CnWall,CtWall;			// Damping coefficient in the Normal and Tangential direction for a wall
double coeff_friction;			// Coefficient of friction
double coeff_restitution;		// Coefficient of restitution
double OLD_DELTAN[MAXNUM][MAXNUM+7];	// Time history of Tangential Forces
double OLD_DELTAT[MAXNUM][MAXNUM+7];	// Time history of Normal Forces
double mark;				// Tolerance margin for collision detection
double STATICtime;			// Standard time step
double sim_time;			// Total simulation time for DEM
int nstep;				// Number of timesteps in the DEM simulation
double rp; 				// Radius of particle.
int packing_option; 			// Packing Options
double limitH;				// Particle Generation height

double x1_t, x2_t, x3_t, x4_t;
double y1_t, y2_t, y3_t, y4_t;
double ke_system;
// ------------------CORRELATION-BASED METHOD VARIABLES
int condition;				// Indicates if the velocity is increasing or decreasing
double raceway_radius;			// Radius of the circular raceway predicted
double Raceway_Center_X; 		// X co-ordinate of the raceway center
double Raceway_Center_Y;		// Y co-ordinate of the raceway center
//---------------------------------------------------------------------------------
// LIQUID FLOW COMPUTATION
double ul[ncx+1][ncy+1];		// Liquid velocity in the X direction at each grid point
double vl[ncx+1][ncy+1];		// Liquid velocity in the Y direction at each grid point
double Area_SL[ncx+1][ncy+1];		// Liquid velocity in the Y direction at each grid point
double Area_GL[ncx+1][ncy+1];		// Liquid velocity in the Y direction at each grid point
//---------------------------------------------------------------------------------
// FINES FLOW COMPUTATION
double uf[ncx+1][ncy+1];		// Corrected Fines velocity in the X direction at each grid point
double vf[ncx+1][ncy+1];		// Corrected Fines velocity in the Y direction at each grid point
double ff[ncx+1][ncy+1];		// Projected Fines velocity in the X direction at each grid point
double gf[ncx+1][ncy+1];		// Projected Fines velocity in the Y direction at each grid point
double vorticity_f[ncx+1][ncy+1];	// Fines vorticity at each grid point
double epfd[ncx+1][ncy+1];		// 
double epfd1[ncx+1][ncy+1];		// 
double epfd2[ncx+1][ncy+1];		// 
double epfs[ncx+1][ncy+1];		// 
double Force_on_Par_X[ncx+1][ncy+1];		// 
double Force_on_Par_Y[ncx+1][ncy+1];		// 
double epld[ncx+1][ncy+1];		// 
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// TURBULENCE MODELLING (k-epsilon)
double lm;				// Mixing length for eddy viscosity model
double c_mu;				// Constant for determination of eddy viscosity
double c_1;				// Constant for dissipation of turbulent kinetic energy
double c_2;				// Constant for dissipation of turbulent kinetic energy
double sigma_k;				// Constant of turbulence modelling used for computation of turbulent kinetic energy
double sigma_e;				// Constant of turbulence modelling used for computation of turbulent kinetic energy dissipation
//---------------------------------------------------------------------------------
// PARAMETERS FOR PHASE COUPLING
double Forcefield_GP[ncx+1][ncy+1];	// Force between the GAS & PARTICLE phases at every grid point
double Forcefield_GL[ncx+1][ncy+1];	// Force between the GAS & LIQUID phases at every grid point
double Forcefield_GF[ncx+1][ncy+1];	// Force between the GAS & FINES phases at every grid point
double Forcefield_PL[ncx+1][ncy+1];	// Force between the PARTICLE & LIQUID phases at every grid point
double Forcefield_PF[ncx+1][ncy+1];	// Force between the PARTICLE & FINES phases at every grid point
double Forcefield_LF[ncx+1][ncy+1];	// Force between the LIQUID & FINES phases at every grid point
//---------------------------------------------------------------------------------
// COMPUTATIONAL CONTROL PARAMETERS
double tol;				// Tolerance for concluding iterative convergence
double mass_bal;
double relaxu;				// Relaxation parameter for gas X velocity update
double relaxv;				// Relaxation parameter for gas Y velocity update
double relaxp;				// Relaxation parameter for gas Pressure update
double relaxke;				// Relaxation parameter for gas Turbulent kinetic energy update
double relaxde;				// Relaxation parameter for gas Turbulent kinetic energy dissipation update
double relaxepfd;			
double omega_u;				// Under relaxation parameter in determining the gas X velocity
double omega_v;				// Under relaxation parameter in determining the gas Y velocity
double omega_p;				// Under relaxation parameter in determining the gas pressure
double omega_ke;			// Under relaxation parameter in determining the gas Turbulent kinetic energy
double omega_de;			// Under relaxation parameter in determining the gas Turbulent kinetic energy dissipation
double omega_uf;			// Under relaxation parameter in determining the fine X velocity
double omega_vf;			// Under relaxation parameter in determining the fine Y velocity
//---------------------------------------------------------------------------------
// COMPUTATION OF RESIDUE
double mass_out_fine;			// Mass of gas exiting the domain through the top
double mass_in_fine;			// Mass of gas entering the domain through the tuyere
double mass_fine;			// Mass of gas retained in the domain
double mass_out_fine_top;
double mass_out_fine_bot;
double Residue_V;			// Difference between previous and current velocity fields
double Residue_P;			// Difference between previous and current pressure fields
double Residue;				// Maximum residue of velocity and pressure fields
//---------------------------------------------------------------------------------
// INTERNAL GEOMETRY DEFINITION VARIABLES
double u_bottom, v_bottom;		// X and Y Velocity boundary condition of the bottom wall
double u_top, v_top;			// X and Y Velocity boundary condition of the top wall
double u_left, v_left; 			// X and Y Velocity boundary condition of the left wall
double u_right, v_right;		// X and Y Velocity boundary condition of the right wall
int n_tuyeres;				// Number of tuyeres (inlets) in the domain
double Tuyere_protrusion[100];		// Protrusion of the tuyere inside the domain
double Tuyere_opening[100];		// Opening of the tuyere. Actual inlet cross section. 
double Tuyere_height[100];		// Height of the tuyere from bottom of the bed.
double Tuyere_angle[100];		// Angle of the tuyere.
int ii01,ii02,jj01,jj02;		// Integer indices of the tuyere			// 
int n_structures;			// Number of interma structures inside the domain
int structure_type[100];		// Type of the structure: (Block, Circular, Triangular)
double A_x[100], A_y[100];		// Co-ordinates of the point A
double B_x[100], B_y[100];		// Co-ordinates of the point B
double C_x[100], C_y[100];		// Co-ordinates of the point C
double radius[100];			// Radius (for circular structure only)
double L_block[100], H_block[100];	// Length and Height (for block structure only)
double voidage_Geom[100];		// Porosity of the structure
char testname[50];			// Filename of text file containing the structure info
//=================================================================================
//			BEGINNING OF SOLID PHASE MODELLING
//=================================================================================
//---------------------------------------------------------------------------------


// Computation


double D_T;
double r;
double h;

int update_count=0;

int GFS_Max, DEM_Max;
double Area_RC_Previous;

double epfs_frac,epfs_max,epfs_limit;

double limit_m;
// BED DROP DETAILS
double L_down;						// Bottom of Bed for removal
double Y_bottom_wall;
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
int Flag_Raceway_Computation;

int Bin0, Bin1, Bin2, Bin3, Bin4, Bin5, Bin6, Bin7, Bin8, Bin9 ,Bin10;

double X1,X2;
double Y1,Y2;
double theta;
double time_step;
					// Radius of particle.

// UNKNOWN
int nt;
double Bed_height;
int Bed_height_idx ;

/*-------------------------------------------*/
/*-------------------------------------------*/
/*-----------Array Based Variables-----------*/
/*-------------------------------------------*/
/*-----------Parameters----------------------*/

/*-------------------------------------------*/

double epfs_total;
double Previous[ncx+1][ncy+1];