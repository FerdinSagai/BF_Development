int SEARCH_DOMAIN[MAXNUM][5];
int	INV_P_CELL[6000][MaxBinPack];
#include "D.1.1.Void_Calc.h"
void initialize_variables(void);
void init_particle_position_read(void);

void add_particles(double ,double ,double ,double ,int ,double);
void delete_particles(double,double,double,double,int);

void Nearest_Neighbour_Algorithm_Brute_Force(int);
void Nearest_Neighbour_Algorithm_L_Search(void);

void Combined_Algorithm(void);

double compute_collision_time(int, int);
void calculation_force_particles(int, int);

double Distance_between_Line_and_Point(int, double, double);
void calculation_force_wall(int, int);

void calculation_force_obstacles(int, double, double, double, double);
void calculation_force_gasdrag(int);
void update_particle_info(int, double, double, double, double);

void save_particle_position(void);
void save_final_particle_position(void);
void reset_neighbours(void);

void DEM_void_fraction(void);
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void DEM_Module()
{
	int i, j;
	int ini, fin;
	int Particles_Added, Particles_Removed;
	double generation_region;
	int pack_actual, pack;
	
	/*------Noting the cpu start time --------------------------------------------------------*/
	time_DEMb	= clock();
	/*------Variable Initializations ---------------------------------------------------------*/	
	initialize_variables();
	Particles_Added		= 0;
	Particles_Removed	= 0;
	ke_system			= 0.00;	
			
	if (update_count == 1)		// Initial packing of the bed structure
	{
		int mode = 1;			// Actually pack the bed OR use a file to read in the particle positions
		
		if(mode == 0)			// Actually pack the bed
		{
			pack_actual			= NUM;
			NUM					= 0;
			generation_region	= (insertion_max_X-insertion_min_X) * (insertion_max_Y-insertion_min_Y); 
			pack				= 0.5*generation_region/(dp*dp);
			if(pack > pack_actual) pack = pack_actual;
	 
			nstep		= (int) ceil(sim_time/STATICtime);
			ini			= 1;
			fin			= nstep;
		}
		else					// Use a file to read in the particle positions	
		{
			FILE *fp1;
			fp1 = fopen("1.5.Particle_Positions_Post.inp","r");
			for(i = 0;i < NUM;i++)
			{
				fscanf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf\n",&xp[i],&yp[i],&d[i],&up[i],&vp[i]);
			}
			fclose(fp1);  
			
			for(i = 0; i < NUM; i++)
			{
				d[i] = dp;
				printf("%d %lf %lf %lf %lf %lf\n",i, xp[i], yp[i], d[i], up[i], vp[i]);
			}
			
			pack		= 0;
			nstep		= 20000;
			ini			= 1;
			fin			= nstep;
		}		
	}
	else
	{
		sim_time	=	0.2;
		STATICtime	=	time_step;
		ini			=	nstep;
		nstep		=	nstep + (int) ceil(sim_time/STATICtime);
		fin			=	nstep;
	}
	
	ke_system	= 0.00; 
	for(i = 0; i < NUM; i++)
	{
		ke_system	=	ke_system + 0.5 * mass[i] * ( up[i]*up[i] + vp[i]*vp[i] );
	}
	/*------Time Marching for the particles-----------------------------------------------*/
   	for (nt = ini; nt <= fin; nt++)
	{
		if(nt == 1)
		{
			printf("nt=%d, time=%f\n",nt,nt*STATICtime);
			cpu_time_NN		= 0.00;
			max_deformation		= 0.00;
			max_gasforce		= 0.00;
			Particles_Added		= 0;
			Particles_Removed	= 0;
			ke_system		= 0.00;				
		}
		if( (nt == 1) || (nt%10000 == 0) )
		{
			if ( (NUM < pack_actual) && (update_count == 1) )
			{
				add_particles(insertion_min_X,insertion_max_X,insertion_min_Y,insertion_max_Y,pack,dp);
				printf("%d %d\n", pack, NUM );
				if( (pack_actual - NUM) < pack ) 
				{
					pack = pack_actual - NUM;
				}
			}
	
			Particles_Removed = NUM;
			delete_particles(deletion_min_X,deletion_max_X,deletion_min_Y,deletion_max_Y,discharge_rate);	
			Particles_Removed = Particles_Removed - NUM;
			
			printf("Update Count = %d :: nt=%d out of nstep=%d\n",update_count,nt,fin);
			printf("\t Number of particles....................................%d \n",NUM);
			printf("\t Particles Added........................................%d \n",Particles_Added);
			printf("\t Particles Removed......................................%d \n",Particles_Removed);
			printf("\t time for DEM Process ..................................%lf seconds\n",cpu_time_NN);
			printf("\t Kinetic Energy of the system ..........................%lf Joules \n",ke_system);
			printf("\t Maximum deformation of a particle......................%lf m \n",max_deformation);
			printf("\t Maximum gas force of a particle........................%lf N \n",max_gasforce);
						
			save_particle_position();
			if ( (update_count == 1) && (nt > 100) )
			{
				if (ke_system < 0.050)
				{
					break;
				}				
			}			
			cpu_time_NN		= 0.00;
			max_deformation		= 0.00;
			max_gasforce		= 0.00;
			Particles_Added		= 0;
			Particles_Removed	= 0;
			ke_system		= 0.00;
		}
		//Nearest_Neighbour_Algorithm_Brute_Force();
		Nearest_Neighbour_Algorithm_L_Search();
		Combined_Algorithm();			
		
		if ( (xp[i] != xp[i]) || (yp[i] != yp[i]) )
		{
			printf("Particle ID:: %d\n", i);
			printf("Velocity of the particle :: %lf %lf \n", up[i], vp[i]);
			int pause;
			scanf("%d", &pause);
		}
		
	}  

	if( update_count == 1 ) 
	{
		FILE *fp1;
		fp1 = fopen("9.0.Particle_Positions_Post_Packing.nxt","w");
		for(i = 0;i < NUM;i++)
		{
			fprintf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf\n",xp[i],yp[i],d[i],up[i],vp[i]);
		}
		fclose(fp1);  
	}
	
	DEM_void_fraction();
	save_particle_position();
	Particle_position_count++;
	/*------Saving Final position of Particles -------------------------------------------*/
	save_final_particle_position();
	/*------Noting the cpu end time ------------------------------------------------------*/
	time_DEMa	= clock();	
	/*------Noting the total cpu time for this computation--------------------------------*/
	cpu_time_DEM	= cpu_time_DEM + ((double)(time_DEMa - time_DEMb) )/(nthreads*CLOCKS_PER_SEC);	
	printf("\t Total time for the system ..........................%lf seconds \n",cpu_time_DEM);
	printf("\t FUTURE TEST PURPOSE ....... Running with Gas\n");
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void initialize_variables()
{
	int i,j, k;	
	// Left Surface
	Surf_Xi[0] = 0.00;
	Surf_Yi[0] = 0.00;
	Surf_Xf[0] = 0.00;
	Surf_Yf[0] = Ly;
	// Bottom Surface
	Surf_Xi[1] = Lx;
	Surf_Yi[1] = 0.00;
	Surf_Xf[1] = 0.00;
	Surf_Yf[1] = 0.00;
	// Right Surface
	Surf_Xi[2] = Lx;
	Surf_Yi[2] = Ly;
	Surf_Xf[2] = Lx;
	Surf_Yf[2] = 0.00;

	// Surface Properties
	for (int i = 0; i < NWALLS; i++)
	{
		Surf_dx[i]		= Surf_Xf[i] - Surf_Xi[i];
		Surf_dy[i]		= Surf_Yf[i] - Surf_Yi[i];
		
		if (Surf_dx[i] == 0.00 )
		{
			if (Surf_dy[i] > 0.00 )		Surf_Angle[i]	= (PI/2); //Check
			else if (Surf_dy[i] < 0.00 )Surf_Angle[i]	= 3*PI/2; //Check
		}
		else if (Surf_dy[i] == 0.00 )
		{
			if (Surf_dx[i] > 0.00 )		Surf_Angle[i]	= 2*PI;
			else if (Surf_dx[i] < 0.00 )Surf_Angle[i]	= PI;
		}		
		else
		{
			Surf_Angle[i]	= atan(Surf_dy[i]/Surf_dx[i]);
		}
		Surf_Normal_Angle[i]= Surf_Angle[i] - (PI/2);
		if( Surf_Normal_Angle[i]<0.00) Surf_Normal_Angle[i]= Surf_Normal_Angle[i] + 2*PI;
		
		Surf_n[i][0] 		= cos(Surf_Normal_Angle[i]);
		Surf_n[i][1] 		= sin(Surf_Normal_Angle[i]);		
		Surf_n[i][2] 		= 0.00;		
		Surf_len[i] = sqrt( sqr(Surf_dx[i]) + sqr (Surf_dy[i]) );
		
		printf("%d %lf %lf \n", i, Surf_Normal_Angle[i]*180/PI, Surf_Angle[i]*180/PI);
	}
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void add_particles(double minX,double maxX,double minY,double maxY,int N,double diameter)
{
	double mass_min = 10000;	
	double range_X = (maxX - minX);
	double range_Y = (maxY - minY);
	
	srand(time(0));
	for(int i = NUM; i < (NUM + N); i++)
	{
		xp[i]		= minX + range_X * rand()/RAND_MAX;
		yp[i]		= minY + range_Y * rand()/RAND_MAX;
		up[i]		= 0.00;   
		vp[i]		= 0.00;  
		omega[i]	= 0.00;
		acc_x[i]	= 0.00;
		acc_y[i]	= 0.00;
		acc_ang[i]	= 0.00;
		d[i]  		= diameter;
		mass[i]		= rho_s*(4.0/3.0)*PI*pow(d[i]/2.0,3.0);
		MMOI[i]		= (2.0/5.0)*mass[i]*pow(d[i]/2.0,2.0); 
		mass_min	=	min(mass_min, mass[i]);		
		printf("Adding Particle :: %d %le %le %lf %lf %lf %lf\n",i, mass[i], MMOI[i], xp[i], yp[i], up[i], vp[i]);
	}		
	NUM = NUM + N;
	
	if (time_step>0.1*sqrt(mass_min/Kn))
	{
		time_step=0.1*sqrt(mass_min/Kn);
	}
	STATICtime	=	time_step;
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void delete_particles(double minX,double maxX,double minY,double maxY,int N)
{
	int deleted_N = 0;
	double range_X = (maxX - minX);
	double range_Y = (maxY - minY);
	
	if (N > 0)
	{
		srand(time(0));	
		for(int i = 0; i < NUM; i++)
		{
			if ( (xp[i] < maxX) && (xp[i] > minX) && (yp[i] < maxY) && (yp[i] > minY))
			{
				//double r = float(rand())/float(RAND_MAX);
				//if( r  > 0.50 )	
				{
					xp[i]		= xp[NUM-1];
					yp[i]		= yp[NUM-1];
					up[i]		= up[NUM-1];
					vp[i]		= vp[NUM-1];
					omega[i]	= omega[NUM-1];
					acc_x[i]	= acc_x[NUM - 1];
					acc_y[i]	= acc_y[NUM - 1];
					acc_ang[i]	= acc_ang[NUM - 1];
					NUM			= NUM - 1;
					deleted_N++;		
					if ( deleted_N >= N ) break;
				}
			}						
		}
	}
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void Combined_Algorithm()
{
	//------Noting the cpu start time ------------------------------------------------------
	time_NNb	= clock();
	//------Determination of neighbours for particle i by searching its search domain-------
	#pragma omp parallel private(tid)
	{		
		tid			= omp_get_thread_num();
		nthreads	= omp_get_max_threads();
		int start	= 1 + tid*int(NUM/nthreads);
		int finish	= start + NUM/nthreads - 1;
		if(tid == 0) start = 0;
		if(tid == nthreads-1)	finish = NUM-1;
		for(int i = start; i <= finish; i++)
		{	
			for(int N = 0; N < 5; N++)
			{
				int MC	=	SEARCH_DOMAIN[i][N];
				if (MC >= 0)
				{
					for(int idx = 0; idx < MaxBinPack; idx++)
					{		
						int j = INV_P_CELL[MC][idx];
						if (j >= 0 && j != i) 
						{			
							double d_avg = 0.50 * (d[i] + d[j]);					
							double R	= ( sqr(xp[i] - xp[j]) +  sqr(yp[i] - yp[j]) );
							if ( R <= sqr(mark*d_avg) )
							{
								calculation_force_particles(i,j);
								if(SEARCH_DOMAIN[i][N] != SEARCH_DOMAIN[j][N])
								{
									calculation_force_particles(j,i);
								}
							}	
						}
					}
				}
			}
			
			for(int k = 0; k < 3; k++)
			{	
				calculation_force_wall(i,k);
			}	
/*			
			double top,  bot,  left, right;
			top	= h + (0.50 * D_T);
			bot	= h - (0.50 * D_T);
			left	= 0.00; 
			right	= r;  
			calculation_force_obstacles(i, top, bot, left, right);
			
			for(int N = 1; N<= n_structures; N++)
			{
				top	= C_y[N] + 0.50 * H_block[N];
				bot	= C_y[N] - 0.50 * H_block[N];
				left	= C_x[N] - 0.50 * L_block[N];
				right	= C_x[N] + 0.50 * L_block[N];
				calculation_force_obstacles(i, top, bot, left, right);
			}
*/
			calculation_force_gasdrag(i);
			
			acc_y[i] = acc_y[i] - gravity;
			xp[i] = xp[i] + up[i] * STATICtime + 0.5 * acc_x[i] * pow(STATICtime, 2.0);
			yp[i] = yp[i] + vp[i] * STATICtime + 0.5 * acc_y[i] * pow(STATICtime, 2.0);
			up[i] = up[i] + acc_x[i] * STATICtime;
			vp[i] = vp[i] + acc_y[i] * STATICtime;
			omega[i] = omega[i];// +acc_ang[i] * STATICtime;
			ke_system	=	ke_system + 0.5 * mass[i] * ( up[i]*up[i] + vp[i]*vp[i] );
			
			acc_x[i] = 0.00;
			acc_y[i] = 0.00;
		}	
	}
	/*------Thread synchronization is performed at this point--------------------*/	
	#pragma omp barrier	
	/*------Noting the cpu end time ------------------------------------------------------*/
	time_NNa	= clock();			
	/*------Noting the total cpu time for this computation--------------------------------*/
	cpu_time_NN	= cpu_time_NN + ((double)(time_NNa - time_NNb) )/(nthreads*CLOCKS_PER_SEC);
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void update_particle_info(int i, double ForceX, double ForceY, double Torque, double Contact_time)
{
	acc_x[i] = acc_x[i] + (ForceX / mass[i]);
	acc_y[i] = acc_y[i] + (ForceY / mass[i]);
	acc_ang[i] = (Torque / MMOI[i]);
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
double compute_collision_time(int i, int j)
{
	double contact_time, travel_time;
	double collision_time;
	double xij,yij;
	double xfij,yfij;
	double rij,rfij;
	double vij;
	double rvij;
	double vxij,vyij;
	double det;
	double coeff_a, coeff_b, coeff_c;
	double root1, root2, discr;

	double d_avg	= 0.50*(d[i] + d[j]); 
	vxij	= up[i] - up[j];							// Relative Velocity between particle i and j (X direction)
	vyij	= vp[i] - vp[j];							// Relative Velocity between particle i and j (Y direction)
	xij	= xp[i] - xp[j];							// Distance between particle i and j (X direction)
	yij	= yp[i] - yp[j];							// Distance between particle i and j (Y direction)
	xfij	= xij + STATICtime*vxij;						// Predicated distance between particle i and j (X direction)
	yfij	= yij + STATICtime*vyij;						// Predicated distance between particle i and j (Y direction)

	coeff_a	= vxij* vxij + vyij* vyij;	 					// Square of Relative velocity  (Coefficient a) 
	coeff_b	= 2*(xij * vxij + yij * vyij); 						// Approach vector (Coefficient b) 
	coeff_c	= ( (xij * xij  + yij * yij) - (xij * xij  + yij * yij) ) - (d_avg * d_avg);	// Difference between initial and contact position (Coefficient c) 

	rvij	= xij*vxij + yij*vyij; 							// Approach vector
	vij	= pow(vxij,2.0) + pow(vyij,2.0); 					// Square of Relative velocity between particles i and j
	rij	= pow(xij,2.0) + pow(yij,2.0); 						// Square of distance between particles i and j
	rfij	= pow(xfij,2.0) + pow(yfij,2.0);					// Expected New position of particle j wrt particle i after STATICtime
	discr	= pow(coeff_b,2.0) - 4.0 * coeff_a * coeff_c;				// Discriminant (b^2- 4ac)

	// --------------
	if ( rij > pow(d_avg,2.0) ) 
	{
		if ( rfij < pow(d_avg,2.0) )
		{	// Option A
			if (discr >= 0.0)
			{
				root1	= (-coeff_b + sqrt(discr)) / (2.0*coeff_a);
				root2	= (-coeff_b - sqrt(discr)) / (2.0*coeff_a);
				if ((STATICtime - root1) > 0.0)
				{
					travel_time	= root1;
				}
				else
				{
					travel_time	= root2;
				}
				contact_time	= STATICtime - travel_time;
			}
			else if (discr < 0.0) 
			{
				contact_time	= 0.0;					// Anamoly condition ?
				travel_time	= STATICtime;
			}
		}
		else
		{	// Option B
			travel_time	= STATICtime;
			contact_time	= 0.0;
		}
	}
	else
	{
		if ( rfij > pow(d_avg,2.0) ) 
		{	// Option D		
			if (discr >= 0.0) 
			{
				root1	= (-coeff_b + sqrt(discr)) / (2.0*coeff_a);
				root2	= (-coeff_b - sqrt(discr)) / (2.0*coeff_a);
				if ( (STATICtime - root1) >0.0)
				{
					contact_time	= root1;
				}
				else
				{
					contact_time	= root2;
				}
				travel_time	= STATICtime - contact_time;
			}
			else if (discr < 0.0) 
			{								// Anamoly condition ?
				contact_time	= STATICtime;
				travel_time	= 0.00;
			}
		}
		else
		{	// Option C
				contact_time	= STATICtime;
				travel_time	= 0.00;
		} 
	}
	
	//Handling anamolies due to round off errors
	if (contact_time < 0.0)
	{
		contact_time	= 0.0;
	}

	if (contact_time > STATICtime)
	{
		contact_time	= STATICtime;
	}
	
	return contact_time;
}
//====================================================================================================================
void calculation_force_particles(int i, int j)
{
	double delta_N, delta_T;
	double velocity_N, velocity_T;
	double sliding;
	double Contact_time, travel_time;
	double angle;

	double Cn = 2.0 * (-log(coeff_restitution) / (sqrt(pow(PI, 2.0) + pow(log(coeff_restitution), 2.0)))) * sqrt(mass[i] * Kn / 2.0);
	double Ct = Cn;
	// =======================================================================================================
	double vxij = up[i] - up[j];    				// Relative Velocity between particle i and j (X direction)
	double vyij = vp[i] - vp[j];    				// Relative Velocity between particle i and j (Y direction)
	double xij = xp[i] - xp[j];    					// Distance between particle i and j (X direction)
	double yij = yp[i] - yp[j];    					// Distance between particle i and j (Y direction)
	double rij = pow(xij, 2.0) + pow(yij, 2.0);		// Relative distance between the particles
	double d_avg = 0.5 * (d[i] + d[j]);
	if (rij < pow(mark * d_avg, 2.0))
	{
		max_deformation = 0.00;
		//  Angle computation.
		angle = (180.0 / PI) * atan2((yp[j] - yp[i]), (xp[j] - xp[i]));
		if (angle < 0.0)
		{
			angle = 360.0 + angle;
		}
		velocity_N = vxij * cos(angle * PI / 180.0) + vyij * sin(angle * PI / 180.0);
		velocity_T = vxij * sin(angle * PI / 180.0) - vyij * cos(angle * PI / 180.0);
		//  Contact time computation		
		Contact_time = compute_collision_time(i, j);
		//  Normal incremental deformation calculation
		delta_N = fabs(velocity_N * Contact_time);
		if (Contact_time == STATICtime)
		{
			delta_N = 0.5 * fabs(d_avg - sqrt(rij));
		}	
		//  Tangential incremental deformation calculation
		delta_T = (velocity_T + 0.5 * d_avg * (omega[i] - omega[j])) * Contact_time;
	
		//  Columb Criteria for Sliding. If particles slide, then only frictional forces are in effect. If sliding occurs then the sign of the frictional forces will be opposite to the deformation.															
		if (fabs(Kt * delta_T) >= coeff_friction * fabs(Kn * delta_N))			sliding = 1.00;
		//  Columb Criteria for Sliding. If particles roll, then, dashpot forces are in effect. In the absence of sliding, the tangential spring and damper are activated. Ct will always oppose motion. 
		if (fabs(Kt * delta_T) < coeff_friction * fabs(Kn * delta_N))			sliding = 0.00;		
		//  Determination of Cn and Ct in the X direction			
		double Force_Normal		= (-fabs(Kn * delta_N) + (1.0 - sliding) * Cn * fabs(velocity_N));
		double Force_Tangential	= (1.0 - sliding) * (-(Kt * delta_T) - Ct * velocity_T);
		double Force_Frictional	= sliding * -sign(velocity_T) * fabs(coeff_friction * Kn * delta_N);
		double ForceX			= Force_Normal * cos(angle * PI / 180.0) + Force_Tangential * sin(angle * PI / 180.0) + Force_Frictional * sin((angle) * PI / 180.0);
		double ForceY			= Force_Normal * sin(angle * PI / 180.0) - Force_Tangential * cos(angle * PI / 180.0) + Force_Frictional * cos((angle) * PI / 180.0);
		double Torque 			= Kt * delta_T * 0.5 * d_avg;

		//  Determination of Aggregate Displacement	
		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
	}
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
						// FORCE DUE TO WALLS
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
						// FORCE DUE TO OBSTACLES
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void calculation_force_obstacles(int i, double top, double bot, double left, double right)
{	
	double Cn, Ct;
	double delta_N, delta_T;
	double velocity_N, velocity_T;
	double sliding, friction, frictionx, frictiony;
	double Contact_time, travel_time;
	double ForceX, ForceY, Torque;
	double Particle_MinX, Particle_MinY, Particle_MaxX, Particle_MaxY;
	double Margin_MinX, Margin_MinY, Margin_MaxX, Margin_MaxY;
	double distance_from_wall;
	
	CnVAL	=	2.0*(-log(coeff_restitution)/( sqrt(pow(PI,2.0) + pow(log(coeff_restitution),2.0))))*sqrt(mass[i]*Kn/2.0);
	CtVAL	=	CnVAL; 
	CnWall	=	CnVAL;
	CtWall	=	CtVAL;

	Particle_MinX	=	xp[i] - (0.50*d[i]);
	Particle_MinY	=	yp[i] - (0.50*d[i]);
	Particle_MaxX	=	xp[i] + (0.50*d[i]);
	Particle_MaxY	=	yp[i] + (0.50*d[i]);

	Margin_MinX	=	xp[i] - (0.50*mark*d[i]);
	Margin_MinY	=	yp[i] - (0.50*mark*d[i]);
	Margin_MaxX	=	xp[i] + (0.50*mark*d[i]);
	Margin_MaxY	=	yp[i] + (0.50*mark*d[i]);
	/*-------------------------------------------LEFT SIDE OF OBSTACLE-----------------------------------------------------------*/		
	if( ( Margin_MaxX > left ) && ( Margin_MaxX < right ) && (Particle_MinY >= bot) && (Particle_MaxY <= top) )
	{
		//  Angle computation -- X not applicable for walls.. 		
		velocity_N		= up[i];
		velocity_T		= vp[i];				
		//  Contact time computation
		if (Particle_MaxX > left) 					// In contact with the wall
		{
			travel_time		=	0.00;
			Contact_time		=	STATICtime;
		}			
		else if (Particle_MaxX <= left) 				// Close enough to consider a hit
		{
			if (velocity_N != 0.0)
			{
				distance_from_wall	=	left - Particle_MaxX;
				travel_time		=	fabs(distance_from_wall/velocity_N);
				
				if ( (travel_time < STATICtime)  &&  (velocity_N > 0.0) )
				{
					Contact_time		=	STATICtime - travel_time;
					
				}
				else 
				{
					Contact_time		=	0.0;
				}
			}
			else
			{
				Contact_time			=	0.0;
			}
		}
		//  Normal deformation calculation
		delta_N			=	fabs(velocity_N * Contact_time);
		if (Contact_time == STATICtime) 
		{
			delta_N		=	fabs(Particle_MaxX - left);          
		}	
		
		//  Tangential deformation calculation
		delta_T		=	( velocity_T + (0.5 * d[i] * omega[i]) ) * Contact_time;
		
		//  Columb Criteria for Sliding. If particles slide, then, only frictional forces are in effect.
		if (fabs(KtWall*delta_T) >= coeff_friction*fabs(KnWall*delta_N))
		{
			if (delta_T<0.0)
			{
				friction	=	fabs(KnWall*delta_N*coeff_friction);
			}
			else
			{
				friction	=	-fabs(KnWall*delta_N*coeff_friction);
			}
			sliding		=	1.0;
		}
		
		//  Columb Criteria for Sliding. If particles roll, then, frictional and tangential forces are in effect.
		if(fabs(KtWall*delta_T) < coeff_friction*fabs(KnWall*delta_N))
		{
			friction	=	-KtWall * delta_T - CtWall * velocity_T;
			sliding		=	0.0;
		}
					
		ForceX		=	fabs(KnWall * delta_N) - (1.0-sliding) * CnWall * velocity_N;
		ForceY		=	friction;
		Torque		=	KtWall * delta_T * 0.50 * d[i]; 
		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
	}	
	/*-------------------------------------------RIGHT SIDE OF OBSTACLE-----------------------------------------------------------*/		
	if( ( Margin_MinX > left ) && ( Margin_MinX < right ) && (Particle_MinY >= bot) && (Particle_MaxY <= top) )
	{
		//  Angle computation -- X not applicable for walls.. 		
		velocity_N		= up[i];
		velocity_T		= vp[i];		
		//  Contact time computation
		if (Particle_MinX < right) 					// In contact with the wall
		{
			travel_time		=	0.00;
			Contact_time		=	STATICtime;
		}			
		else if (Particle_MinX >= right) 				// Close enough to consider a hit
		{
			if (velocity_N != 0.0)
			{
				distance_from_wall	=	right - Particle_MinX;
				travel_time		=	fabs(distance_from_wall/velocity_N);
				
				if ( (travel_time < STATICtime)  &&  (velocity_N < 0.0) )
				{
					Contact_time		=	STATICtime - travel_time;
					
				}
				else 
				{
					Contact_time		=	0.0;
				}
			}
			else
			{
				Contact_time			=	0.0;
			}
		}
		//  Normal deformation calculation
		delta_N			=	fabs(velocity_N * Contact_time);
		if (Contact_time == STATICtime) 
		{
			delta_N		=	fabs(right - Particle_MinX);          
		}	
		
		//  Tangential deformation calculation
		delta_T		=	( velocity_T + (0.5 * d[i] * omega[i]) ) * Contact_time;
		
		//  Columb Criteria for Sliding. If particles slide, then, only frictional forces are in effect.
		if (fabs(KtWall*delta_T) >= coeff_friction*fabs(KnWall*delta_N))
		{
			if (delta_T<0.0)
			{
				friction	=	fabs(KnWall*delta_N*coeff_friction);
			}
			else
			{
				friction	=	-fabs(KnWall*delta_N*coeff_friction);
			}
			sliding		=	1.0;
		}
		
		//  Columb Criteria for Sliding. If particles roll, then, frictional and tangential forces are in effect.
		if(fabs(KtWall*delta_T) < coeff_friction*fabs(KnWall*delta_N))
		{
			friction	=	-KtWall * delta_T - CtWall * velocity_T;
			sliding		=	0.0;
		}
					
		ForceX		=	fabs(KnWall * delta_N) - (1.0-sliding) * CnWall * velocity_N;
		ForceY		=	friction;
		Torque		=	KtWall * delta_T * 0.50 * d[i]; 
		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
	}	
	/*-------------------------------------------TOP_TUYERE BOUNDARY-------------------------------------------------------------*/
	if( ( Margin_MinY < top) && (Margin_MinY > bot)  && (Particle_MinX >= left) && (Particle_MaxX <= right) )
	{
		//  Angle computation -- X not applicable for walls.. 	
		velocity_N		= vp[i];
		velocity_T		= up[i];
		//  Contact time computation
		if (Particle_MinY < top)
		{
			Contact_time	= STATICtime;
			travel_time	= 0.00;
		}
		else
		{			
			if (velocity_N != 0.00)
			{
				distance_from_wall	=	Particle_MinY - top;
				travel_time		=	fabs( distance_from_wall/velocity_N );
				if ( (travel_time < STATICtime)  &&  (vp[i] < 0.0))
				{
					Contact_time		=	STATICtime - travel_time;
				}
				else
				{
					Contact_time		=	0.0;
				}
			}
			else
			{
				Contact_time		=	0.0;
			}
		}
		//  Normal deformation calculation
		delta_N			=	fabs(velocity_N * Contact_time);
		if (Contact_time == STATICtime) 
		{
			delta_N	=	fabs( Particle_MinY - top );       
		}	

		//  Tangential deformation calculation
		delta_T		=	( velocity_T + (0.5 * d[i] * omega[i]) ) * Contact_time;
		
		//  Columb Criteria for Sliding. If particles slide, then, only frictional forces are in effect.
		if (fabs(KtWall*delta_T) >= coeff_friction*fabs(KnWall*delta_N))
		{
			if (delta_T < 0.0)
			{
				friction	=	fabs(KnWall*delta_N*coeff_friction);
			}
			else
			{
				friction	=	-fabs(KnWall*delta_N*coeff_friction);
			}
			sliding		=	1.0;
		}
		else
		{
			friction	=	-KtWall * delta_T - CtWall * velocity_T;
			sliding		=	0.0;
		}		
		
		ForceX	=	friction;
		ForceY	=	fabs(KnWall * delta_N) - (1.0-sliding) * CnWall * velocity_N;
		Torque	=	KtWall * delta_T * d[i] * 0.50;
		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
	}
	/*-------------------------------------------BOTTOM_TUYERE BOUNDARY----------------------------------------------------------*/
	if  ( (Margin_MaxY > bot) && ( Margin_MaxY < top)  && ( Particle_MinX >= left) && (Particle_MaxX <= right) )
	{
		//  Angle computation -- X not applicable for walls.. 		
		velocity_N		= vp[i];
		velocity_T		= up[i];		
		//  Contact time computation
		if (Particle_MaxY > bot)
		{
			Contact_time	= STATICtime;
			travel_time	= 0.00;
		}
		else
		{			
			if (velocity_N != 0.00)
			{
				distance_from_wall	=	bot - Particle_MaxY;
				travel_time		=	fabs( distance_from_wall / velocity_N);
				if ( (travel_time < STATICtime)  &&  (velocity_N > 0.0))
				{
					Contact_time		=	STATICtime - travel_time;
				}			
				else
				{
					Contact_time		=	0.0;
				}
			}
			else
			{
				Contact_time		=	0.0;
			}
		}

		//  Normal deformation calculation
		delta_N			=	fabs(velocity_N * Contact_time);
		if (Contact_time == STATICtime) 
		{
			delta_N	=	fabs( Particle_MaxY - bot );       
		}	

		//  Tangential deformation calculation
		delta_T		=	( velocity_T + (0.5 * d[i] * omega[i]) ) * Contact_time;
		
		//  Columb Criteria for Sliding. If particles slide, then, only frictional forces are in effect.
		if (fabs(KtWall*delta_T) >= coeff_friction*fabs(KnWall*delta_N))
		{
			if (delta_T < 0.0)
			{
				friction	=	fabs(KnWall*delta_N*coeff_friction);
			}
			else
			{
				friction	=	-fabs(KnWall*delta_N*coeff_friction);
			}
			sliding		=	1.0;
		}
		
		//  Columb Criteria for Sliding. If particles roll, then, frictional and tangential forces are in effect.
		if (fabs(KtWall*delta_T) < coeff_friction*fabs(KnWall*delta_N))
		{
			friction	=	-KtWall * delta_T - CtWall * velocity_T;
			sliding		=	0.0;
		}		

		ForceX	=	friction;
		ForceY	=	-fabs(KnWall * delta_N) - (1.0-sliding) * CnWall * velocity_N;
		Torque	=	KtWall * delta_T * d[i] * 0.50;
		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
	}		
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
						// FORCE DUE TO GAS
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void calculation_force_gasdrag(int i)
{
	int j, k;
	int nx, ny;
	double ug, vg;
	double ForceX, ForceY, Torque, Contact_time;

 	nx	=	floor(xp[i]/dx[1]);
	ny	=	floor(yp[i]/dy[1]);
	
	ug	=	0.50*( u[nx-1][ny] + u[nx][ny]);
	vg	=	0.50*( v[nx][ny-1] + v[nx][ny]);	

 	if (vfrac[nx][ny] != 0.00)
	{
		ForceX			= 0.5 * CD * PI * pow(d[i],2.0) * rho_g * fabs(ug - up[i])*(ug - up[i]) * pow(vfrac[nx][ny],-2.7); 
		ForceY			= 0.5 * CD * PI * pow(d[i],2.0) * rho_g * fabs(vg - vp[i])*(vg - vp[i]) * pow(vfrac[nx][ny],-2.7); 
		Torque			= 0.00;
		Contact_time	= STATICtime;
		
		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
		
		Force_on_Par_X[nx][ny] = ForceX;
		Force_on_Par_Y[nx][ny] = ForceY;
		
		max_gasforce	=	 max( max_gasforce, ForceX);
		max_gasforce	=	 max( max_gasforce, ForceY);
	}
	else
	{
		ForceX		= 0.00;
		ForceY		= 0.00;
		Torque		= 0.00;
		Contact_time	= 0.00;
	}
 }
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void DEM_void_fraction()
{
	int i,j,k;
	double phvo, delx, dely, s;
	
	for(i = 0; i < IMAX-1; i++)
	{
		for(j = 0;j < JMAX-1; j++)
		{
			vfrac[i][j]  = 1.00;
			double s=0;		
			for(int k = 0;k < NUM; k++)
			{
				if(((( (x[i] <= (xp[k]+(dp/2))) && (x[i]>(xp[k]-(dp/2)))) 
				||     ((x[i+1]<(xp[k]+(dp/2))) && (x[i+1]>=(xp[k]-(dp/2))))) 
				&& (   ((y[j]<=(yp[k]+(dp/2))) && (y[j]>(yp[k]-(dp/2)))) 
				||     ((y[j+1]<(yp[k]+(dp/2))) && (y[j+1]>=(yp[k]-(dp/2)))))) 
				|| (   ((x[i+1]>=(xp[k]+dp/2)) && (x[i]<=(xp[k]-dp/2))) 
				&& (   ((y[j+1]>=(yp[k]+dp/2)) && (y[j]<=(yp[k]-dp/2)) )) ))
				{									
					delx=(min(x[i+1],(xp[k]+(dp/2)))-max(x[i],(xp[k]-(dp/2))));
					dely=(min(y[j+1],(yp[k]+(dp/2)))-max(y[j],(yp[k]-(dp/2))));
					phvo=(PI*delx*dely/4);
					s=s+phvo;
				}
				
				vfrac[i][j] = 1 - (s/(dx[i]*dy[j]));
				
				if (vfrac[i][j] <= 0.3)
				{
					vfrac[i][j] = 0.3;
				}
				else if (vfrac[i][j] == 1)
				{
					vfrac[i][j] = 0.95;
				}
				else 
				{
					vfrac[i][j] = vfrac[i][j];
				}
			}
		}
	}	
	
	FILE *gp;			
	gp=fopen("3.0.DEM_Void_Fraction.dat","a");
	if(Particle_position_count == 1)
	{
		fprintf(gp,"TITLE = States of the domain\n");
		fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"State\",\"Void Fraction\"\n");
	}
	fprintf(gp, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",Particle_position_count, IMAX, JMAX);
	
	double yy = 0.5*dy[1];
	for (j=1;j<=JMAX;j++)
	{
		double xx  =0.5*dx[1];
		for (i=1;i<=IMAX;i++)
		{
			fprintf(gp,"%e %e %d %e\n",xx, yy,State[i][j],vfrac[i][j]);
			// UPDATE X POSITION
			xx=xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy=yy + 0.5*(dy[j]+dy[j+1]);
	}
	fclose(gp);
}
//=============================================================================================================================================
void Nearest_Neighbour_Algorithm_Brute_Force(int i)
{
	printf("The particle %d has coordinates %lf %lf\n", i, xp[i], yp[i]);
	int nnb = 0;
	for (int j = 0; j < NUM; j++)
	{
		if (i != j)
		{
			double rij = sqrt(sqr(xp[i] - xp[j]) + sqr(yp[i] - yp[j]));
			double dij = 0.5 * (d[i] + d[j]);
			if (rij <= (mark * dij))
			{
				if (nnb < MAX_NEIGHBOUR)
				{
					CNUM[i][nnb] = j;
					nnb++;
				}
				else
				{
					for (int k = 0; k < nnb; k++)
					{
						int ID_1 = CNUM[i][k];
						double rij_n = sqrt(sqr(xp[i] - xp[ID_1]) + sqr(yp[i] - yp[ID_1]));
						double rij_e = sqrt(sqr(xp[i] - xp[ID_1]) + sqr(yp[i] - yp[ID_1]));
						if (rij_n <= rij_e)
						{
							CNUM[i][k] = j;
							break;
						}
					}
				}
			}
		}
	}
}
//=============================================================================================================================================
void Nearest_Neighbour_Algorithm_L_Search()
{
	//------Noting the cpu start time --------------------------------------------------------
	time_NNb = clock();
	//------Local variables declaration--------------------------------------------------------
	int i, j, k;
	int nbin_x, nbin_y;
	int bx, by;
	double BIN_S;
	double rv_ij;
	int index;
	int N_Neigh[MAXNUM];
	int NeighCount;
	int TBINS;
	/*------Determine the number of bins in the domain----------------------------------------*/
	BIN_S = (mark * d_max);
	nbin_x = ceil(Lx / BIN_S);
	nbin_y = ceil((3.0 * Ly) / BIN_S);
	TBINS = (nbin_x + 1) * (nbin_y + 1);	
	//------Initialize array for particles in bin k------------------------------------------
	for (k = 0; k < TBINS; k++)
	{
		for (int idx = 0; idx < MaxBinPack; idx++)
		{
			INV_P_CELL[k][idx] = -1;
		}
	}
	//------Determination of search domains--------------------------------------------------
	for (i = 0; i < NUM; i++)
	{
		//------Initialize search domain for particle i----------------------------------
		SEARCH_DOMAIN[i][0] = -1;
		SEARCH_DOMAIN[i][1] = -1;
		SEARCH_DOMAIN[i][2] = -1;
		SEARCH_DOMAIN[i][3] = -1;
		SEARCH_DOMAIN[i][4] = -1;
		//------Determine the cell indices for particle i--------------------------------
		bx = ((xp[i]) / (BIN_S));
		by = ((yp[i]) / (BIN_S));
		if (bx < 0)
		{
			xp[i] = 0.5 * d[i];
			bx = ((xp[i]) / (BIN_S));
		}
		if (by < 0)
		{
			yp[i] = 0.5 * d[i];
			by = ((yp[i]) / (BIN_S));
		}
		if (bx > nbin_x)
		{
			xp[i] = Lx - (0.5 * d[i]);
			bx = ((xp[i]) / (BIN_S));
		}
		if (by > nbin_y)
		{
			yp[i] = Ly - (0.5 * d[i]);
			by = ((yp[i]) / (BIN_S));
		}
		int MC = ((by)*nbin_x) + bx;

		//------Populate the array for particles in bin k--------------------------------
		index = 0;
		for (k = 0; k < MaxBinPack; k++)
		{
			if (INV_P_CELL[MC][k] >= 0)
			{
				index = index + 1;
			}
			else
			{
				break;
			}
		}
		INV_P_CELL[MC][index] = i;

		//------Populate the search domain for particle i--------------------------------
		SEARCH_DOMAIN[i][0] = MC;
		SEARCH_DOMAIN[i][1] = (by)*nbin_x + (bx - 1); 	// WEST (SAME PLANE)
		SEARCH_DOMAIN[i][2] = (by + 1) * nbin_x + (bx - 1); 	// NORTH-WEST(SAME PLANE)
		SEARCH_DOMAIN[i][3] = (by + 1) * nbin_x + bx; 	// NORTH(SAME PLANE)
		SEARCH_DOMAIN[i][4] = (by + 1) * nbin_x + (bx + 1);  	// NORTH-EAST(SAME PLANE)	
		//------Treatment for particles near bounding walls------------------------------
		// Treatment of Western wall
		if (bx == 0)
		{
			SEARCH_DOMAIN[i][1] = -1; 			// WEST (SAME PLANE)
			SEARCH_DOMAIN[i][2] = -1; 			// NORTH-WEST(SAME PLANE)
		}

		// Treatment of Eastern wall		
		if (bx == nbin_x - 1)
		{
			SEARCH_DOMAIN[i][4] = -1; 			// NORTH-EAST(SAME PLANE)	
		}

		// Treatment of Northern wall	
		if (by == nbin_y - 1)
		{
			SEARCH_DOMAIN[i][2] = -1; 			// NORTH-WEST(SAME PLANE)
			SEARCH_DOMAIN[i][3] = -1; 			// NORTH(SAME PLANE)
			SEARCH_DOMAIN[i][4] = -1; 			// NORTH-EAST(SAME PLANE)
		}
	}
}
//====================================================================================================================
double Distance_between_Line_and_Point(int Surf_id, double Px, double Py)
{
	double Parameter_1x  = Px - Surf_Xi[Surf_id];
	double Parameter_1y  = Py - Surf_Yi[Surf_id];
	
	double VdotPA 	= Parameter_1x * Surf_dx[Surf_id] + Parameter_1y * Surf_dy[Surf_id];
	double Vmag	  	= Surf_len[Surf_id];
	double t1		= VdotPA/(Vmag*Vmag);
	double X_int 	= Surf_dx[Surf_id] * t1 + Surf_Xi[Surf_id];
	double Y_int 	= Surf_dy[Surf_id] * t1 + Surf_Yi[Surf_id];
	double distance = sqrt( sqr(Px - X_int) + sqr(Py - Y_int) );
	
	// FOR DEBUGGING PURPOSES: DO NOT DELETE
	if(Surf_id == 3)
	{
		printf( "\t Surface :: %d \n", Surf_id);
		printf( "\t Surface :: %lf %lf \n", Surf_Xi[Surf_id], Surf_Yi[Surf_id]);
		printf( "\t Surface :: %lf %lf \n", Surf_Xf[Surf_id], Surf_Yf[Surf_id]);
		printf( "\t Surface :: %lf \n", Vmag);
		printf( "\t Direction Vec :: %lf %lf\n", Surf_dx[Surf_id], Surf_dy[Surf_id]);
		printf( "\t Particle :: %lf %lf\n", Px, Py);
		printf( "\t Numerator :: %lf\n", VdotPA);
		printf( "\t t  :: %lf\n", t1);
		printf( "\t Intersection :: %lf %lf\n", X_int, Y_int);
		printf( "\tDistance :: %lf\n", distance);
	}	
	return distance;
}
//====================================================================================================================
void calculation_force_wall(int i, int Surf_id)
{ 
	double delta_N, delta_T;
	CnWall = 2.0 * (-log(coeff_restitution) / (sqrt(pow(PI, 2.0) + pow(log(coeff_restitution), 2.0)))) * sqrt(mass[i] * Kn / 2.0); 
	CtWall = CnWall;
	
	double Wall_Angle			= Surf_Normal_Angle[Surf_id]* 180.0/PI;
	double velocity_N			= up[i] * cos(Wall_Angle * PI / 180.0) + vp[i] * sin(Wall_Angle * PI / 180.0);
	double velocity_T			= up[i] * sin(Wall_Angle * PI / 180.0) - vp[i] * cos(Wall_Angle * PI / 180.0);
	double contact_time			= 0.0;
	double travel_time			= 0.0;
	double sliding				= 0.0;
	double distance_from_wall	= Distance_between_Line_and_Point(Surf_id,xp[i],yp[i]);
	if (distance_from_wall < (0.50 * d[i]))							// In contact with the wall
	{
		travel_time = 0.00;
		contact_time = STATICtime - travel_time;
	}
	else if (distance_from_wall < mark * (0.50 * d[i]))				// Close enough to consider a hit
	{
		if (velocity_N != 0.0)
		{
			travel_time = fabs(distance_from_wall / velocity_N);
			if (travel_time < STATICtime)
			{
				contact_time = STATICtime - travel_time;
			}
			else
			{
				contact_time = 0.0;
			}
		}
		else
		{
			contact_time = 0.0;
		}
	}
	
	if (contact_time != 0.00)
	{
		//  Incremental Normal deformation calculation
		delta_N = fabs(velocity_N * contact_time);
		if (contact_time == STATICtime)
		{
			delta_N = fabs(distance_from_wall);
		}
		//  Incremental Tangential deformation calculation
		delta_T = (velocity_T + (0.5 * d[i] * omega[i])) * contact_time;
		//  Columb Criteria for Sliding. If particles slide, then, only frictional forces are in effect.
		if (fabs(KtWall * delta_T) >= coeff_friction * fabs(KnWall * delta_N))			sliding = 1.0;
		//  Columb Criteria for Sliding. If particles roll, then, frictional and tangential forces are in effect.
		if (fabs(KtWall * delta_T) < coeff_friction * fabs(KnWall * delta_N))			sliding = 0.0;
		double Force_Normal		= fabs(KnWall * delta_N) + (1.0 - sliding) * CnWall * fabs(velocity_N);
		double Force_Tangential	= (1.0 - sliding) * (-(KtWall * delta_T) - CtWall * velocity_T);
		double Force_Frictional	= sliding * -sign(velocity_T) * fabs(KnWall * delta_N * coeff_friction);
		double ForceX 			= Force_Normal * cos(Wall_Angle * PI / 180.0) + Force_Tangential * sin(Wall_Angle * PI / 180.0) + Force_Frictional * sin((Wall_Angle) * PI / 180.0);
		double ForceY 			= Force_Normal * sin(Wall_Angle * PI / 180.0) - Force_Tangential * cos(Wall_Angle * PI / 180.0) + Force_Frictional * cos((Wall_Angle) * PI / 180.0);
		double Torque 			= KtWall * delta_T * 0.50 * d[i];
		update_particle_info(i, ForceX, ForceY, Torque, contact_time);
	}
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void save_particle_position()
{
	int i,j;
	char name[80];
	double xx, yy;
	
	FILE *gp1;
	gp1 = fopen("3.1.Transient_Particle_Positions.dat", "a");	
	if(Particle_position_count == 1)
	{
		fprintf(gp1,"TITLE = Properties at transient states\n");
		fprintf(gp1, "VARIABLES=\"X\",\"Y\",\"Dp\"\n");
	}

	fprintf(gp1, "ZONE T = \"%d\"\n",nt);
	for (i=0; i<NUM; i++)
	{
		fprintf(gp1,"%e %e %e\n",fabs(xp[i]), fabs(yp[i]), d[i]);
	} 
	fclose(gp1);
	
	FILE *fp1;
	fp1=fopen("3.7.Results_Particle_Force.dat","wt");
	fprintf(fp1,"TITLE = Properties at Final state\n");
	fprintf(fp1, "VARIABLES=\"X (m)\",\"Y (m)\",\"F_gas-particles X\",\"F_gas-particles Y\"\n");
	fprintf(fp1, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",update_count, IMAX-1, JMAX-1);
	
	yy = 0.5*dy[1];
	for (j=1;j<=JMAX-1;j++)
	{
		xx  =0.5*dx[1];
		for (i=1;i<=IMAX-1;i++)
		{ 
			fprintf(fp1,"%e %e %e %e\n",xx,yy,Force_on_Par_X[i][j],Force_on_Par_Y[i][j]);
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}
	fclose(fp1);	 
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void save_final_particle_position()
{
	FILE *fp;
	fp = fopen("2.3.Particle_Positions.int", "wt");
		for (int i = 0; i < NUM; i++)
		{
			fprintf(fp,"%lf %lf %lf\n",xp[i],yp[i],d[i]);
		}
	fclose(fp);
}