#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable : 4996)
void file_read()
{
	printf("Checking File Fidelity ...\n");
	//---------------------------------------------------------------------
	int i, j;
	/*-----------------------File Presence Tests -----------------------------------------*/   
	FILE *ip1,*ip2,*ip3,*ip4,*ip5;
	ip2=fopen("0.4.Input_Geometry.tmp","r");
	//---------------------------------------------------------------------
	if (ip2 == NULL)
	{
		printf("\tFile '0.4.Input_Geometry.tmp' \t\t Absent\n");	
	}
	else
	{
		printf("\tFile '0.4.Input_Geometry.tmp' \t\t Present\n");	
		fclose(ip2);
	}	
	printf("Checking File Fidelity ......................Complete\n");	
	printf("\n");	
	printf("\n");		
	//---------------------------------------------------------------------
	printf("Starting File Reading ...\n");
	//---------------------------------------------------------------------
		grid_type	= 0;	
		Lx = 0.380;	
		Ly = 0.700;
		Lz = 0.048;

		IMAX  =75;	
		JMAX = 150 ;
		KMAX = 1;
	
		u_bottom = 0.000;
		u_top = 0.000;
		u_left = 0.000;
		u_right = 0.000;
		v_bottom = 0.000;
		v_top = 0.000;
		v_left = 0.000;
		v_right = 0.000;

		Condition				=	-1;
		gas_flowrate			=	500.00;	
		rho_g					=	1.178;
		mu_g					=	1.983e-5;
		
		ROTAMETER_FLOWRATE		=	0.00;
		rho_l					=	1000.0;
		mu_l					=	0.001;
		CONTACT_ANGLE			=	80;
		SURFACE_TENSION_LIQUID	=	0.0723;

		Gf						=	0.40;
		dpf						=	116e-6;
		phi_f					=	1.00;
		rho_f					=	2500.00;
		mu_f					=	0.8;
		
		tol						=	1.0e-1;
		mass_bal				=	1.0e-10;
		relaxu					=	1.00;
		relaxv					=	1.00;
		relaxp					=	1.00;
		relaxke					=	1.00;
		relaxde					=	1.00;
		relaxepfd				=	1.00;

		omega_u					=	0.75;;
		omega_v					=	0.75;;
		omega_p					=	1.40;;
		omega_ke				=	0.75;;
		omega_de				=	0.75;;
		omega_uf				=	0.75;;
		omega_vf				=	0.75;;

		c_1						=	1.44;
		c_2						=	1.92;
		c_mu					=	0.09;
		sigma_k					=	1.00;
		sigma_e					=	1.22;
		
		N_Pressure_ports = 1 ;
		if (N_Pressure_ports > 0)
		{
			Port_X[0]	=	0.06;
			Port_Y[0]	=	0.13;
		}
		
	Flag_Raceway_Computation	=	1;
	N_Particle_type = 1;
	if (N_Particle_type > 0)
	{
		NUM					=	1000;
		dp					=	5e-3;
		insertion_min_X		=	0.00;
		insertion_min_Y		=	0.20;
		insertion_max_X		=	0.05;
		insertion_max_Y		=	0.25;
		insertion_rate		=	100;
		
		deletion_min_X		=	0.00;
		deletion_min_Y		=	0.00;
		deletion_max_X		=	0.10;
		deletion_max_Y		=	0.05;
		discharge_rate		=	0;
	}
	
			
	phi_s				=	1.00;
	rho_s				=	2500.00;
	coeff_friction		=	0.30;
	coeff_restitution	=	0.80;	
	Kn					=	1000.00;			
	sim_time			=	2.00;
	
	mark				=	1.10;
	time_step			=	1.0e-6;		

	Raceway_Voidage		=	0.85;
	Domain_Voidage		=	0.40;
		
	//---------------------------------------------------------------------
	ip2=fopen("0.4.Input_Geometry.tmp","r");
		fscanf(ip2,"%d\n",&n_tuyeres);
		if(n_tuyeres>0)
		{
			for (i=1; i<=n_tuyeres; i++)
			{
				fscanf(ip2,"%lf\n",&Tuyere_height[i]);
				fscanf(ip2,"%lf\n",&Tuyere_opening[i]);
				fscanf(ip2,"%lf\n",&Tuyere_protrusion[i]);
				fscanf(ip2,"%lf\n",&Tuyere_angle[i]);
			}
		}
		
		fscanf(ip2,"%d\n",&NUM_OF_ROTAMETERS);
		if(NUM_OF_ROTAMETERS > 0)
		{
			for (i = 0; i < NUM_OF_ROTAMETERS; i++)
			{				
				fscanf(ip2,"%lf\n",&SOURCE_OF_LIQUID[i][0]);
				fscanf(ip2,"%lf\n",&SOURCE_OF_LIQUID[i][1]);
				fscanf(ip2,"%lf\n",&ROTAMETER_OPENING);
				printf("Rotameter info :: %d %lf %lf\n", i, SOURCE_OF_LIQUID[i][0], SOURCE_OF_LIQUID[i][1]);
			}
		}
		
		fscanf(ip2,"%d\n",&n_structures);
		if(n_structures>0)
		{
			for (i=1; i<=n_structures; i++)
			{
				fscanf(ip2,"%d\n",&structure_type[i]);
				if (structure_type[i] == 1)
				{
					fscanf(ip2,"%lf %lf\n",&C_x[i],&C_y[i]);
					fscanf(ip2,"%lf %lf\n",&B_x[i],&B_y[i]);
					fscanf(ip2,"%lf %lf\n",&A_x[i],&A_y[i]);
					fscanf(ip2,"%lf\n",&voidage_Geom[i]);
				}
				else if (structure_type[i] == 2)
				{
					fscanf(ip2,"%lf %lf\n",&C_x[i],&C_y[i]);
					fscanf(ip2,"%lf\n",&radius[i]);
					fscanf(ip2,"%lf\n",&voidage_Geom[i]);
				}
				else if (structure_type[i] == 3)
				{
					fscanf(ip2,"%lf %lf\n",&L_block[i],&H_block[i]);
					fscanf(ip2,"%lf %lf\n",&C_x[i],&C_y[i]);
					fscanf(ip2,"%lf\n",&voidage_Geom[i]);
				}
			}
		}
		else if (n_structures== -1)
		{
			fscanf(ip2,"%s\n",testname);
		}
	fclose(ip2);
	//---------------------------------------------------------------------
	printf("Input File Reading .........................Complete\n");	
	printf("\n");	
	printf("\n");	
}