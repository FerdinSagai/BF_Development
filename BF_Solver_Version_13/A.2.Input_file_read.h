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
	ip1=fopen("0.3.Input_BEAST.tmp","r");
	ip2=fopen("0.4.Input_Geometry.tmp","r");
	//---------------------------------------------------------------------
	if (ip1 == NULL)
	{
		printf("\tFile '0.3.Input_BEAST.tmp' \t\t Absent\n");	
	}
	else
	{
		printf("\tFile '0.3.Input_BEAST.tmp' \t\t Present\n");	
		fclose(ip1);
	}
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
	ip1=fopen("0.3.Input_BEAST.tmp","r");
		fscanf(ip1,"%d\n",&grid_type);	
		fscanf(ip1,"%lf\n",&Lx);	
		fscanf(ip1,"%lf\n",&Ly);
		fscanf(ip1,"%lf\n",&Lz);

		fscanf(ip1,"%d\n",&IMAX);	
		fscanf(ip1,"%d\n",&JMAX);
		fscanf(ip1,"%d\n",&KMAX);
				
		fscanf(ip1,"%lf\n",&u_bottom);
		fscanf(ip1,"%lf\n",&u_top);	
		fscanf(ip1,"%lf\n",&u_left);	
		fscanf(ip1,"%lf\n",&u_right);	

		fscanf(ip1,"%lf\n",&v_bottom);	
		fscanf(ip1,"%lf\n",&v_top);	
		fscanf(ip1,"%lf\n",&v_left);	
		fscanf(ip1,"%lf\n",&v_right);	
		
		fscanf(ip1,"%d\n", &Condition);
		fscanf(ip1,"%lf\n",&gas_flowrate);	
		fscanf(ip1,"%lf\n",&rho_g);
		fscanf(ip1,"%lf\n",&mu_g);
		
		fscanf(ip1,"%lf\n",&ROTAMETER_FLOWRATE);
		fscanf(ip1,"%lf\n",&rho_l);
		fscanf(ip1,"%lf\n",&mu_l);
		fscanf(ip1,"%lf\n",&CONTACT_ANGLE);
		fscanf(ip1,"%lf\n",&SURFACE_TENSION_LIQUID);
		
		fscanf(ip1,"%lf\n",&Gf);
		fscanf(ip1,"%lf\n",&dpf);
		fscanf(ip1,"%lf\n",&phi_f);
		fscanf(ip1,"%lf\n",&rho_f);
		fscanf(ip1,"%lf\n",&mu_f);
		
		fscanf(ip1,"%lf\n",&tol);
		fscanf(ip1,"%lf\n",&mass_bal);
		fscanf(ip1,"%lf\n",&relaxu);
		fscanf(ip1,"%lf\n",&relaxv);
		fscanf(ip1,"%lf\n",&relaxp);
		fscanf(ip1,"%lf\n",&relaxke);
		fscanf(ip1,"%lf\n",&relaxde);
		fscanf(ip1,"%lf\n",&relaxepfd);
		fscanf(ip1,"%lf\n",&omega_u);
		fscanf(ip1,"%lf\n",&omega_v);
		fscanf(ip1,"%lf\n",&omega_p);
		fscanf(ip1,"%lf\n",&omega_ke);
		fscanf(ip1,"%lf\n",&omega_de);
		fscanf(ip1,"%lf\n",&omega_uf);
		fscanf(ip1,"%lf\n",&omega_vf);
		
		fscanf(ip1,"%lf\n",&c_1);
		fscanf(ip1,"%lf\n",&c_2);
		fscanf(ip1,"%lf\n",&c_mu);
		fscanf(ip1,"%lf\n",&sigma_k);
		fscanf(ip1,"%lf\n",&sigma_e);

		fscanf(ip1,"%d\n",&N_Pressure_ports);
		for (i = 0; i < N_Pressure_ports; i++)
		{
			fscanf(ip1,"%lf\n",&Port_X[i]);
			fscanf(ip1,"%lf\n",&Port_Y[i]);
		}

		fscanf(ip1,"%d\n",&Flag_Raceway_Computation);
		fscanf(ip1,"%d\n",&N_Particle_type);
		for (i = 1; i <= N_Particle_type; i++)
		{
			fscanf(ip1,"%d\n",&NUM);
			fscanf(ip1,"%lf\n",&dp);
			fscanf(ip1,"%lf\n",&insertion_min_X);
			fscanf(ip1,"%lf\n",&insertion_min_Y);
			fscanf(ip1,"%lf\n",&insertion_max_X);
			fscanf(ip1,"%lf\n",&insertion_max_Y);
			fscanf(ip1,"%d\n",&insertion_rate);
			fscanf(ip1,"%lf\n",&deletion_min_X);
			fscanf(ip1,"%lf\n",&deletion_min_Y);
			fscanf(ip1,"%lf\n",&deletion_max_X);
			fscanf(ip1,"%lf\n",&deletion_max_Y);
			fscanf(ip1,"%d\n",&discharge_rate);
		}
		
		fscanf(ip1,"%lf\n",&phi_s);
		fscanf(ip1,"%lf\n",&rho_s);
		fscanf(ip1,"%lf\n",&coeff_friction);
		fscanf(ip1,"%lf\n",&coeff_restitution);	
		fscanf(ip1,"%lf\n",&Kn);			
		fscanf(ip1,"%lf\n",&sim_time);
		
		fscanf(ip1,"%lf\n",&mark);
		fscanf(ip1,"%lf\n",&time_step);		

		fscanf(ip1,"%lf\n",&Raceway_Voidage);
		fscanf(ip1,"%lf\n",&Domain_Voidage);		
	fclose(ip1);
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