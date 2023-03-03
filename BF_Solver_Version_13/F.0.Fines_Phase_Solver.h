#include "F.1.Fines_Velocity_Update_U.h"
#include "F.2.Fines_Velocity_Update_V.h"
#include "F.3.Fines_Velocity_Correction.h"
#include "F.4.Dynamic_holdup_calculations.h"
#include "F.5.Static_holdup_calculations.h"
#include "F.6.Fines_Vorticity_Calculation.h"
#include "F.7.Fines_Convergence.h"
#include "F.8.Fines_Outputs.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void initialise_fines_phase()
{
	init_fines_phase				=	true;	
	Fines_Convergence_count 		= 0;
	Fines_Transient_output_count	=	1;
	loading	= (Gf/(vin*rho_f));
	for (int i = 0; i <= IMAX; i++)
	{
		for (int j = 0; j <= JMAX; j++)
		{ 
			uf[i][j]	=	u[i][j];
			vf[i][j]	=	v[i][j];

			ff[i][j]	=	uf[i][j];
			gf[i][j]	=	vf[i][j];
			
			epfs[i][j]	=	0.0;
			epfd[i][j]	=	0.1*loading;
			epfd1[i][j]	=	0.1*loading;
			epfd2[i][j]	=	0.1*loading;
		}
	}
	steady_state_gas = false;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Fines_Phase_Solver()
{   
	int i, j;
	enforce_domain_physics();
	//-----Computing approximate U velocity as ff-------------------------------------
	uf_update();
	printf("\tX Fines Velocity...................:: Iterations = %d      ||\tCPU time = %lf seconds ||\tMax Value = %lf ms-1\n",Conv_Uf,cpu_time_U_fines,max_Uf);
	//-----Computing approximate V velocity as gf--------------------------------------
	vf_update();
	printf("\tY Fines Velocity...................:: Iterations = %d      ||\tCPU time = %lf seconds ||\tMax Value = %lf ms-1\n",Conv_Vf,cpu_time_V_fines,max_Vf);
	//-----Correcting U and V(from f and g) resp.-------------------------------------
	uvf_correct();					
	enforce_domain_physics();
 	//-----------Dynamic holdup calculations------------------------------------------			
 	dynamic();
	printf("\tDynamic Holdup.....................:: Iterations = %d      ||\tCPU time = %lf seconds ||\tMax Value = %lf m2\n",Conv_epfd,cpu_time_dynamic_holdup,max_epfd);
  	//-----------Static holdup calculations-------------------------------------------			
	Static_holdup_correlations(); 
	//-----Computing vorticity of the fines--------------------------------------------
	vorticity_fines_calc();	
	//-----------Mass conservation of fines verification-------------------------------		
	Fines_Convergence_Check();
	//if ( (GFS_iter%20) == 0)
	{
		printf("------------------------------------------------------------\n");
		printf("Wrtiting Transient and steady-state hydrodynamics for gas-fines\n");
		printf("------------------------------------------------------------\n");
		fines_convergence_data();
		save_output_fines_transient();
		save_output_fines_steady_state();
	}
	
	if (mass_in_fine != 0.0)
	{
		if (mass_out_fine >= 0.90 * mass_in_fine)
		{
			steady_state_fines	=	true;
			save_output_fines_steady_state();
			printf("------------------------\n");
			printf("Steady state achieved - Gas-Fines Flows\n");
			printf("------------------------\n");
		}
	}
}