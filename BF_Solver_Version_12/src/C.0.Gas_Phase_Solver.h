#include "C.1.Gas_Velocity_Update_U.h"
#include "C.2.Gas_Velocity_Update_V.h"
#include "C.3.Pressure.h"
#include "C.4.Gas_Velocity_Correction.h"
#include "C.5.Gas_Turbulence_KE.h"
#include "C.6.Gas_Turbulence_DE.h"
#include "C.7.Gas_Vorticity_Calculation.h"
#include "C.8.Gas_Convergence.h"
#include "C.9.Gas_Outputs.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void initialise_gas_phase()
{
	init_gas_phase				=	true;	
	Gas_Convergence_count		= 0;
	Gas_Transient_output_count	=	1;
	
//----Initializaing hydrodynamic quantities for Gas and fines------------------
	for (int i = 0;i < IMAX;i++)
	{	
		for (int j = 0;j < JMAX; j++)
		{ 
			u[i][j]			= 0.0;
			v[i][j]			= 0.0;
			p[i][j]			= 0.0;
			ke[i][j]		= 1.0;
			de[i][j]		= pow(ke[i][j],1.5);	
			
			f[i][j]			=	u[i][j];
			g[i][j]			=	v[i][j];
			ke1[i][j]		=	ke[i][j];
			de1[i][j]		=	de[i][j];			
			
			vorticity_g[i][j]	= 0.0;			
			un[i][j]		= u[i][j];
			vn[i][j]		= v[i][j];	
		}
	}
//----Initializaing intermediate hydrodynamic quantities for gas and fines------------
	for (int i = 0; i < IMAX; i++)
	{	
		for (int j = 0; j < JMAX; j++)
		{ 	
			uf[i][j]			= 0.0;
			vf[i][j]			= 0.0;
			ff[i][j]			= uf[i][j];
			gf[i][j]			= vf[i][j];
			vorticity_f[i][j]	= 0.0;
		}
	} 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Gas_Phase_Solver()
{   
	enforce_domain_physics();
	//-----------Computing approximate U velocity as f--------------------------------------
	u_update();
	printf("\tX Gas Velocity.....................:: Iterations = %d      ||\tCPU time = %lf seconds ||\tMax Value = %lf ms-1\n",Conv_U,cpu_time_U_gas,max_Ug);
	//-----------Computing approximate V velocity as g--------------------------------------		
	v_update();
	printf("\tY Gas Velocity.....................:: Iterations = %d      ||\tCPU time = %lf seconds ||\tMax Value = %lf ms-1\n",Conv_V,cpu_time_V_gas,max_Vg);
	//-----------Solving the pressure linking equation--------------------------------------		
	p_update();
	printf("\tP Pressure.........................:: Iterations = %d ||\tCPU time = %lf seconds ||\tMax Value = %lf Pa\n",Conv_P,cpu_time_Pressure,max_P);	
	//-----------Correcting U and V from f and g resp.---------------------------------------
	uv_correct();
	//u_correct();
	//v_correct();
	enforce_domain_physics();
  	//-----------Computing turbulent kinetic energy------------------------------------------	
 	ke_update();
	printf("\tGas Kinetic Energy.................:: Iterations = %d      ||\tCPU time = %lf seconds ||\tMax Value = %lf J\n",Conv_KE,cpu_time_KE_gas,max_KE);
	//-----------Computing dissipation of turbulent kinetic energy----------------------------
	de_update();
	printf("\tGas Dissipation of Kinetic Energy..:: Iterations = %d      ||\tCPU time = %lf seconds ||\tMax Value = %lf J\n",Conv_de,cpu_time_de_gas,max_de);
	//-----------Computing vorticity of the gas-----------------------------------------------
	vorticity_gas_calc();	
	//-----------Computing dissipation of turbulent kinetic energy----------------------------
	Gas_Convergence_Check();
	
	if ( ((GFS_iter%20) == 0) || (GFS_iter == 1) )
	{
		printf("------------------------------------------------------------\n");
		printf("Wrtiting Transient and steady-state hydrodynamics for gas\n");
		printf("------------------------------------------------------------\n");
		gas_convergence_data();
		save_output_gas_transient();
	}
	
	if( Residue <= 10*(tol/tol) )
	{	
		steady_state_gas = true;
		save_output_gas_steady_state();
		printf("------------------------\n");
		printf("Steady state achieved - Pure Gas Flows\n");
		printf("------------------------\n");
	}
}