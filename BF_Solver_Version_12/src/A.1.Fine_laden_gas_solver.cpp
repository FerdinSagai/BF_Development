// System Header files
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
// User-defined Header files
#include "Z.0.1.Global_Variables.h" 
#include "Z.0.2.Functions.h"
#include "A.2.Input_file_read.h" 
#include "A.3.0.Domain_Setup.h"
#include "A.4.Initialization.h"
#include "A.5.Report.h"
#include "A.6.Monitoring.h"

#include "B.8.Part_force.h"
#include "C.0.Gas_Phase_Solver.h"
#include "D.0.Solid_Phase_Solver.h"
#include "E.0.Liquid_Phase_Solver.h"
#include "F.0.Fines_Phase_Solver.h"
// Function Prototypes
void Print_Status_Label(void);
void Determine_Phase_Status(void); 
void Read_Restart_File(void);
void Write_Restart_File(void);
//--------------------------------------------------------------------------------------------------------------*/
//--------------------------------------------------MAIN PROGRAM------------------------------------------------*/
//--------------------------------------------------------------------------------------------------------------*/
int main(int argc, char** argv) 
{   
	GFS_Max	= 500;
	restart_flag = false;
	
//-----------------------Program Title---------------------------------------------------------------------------
	printf("\n\n");
	printf("===========================================================================\n");
	printf("\tBEAST : Blast-furnace Evaluator Analyser & Simulation Toolkit\n");
	printf("===========================================================================\n");
	printf("\n\n");
//-----------------------Input File Reading ---------------------------------------------------------------------
	printf("\t\t---------------------------------------\n");
	printf("\t\t\tINPUT FILE READING\n");
	printf("\t\t---------------------------------------\n");
	file_read();	
//-----------------------Domain Setup ---------------------------------------------------------------------------
	printf("\t\t---------------------------------------\n");
	printf("\t\t\tDOMAIN SETUP\n");
	printf("\t\t---------------------------------------\n");	
	Domain_Setup();
//-----------------------Initialiization of Hydrodynamic props---------------------------------------------------
	printf("\t\t---------------------------------------\n");
	printf("\t\t\tINITIALIZATION \n");
	printf("\t\t---------------------------------------\n");
	initialise();   	
//-----------------------Writing a Report file for given Inputs--------------------------------------------------
	printf("\t\t---------------------------------------\n");
	printf("\t\t\tREPORT WRITING & PRE-RUN CHECKS\n");
	printf("\t\t---------------------------------------\n");
	report_data_GFS();
	prerun_checks();
/*================================================================================================================*/
	printf("------------------------------------------------------------------------------------------------\n");
	printf("+++++++++++++++++++++++++++++++++++++ B.E.A.S.T. SOLVER ++++++++++++++++++++++++++++++++++++++++\n");
	printf("------------------------------------------------------------------------------------------------\n");
	//--------------------GAS-FINES COMPUTATION LOOP---------------------------------------------------------------
	restart_read = false;
	for (GFS_iter = 0; GFS_iter <= GFS_Max; GFS_iter++) 
	{
		//--------------------READ RESTART FILE--------------------------------------------------------------------
		if ( (restart_flag) && (!restart_read) )
		{
			Read_Restart_File();
			restart_read = true;
		}
		//--------------------PRINT STATUS LABEL--------------------------------------------------------------------
		Print_Status_Label();
		//--------------------TOTAL PHASE CONTROL CENTER------------------------------------------------------------
		printf("---------------------------------------------------------------------------\n");
		printf("\t\t TOTAL PHASE CONTROL CENTER OF B.E.A.S.T.\n");
		printf("---------------------------------------------------------------------------\n");		
		Determine_Phase_Status();		
		/*---------------------------------------------------------------------------------------------------------*/
		/*--------------------------------------SOLVERS FOR EACH PHASES--------------------------------------------*/
		/*---------------------------------------------------------------------------------------------------------*/
		
 		//--------------------SOLID PHASE COMPUTATION---------------------------------------------------------------
		if (status_solid_phase)
		{
			printf("---------------------------------------------------------------------------\n");
			printf("\t\t SOLID COMPUTATION\n");
			printf("---------------------------------------------------------------------------\n");
			if(!init_solid_phase)
			{
				initialise_solid_phase();
				enforce_domain_physics();
			}
			Solid_Phase_Solver();	
			enforce_domain_physics();
			status_solid_phase	=	false;
		}
		//--------------------GAS PHASE COMPUTATION-----------------------------------------------------------------
		if (status_gas_phase)
		{
			printf("---------------------------------------------------------------------------\n");
			printf("\t\t GAS COMPUTATION\n");
			printf("---------------------------------------------------------------------------\n");
			if(!init_gas_phase)
			{
				initialise_gas_phase();
				enforce_domain_physics();
			}
			Gas_Phase_Solver();
			enforce_domain_physics();
		}
		//-------------------LIQUID PHASE COMPUTATION----------------------------------------------------------------
  		if (status_liquid_phase)
		{		
			printf("---------------------------------------------------------------------------\n");
			printf("\t\t LIQUID COMPUTATION\n");
			printf("---------------------------------------------------------------------------\n");
			if(!init_liquid_phase)
			{	
				initialise_liquid_phase();
				enforce_domain_physics();
			}
			printf("Void Definition Module.........Initiated\n");	
			Void_Definition_Module();	
			printf("Void Definition Module.........Complete\n");		
			Liquid_Phase_Solver();
			enforce_domain_physics();
			status_liquid_phase	=	false;	
		}
		//--------------------FINES PHASE COMPUTATION-----------------------------------------------------------------
 		if (status_fines_phase)
		{
			printf("---------------------------------------------------------------------------\n");
			printf("\t\t FINES COMPUTATION\n");
			printf("---------------------------------------------------------------------------\n");
			if(!init_fines_phase)
			{	
				steady_state_gas	=	false;	
				initialise_fines_phase();
				enforce_domain_physics();
			}
			Fines_Phase_Solver();			
			enforce_domain_physics();
		}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if(GFS_iter%10 == 0)
		{
			printf("---------------------------------------------------------------------------\n");
			printf("State of the simulation: \n");
			printf("---------------------------------------------------------------------------\n");
			printf("\tGas Flow Pressure Residual............:: Value = %lf \n",Residue_P);
			printf("\tGas Flow Velocity Residual............:: Value = %lf \n",Residue_V);
			printf("\tGas Flow Residual.....................:: Value = %lf \n",Residue);
			printf("\tMass_in...............................:: Value = %e kg/s\n",mass_in_fine);
			printf("\tMass_out..............................:: Value = %e kg/s\n",mass_out_fine);
			printf("\tMass_balance_fine.....................:: Value = %e \n",mass_fine);
			printf("\tPrevious Raceway Size.................:: Value = %e m2 \n",Area_RC_Previous);
			printf("\tCurrent Raceway Size..................:: Value = %e m2 \n",Area_RC_display);
			printf("\tRaceway Growth........................:: Value = %lf percent \n",Growth_RC_Area);
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			printf("------------------------------------------------------------\n");
			printf("Recording convergence and performance\n");
			save_convergence_data();			
			printf("------------------------------------------------------------\n");
			//-----------------------Writing a Restart File-----------------------------------------------------------
			Write_Restart_File();
		}
		//-----------------------Break Loop if COnvergence achieved---------------------------------------------------
		if (steady_state_fines && steady_state_gas) goto AA;
	} 
/*===================================================================================================================*/
/*--------------Writing Output Files--------------------*/	
	AA:	printf("End of run: Reason for run::\n");
	for (int i = 0; i<=argc; i++)
	{
		printf("%s\t",argv[i]);
	}
	return 1;
} 
// ======================================================================================================================
// ======================================================================================================================
// ======================================================================================================================
void Print_Status_Label()
{
	solution_time = solution_time + GFS_iter*dt1;
	printf("\n\n");
	printf("===========================================================================\n");
	printf("Iteration of Solver::%d || Iteration of Solid Update::%d\n",GFS_iter, update_count);
	printf("Solution time of the Solver=%lf || \tTime step = %e\n",solution_time, dt1);
	printf("===========================================================================\n");
	printf("Case Info:\n");
	printf("\tDensity = %lf\n",rho_s);
	printf("\tGas flow rate = %lf lpm |Fines Flux = %lf kg/m2s\n",gas_flowrate, Gf);			
	printf("\tParticle Size = %e |Fines Size = %e\n", dp, dpf);
	printf("\tFines loading = %e |\n", loading);
	printf("===========================================================================\n");
	//--------------------COMPUTATION OF COUPLING FORCES-------------------------------------
	printf("---------------------------------------------------------------------------\n");
	printf("\t\t COUPLING FORCES COMPUTATION\n");
	printf("---------------------------------------------------------------------------\n");
	Forcefield_Computation();
	printf("\tMaximum Gas-Particle Force............:: Value = %lf N\n",max_Fgp);
	printf("\tMaximum Gas-Liquid Force..............:: Value = %lf N\n",max_Fgl);
	printf("\tMaximum Gas-Fines Force...............:: Value = %lf N\n",max_Fgf);
	printf("\tMaximum Particle-Liquid Force.........:: Value = %lf N\n",max_Fpl);
	printf("\tMaximum Particle-Fines Force..........:: Value = %lf N\n",max_Fpf);
	printf("\tMaximum Liquid-Fines Force............:: Value = %lf N\n",max_Flf);
}
// ======================================================================================================================
void Determine_Phase_Status()
{
	if (GFS_iter > 0)
	{
		if( Residue <= (100 * (tol/tol))															)	status_solid_phase	= true;	
		if( gas_flowrate != 0.00																	)	status_gas_phase	= true;	
		if( Flag_Raceway_Computation == 2 && Residue <= (100 * (tol/tol)) && NUM_OF_ROTAMETERS > 0	)	status_liquid_phase	= true;	
		if( GFS_iter > 10				  && steady_state_gas 	          && Gf > 0.00				)	status_fines_phase	= true;			
	}		
	// 	-----------------------------------------------------------------------------------------------------------------
	if (status_solid_phase)	
	{
		printf("\tSOLID PHASE ........................ACTIVE\n");
	}
	else
	{
		printf("\tSOLID PHASE ......................INACTIVE\n");
	}
	// 	-----------------------------------------------------------------------------------------------------------------
	if (status_gas_phase)
	{
		printf("\tGAS PHASE ..........................ACTIVE\n");
	}
	else
	{
		printf("\tGAS PHASE ........................INACTIVE\n");
	}
	// 	-----------------------------------------------------------------------------------------------------------------
	if (status_liquid_phase)	
	{
		printf("\tLIQUID PHASE ........................ACTIVE\n");
	}
	else
	{
		printf("\tLIQUID PHASE .....................INACTIVE\n");
	}
	// 	-----------------------------------------------------------------------------------------------------------------
	if (status_fines_phase)
	{
		printf("\tFINES PHASE ........................ACTIVE\n");
	}
	else
	{
		printf("\tFINES PHASE ......................INACTIVE\n");
	}
}
// ======================================================================================================================
void Read_Restart_File()
{
	FILE *restart;
	restart = fopen("1.0.Binary_Restart_File.rst","r");
	int ssp, slp, sgp, sfp;
	int isp, ilp, igp, ifp;
	fscanf(restart,"%d\n",&GFS_iter);
	fscanf(restart,"%d\n",&ssp);
	fscanf(restart,"%d\n",&slp);
	fscanf(restart,"%d\n",&sgp);
	fscanf(restart,"%d\n",&sfp);
	
	fscanf(restart,"%d\n",&isp);
	fscanf(restart,"%d\n",&ilp);
	fscanf(restart,"%d\n",&igp);
	fscanf(restart,"%d\n",&ifp);
	fscanf(restart,"%d\n",&update_count);
	fscanf(restart,"%d\n",&nstep);
	if (ssp == 1)	status_solid_phase = true;
	if (slp == 1)	status_liquid_phase = true;
	if (sgp == 1)	status_gas_phase = true;
	if (sfp == 1)	status_fines_phase = true;
	
	if (isp == 1)	init_solid_phase = true;
	if (ilp == 1)	init_liquid_phase = true;
	if (igp == 1)	init_gas_phase = true;
	if (ifp == 1)	init_fines_phase = true;

	for (int j = 1; j <= JMAX-1; j++)
	{
		for (int i = 0; i <= IMAX-1; i++)
		{	
			fscanf(restart,"%lf %lf\n",&u[i][j],&v[i][j]);
		}
	}
		
	for (int j = 1; j <= JMAX-1; j++)
	{
		for (int i = 0; i <= IMAX-1; i++)
		{	
			fscanf(restart,"%lf %lf\n",&uf[i][j],&vf[i][j]);
		}
	}	
	
	for (int j = 1; j <= JMAX-1; j++)
	{
		for (int i = 1; i <= IMAX-1; i++)
		{
			fscanf(restart,"%lf %lf %lf %lf %lf %lf\n",&p[i][j], &ke[i][j], &de[i][j], &vfrac[i][j],&epfd[i][j],&epfs[i][j]);
		}
	}

	fscanf(restart,"%lf %lf\n",&d_max, &STATICtime);
	fscanf(restart,"%lf %lf %lf %lf\n", &Kn, &Kt, &KnWall, &KtWall);		
	fscanf(restart,"%lf %lf\n",&coeff_restitution,&coeff_friction);
	for (int i = 0; i < NUM; i++)
	{
		fscanf(restart,"%lf %le %lf %lf %lf %lf %lf %lf\n",&mass[i], &MMOI[i], &xp[i], &yp[i], &d[i], &up[i], &vp[i], &omega[i]);
	}
	fclose(restart);		
}
// ======================================================================================================================
void Write_Restart_File()
{
	FILE *restart;	
	restart = fopen("1.0.Binary_Restart_File.rst","w");
	
	fprintf(restart,"%d\n",GFS_iter);
	fprintf(restart,"%d\n",status_solid_phase);
	fprintf(restart,"%d\n",status_liquid_phase);
	fprintf(restart,"%d\n",status_gas_phase);
	fprintf(restart,"%d\n",status_fines_phase);
	
	fprintf(restart,"%d\n",init_solid_phase);
	fprintf(restart,"%d\n",init_liquid_phase);
	fprintf(restart,"%d\n",init_gas_phase);
	fprintf(restart,"%d\n",init_fines_phase);
	fprintf(restart,"%d\n",update_count);
	fprintf(restart,"%d\n",nstep);
	for (int j = 1; j <= JMAX-1; j++)
	{
		for (int i = 0; i <= IMAX-1; i++)
		{	
			fprintf(restart,"%e %e\n",u[i][j],v[i][j]);
		}
	}
		
	for (int j = 1; j <= JMAX-1; j++)
	{
		for (int i = 0; i <= IMAX-1; i++)
		{	
			fprintf(restart,"%e %e\n",uf[i][j],vf[i][j]);
		}
	}	

	for (int j = 1; j <= JMAX-1; j++)
	{
		for (int i = 1; i <= IMAX-1; i++)
		{
			fprintf(restart,"%e %e %e %e %e %e\n",p[i][j], ke[i][j], de[i][j], vfrac[i][j],epfd[i][j],epfs[i][j]);
		}
	}
	
	fprintf(restart,"%lf %lf\n",d_max, STATICtime);
	fprintf(restart,"%lf %lf %lf %lf\n", Kn, Kt, KnWall, KtWall);		
	fprintf(restart,"%lf %lf\n",coeff_restitution,coeff_friction);
	for (int i = 0; i < NUM; i++)
	{
		fprintf(restart,"%lf %e %lf %lf %lf %lf %lf %lf\n",mass[i], MMOI[i], xp[i], yp[i], d[i], up[i], vp[i], omega[i]);
	}
	fclose(restart);
}