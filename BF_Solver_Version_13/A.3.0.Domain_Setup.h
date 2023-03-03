#include "A.3.1.Grid.h"
#include "A.3.2.State_Designation.h"
#include "A.3.3.Void_Fraction.h"
#include "A.3.4.Stress.h"

void Domain_Setup()
{   
	double dtt1, dtt2;
	printf("Domain Setup Initiated ...\n");
	//-----------------------Initial Calculations-------------------------------------------
	printf("\t Initial Calculations ... In Progress\n");
	D_T					=	Tuyere_opening[1];
	h					=	Tuyere_height[1];
	r					=	Tuyere_protrusion[1];
	vin					=	gas_flowrate/(1000*D_T*Lz*60);
	v_super				=	gas_flowrate/(1000*Lx*Lz*60);
	Gg					=	vin*rho_g;
	loading				=	(Gf/(vin*rho_f));
	lm					= 	Tuyere_opening[1] * 0.5; 					// Ref.1
	//------------------Grid Generation------------------------------------------------------------------------------  
	printf("\t Grid Generation ... In Progress\n");
	if(grid_type==0)
	{
		printf("\t \t Grid Resolution :: IMAX=%d x JMAX=%d \n",IMAX, JMAX);
		printf("\t \t Grid Type :: UNIFORM GRID\n");
		grid_uniform();
	}
	else if(grid_type==1)
	{
		grid_non_uniform();	
	}
	//-----------------------Determining the state of the grid points-------------------------------------------------   
	printf("\t Cell State Determination ... In Progress\n");
	state_calc();
	//-----------------------Determining the internal stress of the bed-----------------------------------------------   
	printf("\t Bed Stress Determination ... In Progress\n");
	stress_calc();
	interp_printout();
	FILE* ip5;
	ip5=fopen("1.4.Stress_Distribution.inp","r");
		for (int j = 1; j < JMAX; j++)
		{ 
			for (int i = 1; i < IMAX; i++)
			{
				fscanf(ip5,"%lf %lf\n",&sigma[i][j],&sai[i][j]);
			}
		}
	fclose(ip5);	
	//-----------------------Determining the void fraction of the bed-----------------------------------------------  
	printf("\t Void Fraction Determination ... In Progress\n");
	void_fraction();
	//-----------------------Time Step & CFL------------------------------------------------------------------------	
	printf("\t Time Step & CFL Computation ... In Progress\n");
	Re					=	rho_g*v_super*dp/mu_g;
	dtt1				=	0.8*(0.5*sqr(dx[1]*dy[1]/(Ly*Ly))*Re)/(sqr(dx[1]/Ly) + sqr(dy[1]/Ly));
	dtt2				=	0.8*min(dx[1]/Ly,dy[1]/Ly);
	dt1					=	0.025*min(dtt1,dtt2)*Ly/vin;

	Bed_height			=	(NUM*dp)/(Lx);
	min_x				=	dx[0];
	max_x				=	dx[0];
	for(int i = 1; i < IMAX; i++)
	{
		if(min_x > dx[i])
		{
			min_x	=	dx[i];  
		}			
		if(max_x<dx[i])
		{
			max_x	=	dx[i];       
		}
	}
	CFL_x		=	vin*dt1/max_x;
	
	min_y		=	dy[0];
	max_y		=	dy[0];
	for(int j = 1; j < JMAX; j++)
	{
		if(min_y > dy[j])
		{
			min_y = dy[j];   
		}
		if(max_y < dy[j])
		{
			max_y = dy[j];       
		}
	}
	CFL_y		=	vin*dt1/max_y;	
	//---------------------------------------------------------------------
	printf("Domain Setup ...............................Complete\n");	
	printf("\n");	
	printf("\n");	
}