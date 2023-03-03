#include "D.1.DEM_Module.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void initialise_solid_phase()
{
	init_solid_phase	=	true;
	int i, j;
	double fai=28.213*3.14/180;
	double delta=fai*(PI/180.0);
	double eps_c=0.5*(0.5*PI - delta);
	double mass_min;
	
	for (i=0;i<IMAX;i++)
	{
		for (j=0;j<JMAX;j++)
		{
			sigma_xx[i][j]=sigma[i][j]*(1.0+(sin(fai)*cos(2.0*sai[i][j])));
			sigma_yy[i][j]=sigma[i][j]*(1.0-(sin(fai)*cos(2.0*sai[i][j])));
			tau_xy[i][j]=sigma[i][j]*sin(fai)*sin(2.0*sai[i][j]);

			sigma_alpha[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*(sai[i][j]-eps_c))+tau_xy[i][j]*sin(2.0*(sai[i][j]-eps_c));
			tau_alpha[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*(sai[i][j]-eps_c))+tau_xy[i][j]*cos(2.0*(sai[i][j]-eps_c));
			
			sigma_beta[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*(sai[i][j]-eps_c))+tau_xy[i][j]*sin(2.0*(sai[i][j]-eps_c));
			tau_beta[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*(sai[i][j]-eps_c))+tau_xy[i][j]*cos(2.0*(sai[i][j]-eps_c));
		} 
	}
	
	for (j=1;j<=JMAX-1;j++)
	{
		for (i=1;i<=IMAX-1;i++)
		{
			sigma_eff[i][j] = sigma_alpha[i][j];
		}
	} 

	for (i = 0; i <= IMAX; i++)
	{
		for (j = 0; j <= JMAX; j++)
		{ 
			us[i][j]	=	0.00;
			vs[i][j]	=	0.00;
		}
	}
	
	// Initializing the key variables------------------------------------------------------------------------
	for (i = 0; i < NUM; i++)
	{
		xp[i]		=	0.00;
		yp[i]		=	0.00;
		up[i]		=	0.00;   
		vp[i]		=	0.00;  
		omega[i]	=	0.00;
	}
	
	 for (i = 0; i < NUM; i++)
	{
		d[i] = dp;
	}
/* 	for (i = (0.5*NUM); i < NUM; i++)
	{
		d[i] = 0.75 * dp;
	} */
	
	for (i = 0; i < NUM; i++)
	{	
		mass[i]	=	rho_s*(4.0/3.0)*PI*pow(d[i]/2.0,3.0);
		MMOI[i]	=	(2.0/5.0)*mass[i]*pow(d[i]/2.0,2.0); 
	} 
	
	Kt		=	Kn;
	KnWall	=	Kn;
	KtWall	=	Kt;

	
	mass_min= -1000.00;
	mass_max= 0.00;
	d_max	= 0.00;
	for (i = 0; i < NUM; i++)
	{
		mass_min	=	min(mass_min, mass[i]);
		mass_max	=	max(mass_max, mass[i]);
		d_max		=	max(d_max, d[i]);
	}
	
	CnVAL	=	0.00;
	if (coeff_restitution != 0.00)
	{
		CnVAL	=	2.0*(-log(coeff_restitution)/( sqrt(pow(PI,2.0) + pow(log(coeff_restitution),2.0))))*sqrt(mass_max*Kn/2.0);
	}
	CtVAL	=	CnVAL; 
	CnWall	=	CnVAL;
	CtWall	=	CtVAL;

	if (time_step>0.1*sqrt(mass_min/Kn))
	{
		time_step=0.1*sqrt(mass_min/Kn);
	}
	STATICtime	=	time_step;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Solid_Phase_Solver()
{
	/*------Local variables declaration--------------------------------------------------------*/
	int i, j; 
	double search_radius, radi;	
	/*------Increment the update counter-------------------------------------------------------*/
	update_count++;
	/*------Determine the reference raceway radius based on correlation------------------------*/
 	if (Condition > 0)
	{
		raceway_radius = 0.5*increase_corr(Lx,Ly-h,vin,D_T,dp,phi_s,rho_s,Domain_Voidage);
	}
	if (Condition < 0)
	{
		raceway_radius = 0.5*decrease_corr(Lx,Ly-h,vin,D_T,dp,phi_s,rho_s,Domain_Voidage);
	}		
	/*------Neglecting the presence and effect of the raceway----------------------------------*/
	if (Flag_Raceway_Computation == -1)
	{
		printf("No Raceway modelling in effect\n");
		printf("---------------------------------------------------\n\n");
		for(i = 0; i < IMAX; i++)
		{
			for(j = 0; j < JMAX; j++)
			{
				vfrac[i][j]=Domain_Voidage;
			}
		}
		raceway_radius = 0.00;
	}
	/*------Using the correlation to model the raceway-----------------------------------------*/
	if (Flag_Raceway_Computation == 0)
	{
		printf("\tRaceway Evolution through Correlation\n");
		if (update_count > 10)
		{
			update_count	= 10;	
		}
		
		/* search_radius = (0.1 * update_count) * raceway_radius; 
		Raceway_Center_Y = Tuyere_height[1];
		Raceway_Center_X = Tuyere_protrusion[1] + (0.1 * update_count) * raceway_radius;		 */
		
		search_radius = raceway_radius; 
		Raceway_Center_Y = Tuyere_height[1];
		Raceway_Center_X = Tuyere_protrusion[1] + raceway_radius;		
		
		printf("\tUpdate Count.......................:: Value = %d\n", update_count);
		printf("\tVelocity Condition.................:: Value = %d\n", Condition);
		printf("\tSearch Radius......................:: Value = %lf m\n", search_radius);
		printf("\tRaceway location X.................:: Value = %lf m\n", Raceway_Center_X);
		printf("\tRaceway location Y.................:: Value = %lf m\n", Raceway_Center_Y);
		
		for(i = 0; i < IMAX; i++)
		{
			for(j = 0; j < JMAX; j++)
			{
				radi = sqrt(sqr(x[i] - Raceway_Center_X) + sqr(y[j] - Raceway_Center_Y));
				if (radi < search_radius)
				{
					vfrac[i][j]	=	Raceway_Voidage;
				}
			}
		}		
	}
	/*------Using the isostress and stress profile to model the raceway------------------------*/
	if (Flag_Raceway_Computation == 1)
	{
		printf("\tRaceway Evolution through Iso-stress modeling\n");
		search_radius		= 0.1 * update_count*raceway_radius;
		if (update_count > 10)
		{
			update_count	= 10;	
		}
		Raceway_Center_Y	= Tuyere_height[1];
		Raceway_Center_X	= Tuyere_protrusion[1] + (0.1 * update_count) * raceway_radius;		
		
		printf("\tUpdate Count.......................:: Value = %d\n", update_count);
		printf("\tIsostress Value....................:: Value = %lf N/m^2\n", iso_stress);
		printf("\tSearch Radius......................:: Value = %lf m\n", search_radius);
		printf("\tRaceway location X.................:: Value = %lf m\n", Raceway_Center_X);
		printf("\tRaceway location Y.................:: Value = %lf m\n", Raceway_Center_Y);
		
		for(i = 1; i < IMAX; i++)
		{
			for(j = 1; j < JMAX; j++)
			{ 
				if (State[i][j] == 1) 
				{
					sigma_eff[i][j]	=	sigma_alpha[i][j];
				}
				else if (State[i][j] == -1) 
				{
					sigma_eff[i][j]	=	iso_stress;
				}
				else
				{
					sigma_eff[i][j] =	sigma_alpha[i][j] - p[i][j];
				}
			}
		}	
		for(i=1;i<IMAX;i++)
		{
			for(j=1;j<JMAX;j++)
			{
				if (State[i][j]!=-1) 
				{
					radi = sqrt(sqr(x[i]-Raceway_Center_X) + sqr(y[j]-Raceway_Center_Y));
					if (radi < search_radius)
					{
						if(sigma_eff[i][j]<iso_stress)
						{ 
 							if( (vfrac[i-1][j]   != Domain_Voidage)		// WEST
							||  (vfrac[i-1][j-1] != Domain_Voidage)		// SOUTH-WEST
							||  (vfrac[i][j-1]   != Domain_Voidage)		// SOUTH
							||  (vfrac[i][j-1]   != Domain_Voidage)		// SOUTH
							||  (vfrac[i+1][j-1] != Domain_Voidage)		// SOUTH-EAST
							||  (vfrac[i+1][j]   != Domain_Voidage)		// EAST
							||  (vfrac[i+1][j+1] != Domain_Voidage)		// NORTH-EAST
							||  (vfrac[i][j+1]   != Domain_Voidage)		// NORTH
							||  (vfrac[i-1][j+1] != Domain_Voidage) )	// NORTH-WEST
							{
								vfrac[i][j] = Raceway_Voidage;
								//vfrac[i][j] = 1.0 - 0.5 * (sigma_eff[i][j]/iso_stress);		// Should be the implementation but it is bugged..
							} 
						}
					}
					else
					{
						vfrac[i][j] = vfrac[i][j];
					}
				}				
			}	
		}	
	}
	
	/*------Using DEM to model the solid phase and determine the raceway-----------------------*/
	if(Flag_Raceway_Computation == 2)
	{
		printf("No Raceway modelling in effect. DEM determination.\n");
		printf("---------------------------------------------------\n\n");
		if (update_count ==1)
		{
			DEM_Module();
		}
		else
		{
			DEM_Module();
		}
	}
	
	/*------Using an experimental tracing of the raceway for the simulation--------------------*/
	if (Flag_Raceway_Computation == 3)
	{
		printf("Raceway modelling through experimental tracing\n");
		printf("---------------------------------------------------\n\n");
		FILE* gp;	
		gp=fopen("1.3.Raceway_Map.inp","r");
		
		for (i = 1; i <= IMAX; i++)
		{
			for (j = 1; j <= JMAX; j++)
			{
				fscanf(gp,"%lf ",&vfrac[i][j]);
			}
			fscanf(gp,"\n");
		}
		fclose(gp);
	}
	/*------Measuring the growth of the raceway------------------------------------------------*/
	Area_RC = 0.00;
	for (i = 0; i <= IMAX; i++)
	{
		for (j = 0;j <= JMAX;j++)
		{
			if (vfrac[i][j] == Raceway_Voidage)
			{				
				Area_RC = Area_RC + (dx[i]*dy[j]);
			}
		}
	}
	Area_RC_display		= Area_RC;
	Growth_RC_Area		= (Area_RC - Area_RC_Previous) * 100 /Area_RC;
	Area_RC_Previous	= Area_RC;
}