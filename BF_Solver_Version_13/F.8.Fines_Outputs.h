void fines_convergence_data()
{
	FILE *fp;
	fp = fopen("5.2.Fines_Convergence.dat","a");
	if(Fines_Convergence_count==1){

		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Residue\",\"Mass of fines IN\",\"Mass of fines OUT\",\"Mass balance of fines\"\n");
	}
		fprintf(fp,"%d %e %e %e %e\n",Fines_Convergence_count,Residue,mass_in_fine,mass_out_fine,mass_fine);
	fclose(fp);
	Fines_Convergence_count++;
}
//------------------------------------------------------------------------------
void save_output_fines_transient()
{
	double xx, yy;
	double Ug_interstitial,Vg_interstitial;
	double Uf_interstitial,Vf_interstitial, speed, vorticity_fines;
	double Void_Fraction;
	double Fgf_x, Fgf_y, Fgf;
	double Fpf_x, Fpf_y, Fpf;
	double F_gravity;
	double Static_Holdup, Dynamic_Holdup, Total_Holdup;
	double Powder_Pressure;
	FILE *fp;
//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	fp=fopen("3.1.Fines_Flow_Transient.dat","a");
	if(Fines_Transient_output_count == 1)
	{
		fprintf(fp,"TITLE = Properties at Final state\n");
		fprintf(fp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Uf(m/s)\",\"Vf(m/s)\",\"Speed_f(m/s)\",\"Vorticity_f\",\"Powder Pressure (Pa)\",\"Void Fraction\", \"F_gas-fines X\",\"F_gas-fines Y\",\"F_gas-fines\",\"F_particle-fines X\",\"F_particle-fines Y\",\"F_particle-fines\",\"F_gravity\",\"Static Holdup\",\"Dynamic Holdup\",\"Total Holdup\",\"State\"\n");
	}
	fprintf(fp, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",GFS_iter, IMAX-1, JMAX-1);	
	yy = 0.5*dy[1];
	for (int j = 1; j <= JMAX-1; j++)
	{
		xx  =0.5*dx[1];
		for (int i = 1; i <= IMAX-1; i++)
		{
			// INTERSTITIAL VELOCITY COMPUTATION		
			Ug_interstitial 		= 0.5*(u[i][j]+u[i-1][j]);
			Vg_interstitial 		= 0.5*(v[i][j]+v[i][j-1]);
			Uf_interstitial	= 0.5*(uf[i][j] + uf[i-1][j]);
			Vf_interstitial	= 0.5*(vf[i][j] + vf[i][j-1]);
			speed			= sqrt( Uf_interstitial * Uf_interstitial + Vf_interstitial * Vf_interstitial);
			vorticity_fines	= vorticity_f[i][j];
			// SCALAR COMPUTATIONS
			Void_Fraction	=	void_frac(i,j);		
			Static_Holdup	=	epfs[i][j];
			Dynamic_Holdup	=	epfd[i][j];
			Total_Holdup	=	Static_Holdup + Dynamic_Holdup;
			Powder_Pressure	=	(1/150)*exp(150*(Total_Holdup - 0.6));
			// FORCE COMPUTATION			
			Fgf_x			=	Fg_f(i,j) * (Ug_interstitial - Uf_interstitial);
			Fgf_y			=	Fg_f(i,j) * (Vg_interstitial - Vf_interstitial);
			Fgf				=	sqrt(sqr(Fgf_x) + sqr(Fgf_y));
			
			Fpf_x 			= 	Fp_f(i,j) * Uf_interstitial;
			Fpf_y 			= 	Fp_f(i,j) * Vf_interstitial;
			Fpf				=	sqrt(sqr(Fpf_x) + sqr(Fpf_y));
			
			F_gravity		=	rho_f*gravity*epfd[i][j];
			// CONSIDERATION OF GEOMETRY
			if (State[i][j] == -1)
			{			
				Uf_interstitial 	= 0.00;
				Vf_interstitial 	= 0.00;
				speed		 		= 0.00;
				vorticity_f[i][j] 	= 0.00;
				
				Fgf_x				= 0.00;
				Fgf_y				= 0.00;
				Fgf					= 0.00;
				
				Fpf_x				= 0.00;
				Fpf_y				= 0.00;
				Fpf					= 0.00;
				
				F_gravity			= 0.00;
			}

			fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",
			xx, 
			yy, 
			Uf_interstitial, 
			Vf_interstitial, 
			speed, 
			vorticity_fines,
			Powder_Pressure, 
			Void_Fraction, 
			Fgf_x,
			Fgf_y,
			Fgf, 
			Fpf_x, 
			Fpf_y, 
			Fpf, 
			F_gravity,  
			Static_Holdup, 
			Dynamic_Holdup, 
			Total_Holdup, 
			State[i][j]);
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}	
	fclose(fp);
	// CRITICAL UPDATE OF count	
	Fines_Transient_output_count++;
}
//------------------------------------------------------------------------------
void save_output_fines_steady_state()
{
	double xx, yy;
	double Ug_interstitial,Vg_interstitial;
	double Uf_interstitial,Vf_interstitial, speed, vorticity_fines;
	double Void_Fraction;
	double Fgf_x, Fgf_y, Fgf;
	double Fpf_x, Fpf_y, Fpf;
	double F_gravity;
	double Static_Holdup, Dynamic_Holdup, Total_Holdup;
	double Powder_Pressure;
	FILE *fp;
	//-------------------------HYDRODYNAMIC PROPERTIES AT STEADY STATE
	fp=fopen("3.2.Fines_Flow_Steady_State.dat","wt");
	fprintf(fp,"TITLE = Properties at Final state\n");
	fprintf(fp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Uf(m/s)\",\"Vf(m/s)\",\"Speed_f(m/s)\",\"Vorticity_f\",\"Powder Pressure (Pa)\",\"Void Fraction\", \"F_gas-fines X\",\"F_gas-fines Y\",\"F_gas-fines\",\"F_particle-fines X\",\"F_particle-fines Y\",\"F_particle-fines\",\"F_gravity\",\"Static Holdup\",\"Dynamic Holdup\",\"Total Holdup\",\"State\"\n");
	fprintf(fp, "ZONE T = \"%lf\", I=%d, J=%d F=POINT\n",solution_time, IMAX-1, JMAX-1);	
	yy = 0.5*dy[1];
	for (int j = 1; j <= JMAX-1; j++)
	{
		xx  =0.5*dx[1];
		for (int i = 1; i <= IMAX-1; i++)
		{
			// INTERSTITIAL VELOCITY COMPUTATION		
			Ug_interstitial = 0.5*(u[i][j]+u[i-1][j]);
			Vg_interstitial = 0.5*(v[i][j]+v[i][j-1]);
			Uf_interstitial	= 0.5*(uf[i][j] + uf[i-1][j]);
			Vf_interstitial	= 0.5*(vf[i][j] + vf[i][j-1]);
			speed			= sqrt( Uf_interstitial * Uf_interstitial + Vf_interstitial * Vf_interstitial);
			vorticity_fines	= vorticity_f[i][j];
			// SCALAR COMPUTATIONS
			Void_Fraction	=	void_frac(i,j);		
			Static_Holdup	=	epfs[i][j];
			Dynamic_Holdup	=	epfd[i][j];
			Total_Holdup	=	Static_Holdup + Dynamic_Holdup;
			Powder_Pressure	=	(1/150)*exp(150*(Total_Holdup - 0.6));
			// FORCE COMPUTATION			
			Fgf_x			=	Fg_f(i,j) * (Ug_interstitial - Uf_interstitial);
			Fgf_y			=	Fg_f(i,j) * (Vg_interstitial - Vf_interstitial);
			Fgf				=	sqrt(sqr(Fgf_x) + sqr(Fgf_y));
			
			Fpf_x 			= 	Fp_f(i,j) * Uf_interstitial;
			Fpf_y 			= 	Fp_f(i,j) * Vf_interstitial;
			Fpf				=	sqrt(sqr(Fpf_x) + sqr(Fpf_y));
			
			F_gravity		=	rho_f*gravity*epfd[i][j];
			// CONSIDERATION OF GEOMETRY
			if (State[i][j] == -1)
			{			
				Uf_interstitial 	= 0.00;
				Vf_interstitial 	= 0.00;
				speed		 		= 0.00;
				vorticity_f[i][j] 	= 0.00;
				
				Fgf_x				= 0.00;
				Fgf_y				= 0.00;
				Fgf					= 0.00;
				
				Fpf_x				= 0.00;
				Fpf_y				= 0.00;
				Fpf					= 0.00;
				
				F_gravity			= 0.00;
			}

			fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",
			xx, 
			yy, 
			Uf_interstitial, 
			Vf_interstitial, 
			speed, 
			vorticity_fines,
			Powder_Pressure, 
			Void_Fraction, 
			Fgf_x,
			Fgf_y,
			Fgf, 
			Fpf_x, 
			Fpf_y, 
			Fpf, 
			F_gravity,  
			Static_Holdup, 
			Dynamic_Holdup, 
			Total_Holdup, 
			State[i][j]);
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}		
	fclose(fp);
	
	//-------------------------FINESS VELOCITY DETAILS FOR CFD SOLVER
	fp=fopen("2.4.Fines_Velocity_details.int","wt");
	for (int j = 0; j <= JMAX; j++)
	{
		for (int i = 0; i <= IMAX; i++)
		{
			if(i == 0)
			{
				Uf_interstitial = u_left;
			}
			else if(i == IMAX)
			{
				Uf_interstitial = u_right;
			}
			else if(j == 0)
			{
				Vf_interstitial = v_bottom;
			}
			else if(j == JMAX)
			{
				Vf_interstitial = v_top;
			}
			else
			{
				Uf_interstitial 		= 0.5*(uf[i][j] + uf[i-1][j]);
				Vf_interstitial 		= 0.5*(vf[i][j] + vf[i][j-1]);
			}
			fprintf(fp,"%e %e\n",	Uf_interstitial, Vf_interstitial);	
		}
	}		
	fclose(fp);	
}