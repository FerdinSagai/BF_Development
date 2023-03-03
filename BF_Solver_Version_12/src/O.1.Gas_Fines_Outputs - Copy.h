//------------------------------------------------------------------------------
void save_output_gas_transient()
{
	double xx, yy;
	double Ug_interstitial,Vg_interstitial, speed, vorticity_gas;
	FILE *fp;
//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	fp=fopen("3.1.Gas_Flow_Transient.dat","a");
	if(Transient_output_count == 1)
	{
		fprintf(fp,"TITLE = Properties at Final state\n");
		fprintf(fp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Ug\",\"Vg\",\"Speed\",\"Vorticity\",\"Pressure (Pa)\",\"Void Fraction\",\"Uf\",\"Vf\",\"Speed_f\",\"Vorticity_f\",\"Powder Pressure (Pa)\",\"F_gas-particles X\",\"F_gas-particles Y\",\"F_gas-particles\", \"F_gas-fines X\",\"F_gas-fines Y\",\"F_gas-fines\",\"F_particle-fines X\",\"F_particle-fines Y\",\"F_particle-fines\",\"F_gravity\",\"KE\",\"DE\",\"Static Holdup\",\"Dynamic Holdup\",\"Total Holdup\",\"Stress\",\"State\"\n");
	}
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
			speed		= sqrt(Ug_interstitial*Ug_interstitial + Vg_interstitial*Vg_interstitial);
			vorticity_gas	= vorticity_g[i][j];
			
			Uf_interstitial	= 0.5*(uf[i][j] + uf[i-1][j]);
			Vf_interstitial	= 0.5*(vf[i][j] + vf[i][j-1]);
			speed_f		= sqrt( Uf_interstitial * Uf_interstitial + Vf_interstitial * Vf_interstitial);
			vorticity_fines	= vorticity_f[i][j];
							
			// SCALAR COMPUTATIONS
			Pressure	=	p[i][j];							// UNEXPLAINABLE CORRECTION
			Void_Fraction	=	void_frac(i,j);		
			Kinetic_Energy	=	ke[i][j];
			Dissipation_KE	=	de[i][j];
			Stress		=	sigma_eff[i][j];
			Static_Holdup	=	epfs[i][j];
			Dynamic_Holdup	=	epfd[i][j];
			Total_Holdup	=	Static_Holdup + Dynamic_Holdup;
			Powder_Pressure	=	(1/150)*exp(150*(Total_Holdup - 0.6));
			// FORCE COMPUTATION
			Fpg_x		=	Fg_p(i,j)* 0.5 * (u[i][j] + u[i-1][j]);
			Fpg_y		=	Fg_p(i,j)* 0.5 * (v[i][j] + v[i][j-1]);
			Fpg		=	sqrt(sqr(Fpg_x) +sqr(Fpg_y));
			
			Fgf_x		=	Fg_f(i,j) * (Ug_interstitial - Uf_interstitial);
			Fgf_y		=	Fg_f(i,j) * (Vg_interstitial - Vf_interstitial);
			Fgf		=	sqrt(sqr(Fgf_x) + sqr(Fgf_y));
			
			Fpf_x 		= 	Fp_f(i,j) * Uf_interstitial;
			Fpf_y 		= 	Fp_f(i,j) * Vf_interstitial;
			Fpf		=	sqrt(sqr(Fpf_x) + sqr(Fpf_y));
			
			F_gravity	=	rho_f*gravity*epfd[i][j];
			
			// CONSIDERATION OF GEOMETRY
			if (State[i][j] == -1)
			{
				Ug_interstitial 	= 0.00;
				Vg_interstitial 	= 0.00;
				speed 			= 0.00;
				vorticity_g[i][j] 	= 0.00;
				
				Uf_interstitial 	= 0.00;
				Vf_interstitial 	= 0.00;
				speed_f 		= 0.00;
				vorticity_f[i][j] 	= 0.00;
				
				Pressure		= 0.00;
				Kinetic_Energy		= 0.00;
				Dissipation_KE		= 0.00;

				Fpg_x			= 0.00;
				Fpg_y			= 0.00;
				Fpg			= 0.00;
				
				Fgf_x			= 0.00;
				Fgf_y			= 0.00;
				Fgf			= 0.00;
				
				Fpf_x			= 0.00;
				Fpf_y			= 0.00;
				Fpf			= 0.00;
				
				F_gravity		= 0.00;
			}

			fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",
			xx, 
			yy, 
			Ug_interstitial, 
			Vg_interstitial, 
			speed, 
			vorticity_gas,
			Pressure, 
			Void_Fraction, 
			Uf_interstitial, 
			Vf_interstitial, 
			speed_f, 
			vorticity_fines,
			Powder_Pressure, 
			Fpg_x,
			Fpg_y,
			Fpg, 
			Fgf_x,
			Fgf_y,
			Fgf, 
			Fpf_x, 
			Fpf_y, 
			Fpf, 
			F_gravity, 
			Kinetic_Energy, 
			Dissipation_KE,  
			Static_Holdup, 
			Dynamic_Holdup, 
			Total_Holdup, 
			Stress, 
			State[i][j]);			
	fclose(fp);
}
//------------------------------------------------------------------------------
void save_output_gas_steady_state()
{
	int i,j;
	double xx, yy;
	
	double Ug_interstitial,Vg_interstitial, speed, vorticity_gas;
	double Uf_interstitial,Vf_interstitial, speed_f, vorticity_fines;

	double Pressure, Void_Fraction, Kinetic_Energy, Dissipation_KE, Stress;
	
	double Fpg_x, Fpg_y, Fpg;
	double Fgf_x, Fgf_y, Fgf;
	double Fpf_x, Fpf_y, Fpf;
	double F_gravity;

	double Static_Holdup, Dynamic_Holdup, Total_Holdup;
	double Powder_Pressure;

	FILE *gp;
//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS

	gp=fopen("3.2.Results_Steady_State.dat","wt");
	fprintf(gp,"TITLE = Properties at Final state\n");
	fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Ug\",\"Vg\",\"Speed\",\"Vorticity\",\"Pressure (Pa)\",\"Void Fraction\",\"Uf\",\"Vf\",\"Speed_f\",\"Vorticity_f\",\"Powder Pressure (Pa)\",\"F_gas-particles X\",\"F_gas-particles Y\",\"F_gas-particles\", \"F_gas-fines X\",\"F_gas-fines Y\",\"F_gas-fines\",\"F_particle-fines X\",\"F_particle-fines Y\",\"F_particle-fines\",\"F_gravity\",\"KE\",\"DE\",\"Static Holdup\",\"Dynamic Holdup\",\"Total Holdup\",\"Stress\",\"State\"\n");
	fprintf(gp, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",Transient_output_count, IMAX-1, JMAX-1);
	yy = 0.5*dy[1];
	for (j=1;j<=JMAX-1;j++)
	{
		xx  =0.5*dx[1];
		for (i=1;i<=IMAX-1;i++)
		{
			// INTERSTITIAL VELOCITY COMPUTATION
			Ug_interstitial = 0.5*(u[i][j]+u[i-1][j]);
			Vg_interstitial = 0.5*(v[i][j]+v[i][j-1]);
			speed		= sqrt(Ug_interstitial*Ug_interstitial + Vg_interstitial*Vg_interstitial);
			vorticity_gas	= vorticity_g[i][j];
			
			Uf_interstitial	= 0.5*(uf[i][j] + uf[i-1][j]);
			Vf_interstitial	= 0.5*(vf[i][j] + vf[i][j-1]);
			speed_f		= sqrt( Uf_interstitial * Uf_interstitial + Vf_interstitial * Vf_interstitial);
			vorticity_fines	= vorticity_f[i][j];
							
			// SCALAR COMPUTATIONS
			Pressure	=	p[i][j];							// UNEXPLAINABLE CORRECTION
			Void_Fraction	=	void_frac(i,j);		
			Kinetic_Energy	=	ke[i][j];
			Dissipation_KE	=	de[i][j];
			Stress		=	sigma_eff[i][j];
			Static_Holdup	=	epfs[i][j];
			Dynamic_Holdup	=	epfd[i][j];
			Total_Holdup	=	Static_Holdup + Dynamic_Holdup;
			Powder_Pressure	=	(1/150)*exp(150*(Total_Holdup - 0.6));
			// FORCE COMPUTATION
			Fpg_x		=	Fg_p(i,j)* 0.5 * (u[i][j] + u[i-1][j]);
			Fpg_y		=	Fg_p(i,j)* 0.5 * (v[i][j] + v[i][j-1]);
			Fpg		=	sqrt(sqr(Fpg_x) +sqr(Fpg_y));
			
			Fgf_x		=	Fg_f(i,j) * (Ug_interstitial - Uf_interstitial);
			Fgf_y		=	Fg_f(i,j) * (Vg_interstitial - Vf_interstitial);
			Fgf		=	sqrt(sqr(Fgf_x) + sqr(Fgf_y));
			
			Fpf_x 		= 	Fp_f(i,j) * Uf_interstitial;
			Fpf_y 		= 	Fp_f(i,j) * Vf_interstitial;
			Fpf		=	sqrt(sqr(Fpf_x) + sqr(Fpf_y));
			
			F_gravity	=	rho_f*gravity*epfd[i][j];
			
			// CONSIDERATION OF GEOMETRY
			if (State[i][j] == -1)
			{
				Ug_interstitial 	= 0.00;
				Vg_interstitial 	= 0.00;
				speed 			= 0.00;
				vorticity_g[i][j] 	= 0.00;
				
				Uf_interstitial 	= 0.00;
				Vf_interstitial 	= 0.00;
				speed_f 		= 0.00;
				vorticity_f[i][j] 	= 0.00;
				
				Pressure		= 0.00;
				Kinetic_Energy		= 0.00;
				Dissipation_KE		= 0.00;

				Fpg_x			= 0.00;
				Fpg_y			= 0.00;
				Fpg			= 0.00;
				
				Fgf_x			= 0.00;
				Fgf_y			= 0.00;
				Fgf			= 0.00;
				
				Fpf_x			= 0.00;
				Fpf_y			= 0.00;
				Fpf			= 0.00;
				
				F_gravity		= 0.00;
			}
			
			fprintf(gp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",
			xx, 
			yy, 
			Ug_interstitial, 
			Vg_interstitial, 
			speed, 
			vorticity_gas,
			Pressure, 
			Void_Fraction, 
			Uf_interstitial, 
			Vf_interstitial, 
			speed_f, 
			vorticity_fines,
			Powder_Pressure, 
			Fpg_x,
			Fpg_y,
			Fpg, 
			Fgf_x,
			Fgf_y,
			Fgf, 
			Fpf_x, 
			Fpf_y, 
			Fpf, 
			F_gravity, 
			Kinetic_Energy, 
			Dissipation_KE,  
			Static_Holdup, 
			Dynamic_Holdup, 
			Total_Holdup, 
			Stress, 
			State[i][j]);
			
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}
	fclose(gp);
	// CRITICAL UPDATE OF count	
	Transient_output_count++;
	//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS ON LOCATIONS OF INTEREST
	int ii = 1;
	double xx0 = 0.0;
	while (xx0 < 0.06)
	{
		xx0	=	xx0 + dx[ii];
		ii++;
	}
	
	fp=fopen("3.3.Pure_Gas_Pressure_Profiles.dat","wt");
		fprintf(fp,"TITLE = Properties at X Station\n");
		fprintf(fp, "VARIABLES=\"X\",\"Y\",\"P\"\n");
		fprintf(fp, "ZONE T = \"Height Variation\"\n");
		
		yy = 0.5*dy[1];
		for (j = 1; j <= JMAX-1; j++)
		{
			// PRINT INTO FILE STEADY STATE FILE
			fprintf(fp,"%e %e %e\n",xx0, yy, p[ii][j]);
			// UPDATE Y POSITION
			yy = yy + 0.5*(dy[j]+dy[j+1]);
		}
	fclose(fp);  
}
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void save_convergence_data()
{
	FILE *fp;
	fp = fopen("5.1.Convergence.dat","a");
	if(Convergence_output_count==1){

		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Residue\",\"Mass of fines IN\",\"Mass of fines OUT\",\"Mass balance of fines\"\n");
	}
		fprintf(fp,"%d %e %e %e %e\n",Convergence_output_count,Residue,mass_in_fine,mass_out_fine,mass_fine);
	fclose(fp);
	Convergence_output_count++;
	//========================================================================================================
	fp = fopen("5.2a.Performance_Scheme.dat","a");
	if(Performance_scheme_count==1)
	{
		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Velocity X\",\"Gas Velocity Y\",\"Gas Pressure\",\"Gas Turbulent KE\",\"Gas Turbulent DE\",\"Fines Velocity X\",\"Fines Velocity Y\",,\"Dynamic Holdup\"\n");
	}
		fprintf(fp,"%d %d %d %d %d %d %d %d %d\n",Performance_scheme_count,Conv_U,Conv_V,Conv_P,Conv_KE,Conv_de,Conv_Uf,Conv_Vf,Conv_epfd);
	fclose(fp);
	Performance_scheme_count++;
	//========================================================================================================
	fp = fopen("5.2b.Performance_Implementation.dat","a");
	if(Performance_method_count == 1)
	{
		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Velocity X\",\"Gas Velocity Y\",\"Gas Pressure\",\"Gas Turbulent KE\",\"Gas Turbulent DE\",\"Fines Velocity X\",\"Fines Velocity Y\",\"Dynamic Holdup\"\n");
	}
		fprintf(fp,"%d %e %e %e %e %e %e %e %e\n",Performance_method_count,cpu_time_U_gas,cpu_time_V_gas,cpu_time_Pressure,cpu_time_KE_gas,cpu_time_de_gas,cpu_time_U_fines,cpu_time_V_fines,cpu_time_dynamic_holdup);
	fclose(fp);
	Performance_method_count++;
	//========================================================================================================
	Conv_U = 0;
	Conv_V = 0;
	Conv_P = 0;
	Conv_KE = 0;
	Conv_de = 0;
	Conv_Uf = 0;
	Conv_Vf = 0;
	Conv_epfd =0;
}
