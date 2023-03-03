//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void gas_convergence_data()
{
	FILE *fp;
	fp = fopen("3.0.Gas_Convergence.dat","a");
	if(Gas_Convergence_count==1)
	{
		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Residue\",\"Mass of fines IN\",\"Mass of fines OUT\",\"Mass balance of fines\"\n");
	}
		fprintf(fp,"%d %e %e\n",Gas_Convergence_count,Residue_V,Residue_P);
	fclose(fp);
	Gas_Convergence_count++;
}
//------------------------------------------------------------------------------
void save_output_gas_transient()
{
	double xx, yy;
	double Ug_interstitial,Vg_interstitial, speed, vorticity_gas;
	double Pressure, Void_Fraction, Kinetic_Energy, Dissipation_KE, Stress;
	double Fpg_x, Fpg_y, Fpg;
	FILE *fp;
//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	fp=fopen("3.1.Gas_Flow_Transient.dat","a");
	if(Gas_Transient_output_count == 1)
	{
		fprintf(fp,"TITLE = Properties at Final state\n");
		fprintf(fp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Ug(m/s)\",\"Vg(m/s)\",\"Speed(m/s)\",\"Vorticity\",\"Pressure (Pa)\",\"Void Fraction\",\"F_gas-particles X\",\"F_gas-particles Y\",\"F_gas-particles\",\"KE\",\"DE\",\"Stress\",\"State\"\n");
	}
	fprintf(fp, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",GFS_iter, IMAX-1, JMAX-1);	
	
	yy = 0.5*dy[1];
	for (int j = 1; j <= JMAX-1; j++)
	{
		xx  = 0.5*dx[1];
		for (int i = 1; i <= IMAX-1; i++)
		{
			// INTERSTITIAL VELOCITY COMPUTATION
			Ug_interstitial 		= 0.5*(u[i][j]+u[i-1][j]);
			Vg_interstitial 		= 0.5*(v[i][j]+v[i][j-1]);
			speed					= sqrt(Ug_interstitial*Ug_interstitial + Vg_interstitial*Vg_interstitial);
			vorticity_gas			= vorticity_g[i][j];
			// SCALAR COMPUTATIONS
			Pressure				=	p[i][j];
			Void_Fraction			=	void_frac(i,j);		
			Kinetic_Energy			=	ke[i][j];
			Dissipation_KE			=	de[i][j];
			Stress					=	sigma_eff[i][j];
			// FORCE COMPUTATION
			Fpg_x					=	Fg_p(i,j)* 0.5 * (u[i][j] + u[i-1][j]);
			Fpg_y					=	Fg_p(i,j)* 0.5 * (v[i][j] + v[i][j-1]);
			Fpg						=	sqrt(sqr(Fpg_x) +sqr(Fpg_y));
			// CONSIDERATION OF GEOMETRY
			if (State[i][j] == -1)
			{
				Ug_interstitial 	= 0.00;
				Vg_interstitial 	= 0.00;
				speed 				= 0.00;
				vorticity_g[i][j] 	= 0.00;
							
				Pressure			= 0.00;
				Kinetic_Energy		= 0.00;
				Dissipation_KE		= 0.00;

				Fpg_x				= 0.00;
				Fpg_y				= 0.00;
				Fpg					= 0.00;
			}

			fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",
			xx, 
			yy, 
			Ug_interstitial, 
			Vg_interstitial, 
			speed, 
			vorticity_gas,
			Pressure, 
			Void_Fraction, 
			Fpg_x,
			Fpg_y,
			Fpg, 
			Kinetic_Energy, 
			Dissipation_KE,  
			Stress, 
			State[i][j]);	
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}		
	fclose(fp);
//-------------------------HYDRODYNAMIC PRESSURE AT PRESSURE PORTS
	fp=fopen("3.3.Gas_Pressure_Port.dat","a");
	if(Gas_Transient_output_count == 1)
	{
		fprintf(fp,"TITLE = Pressure at ports\n");
		fprintf(fp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Pressure (Pa)\"\n");
	}
	fprintf(fp, "ZONE T = \"%d\"\n",GFS_iter);	
	for (int i = 0; i < N_Pressure_ports; i++)
	{	
 		int x_id	=	int(2 * Port_X[i]/dx[1]); 
		int y_id	=	int(2 * Port_Y[i]/dy[1]);
		Port[i]		=	p[x_id][y_id];
		fprintf(fp,"%e %e %e\n", Port_X[i], Port_Y[i], Port[i]); 
	}	
	fclose(fp);
	
	// CRITICAL UPDATE OF count	
	Gas_Transient_output_count++;
}
//------------------------------------------------------------------------------
void save_output_gas_steady_state()
{
	double xx, yy;
	double Ug_interstitial,Vg_interstitial, speed, vorticity_gas;
	double Pressure, Void_Fraction, Kinetic_Energy, Dissipation_KE, Stress;
	double Fpg_x, Fpg_y, Fpg;
	FILE *fp;
//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	fp=fopen("3.2.Gas_Flow_Steady_State.dat","wt");
	fprintf(fp,"TITLE = Properties at Final state\n");
	fprintf(fp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Ug(m/s)\",\"Vg(m/s)\",\"Speed(m/s)\",\"Vorticity\",\"Pressure (Pa)\",\"Void Fraction\",\"F_gas-particles X\",\"F_gas-particles Y\",\"F_gas-particles\",\"KE\",\"DE\",\"Stress\",\"State\"\n");
	fprintf(fp, "ZONE T = \"%lf\", I=%d, J=%d F=POINT\n",solution_time, IMAX-1, JMAX-1);	
	
	yy = 0.5*dy[1];
	for (int j = 1; j <= JMAX-1; j++)
	{
		xx  = 0.5*dx[1];
		for (int i = 1; i <= IMAX-1; i++)
		{
			// INTERSTITIAL VELOCITY COMPUTATION
			Ug_interstitial 		= 0.5*(u[i][j]+u[i-1][j]);
			Vg_interstitial 		= 0.5*(v[i][j]+v[i][j-1]);
			speed					= sqrt(Ug_interstitial*Ug_interstitial + Vg_interstitial*Vg_interstitial);
			vorticity_gas			= vorticity_g[i][j];
			// SCALAR COMPUTATIONS
			Pressure				=	p[i][j];
			Void_Fraction			=	void_frac(i,j);		
			Kinetic_Energy			=	ke[i][j];
			Dissipation_KE			=	de[i][j];
			Stress					=	sigma_eff[i][j];
			// FORCE COMPUTATION
			Fpg_x					=	Fg_p(i,j)* 0.5 * (u[i][j] + u[i-1][j]);
			Fpg_y					=	Fg_p(i,j)* 0.5 * (v[i][j] + v[i][j-1]);
			Fpg						=	sqrt(sqr(Fpg_x) +sqr(Fpg_y));
			// CONSIDERATION OF GEOMETRY
			if (State[i][j] == -1)
			{
				Ug_interstitial 	= 0.00;
				Vg_interstitial 	= 0.00;
				speed 				= 0.00;
				vorticity_g[i][j] 	= 0.00;
							
				Pressure			= 0.00;
				Kinetic_Energy		= 0.00;
				Dissipation_KE		= 0.00;

				Fpg_x				= 0.00;
				Fpg_y				= 0.00;
				Fpg					= 0.00;
			}

			fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",
			xx, 
			yy, 
			Ug_interstitial, 
			Vg_interstitial, 
			speed, 
			vorticity_gas,
			Pressure, 
			Void_Fraction, 
			Fpg_x,
			Fpg_y,
			Fpg, 
			Kinetic_Energy, 
			Dissipation_KE,  
			Stress, 
			State[i][j]);	
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}		
	fclose(fp);
	
	//-------------------------GAS VELOCITY DETAILS FOR CFD SOLVER
	fp=fopen("2.2.Gas_Velocity_details.int","wt");
	for (int j = 0; j <= JMAX; j++)
	{
		for (int i = 0; i <= IMAX; i++)
		{
			if(i == 0)
			{
				Ug_interstitial = u_left;
			}
			else if(i == IMAX)
			{
				Ug_interstitial = u_right;
			}
			else if(j == 0)
			{
				Vg_interstitial = v_bottom;
			}
			else if(j == JMAX)
			{
				Vg_interstitial = v_top;
			}
			else
			{
				Ug_interstitial 		= 0.5*(u[i][j]+u[i-1][j]);
				Vg_interstitial 		= 0.5*(v[i][j]+v[i][j-1]);
			}
			fprintf(fp,"%e %e\n",	Ug_interstitial, Vg_interstitial);	
		}
	}		
	fclose(fp);
	
//-------------------------HYDRODYNAMIC PRESSURE AT PRESSURE PORTS
	fp=fopen("3.5.Gas_Pressure_Port_Steady.dat","w");
		fprintf(fp,"TITLE = Pressure at ports\n");
		for (int i = 0; i < N_Pressure_ports; i++)
		{	
			int x_id	=	int(2 * Port_X[i]/dx[1]); 
			int y_id	=	int(2 * Port_Y[i]/dy[1]);
			Port[i]		=	p[x_id][y_id];	
		}		
		fprintf(fp,"PRESSURE PORTS\n");
		fprintf(fp,"\tA1 ...........................................................%e\n",Port[0]);
		fprintf(fp,"\tA2 ...........................................................%e\n",Port[1]);
		fprintf(fp,"\tA3 ...........................................................%e\n",Port[2]);
		fprintf(fp,"\tA4 ...........................................................%e\n",Port[3]);
		fprintf(fp,"\tA5 ...........................................................%e\n",Port[4]);
		fprintf(fp,"\tA6 ...........................................................%e\n",Port[5]);
		fprintf(fp,"\tA7 ...........................................................%e\n",Port[6]);
		fprintf(fp,"\tA8 ...........................................................%e\n",Port[7]);
		fprintf(fp,"\tA9 ...........................................................%e\n",Port[8]);
		fprintf(fp,"\tA10...........................................................%e\n",Port[9]);
		
		fprintf(fp,"\n");
		fclose(fp);
}
