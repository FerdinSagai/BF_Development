void v_update()
{
	//------Local variables declaration--------------------------------------------------------
	int i, j;
	double bound = 100.0;
	double Fe[IMAX+1][JMAX+1],Fw[IMAX+1][JMAX+1],Fn[IMAX+1][JMAX+1],Fs[IMAX+1][JMAX+1];
	double De[IMAX+1][JMAX+1],Dw[IMAX+1][JMAX+1],Dn[IMAX+1][JMAX+1],Ds[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1],Peclet_y[IMAX+1][JMAX+1];
	double Source[IMAX+1][JMAX+1];
	double Previous_y[IMAX+1][JMAX+1];
	//------Noting the cpu start time --------------------------------------------------------
	double time_before = clock();
	//=============================================== FLUX COMPUTATION ==============================================
	for (j = 1; j <= JMAX-2; j++)
	{
		for (i = 1; i <= IMAX-1; i++)
		{
			//------Define cell size and distances----------
			double dx_p	= dx[i];
			double dy_p	= 0.5*(dy[j] + dy[j+1]);
			double dx_e	= 0.5*(dx[i] + dx[i+1]);
			double dx_w;
			if (i == 1) 	dx_w	=	dx[i]; 
			else 		dx_w	=	0.5*(dx[i-1] + dx[i]);
			double dy_n	= dy[j+1];
			double dy_s	= dy[j] ;
			//------Computation of values at face-centers--
			double mue_n = mue(i,j+1);
			double mue_s = mue(i,j);
			double mue_e = 0.25*(mue(i,j)+mue(i+1,j)+mue(i,j+1)+mue(i+1,j+1));
			double mue_w = 0.25*(mue(i,j)+mue(i-1,j)+mue(i,j+1)+mue(i-1,j+1));

			double vf_n = void_frac(i,j+1);
			double vf_s = void_frac(i,j);
			double vf_e = 0.25*(void_frac(i,j)+void_frac(i+1,j)+void_frac(i,j+1)+void_frac(i+1,j+1));
			double vf_w = 0.25*(void_frac(i,j)+void_frac(i-1,j)+void_frac(i,j+1)+void_frac(i-1,j+1));
			//------Adjustment of values for boundaries--
			if (i == 1)				// Left Wall
			{
				mue_w	= 0.25 * (mue(i,j)+mue(i,j)+mue(i,j+1)+mue(i,j+1));
				vf_w	= 0.25 * (void_frac(i,j)+void_frac(i,j)+void_frac(i,j+1)+void_frac(i,j+1));
			}
			else if(i == IMAX-1)	// Right Wall
			{
				mue_e	= 0.25 * (mue(i,j) + mue(i,j) + mue(i,j+1) + mue(i,j+1));
				vf_e	= 0.25 * (void_frac(i,j) + void_frac(i,j) + void_frac(i,j+1) + void_frac(i,j+1));
			}
			//------Y Staggered velocity at cell faces-----
			double V_p	= v[i][j];
			double V_pn	= v[i][j+1];
			double V_ps	= v[i][j-1];
			double V_n	= 0.5*(V_p + V_pn);
			double V_s	= 0.5*(V_p + V_ps);
			//------X Staggered velocity at cell faces-----
			double U_nw	= u[i-1][j+1];
			double U_ne	= u[i][j+1];
			double U_sw	= u[i-1][j];
			double U_se	= u[i][j];	
			double U_n	= 0.5*(U_nw + U_ne);
			double U_s	= 0.5*(U_sw + U_se);
			double U_e	= 0.5*(U_ne + U_se);
			double U_w	= 0.5*(U_nw + U_sw);
			//------ Convective Fluxes in the Y direction--	
			Fe[i][j]	= (1.00/dx_p) * 0.50 * rho_g * vf_e * U_e;
			Fw[i][j]	= (1.00/dx_p) * 0.50 * rho_g * vf_w * U_w;
			Fn[i][j]	= (1.00/dy_p) * 0.50 * rho_g * vf_n * V_n;
			Fs[i][j]	= (1.00/dy_p) * 0.50 * rho_g * vf_s * V_s;
			//------ Diffusive Fluxes in the Y direction--
			De[i][j]	= (vf_e*mue_e/(dx_p*dx_e));
			Dw[i][j]	= (vf_w*mue_w/(dx_p*dx_w));
			Dn[i][j]	= (vf_n*mue_n/(dy_p*dy_n));
			Ds[i][j]	= (vf_s*mue_s/(dy_p*dy_s));
			//------ Sources and Sinks in the Y direction--
			//================================================TURBULENT KINETIC ENERGY COMPUTATION====================//
			double KE_N 			= ke[i][j+1];
			double KE_P 			= ke[i][j];
			double Turbulent_KE_term= (2.0/3.0) * rho_g *  (vf_n * KE_N - vf_s * KE_P)/dy_p;	
			double Fgp				= interp1D( 0.5*dy[j], Forcefield_GP[i][j], 0.5*dy[j+1], Forcefield_GP[i][j+1] );
			double Fgf				= interp1D( 0.5*dy[j], Forcefield_GF[i][j], 0.5*dy[j+1], Forcefield_GF[i][j+1] );
			double Fgl				= interp1D( 0.5*dy[j], Forcefield_GL[i][j], 0.5*dy[j+1], Forcefield_GL[i][j+1] );			
			double Ext_forces 		= - Fgp * (v[i][j] - vs[i][j]) - Fgf * (v[i][j] - vf[i][j]) - Fgl * (v[i][j] - vl[i][j]);
			Source[i][j]			= Ext_forces - Turbulent_KE_term;
			// Determination of Peclet number
			Peclet_x[i][j]	= (U_e + U_w)* dx_p/(mue_e + mue_w);
			Peclet_y[i][j]	= (V_n + V_s)* dy_p/(mue_n + mue_s);
			// Stores value of the previous iteration
			Previous_y[i][j]	= 0.00;		
		}
	}
	int explicit_convection_flag = 0;
	int explicit_diffusion_flag = 0;	
	/*------Loop for the determination of v* (g)---------------------------------------------*/
	while (bound>tol)
	{
		Conv_V++;
		/*------Begins thread level parallelization--------------------------------------*/
		#pragma omp parallel private(i,j,tid)
		{	
			/*------Domain decomposition for each thread is performed----------------*/
			tid			= omp_get_thread_num();
			nthreads	= omp_get_max_threads();
			int start_j = 1 + tid*int(JMAX/nthreads);
			int finish_j= start_j + JMAX/nthreads - 1;
			if(tid == 0)
			{
				start_j = 0;
			}
			if(tid == nthreads-1)
			{
				finish_j = JMAX-1;
			}
			/*------Traversal of the grid and computation loop begins---------------*/
			for (j = start_j; j <= finish_j; j++)
			{
				for (i = 1; i <= IMAX-1; i++)
				{
					if(j == 0)					// Bottom Wall
					{
						g[i][0]			=	v_bottom;				
					}
					else if(j == JMAX-1)				// Top Wall
					{
						g[i][JMAX-1]	=	g[i][JMAX-2];
					}
					else
					{		
						//------Computation of values at cell-center---
						double vf_p	= 	0.50*(void_frac(i,j+1) + void_frac(i,j));
						//------Define the velocity variables----------
						double V_p	= v[i][j];
						double V_pn	= v[i][j+1];
						double V_ps	= v[i][j-1];
						double V_pe	= v[i+1][j];
						double V_pw	= v[i-1][j];
						//------Define the intermediate velocity variables-----------
						double g_p	= g[i][j];
						double g_pn	= g[i][j+1];
						double g_ps	= g[i][j-1];
						double g_pe	= g[i+1][j];
						double g_pw	= g[i-1][j];			
						//------Adjustment of values for boundaries--
						if (i == 1)						// Left Wall
						{
							g_pw	= 2.00 * v_left - g[i][j];
						}
						else if(i == IMAX-1)			// Right Wall
						{
							g_pe	= 2.00 * v_right - g[i][j];
						}
						//=============================================== CONVECTION COMPUTATION ========================================//						
						//------Convection term(Implicit)
						double DUVDX, DVVDY;
						double alpha_x = 1.00;
						if (Peclet_x[i][j] < - 2) 		alpha_x = -1;
						else if (Peclet_x[i][j] > 2)	alpha_x = 1;
						else							alpha_x = 0;					
						double alpha_y = 1.00;
						if (Peclet_y[i][j] < - 2)		alpha_y = -1;
						else if (Peclet_y[i][j] > 2)	alpha_y = 1;
						else							alpha_y = 0;
						if( explicit_convection_flag == 1)
						{
							DUVDX	= ( Fe[i][j] * ((V_p + V_pe) + alpha_x * (V_p - V_pe)) - Fw[i][j] * ((V_pw + V_p) + alpha_x * (V_pw - V_p)) );
							DVVDY	= ( Fn[i][j] * ((V_p + V_pn) + alpha_y * (V_p - V_pn)) - Fs[i][j] * ((V_ps + V_p) + alpha_y * (V_ps - V_p)) );
						}
						else
						{
							DUVDX	= ( Fe[i][j] * ((g_p + g_pe) + alpha_x * (g_p - g_pe)) - Fw[i][j] * ((g_pw + g_p) + alpha_x * (g_pw - g_p)) );
							DVVDY	= ( Fn[i][j] * ((g_p + g_pn) + alpha_y * (g_p - g_pn)) - Fs[i][j] * ((g_ps + g_p) + alpha_y * (g_ps - g_p)) );
						}
						double Convec_term 	= DUVDX + DVVDY;						
						//================================================DIFFUSION COMPUTATION ========================================//
						/*------Diffusion term(Implicit for stability)*/
						double D2VDX2, D2VDY2;
						if( explicit_diffusion_flag == 1)
						{
							D2VDX2		= (Dw[i][j] * V_pw) + (De[i][j] * V_pe);
							D2VDY2		= (Ds[i][j] * V_ps) + (Dn[i][j] * V_pn);
						}
						else
						{
							D2VDX2		= (Dw[i][j] * g_pw) + (De[i][j] * g_pe);
							D2VDY2		= (Ds[i][j] * g_ps) + (Dn[i][j] * g_pn);							
						}
						double Diffusive_term	= D2VDX2 + D2VDY2;
						//================================================SOURCE COMPUTATION =====================================//
						double Source_term		= Source[i][j];
						/*------Previous time term computation--------*/		
						double Previous_Time_term 	= (vf_p*rho_g/dt1) * V_p;
						/*------Computation of v*---------------------*/
						double aP;
						if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) )
						{
							aP	= 	(vf_p*rho_g/dt1); 	
						}
						else if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 0) )
						{
							aP	= 	(vf_p*rho_g/dt1) + (De[i][j] + Dw[i][j] + Dn[i][j] + Ds[i][j]); 
						}
						else if( ( explicit_convection_flag == 0) && ( explicit_diffusion_flag == 1) )
						{
							aP	= 	(vf_p*rho_g/dt1) + (1.00 + alpha_x)*Fe[i][j] - (1.00 - alpha_x)*Fw[i][j] + (1.00 + alpha_y)*Fn[i][j] - (1.00 - alpha_y)*Fs[i][j]; 
						}
						else
						{
							aP	= 	(vf_p*rho_g/dt1) + (1.00 + alpha_x)*Fe[i][j] - (1.00 - alpha_x)*Fw[i][j] + (1.00 + alpha_y)*Fn[i][j] - (1.00 - alpha_y)*Fs[i][j] + (De[i][j] + Dw[i][j] + Dn[i][j] + Ds[i][j]); 
						}
						if (aP == 0)
						{
							g[i][j]	= 0.0;	// If the grid point is in the solid phase
						}
						else
						{	
							g[i][j]	= (1.00 - omega_v) * g_p + (omega_v/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
						}
					}
				}
			}
		}
		/*------Thread synchronization is performed at this point--------------------*/		
		#pragma omp barrier	
		/*------Enforce a state-based boundary condition-----------------------------*/		
		enforce_domain_physics();
		/*------Computation of residue to determine convergence----------------------*/	
		bound	=	0.0;	        
		for (i = 1; i <= IMAX-1; i++)
		{
			for (j = 0; j <= JMAX-1; j++)
			{
				bound		= bound + fabs(g[i][j] - Previous_y[i][j]);
				Previous_y[i][j]	= g[i][j];
			}
		}
		if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) ) bound = -tol;
	}
	/*------Convergence met. v* has been determined--------------------------------------*/		
	/*------Compute the maximum value of v* for monitoring and debugging-----------------*/		
	max_Vg		= 0.00;
	for (i = 0; i <= IMAX; i++)
	{
		for (j = 0; j <= JMAX; j++)
		{
			if (abs(g[i][j]) > max_Vg)
			{
				max_Vg = abs(g[i][j]);
			}
		}
	}
	/*------Noting the cpu end time ------------------------------------------------------*/
	double time_after  = clock();
	/*------Enforce a state-based boundary condition--------------------------------------*/
	enforce_domain_physics();
	/*------Noting the total cpu time for this computation--------------------------------*/
	cpu_time_V_gas = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);
}