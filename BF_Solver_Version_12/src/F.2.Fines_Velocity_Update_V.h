void vf_update() 
{
	//------Local variables declaration--------------------------------------------------------
	int i, j;
	double bound =	100.0;
	double Fe[IMAX+1][JMAX+1],Fw[IMAX+1][JMAX+1],Fn[IMAX+1][JMAX+1],Fs[IMAX+1][JMAX+1];
	double De[IMAX+1][JMAX+1],Dw[IMAX+1][JMAX+1],Dn[IMAX+1][JMAX+1],Ds[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1],Peclet_y[IMAX+1][JMAX+1];
	double Source[IMAX+1][JMAX+1];
	double Previous_y[IMAX+1][JMAX+1];
	//------Noting the cpu start time ---------------------------------------------------------
	double time_before = clock();
	//------Initialization of the values to be determined--------------------------------------
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
			double mue_n			= mu_f;
			double mue_s			= mu_f;
			double mue_e			= mu_f;
			double mue_w			= mu_f;
			
			double vfd_n	= epfd[i][j+1];
			double vfd_s	= epfd[i][j];
			double vfd_e	= 0.25*(epfd[i][j] + epfd[i+1][j] + epfd[i][j+1] + epfd[i+1][j+1]);
			double vfd_w	= 0.25*(epfd[i][j] + epfd[i-1][j] + epfd[i][j+1] + epfd[i-1][j+1]);
			//------Adjustment of values for boundaries----
			if(i == 1)				// Left Wall
			{
				vfd_w	= 0.25 * (epfd[i][j] + epfd[i][j] + epfd[i][j+1] + epfd[i][j+1]);
			}
			else if(i == IMAX-1)			// Right Wall
			{
				vfd_e	= 0.25 * (epfd[i][j] + epfd[i][j] + epfd[i][j+1] + epfd[i][j+1]);
			}
			//------Y Staggered velocity at cell faces-----
			double V_P	= vf[i][j];
			double V_N	= vf[i][j+1];
			double V_S	= vf[i][j-1];
			double V_n	= 0.50 * (V_P + V_N);
			double V_s	= 0.50 * (V_P + V_S);
			//------X Staggered velocity at cell faces-----
			double U_nw = uf[i-1][j+1];
			double U_ne = uf[i][j+1];
			double U_sw = uf[i-1][j];
			double U_se = uf[i][j];
			double U_n	= 0.50 * (U_nw + U_ne);
			double U_s	= 0.50 * (U_sw + U_se);
			double U_e	= 0.50 * (U_ne + U_se);
			double U_w	= 0.50 * (U_nw + U_sw);
			//------ Convective Fluxes in the Y direction--	
			Fe[i][j]	= (1.00/dx_p) * 0.50 * rho_f * vfd_e * U_e;
			Fw[i][j]	= (1.00/dx_p) * 0.50 * rho_f * vfd_w * U_w;
			Fn[i][j]	= (1.00/dy_p) * 0.50 * rho_f * vfd_n * V_n;
			Fs[i][j]	= (1.00/dy_p) * 0.50 * rho_f * vfd_s * V_s;
			//------ Diffusive Fluxes in the Y direction--
			De[i][j]	= (vfd_e*mue_e/(dx_p*dx_e));
			Dw[i][j]	= (vfd_w*mue_w/(dx_p*dx_w));
			Dn[i][j]	= (vfd_n*mue_n/(dy_p*dy_n));
			Ds[i][j]	= (vfd_s*mue_s/(dy_p*dy_s));
			//------ Sources and Sinks in the X direction--
			//=============================================== PRESSURE ENERGY COMPUTATION====================		
			//------Define the shared pressure--------------
			double P_N	= p[i][j+1];
			double P_P	= p[i][j];
			double P1	= (vfd_n*P_N - vfd_s*P_P)/dy_p;
			//------Define the dynamic holdups--------------
			double vff_n= epfd[i][j+1];
			double vff_s= epfd[i][j];
			//------Powder pressure term computation--------
			double Pf1	= (1/150)*exp(150*(vff_n - 0.6));			// 0.6 is arbitrary
			double Pf2	= (1/150)*exp(150*(vff_s - 0.6));			// 0.6 is arbitrary
			double P2	= (Pf1 - Pf2)/dy_p;		
			//------Total pressure term computation--------
			double Pressure_term		= P1 + P2;			
			double Ffg			= interp1D( 0.5*dy[j], Forcefield_GF[i][j], 0.5*dy[j+1], Forcefield_GF[i][j+1] );
			double Fpf			= interp1D( 0.5*dy[j], Forcefield_PF[i][j], 0.5*dy[j+1], Forcefield_PF[i][j+1] );
			double Fgr			= 0.50* (epfd[i][j+1] + epfd[i][j]) * rho_f * gravity;
			double Ext_forces	= Ffg * (v[i][j] - vf[i][j]) - Fgr - Fpf * (vf[i][j] - vs[i][j]);
			Source[i][j]		= Ext_forces - Pressure_term; 
			// Determination of Peclet number
			Peclet_x[i][j]	= (U_e + U_w)* dy_p/(mue_e + mue_w);
			Peclet_y[i][j]	= (V_n + V_s)* dx_p/(mue_n + mue_s);
			// Stores value of the previous iteration
			Previous_y[i][j]	= 0.00;
		}
	}
	int explicit_convection_flag = 1;
	int explicit_diffusion_flag = 0;
	//------Loop for the determination of vf---------------------------------------------------
	while (bound>tol)
	{
		Conv_Vf++;
		//------Begins thread level parallelization----------------------------------------
		#pragma omp parallel private(i,j,tid)
		{		
			//------Domain decomposition for each thread is performed------------------
			tid			= omp_get_thread_num();
			nthreads	= omp_get_max_threads();
			int start_j	= 1 + tid*int(JMAX/nthreads);
			int finish_j= start_j + JMAX/nthreads - 1;
			if(tid == 0)
			{
				start_j = 0;
			}
			if(tid == nthreads-1)
			{
				finish_j = JMAX-1;
			}
			//------Traversal of the grid and computation loop begins------------------
			for (i = 1; i <= IMAX-1; i++)
			{
				for (j = start_j; j <= finish_j; j++)
				{
					if(j == 0)				// Bottom Wall
					{
						gf[i][0]		= v_bottom;				
					}
					else if(j == JMAX-1)	// Top Wall
					{
						gf[i][JMAX-1]	= gf[i][JMAX-2];
					}
					else
					{	
						//------Computation of values at cell-center---
						double vfd_p			= 0.50 * (epfd[i][j+1] + epfd[i][j]);
						//------Define the velocity variables----------
						double V_P	= vf[i][j];
						double V_N	= vf[i][j+1];
						double V_S	= vf[i][j-1];
						double V_E	= vf[i+1][j];
						double V_W	= vf[i-1][j];
						//------Define the intermediate velocity variables-----------
						double g_P	= gf[i][j];
						double g_N	= gf[i][j+1];
						double g_S	= gf[i][j-1];
						double g_E	= gf[i+1][j];
						double g_W	= gf[i-1][j];
						//------Adjustment of values for boundaries--
						if (i == 1)						// Left Wall
						{
							g_W	= 2.00 * v_left - gf[i][j];
						}
						else if(i == IMAX-1)			// Right Wall
						{
							g_E	= 2.00 * v_right - gf[i][j];
						}
						//=============================================== CONVECTION COMPUTATION ========================================//						
						//------Convection term(Implicit)
						double DUVDX, DVVDY;
						double alpha_x= 1.00;
						if (Peclet_x[i][j] < - 2) 		alpha_x = -1;
						else if (Peclet_x[i][j] > 2)	alpha_x = 1;
						else							alpha_x = 0;					
						double alpha_y= 1.00;
						if (Peclet_y[i][j] < - 2)		alpha_y = -1;
						else if (Peclet_y[i][j] > 2)	alpha_y = 1;
						else							alpha_y = 0;			
						if( explicit_convection_flag == 1)
						{
							DUVDX			= ( Fe[i][j] * ((V_P + V_E) + alpha_x*(V_P - V_E)) - Fw[i][j] * ((V_W + V_P) + alpha_x*(V_W - V_P)) );
							DVVDY			= ( Fn[i][j] * ((V_P + V_N) + alpha_y*(V_P - V_N)) - Fs[i][j] * ((V_S + V_P) + alpha_y*(V_S - V_P)) );
						}
						else
						{
							DUVDX			= ( Fe[i][j] * ((g_P + g_E) + alpha_x*(g_P - g_E)) - Fw[i][j] * ((g_W + g_P) + alpha_x*(g_W - g_P)) );
							DVVDY			= ( Fn[i][j] * ((g_P + g_N) + alpha_y*(g_P - g_N)) - Fs[i][j] * ((g_S + g_P) + alpha_y*(g_S - g_P)) );
						}
						double Convec_term 	= DUVDX + DVVDY;
						//================================================DIFFUSION COMPUTATION ========================================//
						//------Diffusion term(Implicit for stability)
						double D2VDX2, D2VDY2;
						if( explicit_diffusion_flag == 1)
						{
							D2VDX2			= (Dw[i][j] * V_W) + (De[i][j] * V_E);
							D2VDY2			= (Ds[i][j] * V_S) + (Dn[i][j] * V_N);
						}
						else 
						{
							D2VDX2			= (Dw[i][j] * g_W) + (De[i][j] * g_E);
							D2VDY2			= (Ds[i][j] * g_S) + (Dn[i][j] * g_N);
						}
						double Diffusive_term	= D2VDX2 + D2VDY2;
						//================================================SOURCE COMPUTATION =====================================//
						double Source_term		= Source[i][j];
						//================================================PREVIOUS TIME TERM COMPUTATION ============================//		
						double Previous_Time_term 	= (vfd_p*rho_f/dt1) * V_P;
						//------Computation of vf---------------------
						double aP;
						if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) )
						{
							aP	= 	(vfd_p*rho_f/dt1); 	
						}
						else if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 0) )
						{
							aP	= 	(vfd_p*rho_f/dt1) + (De[i][j] + Dw[i][j] + Dn[i][j] + Ds[i][j]); 
						}
						else if( ( explicit_convection_flag == 0) && ( explicit_diffusion_flag == 1) )
						{
							aP	= 	(vfd_p*rho_f/dt1) + (1.00 + alpha_x)*Fe[i][j] - (1.00 - alpha_x)*Fw[i][j] + (1.00 + alpha_y)*Fn[i][j] - (1.00 - alpha_y)*Fs[i][j]; 
						}
						else
						{
							aP	= 	(vfd_p*rho_f/dt1) + (1.00 + alpha_x)*Fe[i][j] - (1.00 - alpha_x)*Fw[i][j] + (1.00 + alpha_y)*Fn[i][j] - (1.00 - alpha_y)*Fs[i][j] + (De[i][j] + Dw[i][j] + Dn[i][j] + Ds[i][j]); 
						}
						if (aP == 0.00)							
						{
							gf[i][j]	= 0.0;	// If the grid point is in the solid phase
						}
						else		
						{
							gf[i][j]	= (1.0-omega_vf)*g_P +(omega_vf/aP)*(Previous_Time_term  - Convec_term + Diffusive_term + Source_term);		
						}		
					}
				}
			}
		}
		//------Thread synchronization is performed at this point--------------------
		#pragma omp barrier
		//------Enforce a state-based boundary condition-----------------------------
		enforce_domain_physics();
		//------Computation of residue to determine convergence----------------------
		bound	=	0.0;								
		for (i = 1; i <= IMAX-1; i++)
		{
			for (j = 0;j <= JMAX-1; j++)
			{
				bound		= bound+fabs(gf[i][j] - Previous_y[i][j]);
				Previous_y[i][j]	= gf[i][j];
			}
		}
	}
	//------Convergence met. vf has been determined--------------------------------------
	//------Compute the maximum value of vf for monitoring and debugging-----------------
	max_Vf		= 0.00;
	for (i=0; i<=IMAX; i++)
	{
		for (j=0;j<=JMAX;j++)
		{
			if (abs(gf[i][j]) > max_Vf)
			{
				max_Vf = abs(gf[i][j]);
			}
		}
	}
	//------Noting the cpu end time ------------------------------------------------------
	double time_after = clock();
	//------Enforce a state-based boundary condition--------------------------------------
	enforce_domain_physics();
	//------Noting the total cpu time for this computation--------------------------------
	cpu_time_V_fines = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);
}