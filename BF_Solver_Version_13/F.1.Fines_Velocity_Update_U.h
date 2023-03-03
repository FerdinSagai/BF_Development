void uf_update() 
{
	//------Local variables declaration--------------------------------------------------------
	int i, j;
	double bound =	100.0;
	double Fe[IMAX+1][JMAX+1],Fw[IMAX+1][JMAX+1],Fn[IMAX+1][JMAX+1],Fs[IMAX+1][JMAX+1];
	double De[IMAX+1][JMAX+1],Dw[IMAX+1][JMAX+1],Dn[IMAX+1][JMAX+1],Ds[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1],Peclet_y[IMAX+1][JMAX+1];
	double Source[IMAX+1][JMAX+1];
	double Previous_x[IMAX+1][JMAX+1];
	//------Noting the cpu start time ---------------------------------------------------------
	double time_before = clock();
	//------Initialization of the values to be determined--------------------------------------
	for (j = 1; j <= JMAX-1; j++)
	{
		for (i = 1; i <= IMAX-2; i++)
		{	
			//------Define cell size and distances----------
			double dx_p	= 0.5*(dx[i] + dx[i+1]);
			double dy_p	= dy[j];
			double dx_e	= dx[i+1];
			double dx_w	= dx[i];
			double dy_n	= 0.5*(dy[j] + dy[j+1]);
			double dy_s;
			if (j == 1)	dy_s	=	dy[j]; 
			else		dy_s	=	0.5*(dy[j-1] + dy[j]);	
			//------Computation of values at face-centers--
			double mue_n = mu_f;
			double mue_s = mu_f;
			double mue_e = mu_f;
			double mue_w = mu_f;

			double vfd_n = 0.25*(epfd[i][j] + epfd[i+1][j] + epfd[i][j+1] + epfd[i+1][j+1]);
			double vfd_s = 0.25*(epfd[i][j] + epfd[i+1][j] + epfd[i][j-1] + epfd[i+1][j-1]);
			double vfd_e = epfd[i+1][j];
			double vfd_w = epfd[i][j];	
			//------Adjustment of values for boundaries----
			if (j == 1)				// Bottom Wall
			{
				vfd_s	=	0.25 * (epfd[i][j] + epfd[i+1][j] + epfd[i][j] + epfd[i+1][j]);
			}
			else if (j == JMAX-1)	// Top Wall
			{
				vfd_n	=	0.25 * (epfd[i][j] + epfd[i+1][j] + epfd[i][j] + epfd[i+1][j]);
			}
			//------X Staggered velocity at cell faces-----
			double U_P = uf[i][j];
			double U_E = uf[i+1][j];
			double U_W = uf[i-1][j];
			double U_e	= 0.5*(U_P + U_E);
			double U_w	= 0.5*(U_P + U_W);
			//------Y Staggered velocity at cell faces-----
			double V_nw= vf[i][j];
			double V_ne= vf[i+1][j];
			double V_sw= vf[i][j-1];
			double V_se= vf[i+1][j-1];
			double V_n	= 0.5*(V_nw + V_ne);
			double V_s	= 0.5*(V_sw + V_se);
			double V_e	= 0.5*(V_ne + V_se);
			double V_w	= 0.5*(V_nw + V_sw);	
			//------ Convective Fluxes in the X direction--	
			Fe[i][j]	= (1.00/dx_p) * 0.50 * rho_f * vfd_e * U_e;
			Fw[i][j]	= (1.00/dx_p) * 0.50 * rho_f * vfd_w * U_w;
			Fn[i][j]	= (1.00/dy_p) * 0.50 * rho_f * vfd_n * V_n;
			Fs[i][j]	= (1.00/dy_p) * 0.50 * rho_f * vfd_s * V_s;
			//------ Diffusive Fluxes in the X direction--
			De[i][j]	= (vfd_e*mue_e/(dx_p*dx_e));
			Dw[i][j]	= (vfd_w*mue_w/(dx_p*dx_w));
			Dn[i][j]	= (vfd_n*mue_n/(dy_p*dy_n));
			Ds[i][j]	= (vfd_s*mue_s/(dy_p*dy_s));
			//------ Sources and Sinks in the X direction--
			//=============================================== PRESSURE ENERGY COMPUTATION====================
			//------Define the shared pressure--------------
			double P_E	= p[i+1][j];
			double P_P	= p[i][j];
			double P1	= (vfd_e*P_E - vfd_e*P_P)/dx_p;
			//------Define the dynamic holdups--------------
			double vff_e= epfd[i+1][j];
			double vff_w= epfd[i][j];			
			//------Powder pressure term computation--------
			double Pf1	= (1/150)*exp(150*(vff_e - 0.6));			// 0.6 is arbitrary
			double Pf2	= (1/150)*exp(150*(vff_w - 0.6));			// 0.6 is arbitrary
			double P2	= (Pf1 - Pf2)/dx_p;
			//------Total pressure term computation--------
			double Pressure_term= P1 + P2;
			double Ffg			= interp1D( 0.5*dx[i], Forcefield_GF[i][j], 0.5*dx[i+1], Forcefield_GF[i+1][j] );
			double Fpf			= interp1D( 0.5*dx[i], Forcefield_PF[i][j], 0.5*dx[i+1], Forcefield_PF[i+1][j] );
			double Fgr			= 0.00;
			double Ext_forces	= Ffg * (u[i][j] - uf[i][j]) - Fgr - Fpf * (uf[i][j] - us[i][j]);
			Source[i][j]		= Ext_forces - Pressure_term; 
			// Determination of Peclet number
			Peclet_x[i][j]	= (U_e + U_w)* dy_p/(mue_e + mue_w);
			Peclet_y[i][j]	= (V_n + V_s)* dx_p/(mue_n + mue_s);
			// Stores value of the previous iteration
			Previous_x[i][j]	= 0.00;
		}
	}
	int explicit_convection_flag = 1;
	int explicit_diffusion_flag = 0;
	//------Loop for the determination of uf---------------------------------------------------
	while (bound>tol)
	{
		Conv_Uf++;	
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
				start_j = 1;
			}
			if(tid == nthreads-1)
			{
				finish_j = JMAX-1;
			}
			//------Traversal of the grid and computation loop begins------------------
			for (j=start_j; j<=finish_j; j++)
			{	
				for (i=0; i<=IMAX-1; i++)
				{	
					if(i == 0)				// Left Wall
					{
						ff[0][j]	= u_left;
					}
					else if (i == IMAX-1)	// Right Wall
					{
						ff[IMAX-1][j]	= u_right;
					}
					else
					{
						//------Computation of values at cell-center---
						double vfd_p			= 0.50 * (epfd[i+1][j] + epfd[i][j]);
						//------Define the velocity variables----------
						double U_P	= uf[i][j];
						double U_N	= uf[i][j+1];
						double U_S	= uf[i][j-1];
						double U_E	= uf[i+1][j];
						double U_W	= uf[i-1][j];	
						//------Define the intermediate velocity variables-----------
						double f_P	= ff[i][j];
						double f_N	= ff[i][j+1];
						double f_S	= ff[i][j-1];
						double f_E	= ff[i+1][j];
						double f_W	= ff[i-1][j];	
						//------Adjustment of values for boundaries----
						if (j == 1)					// Bottom Wall
						{
							f_S	=	2.00 * u_bottom - ff[i][j];
						}
						else if (j == JMAX-1)		// Top Wall
						{
							f_N	=	2.00 * u_top - ff[i][j];
						}						
						//=============================================== CONVECTION COMPUTATION ========================================//
						//------Convection term(Explicit for simplicity)
						double DUUDX, DUVDY;
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
							DUUDX			= ( Fe[i][j] * ((U_P + U_E) + alpha_x*(U_P - U_E)) - Fw[i][j] * ((U_P + U_W) + alpha_x*(U_W - U_P)) );
							DUVDY			= ( Fn[i][j] * ((U_P + U_N) + alpha_y*(U_P - U_N)) - Fs[i][j] * ((U_P + U_S) + alpha_y*(U_S - U_P)) );							
						}
						else
						{
							DUUDX			= ( Fe[i][j] * ((f_P + f_E) + alpha_x*(f_P - f_E)) - Fw[i][j] * ((f_P + f_W) + alpha_x*(f_W - f_P)) );
							DUVDY			= ( Fn[i][j] * ((f_P + f_N) + alpha_y*(f_P - f_N)) - Fs[i][j] * ((f_P + f_S) + alpha_y*(f_S - f_P)) );							
						}
						double Convec_term		= DUUDX + DUVDY;				
						//================================================DIFFUSION COMPUTATION ========================================//
						//------Diffusion term(Implicit)
						double D2UDX2, D2UDY2;
						if( explicit_diffusion_flag == 1)
						{
							D2UDX2			= (Dw[i][j] * U_W) + (De[i][j] * U_E);
							D2UDY2			= (Ds[i][j] * U_S) + (Dn[i][j] * U_N);
						}
						else
						{
							D2UDX2			= (Dw[i][j] * f_W) + (De[i][j] * f_E);
							D2UDY2			= (Ds[i][j] * f_S) + (Dn[i][j] * f_N);							
						}
						double Diffusive_term	= D2UDX2 + D2UDY2;
						//================================================SOURCE COMPUTATION =====================================//
						double Source_term		= Source[i][j]; 
						//================================================PREVIOUS TIME TERM COMPUTATION ============================//
						double Previous_Time_term	= (vfd_p*rho_f/dt1) * U_P;	
						//------Computation of uf---------------------
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
							ff[i][j]	= 0.0;	// If the grid point is in the solid phase
						}
						else			
						{
							ff[i][j]	= (1.0-omega_uf)*f_P + (omega_uf/aP)*(Previous_Time_term  - Convec_term + Diffusive_term + Source_term);
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
		for (i = 0; i <= IMAX-1; i++)
		{
			for (j = 1;j <= JMAX-1; j++)
			{
				bound			= bound + fabs(ff[i][j] - Previous_x[i][j]);
				Previous_x[i][j]	= ff[i][j];
			}
		}
	}
	//------Convergence met. uf has been determined--------------------------------------
	//------Compute the maximum value of uf for monitoring and debugging-----------------
	max_Uf		= 0.00;
	for (i = 0; i <= IMAX; i++)
	{
		for (j=0;j<=JMAX;j++)
		{
			if (abs(ff[i][j]) > max_Uf)
			{
				max_Uf = abs(ff[i][j]);
			}
		}
	}
	//------Noting the cpu end time ------------------------------------------------------
	double time_after = clock();
	//------Enforce a state-based boundary condition--------------------------------------
	enforce_domain_physics();
	//------Noting the total cpu time for this computation--------------------------------
	cpu_time_U_fines = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);
}