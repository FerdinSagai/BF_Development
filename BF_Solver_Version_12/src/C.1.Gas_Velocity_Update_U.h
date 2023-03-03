void u_update()
{
	//------Local variables declaration--------------------------------------------------------
	int i, j;
	double bound = 100.0;
	double Fe[IMAX+1][JMAX+1],Fw[IMAX+1][JMAX+1],Fn[IMAX+1][JMAX+1],Fs[IMAX+1][JMAX+1];
	double De[IMAX+1][JMAX+1],Dw[IMAX+1][JMAX+1],Dn[IMAX+1][JMAX+1],Ds[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1],Peclet_y[IMAX+1][JMAX+1];
	double Source[IMAX+1][JMAX+1];
	double Previous_x[IMAX+1][JMAX+1];
	//------Noting the cpu start time --------------------------------------------------------
	double time_before = clock();
	//=============================================== FLUX COMPUTATION ==============================================
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
			//------Computation of scalar values at face-centers--
			double mue_n = 0.25*(mue(i,j) + mue(i+1,j) + mue(i,j+1) + mue(i+1,j+1));
			double mue_s = 0.25*(mue(i,j) + mue(i+1,j) + mue(i,j-1) + mue(i+1,j-1));
			double mue_e = mue(i+1,j);
			double mue_w = mue(i,j);
			
			double vf_n = 0.25*(void_frac(i,j)+void_frac(i+1,j)+void_frac(i,j+1)+void_frac(i+1,j+1));
			double vf_s = 0.25*(void_frac(i,j)+void_frac(i+1,j)+void_frac(i,j-1)+void_frac(i+1,j-1));
			double vf_e = void_frac(i+1,j);
			double vf_w = void_frac(i,j);			
			if (j == 1)				// Bottom Wall
			{
				mue_s	=	0.25 * (mue(i,j) + mue(i+1,j) + mue(i,j) + mue(i+1,j));
				vf_s	=	0.25 * (void_frac(i,j) + void_frac(i+1,j) + void_frac(i,j) + void_frac(i+1,j));
			}
			else if (j == JMAX-1)	// Top Wall
			{
				mue_n	=	0.25 * (mue(i,j) + mue(i+1,j) + mue(i,j) + mue(i+1,j));
				vf_n	=	0.25 * (void_frac(i,j) + void_frac(i+1,j) + void_frac(i,j) + void_frac(i+1,j));
			}
			//------X Staggered velocity at cell faces-----
			double U_p	= u[i][j];
			double U_pe	= u[i+1][j];
			double U_pw	= u[i-1][j];		
			double U_e	= 0.5*(U_p + U_pe);
			double U_w	= 0.5*(U_p + U_pw);
			//------Y Staggered velocity at cell faces-----
			double V_nw	= v[i][j];
			double V_ne	= v[i+1][j];
			double V_sw	= v[i][j-1];
			double V_se	= v[i+1][j-1];
			double V_n	= 0.5*(V_nw + V_ne);
			double V_s	= 0.5*(V_sw + V_se);
			double V_e	= 0.5*(V_ne + V_se);
			double V_w	= 0.5*(V_nw + V_sw);
			//------ Convective Fluxes in the X direction--	
			Fe[i][j]	= (1.00/dx_p) * 0.50 * rho_g * vf_e * U_e;
			Fw[i][j]	= (1.00/dx_p) * 0.50 * rho_g * vf_w * U_w;
			Fn[i][j]	= (1.00/dy_p) * 0.50 * rho_g * vf_n * V_n;
			Fs[i][j]	= (1.00/dy_p) * 0.50 * rho_g * vf_s * V_s;
			//------ Diffusive Fluxes in the X direction--
			De[i][j]	= (vf_e*mue_e/(dx_p*dx_e));
			Dw[i][j]	= (vf_w*mue_w/(dx_p*dx_w));
			Dn[i][j]	= (vf_n*mue_n/(dy_p*dy_n));
			Ds[i][j]	= (vf_s*mue_s/(dy_p*dy_s));
			//------ Sources and Sinks in the X direction--
			//================================================TURBULENT KINETIC ENERGY COMPUTATION====================//
			double KE_E				= ke[i+1][j];
			double KE_P 			= ke[i][j];
			double Turbulent_KE_term= (2.0/3.0) * rho_g *  (vf_e *KE_E - vf_w *KE_P)/dx_p; 
			double Fgp				= interp1D( 0.5*dx[i], Forcefield_GP[i][j], 0.5*dx[i+1], Forcefield_GP[i+1][j] );						
			double Fgf				= interp1D( 0.5*dx[i], Forcefield_GF[i][j], 0.5*dx[i+1], Forcefield_GF[i+1][j] );
			double Fgl				= interp1D( 0.5*dx[i], Forcefield_GL[i][j], 0.5*dx[i+1], Forcefield_GL[i+1][j] );
			double Ext_forces 		= - Fgp * (u[i][j] - us[i][j]) - Fgf * (u[i][j] - uf[i][j]) - Fgl * (u[i][j] - ul[i][j]);
			Source[i][j]			= Ext_forces - Turbulent_KE_term; 
			// Determination of Peclet number
			Peclet_x[i][j]			= (U_e + U_w)* dx_p/(mue_e + mue_w);
			Peclet_y[i][j]			= (V_n + V_s)* dy_p/(mue_n + mue_s);
			// Stores value of the previous iteration
			Previous_x[i][j]	= 0.00;		
		}
 	}
	int explicit_convection_flag = 1;
	int explicit_diffusion_flag = 0;	
/*------Loop for the determination of u* (f)---------------------------------------------*/
	while (bound>tol)
	{
		Conv_U++;
		/*------Begins thread level parallelization--------------------------------------*/
		#pragma omp parallel private(i,j,tid)
		{	
			/*------Domain decomposition for each thread is performed----------------*/
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
			/*------Traversal of the grid and computation loop begins---------------*/
			for (j = start_j; j <= finish_j; j++)
			{
				for (i = 0; i <= IMAX-1; i++)
				{	
					if(i == 0)						// Left Wall				
					{
						f[0][j]		=	u_left;
					}
					else if (i == IMAX-1)			// Right Wall	
					{
						f[IMAX-1][j]	=	u_right;
					}
					else
					{
						/*------Computation of values at cell-center---*/
						double vf_p	= 	0.50*(void_frac(i,j) + void_frac(i+1,j));
						/*------Define the velocity variables-----------*/
						double U_p	= u[i][j];
						double U_pn = u[i][j+1];
						double U_ps	= u[i][j-1];
						double U_pe	= u[i+1][j];
						double U_pw	= u[i-1][j];
						/*------Define the intermediate velocity variables-----------*/
						double f_p	= f[i][j];
						double f_pn	= f[i][j+1];
						double f_ps	= f[i][j-1];
						double f_pe	= f[i+1][j];
						double f_pw	= f[i-1][j];			
						//------Adjustment of values for boundaries--
						if (j == 1)						// Bottom Wall
						{
							f_ps	= 2.00 * u_bottom - f[i][j];
						}
						else if (j == JMAX-1)			// Top Wall
						{
							f_pn	= 2.00 * u_top - f[i][j];
						}
						//=============================================== CONVECTION COMPUTATION ========================================//
						//------Convection term(Implicit)
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
							DUUDX		= ( Fe[i][j] * ((U_p + U_pe) + alpha_x * (U_p - U_pe)) - Fw[i][j] * ((U_p + U_pw) + alpha_x * (U_pw - U_p)) );
							DUVDY		= ( Fn[i][j] * ((U_p + U_pn) + alpha_y * (U_p - U_pn)) - Fs[i][j] * ((U_p + U_ps) + alpha_y * (U_ps - U_p)) );					
						}
						else
						{
							DUUDX		= ( Fe[i][j] * ((f_p + f_pe) + alpha_x * (f_p - f_pe)) - Fw[i][j] * ((f_p + f_pw) + alpha_x * (f_pw - f_p)) );
							DUVDY		= ( Fn[i][j] * ((f_p + f_pn) + alpha_y * (f_p - f_pn)) - Fs[i][j] * ((f_p + f_ps) + alpha_y * (f_ps - f_p)) );
						}							
						double Convec_term		= DUUDX + DUVDY;		
						//================================================DIFFUSION COMPUTATION ========================================//
						//------Diffusion term(Implicit)
						double D2UDX2, D2UDY2;
						if( explicit_diffusion_flag == 1)
						{
							D2UDX2		= (Dw[i][j] * U_pw) + (De[i][j] * U_pe);
							D2UDY2		= (Ds[i][j] * U_ps) + (Dn[i][j] * U_pn);
						}
						else
						{
							D2UDX2		= (Dw[i][j] * f_pw) + (De[i][j] * f_pe);
							D2UDY2		= (Ds[i][j] * f_ps) + (Dn[i][j] * f_pn);	
						}
						double Diffusive_term	= D2UDX2 + D2UDY2;
						//================================================SOURCE COMPUTATION =========================================//
						double Source_term		= Source[i][j]; 
						//================================================PREVIOUS TIME TERM COMPUTATION ============================//
						double Previous_Time_term	= (vf_p*rho_g/dt1) * U_p;			
						//================================================COMPUTATION OF DENOMINATOR ================================//
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
						//================================================COMPUTATION of u*  =======================================//
						if (aP == 0)
						{
							f[i][j]	= 0.0;	// If the grid point is in the solid phase
						}
						else
						{			
							f[i][j] = (1.00 - omega_u) * f_p + (omega_u/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
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
		for (i = 0; i <= IMAX-1; i++)
		{
			for (j = 1; j <= JMAX-1; j++)
			{
				bound			= bound + fabs(	f[i][j] - Previous_x[i][j]);
				Previous_x[i][j]= f[i][j];
			}
		}
		if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) ) bound = -tol;
	}	
	/*------Convergence met. u* has been determined--------------------------------------*/		
	/*------Compute the maximum value of u* for monitoring and debugging-----------------*/	
	max_Ug		= 0.00;
	for (i = 0; i <= IMAX; i++)
	{
		for (j = 0;j <= JMAX; j++)
		{
			if (abs(f[i][j]) > max_Ug)
			{
				max_Ug = abs(f[i][j]);
			}
		}
	}	
	/*------Noting the cpu end time ------------------------------------------------------*/
	double time_after = clock();
	/*------Enforce a state-based boundary condition--------------------------------------*/
	enforce_domain_physics();
	/*------Noting the total cpu time for this computation--------------------------------*/
	cpu_time_U_gas = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);
}