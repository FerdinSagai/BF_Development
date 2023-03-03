void ke_update()
{
	//------Local variables declaration-----------------------------------------------------------------------------------
	int i, j;
	double bound = 10.0;
	double Fe[IMAX+1][JMAX+1],Fw[IMAX+1][JMAX+1],Fn[IMAX+1][JMAX+1],Fs[IMAX+1][JMAX+1];
	double De[IMAX+1][JMAX+1],Dw[IMAX+1][JMAX+1],Dn[IMAX+1][JMAX+1],Ds[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1],Peclet_y[IMAX+1][JMAX+1];
	double Source[IMAX+1][JMAX+1];
	double Previous_ke[IMAX+1][JMAX+1];
	//------Noting the cpu start time-------------------------------------------------------------------------------------
	double time_before = clock();
	//------Initialization the variable that stores the previous iteration value------------------------------------------
	for (i = 2; i <= IMAX-2; i++)
	{
		for (j = 2; j <= JMAX-2; j++)
		{
			//------Define cell size and distances-------------------------------------------------------------------------
			double dx_p, dx_e, dx_w;
			double dy_p, dy_n, dy_s;
			dx_p = dx[i];
			dy_p = dy[j];
			if (i==IMAX-1) dx_e=dx[i]; 
			else dx_e=0.5*(dx[i]+dx[i+1]);
			if (i==1) dx_w=dx[i]; 
			else dx_w=0.5*(dx[i-1]+dx[i]);
			if (j==JMAX-1) dy_n=dy[j]; 
			else dy_n=0.5*(dy[j]+dy[j+1]);
			if (j==1) dy_s=dy[j]; 
			else dy_s=0.5*(dy[j-1]+dy[j]);
			/*------Computation of values at face-centers-----------------------------*/			
			double vf_p		= void_frac(i,j);
			double vf_N		= void_frac(i,j+1); 
			double vf_S		= void_frac(i,j-1); 
			double vf_W		= void_frac(i-1,j); 
			double vf_E		= void_frac(i+1,j);

			double mut_P	= mut(i,j)/sigma_k;
			double mut_N	= mut(i,j+1)/sigma_k;
			double mut_S	= mut(i,j-1)/sigma_k;
			double mut_W	= mut(i-1,j)/sigma_k;
			double mut_E	= mut(i+1,j)/sigma_k;

			double vf_n	= 0.50 * (vf_p + vf_N); 
			double vf_s	= 0.50 * (vf_p + vf_S);
			double vf_w	= 0.50 * (vf_p + vf_W); 
			double vf_e	= 0.50 * (vf_p + vf_E);
			
			double mut_n= 0.50 * (mut_P + mut_N); 
			double mut_s= 0.50 * (mut_P + mut_S);
			double mut_w= 0.50 * (mut_P + mut_W); 
			double mut_e= 0.50 * (mut_P + mut_E);	
			//------Define the gas velocity variables---------------------------------
			double U_w	= u[i-1][j]; 
			double U_e	= u[i][j];
			double V_n	= v[i][j]; 
			double V_s	= v[i][j-1];
			//------ Convective Fluxes 
			Fe[i][j]	= (1.00/dx_p) * 0.50 * rho_g * vf_e * U_e;
			Fw[i][j]	= (1.00/dx_p) * 0.50 * rho_g * vf_w * U_w;
			Fn[i][j]	= (1.00/dy_p) * 0.50 * rho_g * vf_n * V_n;
			Fs[i][j]	= (1.00/dy_p) * 0.50 * rho_g * vf_s * V_s;
			//------ Diffusive Fluxes 
			De[i][j]	= (vf_e*mut_e/(dx_p*dx_e));
			Dw[i][j]	= (vf_w*mut_w/(dx_p*dx_w));
			Dn[i][j]	= (vf_n*mut_n/(dy_p*dy_n));
			Ds[i][j]	= (vf_s*mut_s/(dy_p*dy_s));
			//------ Sources and Sinks
			//------Turbulence generation term computation----------------------------------------
			
			double DE_P	= de[i][j];
			double U_nw	= u[i-1][j+1]; 
			double U_ne	= u[i][j+1];
			double U_sw	= u[i-1][j-1]; 
			double U_se	= u[i][j-1];
			
			double V_nw	= v[i-1][j]; 
			double V_ne	= v[i+1][j]; 
			double V_sw	= v[i-1][j-1]; 
			double V_se	= v[i+1][j-1]; 

			double U_n	= 0.25*(U_nw + U_ne + U_w + U_e);
			double U_s	= 0.25*(U_sw + U_se + U_w + U_e); 
			double V_w	= 0.25*(V_nw + V_sw + V_n + V_s);
			double V_e	= 0.25*(V_ne + V_se + V_n + V_s);
			
			double DUDX		= (1.0/dx_p)*(U_e - U_w);
			double DVDY		= (1.0/dy_p)*(V_n - V_s);
			double DUDY		= (1.0/dy_p)*(U_n - U_s);
			double DVDX		= (1.0/dx_p)*(V_e - V_w);
			double Generation_term = vf_p*mut_P*(2.0*(sqr(DUDX)+sqr(DVDY)) + sqr(DUDY + DVDX));
			//------Turbulent dissipation term computation----------------------------------------
			double Dissipation_term = vf_p * rho_g * DE_P;
			Source[i][j]			= Generation_term - Dissipation_term;
			// Determination of Peclet number
			Peclet_x[i][j]			= (U_e + U_w)* dx_p/(mut_e + mut_w);
			Peclet_y[i][j]			= (V_n + V_s)* dy_p/(mut_n + mut_s);
			// Stores value of the previous iteration
			Previous_ke[i][j] = 0.00;
		}
	}
	int explicit_convection_flag = 1;
	int explicit_diffusion_flag = 0;
	//------Loop for the determination of ke------------------------------------------------------------------------
	while (bound>tol)
	{
		Conv_KE++;
		#pragma omp parallel private(i,j,tid)
		{
			//------Domain decomposition for each thread is performed---------------------------------------------
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
			/*------Traversal of the grid and computation loop begins---------------------------------------------*/
			for (j = start_j; j <= finish_j; j++)
			{
				for (i = 1; i <= IMAX-1; i++)
				{
					if(i == 1)							// Left Wall				
					{
						ke1[1][j]			=	ke[i][j];
					}
					else if (i == IMAX-1)				// Right Wall	
					{
						ke1[IMAX-1][j]		=	ke[i][j];
					}
					if(j == 1)							// Bottom Wall
					{
						ke1[i][1]			=	ke[i][j];			
					}
					else if(j == JMAX-1)				// Top Wall
					{
						ke1[i][JMAX-1]		=	ke[i][j];
					}
					else
					{
						//------Computation of values at cell-center---
						double vf_p		= void_frac(i,j);
						//------Define the variables with gas turbulent ke at previous timestep---
						double KE_P		= ke[i][j];
						double KE_N		= ke[i][j+1]; 
						double KE_S		= ke[i][j-1]; 
						double KE_W		= ke[i-1][j]; 
						double KE_E		= ke[i+1][j];
						
						double KE1_P	= ke1[i][j];
						double KE1_N	= ke1[i][j+1]; 
						double KE1_S	= ke1[i][j-1]; 
						double KE1_W	= ke1[i-1][j]; 
						double KE1_E	= ke1[i+1][j];
						//=============================================== CONVECTION COMPUTATION ========================================//
						//------Convection term(Implicit)
						double DUKDX, DVKDY;
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
							DUKDX	= ( Fe[i][j] * ((KE_P + KE_E) + alpha_x*(KE_P - KE_E)) - Fw[i][j] * ((KE_W + KE_P) + alpha_x*(KE_W - KE_P)) );
							DVKDY	= ( Fn[i][j] * ((KE_P + KE_N) + alpha_y*(KE_P - KE_N)) - Fs[i][j] * ((KE_S + KE_P) + alpha_y*(KE_S - KE_P)) );
						}
						else
						{
							DUKDX	= ( Fe[i][j] * ((KE1_P + KE1_E) + alpha_x*(KE1_P - KE1_E)) - Fw[i][j] * ((KE1_W + KE1_P) + alpha_x*(KE1_W - KE1_P)) );
							DVKDY	= ( Fn[i][j] * ((KE1_P + KE1_N) + alpha_y*(KE1_P - KE1_N)) - Fs[i][j] * ((KE1_S + KE1_P) + alpha_y*(KE1_S - KE1_P)) );
						}
						double Convec_term		= DUKDX + DVKDY;
						/*------Diffusion term using upwind (Implicit for stability)--------------*/
						double D2KDX2, D2KDY2;
						if( explicit_diffusion_flag == 1)
						{
							D2KDX2	= (Dw[i][j] * KE_W) + (De[i][j] * KE_E);
							D2KDY2	= (Ds[i][j] * KE_S) + (Dn[i][j] * KE_N);
						}
						else
						{
							D2KDX2	= (Dw[i][j] * KE1_W) + (De[i][j] * KE1_E);
							D2KDY2	= (Ds[i][j] * KE1_S) + (Dn[i][j] * KE1_N);
						}
						double Diffusive_term	= D2KDX2 + D2KDY2;
						//------ Sources and Sinks
						double Source_term = Source[i][j];
						//------Previous time term computation------------------------------------------------
						double Previous_Time_term	= (vf_p * rho_g/dt1) * KE_P;
						//------Computation of gas turbulent ke-----------------------------------
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
						//================================================COMPUTATION of ke1  =======================================
						if (aP == 0)
						{
							ke1[i][j]	= 0.0;	// If the grid point is in the solid phase
						}
						else
						{
							ke1[i][j]	= (1.0-omega_ke)*KE_P + (omega_ke/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
						}
					}
				}
			}
		}
		//------Thread synchronization is performed at this point-----------------------------------------------------
		#pragma omp barrier	
		//------Enforce a state-based boundary condition--------------------------------------------------------------
		enforce_domain_physics();
		//------Computation of residue to determine convergence-------------------------------------------------------
		bound	=	0.0;
		for (i=1;i<=IMAX-1;i++)
		{
			for (j=1;j<=JMAX-1;j++)
			{
				bound			= bound + fabs(ke1[i][j] - Previous_ke[i][j]);
				Previous_ke[i][j]	= relaxke * ke1[i][j];
			}
		}
	}
	//------Convergence met. ke has been determined-----------------------------------------------------------------------
	for (i=1;i<IMAX;i++)
	{
		for (j=1;j<JMAX;j++)
		{
			ke[i][j] = ke1[i][j];
		}
	}
	//------Compute the maximum value of ke for monitoring and debugging--------------------------------------------------
	max_KE		= 0.00;
	for (i=0; i<=IMAX; i++)
	{
		for (j=0;j<=JMAX;j++)
		{
			if (fabs(ke[i][j]) > max_KE)
			{
				max_KE = fabs(ke[i][j]);
			}
		}
	}
	//------Noting the cpu end time --------------------------------------------------------------------------------------
	double time_after = clock();
	//------Enforce a state-based boundary condition----------------------------------------------------------------------
	enforce_domain_physics();
	//------Noting the total cpu time for this computation----------------------------------------------------------------
	cpu_time_KE_gas = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);
}