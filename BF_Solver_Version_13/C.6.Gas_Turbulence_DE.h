void de_update()
{
	//------Local variables declaration-----------------------------------------------------------------------------------
	int i, j;
	double bound = 10.0;
	double Fe[IMAX+1][JMAX+1],Fw[IMAX+1][JMAX+1],Fn[IMAX+1][JMAX+1],Fs[IMAX+1][JMAX+1];
	double De[IMAX+1][JMAX+1],Dw[IMAX+1][JMAX+1],Dn[IMAX+1][JMAX+1],Ds[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1],Peclet_y[IMAX+1][JMAX+1];
	double Source[IMAX+1][JMAX+1];
	double Previous_de[IMAX+1][JMAX+1];
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

			double mut_P	= mut(i,j)/sigma_e;
			double mut_N	= mut(i,j+1)/sigma_e;
			double mut_S	= mut(i,j-1)/sigma_e;
			double mut_W	= mut(i-1,j)/sigma_e;
			double mut_E	= mut(i+1,j)/sigma_e;

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
			double KE_P	= ke[i][j];
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
			double Generation_term	= (DE_P/KE_P)*c_1*(mut_P*vf_p*(2.0*(sqr(DUDX)+sqr(DVDY)) + sqr(DUDY+DVDX)));
			//------Turbulent dissipation term computation----------------------------------------
			double Dissipation_term	= (DE_P/KE_P)*c_2*(rho_g * vf_p * DE_P);
			Source[i][j]			= Generation_term - Dissipation_term;
			// Determination of Peclet number
			Peclet_x[i][j]			= (U_e + U_w)* dx_p/(mut_e + mut_w);
			Peclet_y[i][j]			= (V_n + V_s)* dy_p/(mut_n + mut_s);
			// Stores value of the previous iteration
			Previous_de[i][j] = 0.00;
		}
	}
	int explicit_convection_flag = 1;
	int explicit_diffusion_flag = 0;
	//------Loop for the determination of de------------------------------------------------------------------------
	while (bound>tol)
	{
		Conv_de++;
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
						de1[1][j]			=	de[i][j];
					}
					else if (i == IMAX-1)				// Right Wall	
					{
						de1[IMAX-1][j]		=	de[i][j];
					}
					if(j == 1)							// Bottom Wall
					{
						de1[i][1]			=	de[i][j];			
					}
					else if(j == JMAX-1)				// Top Wall
					{
						de1[i][JMAX-1]		=	de[i][j];
					}
					else
					{
						//------Computation of values at cell-center---
						double vf_p		= void_frac(i,j);
						//------Define the variables with gas turbulent ke dissipation at previous timestep---
						double DE_P		= de[i][j];
						double DE_N		= de[i][j+1]; 
						double DE_S		= de[i][j-1]; 
						double DE_W		= de[i-1][j]; 
						double DE_E		= de[i+1][j];
						
						double DE1_P	= de1[i][j];
						double DE1_N	= de1[i][j+1]; 
						double DE1_S	= de1[i][j-1]; 
						double DE1_W	= de1[i-1][j]; 
						double DE1_E	= de1[i+1][j];
						//=============================================== CONVECTION COMPUTATION ========================================//
						//------Convection term(Implicit)
						double DUeDX, DVeDY;
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
							DUeDX	= ( Fe[i][j] * ((DE_P + DE_E) + alpha_x*(DE_P - DE_E)) - Fw[i][j] * ((DE_W + DE_P) + alpha_x*(DE_W - DE_P)) );
							DVeDY	= ( Fn[i][j] * ((DE_P + DE_N) + alpha_y*(DE_P - DE_N)) - Fs[i][j] * ((DE_S + DE_P) + alpha_y*(DE_S - DE_P)) ); 
						}
						else
						{
							DUeDX	= ( Fe[i][j] * ((DE1_P + DE1_E) + alpha_x*(DE1_P - DE1_E)) - Fw[i][j] * ((DE1_W + DE1_P) + alpha_x*(DE1_W - DE1_P)) );
							DVeDY	= ( Fn[i][j] * ((DE1_P + DE1_N) + alpha_y*(DE1_P - DE1_N)) - Fs[i][j] * ((DE1_S + DE1_P) + alpha_y*(DE1_S - DE1_P)) ); 
						}
						double Convec_term		= DUeDX + DVeDY;
						/*------Diffusion term using upwind (Implicit for stability)--------------*/
						double D2eDX2, D2eDY2;
						if( explicit_diffusion_flag == 1)
						{
							D2eDX2	= (Dw[i][j] * DE_W) + (De[i][j] * DE_E);
							D2eDY2	= (Ds[i][j] * DE_S) + (Dn[i][j] * DE_N);
						}
						else
						{
							D2eDX2	= (Dw[i][j] * DE1_W) + (De[i][j] * DE1_E);
							D2eDY2	= (Ds[i][j] * DE1_S) + (Dn[i][j] * DE1_N);
						}
						double Diffusive_term	= D2eDX2 + D2eDY2;					
						//------ Sources and Sinks
						double Source_term = Source[i][j];
						//------Previous time term computation------------------------------------------------
						double Previous_Time_term	= (vf_p * rho_g/dt1) * DE_P;
						//------Computation of dissipation of gas turbulent ke--------------------------------
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
						//================================================COMPUTATION of de1  =======================================
						if (aP == 0)
						{
							de1[i][j]	= 0.0;	// If the grid point is in the solid phase
						}
						else
						{
							de1[i][j]	= (1.0-omega_de)*DE_P + (omega_de/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
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
				bound			= bound + fabs(de1[i][j] - Previous_de[i][j]);
				Previous_de[i][j]	= relaxde * de1[i][j];
			}
		}
	}
	//------Convergence met. de has been determined-----------------------------------------------------------------------
	for (i=1;i<IMAX;i++)
	{
		for (j=1;j<JMAX;j++)
		{
			de[i][j] = de1[i][j];
		}
	}
	//------Compute the maximum value of de for monitoring and debugging--------------------------------------------------
	max_de		= 0.00;
	for (i=0; i<=IMAX; i++)
	{
		for (j=0;j<=JMAX;j++)
		{
			if (fabs(de[i][j]) > max_de)
			{
				max_de = fabs(de[i][j]);
			}
		}
	}
	//------Noting the cpu end time --------------------------------------------------------------------------------------
	double time_after = clock();
	//------Enforce a state-based boundary condition----------------------------------------------------------------------
	enforce_domain_physics();
	//------Noting the total cpu time for this computation----------------------------------------------------------------
	cpu_time_de_gas = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);
}