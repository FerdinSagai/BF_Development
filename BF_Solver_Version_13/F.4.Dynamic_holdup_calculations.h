void dynamic()
{ 
	double omega_fd = 1.00; 
	double bound=10.0;
	double eps_max;
	
	time_hdb = clock();
	while (bound>tol)
	{
		Conv_epfd++;
		#pragma omp parallel private(tid)
		{
			/*------Private variables for each thread is declared-------------------*/
			double Fe, Fw, Fs, Fn;
			double aE, aW, aS, aN, aP;
			double hd_P, hd_N, hd_S, hd_W, hd_E;
			/*------Domain decomposition for each thread is performed---------------------------------------------*/
			tid			= omp_get_thread_num();
			nthreads	= omp_get_max_threads();
			int start_j = 1 + tid*int(JMAX/nthreads);
			int finish_j= start_j + JMAX/nthreads - 1;
			if(tid == nthreads-1)
			{
				finish_j = JMAX-1;
			}
			/*------Loop for the determination of hd--------------------------------------------------*/
			for (int i = 1; i <= IMAX-1; i++)
			{
				for (int j = start_j; j <= finish_j; j++)
				{	
					hd_P = epfd[i][j];
					hd_N = epfd1[i][j+1];
					hd_S = epfd1[i][j-1];
					hd_W = epfd1[i-1][j];
					hd_E = epfd1[i+1][j];
					/*------Adjustment of values for boundaries--*/			
					if(i == 1)										// Left Wall
					{
						hd_W	= epfd1[i][j]; 
						if(j == 1)									// Left-Bottom Corner
						{
							hd_S	= epfd1[i][j];
						}
						else if(j == JMAX-1)						// Left-Top Corner
						{
							hd_N	= epfd1[i][j];
						}
					}
					else if(i == IMAX-1)							// Right Wall
					{
						hd_E	= epfd1[i][j];
						if(j == 1)									// Right-Bottom Corner
						{
							hd_S	= epfd1[i][j];
						}
						else if(j == JMAX-1)						// Right-Top Corner
						{
							hd_N	= epfd1[i][j];
						}
					}
					else
					{
						if(j == 1)									// Bottom Wall
						{
							hd_S	= epfd1[i][j]; 
						}
						else if(j == JMAX-1)						// Top Wall
						{
							hd_N	= epfd1[i][j];
						}
					}
					/*------Declaration of gas velocity at cell centers-------------------------------------*/
					
 					double U_w = uf[i-1][j]; 
					double U_e = uf[i][j]; 
					double V_n = vf[i][j]; 
					double V_s = vf[i][j-1]; 
					
					Fe = 0.50 * U_e * (1.00/dx[i]);
					Fw = 0.50 * U_w * (1.00/dx[i]);
					Fn = 0.50 * V_n * (1.00/dy[j]);
					Fs = 0.50 * V_s * (1.00/dy[j]);

					aE = max(-Fe,0.0);  
					aW = max( Fw,0.0);  
					aN = max(-Fn,0.0);  
					aS = max( Fs,0.0);  
					
					aP = aE + aW + aN + aS + (Fe - Fw) + (Fn - Fs);

					eps_max=0.6*(1.0-void_frac_solid(i,j));
					
					if(fabs(aP)<1e-6) 
					{
						epfd1[i][j]=1e-10;
					}
					else if(epfs[i][j] == eps_max || fabs(epfs[i][j]-eps_max) < 1e-4)
					{
						epfd1[i][j]=1e-10;
					}
					else
					{
						epfd1[i][j]	=	(1.0 - omega_fd) * hd_P + (omega_fd/aP) * (aE*hd_E + aW*hd_W + aN*hd_N + aS*hd_S); 
					}
				}
			}
		} 
		#pragma omp barrier	
		enforce_domain_physics();
		bound=0.0;							
		for (int i = 1; i <= IMAX-1; i++)
		{
			for (int j = 1; j <= JMAX-1; j++)
			{
				bound		= bound + fabs(epfd1[i][j]-epfd2[i][j]);
				epfd2[i][j]	= relaxepfd*epfd1[i][j];
			}
		}
	}

	for (int i = 1; i <= IMAX-1; i++)
	{
		for (int j = 1; j <= JMAX-1; j++)
		{
			epfd[i][j]=epfd1[i][j];		
		}
	}
	
	max_epfd = 0.00;
	for (int i = 0; i <= IMAX; i++)
	{
		for (int j = 0; j <= JMAX; j++)
		{
			if (epfd[i][j] > max_epfd)
			{				
				max_epfd = epfd[i][j];
			}
		}
	}
	time_hda = clock();	
	enforce_domain_physics();
	cpu_time_dynamic_holdup = ((double) (time_hda - time_hdb)) / (nthreads*CLOCKS_PER_SEC);
}