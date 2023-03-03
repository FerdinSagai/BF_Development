void p_update()
{
	/*------Local variables declaration--------------------------------------------------------*/
	int i, j;
	double bound = 100.0;
	double bound_old = 100.0;
	int loop = 0;
	//double Previous[IMAX+1][JMAX+1];
	/*------Noting the cpu start time ---------------------------------------------------------*/
	double time_before = clock();
	/*------Initialization the variable that stores the previous iteration value---------------*/
	for (i = 1; i <= IMAX-1; i++)
	{
		for (j=1; j<=JMAX-1; j++)
		{				
			Previous[i][j] = p[i][j];
		}
	}
	/*------Loop for the determination of p----------------------------------------------------*/
	while (bound > tol)
	{
		Conv_P++;
		/*------Begins thread level parallelization---------------------------------------*/		
		#pragma omp parallel private(i,j,tid)
		{
			/*------Private variables for each thread is declared---------------------*/
			double DUDX, DVDY;
			double D2PDX2, D2PDY2;
			double RHS, coefficient;
			double dx_p, dx_e, dx_w, dy_p, dy_n ,dy_s;
			double P_P, P_N, P_S, P_W, P_E;
			double f_w, f_e, g_n, g_s;
			/*------Domain decomposition for each thread is performed-----------------------------*/
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
			/*------Traversal of the grid and computation loop begins-----------------------------*/
			for (i = 1; i <= IMAX-1; i++)
			{
				for (j = start_j; j <= finish_j; j++)
				{			
					/*------Define cell size and distances--------------------------------*/
					dx_p=dx[i];
					dy_p=dy[j];
					if (i==IMAX-1) dx_e=dx[i]; 
					else dx_e=0.5*(dx[i]+dx[i+1]);
					
					if (i==1) dx_w=dx[i]; 
					else dx_w=0.5*(dx[i-1]+dx[i]);
					
					if (j==JMAX-1) dy_n=dy[j]; 
					else dy_n=0.5*(dy[j]+dy[j+1]);
					
					if (j==1) dy_s=dy[j]; 
					else dy_s=0.5*(dy[j-1]+dy[j]);
					/*------Define the approximate gas velocity variables-----------------*/
					f_w	= f[i-1][j]; 
					f_e	= f[i][j]; 
					g_n	= g[i][j]; 
					g_s	= g[i][j-1];
					/*------Define the variables with gas pressure at previous timestep---*/
					P_P	= p[i][j];
					P_N	= p[i][j+1]; 
					P_S	= p[i][j-1]; 
					P_W	= p[i-1][j]; 
					P_E	= p[i+1][j];
					/*------Adjustment of values for boundaries--------------------------*/
					if(i == 1)				// Left Wall
					{
						P_W	= p[i][j]; 
						if(j == 1)			// Left-Bottom Corner
						{
							P_S	= p[i][j]; 
						}
						else if(j == JMAX-1)		// Left-Top Corner
						{
							P_N	= 0;  
						}
					}
					else if(i == IMAX-1)			// Right Wall
					{
						P_E	= p[i][j];
						if(j == 1)			// Right-Bottom Corner
						{
							P_S	= p[i][j]; 
						}
						else if(j == JMAX-1)		// Right-Top Corner
						{
							P_N	= 0;
						}
					}
					else
					{
						if(j == 1)			// Bottom Wall
						{
							P_S = p[i][j]; 
						}
						else if(j == JMAX-1)		// Top Wall
						{
							P_N	= 0;
						}
					}	
					/*------Coefficient term computation---------------------------------*/
					double aE= (1.0/(dx_p*dx_e));
					double aW= (1.0/(dx_p*dx_w));
					double aN= (1.0/(dy_p*dy_n));
					double aS= (1.0/(dy_p*dy_s));
					double aP= (aE + aW + aN + aS);
					/*------Computation of pressure derivative terms---------------------*/
					D2PDX2 = aW*P_W + aE*P_E;
					D2PDY2 = aS*P_S + aN*P_N;  
					/*------Computation of divergence of approximate velocity(RHS)-------*/
  					DUDX	=	(f_e - f_w)/dx_p;
					DVDY	=	(g_n - g_s)/dy_p;
					RHS		=	(rho_g/dt1)*(DUDX + DVDY);					
					/*------Computation of pressure term---------------------------------*/
					p[i][j] = ((1.0 - omega_p) * P_P + (omega_p/aP) * (D2PDX2 + D2PDY2 - RHS));
				}
			}
		}
		/*------Thread synchronization is performed at this point------------------------------------*/	
		#pragma omp barrier		
		/*------Computation of residue to determine convergence--------------------------------------*/	
		bound	=	0.00;			
		max_P	=	0.00;
		for (i = 1; i <= IMAX-1; i++)
		{
			for (j = 1; j <= JMAX-1; j++)
			{
				bound = bound + fabs( p[i][j] - Previous[i][j] );
				Previous[i][j] = relaxp * p[i][j];
				if (p[i][j] > max_P)
				{
					max_P = p[i][j];
				}	
			}
		}
		
		double r = float(rand())/float(RAND_MAX);
		if (bound > bound_old)
		{
			omega_p	= 1.00 - r * 0.50;
		}
		else
		{
			omega_p	= 1.00 + r * 0.50;
		}			
		bound_old = bound;
		loop = loop + 1;
		if (loop == 1)			printf("Pressure Evolution :: Iteration\t    Bound\t\tMax P\t    Int.Relaxation\tExt.Relaxation\n");
		if ((loop%10000) == 0)	printf("Pressure Evolution :: %d\t%lf\t    %lf\t%lf\t%lf\n", loop, bound, max_P, omega_p, relaxp);
		
		if (loop == 100000)
		{
			printf("Pressure Evolution :: Forced Exit \n");
			break;
		}
	}
	/*------Convergence met. p has been determined------------------------------------------------------*/		
	/*------Noting the cpu end time ----------------------------------------------------------------------*/
	double time_after = clock();
	/*------Enforce a state-based boundary condition------------------------------------------------------*/
	enforce_domain_physics();
	/*------Noting the total cpu time for this computation------------------------------------------------*/
	cpu_time_Pressure = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);
}