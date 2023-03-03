void uv_correct()
{
	/*------Local variables declaration---------------------------------------------------------------------------*/
	double dx_p, dx_e;
	double dy_p, dy_n;
	double vf_P;
	double f_P;
	double g_P;
	double P_P;
	
	double P_N;
	double P_E;
	double vf_N;
	double vf_E;
	
	double vf_p;
	/*------Loop for the determination of u from approximated velocity u*(f)--------------------------------------*/
	for (int i = 1;i <= IMAX-2; i++)
	{
		for (int j = 1;j <= JMAX-2; j++)
		{
			/*------Define cell size and distances -------------------------------------------------------*/
			dx_p	= dx[i];
			if (i==IMAX-1) dx_e=dx[i]; 
			else dx_e=0.5*(dx[i]+dx[i+1]);
			dy_p	= dy[j];		
			if (j==JMAX-1) dy_n = dy[j]; 
			else dy_n = 0.5*(dy[j] + dy[j+1]);
			/*------Define approximate velocity u*--------------------------------------------------------*/
			f_P 	= f[i][j];
			/*------Define gas pressure-------------------------------------------------------------------*/
			P_P 	= p[i][j];
			P_E 	= p[i+1][j];
			/*------Computation of values at cell-centers-------------------------------------------------*/			
 			vf_P	= void_frac(i,j);
			vf_E	= void_frac(i+1,j);
			/*------Computation of values at face-centers-------------------------------------------------*/			
			vf_p	= 0.50 * (vf_P + vf_E);
			/*------Compute the corrected gas X velocity--------------------------------------------------*/
			u[i][j] = f_P - (dt1/rho_g) * ( P_E - P_P )/dx_e;
			/*------Adjustment of values for boundaries-==================================================-*/
			u[0][j]			= u_left;	
			u[IMAX-1][j]	= u_right;
			
			
			/*------Define approximate velocity v*--------------------------------------------------------*/
			g_P 	= g[i][j];
			/*------Define gas pressure-------------------------------------------------------------------*/
			P_P 	= p[i][j];
			P_N 	= p[i][j+1];
			/*------Computation of values at cell-centers-------------------------------------------------*/			
			vf_P	= void_frac(i,j);
			vf_N	= void_frac(i,j+1); 
			/*------Computation of values at face-centers-------------------------------------------------*/					
			vf_p	= 0.50 * (vf_P + vf_N); 
			/*------Compute the corrected gas X velocity--------------------------------------------------*/
			v[i][j] = g_P - (dt1/rho_g) * ( P_N - P_P )/dy_n;
			/*------Adjustment of values for boundaries-==================================================-*/
			v[i][0]			= v_bottom;
			v[i][JMAX-1]	= v[i][JMAX-2];
		}
	} 
}




void u_correct()
{
	/*------Local variables declaration---------------------------------------------------------------------------*/
	double dx_p, dx_e;
	double f_P;
	double P_P, P_E;
	double vf_P, vf_E;
	double vf_p;
	/*------Loop for the determination of u from approximated velocity u*(f)--------------------------------------*/
	for (int i = 1;i <= IMAX-2; i++)
	{
		for (int j = 1;j <= JMAX-2; j++)
		{
			/*------Define cell size and distances -------------------------------------------------------*/
			dx_p	= dx[i];
			if (i==IMAX-1) dx_e=dx[i]; 
			else dx_e=0.5*(dx[i]+dx[i+1]);
			/*------Define approximate velocity u*--------------------------------------------------------*/
			f_P 	= f[i][j];
			/*------Define gas pressure-------------------------------------------------------------------*/
			P_P 	= p[i][j];
			P_E 	= p[i+1][j];
			/*------Computation of values at cell-centers-------------------------------------------------*/			
 			vf_P	= void_frac(i,j);
			vf_E	= void_frac(i+1,j);
			/*------Computation of values at face-centers-------------------------------------------------*/			
			vf_p	= 0.50 * (vf_P + vf_E);
			/*------Compute the corrected gas X velocity--------------------------------------------------*/
			u[i][j] = f_P - (dt1/rho_g) * ( P_E - P_P )/dx_e;
			/*------Adjustment of values for boundaries-==================================================-*/
			u[0][j]			= u_left;	
			u[IMAX-1][j]	= u_right;
		}
	} 
}
//----------------------------------------------------------------------------------------------------------------------
void v_correct()
{
	/*------Local variables declaration---------------------------------------------------------------------------*/
	int i, j;
	double dy_p, dy_n;
	double g_P;
	double P_P, P_N;
	double vf_P, vf_N;
	double vf_p;
	/*------Loop for the determination of u from approximated velocity u*(f)--------------------------------------*/
	for (i = 1;i <= IMAX-2; i++)
	{
		for (j = 1;j <= JMAX-2; j++)
		{
			/*------Define cell size and distances -------------------------------------------------------*/
			dy_p	= dy[j];		
			if (j==JMAX-1) dy_n = dy[j]; 
			else dy_n = 0.5*(dy[j] + dy[j+1]);
			/*------Define approximate velocity v*--------------------------------------------------------*/
			g_P 	= g[i][j];
			/*------Define gas pressure-------------------------------------------------------------------*/
			P_P 	= p[i][j];
			P_N 	= p[i][j+1];
			/*------Computation of values at cell-centers-------------------------------------------------*/			
			vf_P	= void_frac(i,j);
			vf_N	= void_frac(i,j+1); 
			/*------Computation of values at face-centers-------------------------------------------------*/					
			vf_p	= 0.50 * (vf_P + vf_N); 
			/*------Compute the corrected gas X velocity--------------------------------------------------*/
			v[i][j] = g_P - (dt1/rho_g) * ( P_N - P_P )/dy_n;
			/*------Adjustment of values for boundaries-==================================================-*/
			v[i][0]			= v_bottom;
			v[i][JMAX-1]	= v[i][JMAX-2];
		}
	} 
}