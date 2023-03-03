void vorticity_gas_calc()
{
	int i, j;
	double dx_p, dy_p;	
	double U_n,U_s,V_w, V_e;	
	double DUDY, DVDX;	
	
	for (i=1; i<=IMAX-1; i++)
	{
		for (j=1; j<=JMAX; j++)
		{
			dx_p=dx[i];
			dy_p=dy[j];
						
			U_n = 0.25*(u[i-1][j] + u[i][j] + u[i-1][j+1] + u[i][j+1]);
			U_s = 0.25*(u[i-1][j] + u[i][j] + u[i-1][j-1] + u[i][j-1]);
			V_w = 0.25*(v[i-1][j] + v[i-1][j-1] + v[i][j] + v[i-1][j]);
			V_e = 0.25*(v[i+1][j] + v[i+1][j-1] + v[i][j] + v[i-1][j]);
			
			DUDY = (U_n - U_s)/dy_p;
			DVDX = (V_e - V_w)/dx_p;
			
			vorticity_g[i][j] =  DVDX - DUDY;
		}
	}
}