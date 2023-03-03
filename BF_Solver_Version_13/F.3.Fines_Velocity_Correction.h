/*------------------------------Correction of fines Velocity-----------------------------------------------*/
void uvf_correct ()
{	
	for (int i = 1; i <= IMAX-2; i++)
	{
		for (int j = 1; j <= JMAX-2; j++)
		{
			epfs_max		=	epfs_frac*(1.0-void_frac_solid(i,j));
			epfs_limit		=	0.95*epfs_max;
			if(epfs[i][j]>=epfs_limit)
			{
				uf[i][j]	= 0.0;
				vf[i][j]	= 0.0;
			}
			else
			{
				uf[i][j]	= ff[i][j];
				vf[i][j]	= gf[i][j];
			}
			uf[0][j]		= u_left;	
			uf[IMAX-1][j]	= u_right;	
		} 
		vf[i][0]			= v_bottom;
		vf[i][JMAX-1]		= vf[i][JMAX-2]; //Outflow boundary condition
	}
}