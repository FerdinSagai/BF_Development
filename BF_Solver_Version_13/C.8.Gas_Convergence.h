void Gas_Convergence_Check()
{
	double DIV1;
	Residue_V	= 0.0;
	Residue_P	= 0.0;
	DIV1		= 0.0;
	for(int i = 0; i < IMAX; i++)
	{
		for(int j = 0; j < JMAX; j++)
		{
			DIV1 = fabs(u[i][j] - un[i][j]);
			Residue_V  = max(Residue_V, DIV1);
			DIV1 = fabs(v[i][j] - vn[i][j]);
			Residue_V  = max(Residue_V, DIV1);
			DIV1 = fabs(p[i][j] - pn[i][j]);
			Residue_P  = max(Residue_P, DIV1);
		}
	}
	
	for(int i = 0; i < IMAX; i++)
	{	for(int j = 0; j < JMAX; j++)
		{
			un[i][j]	=	u[i][j];
			vn[i][j]	=	v[i][j];
			pn[i][j]	=	p[i][j];
		}
	}
	Residue  = max(Residue_V, Residue_P);
	Residue_P_old = Residue_P;
}