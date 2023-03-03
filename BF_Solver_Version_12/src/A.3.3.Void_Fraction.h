void void_fraction()
{
	int i, j, k;	
	for(i=0;i<IMAX;i++)
	{
		for(j=0;j<JMAX;j++)
		{
			vfrac[i][j]=Domain_Voidage;
		}
	}			
	// Setting Void Fractions in inlet and outlet regions
	for(i = 0; i <= IMAX; i++)
	{
		for(j = 0; j <= JMAX; j++)
		{
			if(State[i][j] == -1) 
			{
				vfrac[i][j]	= 0.00;
			}
			else if (State[i][j] == 1)
			{
				vfrac[i][j]	= Raceway_Voidage;	
			}	
			else if (State[i][j] == 2)
			{
				vfrac[i][j]	= Raceway_Voidage;
			}
		}
	}
	// Setting void fraction for obstacles
	for(int N = 1; N<= n_structures; N++)
	{
		for( i =0; i < IMAX; i++)
		{
			for(j = 0; j < JMAX; j++)
			{
				if(State[i][j] == -2) 
				{
					vfrac[i][j] = voidage_Geom[N];
				}
			}
		}
	}
}