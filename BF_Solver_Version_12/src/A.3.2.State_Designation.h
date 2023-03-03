	// Setting the bricks states
		/*
			State
			= -1	:- Solid 
			= 0		:- Void Space 
			= 1		:- Inflow 
			= 2		:- Outflow 

		*/
void state_calc()
{
	int ii, jj;
	int ii1, ii2;
	int jj1, jj2;
	int ii01, ii02, jj01, jj02;
	int N;
	double  P_x;
	double  P_y;
	// Initial state of the domain as void
	for (int i = 0; i < IMAX; i++)
	{
		for (int j = 0; j < JMAX; j++)
		{
			State[i][j] = 0;
		}
	}
	
	// State of the exits
	for (int i = 0; i < IMAX; i++)
	{
		State[i][JMAX] = 2;
	}

	double AA_x, AA_y, BB_x, BB_y, CC_x, CC_y;
	if(n_structures>0)
	{
		for(int N = 1; N <= n_structures; N++)
		{
			if(structure_type[N]==1)		
			{
				AA_x = A_x[N];
				AA_y = A_y[N];
				BB_x = B_x[N];
				BB_y = B_y[N];
				CC_x = C_x[N];
				CC_y = C_y[N];
			
				for(int i = 0; i <= IMAX; i++)
				{
					for(int j = 0; j <= JMAX; j++)
					{
						P_x = x[i];
						P_y = y[j];
		
						double D_AB = (P_x - AA_x)*(BB_y - AA_y) - (P_y - AA_y)*(BB_x - AA_x);
						double D_CB = (P_x - CC_x)*(BB_y - CC_y) - (P_y - CC_y)*(BB_x - CC_x);
		
						if ((D_AB>0) &&(D_CB<0)) 
						{
							if (voidage_Geom[N] == 0.00) 
							{
								State[i][j]= -1;
							}
							else
							{
								State[i][j]= -2;	
							}	

						}
					}
				}
			}
			else if(structure_type[N]==2)		
			{
				CC_x = C_x[N];
				CC_y = C_y[N];
						
				for(int i = 0; i <= IMAX; i++)
				{
					for(int j = 0; j <= JMAX; j++)
					{
						P_x = x[i];
						P_y = y[j];

						double D = (CC_x - P_x)*(CC_x - P_x) +  (CC_y - P_y)*(CC_y - P_y);
						double R = radius[N]*radius[N];
						if (D<R) 
						{
							if (voidage_Geom[N] == 0.00) 
							{
								State[i][j]= -1;
							}
							else
							{
								State[i][j]= -2;	
							}	
						}
					}
				}
			}
			else if(structure_type[N]==3)		
			{
				CC_x = C_x[N];
				CC_y = C_y[N];
				double L_brick = L_block[N];
				double H_brick = H_block[N];
				
				ii1	= (CC_x - 0.5*L_brick)/dx[2];
				ii2	= (CC_x + 0.5*L_brick)/dx[2];
				jj1	= (CC_y - 0.5*H_brick)/dy[2];
				jj2	= (CC_y + 0.5*H_brick)/dy[2];
				for(int i = 0; i < IMAX; i++)
				{
					for(int j = 0; j < JMAX; j++)
					{
						if (i >= ii1 && i <= ii2)
						{
							if (j >= jj1 && j <= jj2)	
							{
								if (voidage_Geom[N] == 0.00) 
								{
									State[i][j]= -1;
								}
								else
								{
									State[i][j]= -2;	
								}
							}
						}
					}
				}
			}
		}
	}
	else if(n_structures < 0)
	{
		FILE* gp;
		gp =fopen(testname,"r");
		for(int j = JMAX; j >= 1; j--)
		{
			for(int i = 1; i <= IMAX; i++)
			{
				fscanf(gp,"%d",&State[i][j]);
				if(State[i][j]==1) State[i][j] = -1;
			}	
			fscanf(gp,"\n");			
		}
	}	
	//------------Tuyere---------------------------
	// TUYERES
	for (N = 1; N <= n_tuyeres; N++)
	{
		int j	= 1;
		double xx0=0.0;
		while (xx0 < Tuyere_protrusion[N])
		{
			xx0 = xx0+dx[j];
			j++;
		}
		ii01=1;
		ii02=j;//(int) floor(xx0/dx[1]);
		j=1;
		double yy0=0.0;
		while (yy0<Tuyere_height[N])
		{
			yy0=yy0+dy[j];
			j++;
		}

		jj01=j-1;//(int) floor(h/dy[1]);
		jj02=jj01;//(int) floor(yy0/dy[1]);
	
		// Interior
		//for (ii=ii01-1;ii<=ii02+2;ii++)
		for (ii=ii01-1;ii<=ii02+2;ii++)
		{ 
			for (jj=jj01;jj<=jj02;jj++)
			{
				State[ii][jj] = 1;
			}
		}
		// Walls
		for (ii=ii01-1;ii<=ii02;ii++)
		{ 
			State[ii][jj01-1] = -1;
			State[ii][jj02+1] = -1;
		}
	}	
}