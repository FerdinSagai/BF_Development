void grid_uniform()
{ 
//double x[IMAX+1];			// X co-ordinate of the grid point
//double y[JMAX+1];			// Y co-ordinate of the grid point
	//------------x-direction---------------------------
	x[0] = - 0.50 * Lx/(IMAX);
	for (int i = 1;i<=(IMAX);i++)
	{
		dx[i]=Lx/(IMAX-1);//
		if(i==1)
			x[i]=0.5*dx[i];
		else
			x[i]=x[i-1] + 0.5*(dx[i]+dx[i-1]);
	} 
	x[IMAX+1] = Lx + 0.50 * Lx/(IMAX);	
	dx[0]		=	dx[1];
	dx[IMAX+1]	=	dx[IMAX];
	//------------y-direction---------------------------
	y[0]=0.0;
	for (int  j = 1;j<=(JMAX);j++)
	{
		dy[j]=Ly/(JMAX-1); //
		if(j==1) 
			y[j]=0.5*dy[j];
		else 
			y[j]=y[j-1] + 0.5*(dy[j]+dy[j-1]);
	}
	y[JMAX+1] = Ly + 0.50 * Ly/(JMAX);  
	dy[0]		=	dy[1];
	dy[JMAX+1]	=	dy[JMAX];
	
	FILE *fp;
	//-------------------------GRID DETAILS FOR CFD SOLVER
	fp=fopen("2.1.Grid_details.int","wt");
	fprintf(fp,"%d %d\n",	IMAX, JMAX);	
	for (int j = 0; j <= JMAX; j++)
	{
		for (int i = 0; i <= IMAX; i++)
		{
			fprintf(fp,"%e %e %e %e\n",	x[i], y[j],	dx[i], dy[j] );	
		}
	}		
	fclose(fp);
} 
//------------------------------------------------------------------------------

void grid_non_uniform()       
{ 
/*
double x[IMAX+1];			// X co-ordinate of the grid point
double y[JMAX+1];			// Y co-ordinate of the grid point
	int i,j;
	double mindx=Lx;
	double maxdx=0.0;
	double mindy=Ly;
	double maxdy=0.0,yedge[JMAX];
	double x1[IMAX];
	
	printf("Grid Resolution :: IMAX=%d x JMAX=%d \n",IMAX, JMAX);
	printf("Grid Type :: NON-UNIFORM GRID CONCENTRATION IN RACEWAY\n");
	printf("---------------------------------------------------\n\n");

  for (i=0;i<=(IMAX-1);i++)
    {
      dx[i]=44.7284e-3;//Lx/(IMAX-1);//
      mindx=min(mindx,dx[i]);
      maxdx=max(maxdx,dx[i]); 
    } 

  for (i=IMAX-16;i<=IMAX-1;i++)
    {
      dx[i]=10.0e-3;//Lx/(IMAX-1);//
      mindx=min(mindx,dx[i]);
      maxdx=max(maxdx,dx[i]); 
    } 

  for (i=1;i<=(IMAX-1);i++)
    {
      dx[i]=14.49e-3;//Lx/(IMAX-1);//
      mindx=min(mindx,dx[i]);
      maxdx=max(maxdx,dx[i]); 
    } 

  for (i=1;i<=5;i++)
    {
      dx[i]=7.5e-3;//Lx/(IMAX-1);//
      mindx=min(mindx,dx[i]);
      maxdx=max(maxdx,dx[i]); 
    } 
for (i=6;i<=22;i++)
    {
      dx[i]=5.0e-3;//Lx/(IMAX-1);//
      mindx=min(mindx,dx[i]);
      maxdx=max(maxdx,dx[i]); 
    } 
  for (i=IMAX-5;i<=IMAX-1;i++)
    {
      dx[i]=5.0e-3;//Lx/(IMAX-1);//
      mindx=min(mindx,dx[i]);
      maxdx=max(maxdx,dx[i]); 
    } 

  double Rx=1.19;
  int i01,i0,Nx;
  i01=1;
  for (i=7;i<IMAX;i++)
  {
    if (dx[i]<dx[i-1])
	{ 
		// From coarse to fine grids
	  
	  Nx=int(log(dx[i-1]/dx[i])/log(Rx));//-1;
	  double min=dx[i];
	  i01=i-Nx;
	  if (i01<0) i01=1;
	  for (i0=i01;i0<i;i0++)
	    {
	      dx[i0]=dx[i0-1]/Rx;
              if (dx[i0]<min) dx[i0]=min;
	    }
	}

    if ((dx[i]>dx[i-1]) && (i01<i))
	{ 
		// From fine to coarse grids
	  
	  Nx=int(log(dx[i]/dx[i-1])/log(Rx));
	  double max=dx[i];
	  i01=i+Nx;
	  if (i01>IMAX) i01=IMAX;
	  for (i0=i-1;i0<i01-1;i0++)
	    {
	      dx[i0]=dx[i0-1]*Rx;
	      if (dx[i0]>max) dx[i0]=max;
	    }
	}
  
  }
  
  x[0]=0.0;
      x[1]=dx[1]/2;
      x1[1]=dx[1];
       for(i = 2;i<=IMAX-1;i++)
	   { 
        x[i]=x[i-1]+0.5*(dx[i]+dx[i-1]);
        x1[i]=x1[i-1]+dx[i];
       }
       x[IMAX]=x[IMAX-1]+0.5*dx[IMAX-1];
       
       //dx[IMAX-1]=dx[IMAX-1]-2*(x[IMAX]-Lx);
         x[IMAX]=x[IMAX-1]+0.5*dx[IMAX-1];
  
  for (j=0;j<=(JMAX-1);j++)
    {
      dy[j]=27.37e-3;//Ly/(JMAX-1); //
      mindy=min(mindy,dy[j]);
      maxdy=max(maxdy,dy[j]); 
    }

  for (j=1;j<=50;j++)
    {
      dy[j]=11.0e-3;//Ly/(JMAX-1); //
      mindy=min(mindy,dy[j]);
      maxdy=max(maxdy,dy[j]); 
    }

  for (j=1;j<=35;j++)
    {
      dy[j]=5.0e-3;//Ly/(JMAX-1); 
      mindy=min(mindy,dy[j]);
      maxdy=max(maxdy,dy[j]); 
    }
   
  double Ry=1.17;
  int j01,j0,Ny;
  j01=1;
  for (j=1;j<JMAX;j++)
  {
    if (dy[j]<dy[j-1])
    { 
		// From coarse to fine grids
	  
	  Ny=int(log(dy[j-1]/dy[j])/log(Ry));
	  double min=dy[j];
	  j01=j-Ny;
	  if (j01<0) j01=1;
	  for (j0=j01;j0<j;j0++)
	    {
	      dy[j0]=dy[j0-1]/Ry;
          if (dy[j0]<min) dy[j0]=min;
	    }
	}

    if ((dy[j]>dy[j-1]) && (j01<j))
	{ 
		// From fine to coarse grids
	  
	  Ny=int(log(dy[j]/dy[j-1])/log(Ry));
	  double max=dy[j];
	  j01=j+Ny;
	  if (j01>JMAX) j01=JMAX;
	  for (j0=j-1;j0<j01-1;j0++)
	    {
	      dy[j0]=dy[j0-1]*Ry;
		  if (dy[j0]>max) dy[j0]=max;
	    }
	}
  
  }
  y[0]=0.0;
  yedge[0]=0.0;
  yedge[1]=dy[1];
  
      y[1]=dy[1]/2;
       for(i = 2;i<=JMAX-1;i++)
	   { 
        y[i]=y[i-1]+0.5*(dy[i]+dy[i-1]);
        yedge[i]=yedge[i-1]+dy[i];
       }
       y[JMAX]=y[JMAX-1]+0.5*dy[JMAX-1];
	   */
}  
