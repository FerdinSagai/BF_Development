//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void save_particle_position()
{
	int i,j;
	char name[80];
	double xx, yy;
	
	FILE *gp1;
	gp1 = fopen("3.1.Transient_Particle_Positions.dat", "a");	
	if(Particle_position_count == 1)
	{
		fprintf(gp1,"TITLE = Properties at transient states\n");
		fprintf(gp1, "VARIABLES=\"X\",\"Y\",\"Dp\"\n");
	}

	fprintf(gp1, "ZONE T = \"%d\"\n",nt);
	for (i=0; i<NUM; i++)
	{
		fprintf(gp1,"%e %e %e\n",fabs(xp[i]), fabs(yp[i]), d[i]);
	} 
	fclose(gp1);
	
	FILE *fp1;
	fp1=fopen("3.7.Results_Particle_Force.dat","wt");
	fprintf(fp1,"TITLE = Properties at Final state\n");
	fprintf(fp1, "VARIABLES=\"X (m)\",\"Y (m)\",\"F_gas-particles X\",\"F_gas-particles Y\"\n");
	fprintf(fp1, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",update_count, IMAX-1, JMAX-1);
	
	yy = 0.5*dy[1];
	for (j=1;j<=JMAX-1;j++)
	{
		xx  =0.5*dx[1];
		for (i=1;i<=IMAX-1;i++)
		{ 
			fprintf(fp1,"%e %e %e %e\n",xx,yy,Force_on_Par_X[i][j],Force_on_Par_Y[i][j]);
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}
	fclose(fp1);	 
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void save_final_particle_position()
{
	int i;
	FILE *fp3;
	fp3 = fopen("1.1.Particle_Positions.dat", "wt");
		for (i=0; i<NUM; i++)
		{
			fprintf(fp3,"%lf %lf %lf\n",xp[i],yp[i],d[i]);
		}
	fclose(fp3);
}