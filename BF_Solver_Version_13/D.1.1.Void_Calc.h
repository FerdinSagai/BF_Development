void initialise_void_determination(void);
void nearest_neighbour(void);			
void void_location(void);
void void_identification(int &, int, int[]);
int Num_touches (int, int, int[]);
void Adding_Void(int &, int, int[]);
void particle_and_void(void);
void compute_void_properties(void);
void print_output_dl(void);
bool no_duplicate(double, double,int);
bool isneighbour(int, int);
int pip(int, int[]);
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
void Void_Definition_Module()
{
	printf("Void Definition Module.........Initiated\n");	
	printf("\t Neighbourhood definition.........Initiated\n");		
	initialise_void_determination();
	nearest_neighbour();		
	printf("\t Neighbourhood definition.........Completed\n");		
	
	printf("\t Void determination...............Initiated\n");		
	void_location();
	compute_void_properties();
	particle_and_void();
	printf("\t Void determination...............Completed\n");		
	
	FILE *fp;
	fp = fopen("2.5.Void_Properties.int", "wt");
		fprintf(fp,"%d\n",NUM_OF_VOIDS);
		
		for (int i = 0; i < NUM_OF_VOIDS; i++)
		{
			fprintf(fp,"%lf %lf %le %lf %lf %lf\n",VOID_POS[i][0],VOID_POS[i][1],VOID_POS[i][2],VOID_POS[i][3],VOID_POS[i][4],VOID_POS[i][5]);
			fprintf(fp,"%d\n",NUM_PARTICLE_VOIDS[i]);
			int p = 0;
			while(p < NUM_PARTICLE_VOIDS[i])
			{
				fprintf(fp,"%d\n",VOID_ENCIRCLED_BY[i][p]);	
				fprintf(fp,"%d\n",VOIDP[i][p][0]);	
				fprintf(fp,"%d\n",VOIDP[i][p][1]);	
				fprintf(fp,"%d\n",LEFT_NEIGH_VOIDP[i][p][0]);	
				fprintf(fp,"%d\n",LEFT_NEIGH_VOIDP[i][p][1]);	
				fprintf(fp,"%d\n",RIGHT_NEIGH_VOIDP[i][p][0]);
				fprintf(fp,"%d\n",RIGHT_NEIGH_VOIDP[i][p][1]);
				p++;
			}	
		}
		
	fclose(fp);
	
	printf("\t Print void location..............Initiated\n");		
	print_output_dl();	
	printf("\t Print void location..............Completed\n");		
	printf("Void Definition Module.........Complete\n");	
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
void initialise_void_determination()
{
	Num_suspended_particles = 0;
	
	for (int i = 0; i < NUM; i++)
	{
		for (int j = 0; j < MAX_NEIGHBOUR; j++)
		{
			CNUM[i][j]	=	NO_NEIGHBOUR;			
		}	
	}
	
	for(int i = 0;i < NVOIDS; i++)
	{
		suspended_solid_volume[i] = 0.00;
	}
	
	for(int i = 0; i < NUM; i++)
	{		
		suspended_particles[i] = 0;
	}
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
void nearest_neighbour()
{
	int i, j, nnb;
	time_NNb	= clock();	
	
	for (i = 0; i<NUM; i++)
	{
		nnb = 0;
		for (j = 0; j<NUM; j++)
		{
			if (j != i)
			{
				if ( (pow((xp[i]-xp[j]),2.0) + pow((yp[i]-yp[j]),2.0)) <= pow((mark*(d[i]+d[j])/2.0),2.0) )
				{
					if (nnb < MAX_NEIGHBOUR)
					{
						CNUM[i][nnb] = j;
						nnb++;
					}
				}
			}		
		}
	}
	time_NNa	= clock();			
	cpu_time_NN	= cpu_time_NN + ((double)(time_NNa - time_NNb) )/(nthreads*CLOCKS_PER_SEC);	
}	
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
void void_location()
{	
	int l,m;
	int count;
	int par_void[VOID_SURR_NUM];
	
	count	=	0;	
	double cpu_time = 0.0;

	for(int i = 0; i < NUM; i++)
	{	
		clock_t time_before = clock();
		printf("Particle Id : %d, void count: %d, %lf secs.\n", i, count, cpu_time);
	
		for(int o = 0; o < VOID_SURR_NUM; o++)
		{
			par_void[o] = NO_NEIGHBOUR;									// par_void array tracks the particle surrounding a void
		}
		
		int count_p			= 0; 										// Represents the index position of neighbour particle
		par_void[count_p]	= i;										// The particle i is considered as the first particle forming the void
		
		for(int j = 0; j < MAX_NEIGHBOUR; j++)
		{
			count_p				= 1;									// Represents the index position of neighbour particle	
			par_void[count_p]	= CNUM[par_void[count_p-1]][j];			// The neighbours of particle i are considered as second particle (j) forming the void
			if(par_void[count_p] != NO_NEIGHBOUR)
			{
				for(int k = 0;k < MAX_NEIGHBOUR; k++)
				{
					count_p				= 2;
					par_void[count_p]	= CNUM[par_void[count_p-1]][k];	// The neighbour of particle j is considered as the third particle (k) forming the void
					
					if(par_void[count_p] != par_void[0] && par_void[count_p] != NO_NEIGHBOUR)
					{						
						if( isneighbour (par_void[count_p], i) )		// Are i and k neighbours. If so, then a 3 particle void has been detected
						{	
							Adding_Void(count, count_p, par_void);		
						}
						else
						{	
						 	for(int o = count_p+1; o < VOID_SURR_NUM; o++)
							{
								par_void[o] = NO_NEIGHBOUR;
							}									
							void_identification(count, count_p, par_void);
						}
					}					
				}
			}			
		}			
		clock_t time_after = clock();
		cpu_time = ((double) (time_after - time_before)) / (nthreads*CLOCKS_PER_SEC);	
	}
	NUM_OF_VOIDS = count;
	printf("NUM OF VOIDS DETECTED = %d\n",NUM_OF_VOIDS);
	printf("NUM OF SUSPENDED PARTICLES = %d\n",Num_suspended_particles);
	//==============================================================================================================
	printf("End of Void Location Function\n");
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
void void_identification(int &count, int count_p, int par_void[])
{
	int i,j,k,l;
	bool equal, neighbour;

	//if(count_p < VOID_SURR_NUM)
	if(count_p < 6)
	{		
		count_p++;
		for(l = 0;l < MAX_NEIGHBOUR;l++)
		{		
			par_void[count_p] = CNUM[par_void[count_p-1]][l];						
			if(par_void[count_p] != NO_NEIGHBOUR)
			{				
				equal		=	false;
				neighbour	=	false;
									
				for(i = 1;i < count_p - 1; i++)
				{	
					// Checks if the particle is equal to any in the par void array
					if (par_void[count_p]==par_void[i])			
					{
						equal = true;						
					}										
					// Checks if the particle is a neighbour to any particle in the par void array
					if (isneighbour(par_void[count_p], par_void[i]))
					{
						neighbour = true;												
					}													
				}						
				
				if(!equal && !neighbour)
				{
					if(isneighbour(par_void[count_p], par_void[0])) // Closed loop found: void detected
					{
						Adding_Void(count, count_p, par_void);
					}
					else
					{							
						for(int o = count_p+1; o < VOID_SURR_NUM; o++)
						{
							par_void[o] = NO_NEIGHBOUR;
						}											
						void_identification(count, count_p, par_void);	
					}
				}
			}
		}	
	}
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
void Adding_Void (int &count, int count_p, int par_void[])
{
	double dummyx		=	0.0;
	double dummyy		=	0.0;
	bool valid_void		=  true;
	bool new_suspended_particle = true;
	for(int m = 0;m <= count_p; m++)
	{
		dummyx += xp[par_void[m]];
		dummyy += yp[par_void[m]];
	}
	dummyx /= (count_p + 1);
	dummyy /= (count_p + 1);
	
	
	if( (pip((count_p + 1),par_void) == -1) )
	{
		valid_void		=  true;
	}
	else
	{
		int pid				= pip((count_p + 1),par_void);
		int no_of_touches	= Num_touches(pid, count_p, par_void);
		
		if ( (no_of_touches == 0) || (no_of_touches == 1) )
		{
			valid_void		=  true;	
			for (int i = 0; i < Num_suspended_particles; i++)
			{
				if(pid == suspended_particles[i])
				{
					new_suspended_particle = false;
					break;
				}
			}
	
			if( new_suspended_particle)	
			{
				suspended_solid_volume[count] = suspended_solid_volume[count] + (4/3) * PI * pow( (0.50*d[pid]),3);
				suspended_particles[Num_suspended_particles] = pid;
				Num_suspended_particles++;
			}
		}
		else
		{
			valid_void		=  false;	
		}
	}
	
	if( no_duplicate(dummyx,dummyy,count) && valid_void)
	{			
		VOID_POS[count][0]	=	dummyx;
		VOID_POS[count][1]	=	dummyy;					
		for(int m = 0;m < VOID_SURR_NUM; m++)
		{
			VOID_ENCIRCLED_BY[count][m] = NO_NEIGHBOUR;
		}
		
		for(int m = 0;m <= count_p; m++)			
		{									
			VOID_ENCIRCLED_BY[count][m] = par_void[m];
		}
		count++;										
	}	
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
void compute_void_properties()
{
	int Num_particle_forming_void;
	for(int i = 0;i < NUM_OF_VOIDS; i++)
	{
		Num_particle_forming_void = 0;
		for(int k = 0;k < VOID_SURR_NUM; k++)
		{
			if (VOID_ENCIRCLED_BY[i][k] != NO_NEIGHBOUR)
			{
				Num_particle_forming_void++;
			}
		}
		
		double x_s[3],y_s[3];							
		double vol_s	=	0.0;

		for(int m = 0; m < Num_particle_forming_void; m++)
		{
			x_s[0]		=	xp[VOID_ENCIRCLED_BY[i][m]];
			y_s[0]		=	yp[VOID_ENCIRCLED_BY[i][m]];
			int count_s	=	1;
			
			for(int n = 0;n < Num_particle_forming_void; n++)
			{
				if(isneighbour(VOID_ENCIRCLED_BY[i][m], VOID_ENCIRCLED_BY[i][n]) && m!=n)
				{
					x_s[count_s]=xp[VOID_ENCIRCLED_BY[i][n]];
					y_s[count_s]=yp[VOID_ENCIRCLED_BY[i][n]];
					count_s++;
				}
			}
			double dot		= (x_s[1]-x_s[0])*(x_s[2]-x_s[0])+(y_s[1]-y_s[0])*(y_s[2]-y_s[0]);
			double mod_1	= sqrt(pow((x_s[1]-x_s[0]),2)+pow((y_s[1]-y_s[0]),2));
			double mod_2	= sqrt(pow((x_s[2]-x_s[0]),2)+pow((y_s[2]-y_s[0]),2));
			double theta	= fabs(acos(dot/(mod_1*mod_2)));
			vol_s 			+= theta * pow(d[0],3) * inv_12;	//volume of the sphere s1 inside the polyhedron Reference :: https://en.wikipedia.org/wiki/Spherical_wedge
		}
		
		// Area of Irregular polygon :: Reference  https://www.wikihow.com/Calculate-the-Area-of-a-Polygon
		int p			=	0;
		double pol_area	=	0;
		while(p < Num_particle_forming_void)
		{
			if(p != (Num_particle_forming_void-1) )
			{
				pol_area += (xp[VOID_ENCIRCLED_BY[i][p]]*yp[VOID_ENCIRCLED_BY[i][p+1]])-(xp[VOID_ENCIRCLED_BY[i][p+1]]*yp[VOID_ENCIRCLED_BY[i][p]]);
			}
			else
			{
				pol_area += (xp[VOID_ENCIRCLED_BY[i][p]]*yp[VOID_ENCIRCLED_BY[i][0]])-(xp[VOID_ENCIRCLED_BY[i][0]]*yp[VOID_ENCIRCLED_BY[i][p]]);
			}
			p++;
		}
		pol_area = fabs(0.50 * pol_area);
		
		VOID_POS[i][2]	=	fabs((pol_area) * d[0] - vol_s - suspended_solid_volume[i]);
		VOID_POS[i][3]	=	VOID_POS[i][2]/((pol_area)*d[0]);
		VOID_POS[i][4]	=	2.0 * pow( (VOID_POS[i][2]/(0.75 * PI)),(1.0/3.0) );	
		VOID_POS[i][5]	=	double(Num_particle_forming_void);
		NUM_PARTICLE_VOIDS[i] = Num_particle_forming_void;
	}	
}
//=====================================================================================================================================================
void particle_and_void()
{	
	int i, j, k, count;
	int countP, countV;
	int max_count = 0;
	int PARTICLE_BOUNDS_THE_VOID[NUM][VOID_SURR_NUM];	// PARTICLE_BOUNDS_THE_VOID[i][j]	:: Void number for which Particle i is at location j
	printf("Particle_and_void start\n");
	
	// INITIALIZATION 
	for(i = 0; i < NUM; i++)
	{	
		NUM_VOIDS_PARTICLE[i]	=	0;
		for(j = 0; j < VOID_SURR_NUM; j++)
		{
			PARTICLE_BOUNDS_THE_VOID[i][j]	=	NO_NEIGHBOUR;
		}
	}
	
	for(i = 0;i < NUM_OF_VOIDS; i++)
	{
		for(j = 0;j < VOID_SURR_NUM; j++)
		{
			for(k = 0;k < 3; k++)
			{
				VOIDP[i][j][k]=NO_NEIGHBOUR;				
			}			
		}
	}
	
	// INITIALIZATION 
	for(i = 0;i < NUM; i++)
	{	
		count = 0;
		for(j = 0;j < NUM_OF_VOIDS; j++)
		{	
			for(k = 0;k < VOID_SURR_NUM; k++)
			{
				if(i == VOID_ENCIRCLED_BY[j][k]) 
				{	
					if( count < VOID_SURR_NUM) 
					{	
						PARTICLE_BOUNDS_THE_VOID[i][count] = j;
						count = count+1;
					}
				}
			}
		}
		max_count				=	max(max_count,count);
		NUM_VOIDS_PARTICLE[i]	=	count;
	}
	
	printf("Max Count = %d\n",max_count);
	
	int num_of_voids_in_Q1;
	int num_of_voids_in_Q2;
	int num_of_voids_in_Q3;
	int num_of_voids_in_Q4;
	
	int voids_in_Q1[max_count];
	int voids_in_Q2[max_count];
	int voids_in_Q3[max_count];
	int voids_in_Q4[max_count];
	/*
	Q1 -> Quadrant 1
	Q2 -> Quadrant 2
	Q3 -> Quadrant 3
	Q4 -> Quadrant 4
	*/
	
	for(i = 0 ;i < NUM; i++)
	{		
		for(int k = 0 ;k < max_count; k++)
		{
			voids_in_Q1[k]	=	NO_NEIGHBOUR;
			voids_in_Q2[k]	=	NO_NEIGHBOUR;
			voids_in_Q3[k]	=	NO_NEIGHBOUR;
			voids_in_Q4[k]	=	NO_NEIGHBOUR;
		}		
		num_of_voids_in_Q1	=	0;
		num_of_voids_in_Q2	=	0;
		num_of_voids_in_Q3	=	0;
		num_of_voids_in_Q4	=	0;
		
		for(j = 0;j < VOID_SURR_NUM; j++)
		{
			if(PARTICLE_BOUNDS_THE_VOID[i][j] != NO_NEIGHBOUR)
			{				
				int void_id			= PARTICLE_BOUNDS_THE_VOID[i][j]; 
				double void_x		= VOID_POS[void_id][0];
				double void_y		= VOID_POS[void_id][1];
				double particle_x	= xp[i];
				double particle_y	= yp[i];
				
				if(void_x >= particle_x && void_y >= particle_y)
				{
					//1st quadrant
					voids_in_Q1[num_of_voids_in_Q1]	=	void_id;
					num_of_voids_in_Q1++;
				}
				else if(void_x < particle_x && void_y >= particle_y)
				{
					//2nd quadrant
					voids_in_Q2[num_of_voids_in_Q2]	=	void_id;
					num_of_voids_in_Q2++;
				}
				else if(void_x < particle_x && void_y < particle_y)
				{
					//3rd quadrant
					voids_in_Q3[num_of_voids_in_Q3]	=	void_id;
					num_of_voids_in_Q3++;
				}
				else if(void_x >= particle_x && void_y < particle_y)
				{
					//4th quadrant
					voids_in_Q4[num_of_voids_in_Q4]	=	void_id;
					num_of_voids_in_Q4++;
				}
				else 
				{
					// Should never be realized
				}
			}
		}

		count = 0;
		if(voids_in_Q1[0] != NO_NEIGHBOUR)
		{
		 	if(voids_in_Q1[1] != NO_NEIGHBOUR)
		 	{//1st quadrant
		 		if (VOID_POS[voids_in_Q1[1]][1]<VOID_POS[voids_in_Q1[0]][1])
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q1[1];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q1[0];
		 			count++;
		 		}
		 		else
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q1[0];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q1[1];
		 			count++;
		 		}
		 	}
		 	else
		 	{
		 		PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q1[0];
				count++;
		 	}
		}
		
		if(voids_in_Q2[0]!=NO_NEIGHBOUR)
		{
		 	if(voids_in_Q2[1]!=NO_NEIGHBOUR)
		 	{//2nd quadrant
		 		if (VOID_POS[voids_in_Q2[1]][1]>VOID_POS[voids_in_Q2[0]][1])
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q2[1];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q2[0];
		 			count++;
		 		}
		 		else
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q2[0];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q2[1];
		 			count++;
		 		}
		 	}
		 	else
		 	{
		 		PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q2[0];
				count++;
		 	}
		}
		
		if(voids_in_Q3[0]!=NO_NEIGHBOUR)
		{
		 	if(voids_in_Q3[1]!=NO_NEIGHBOUR)
		 	{//3rd quadrant
		 		if (VOID_POS[voids_in_Q3[1]][1]>VOID_POS[voids_in_Q3[0]][1])
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q3[1];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q3[0];
		 			count++;
		 		}
		 		else
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q3[0];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q3[1];
		 			count++;
		 		}
		 	}
		 	else
		 	{
		 		PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q3[0];
				count++;
		 	}
		}
		
		if(voids_in_Q4[0]!=NO_NEIGHBOUR)
		{
		 	if(voids_in_Q4[1]!=NO_NEIGHBOUR)
		 	{//4th quadrant
		 		if (VOID_POS[voids_in_Q4[1]][1]<VOID_POS[voids_in_Q4[0]][1])
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q4[1];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q4[0];
		 			count++;
		 		}
		 		else
		 		{
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q4[0];
		 			count++;
		 			PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q4[1];
		 			count++;
		 		}
		 	}
		 	else
		 	{
		 		PARTICLE_BOUNDS_THE_VOID[i][count]=voids_in_Q4[0];
				count++;
		 	}
		}
		NUM_VOIDS_PARTICLE[i]=count;
	}
	
	//=============================================================================================================================
	/*
		VOIDP[NVOIDS][Location][0,1,2]=NO_NEIGHBOUR;				
	*/	
	for(int i = 0;i < NUM_OF_VOIDS; i++)
	{
		countP	=	0;
		//for(int j = 0;j < VOID_SURR_NUM; j++)
		for(int j = 0;j < int(VOID_POS[i][5]); j++)
		{			
			int pid					=	VOID_ENCIRCLED_BY[i][j];
			if (pid != NO_NEIGHBOUR)
			{
				countV				=	0;
				for(int k = 0;k < NUM_VOIDS_PARTICLE[pid]; k++)
				{	
					int void_id			= PARTICLE_BOUNDS_THE_VOID[pid][k];
					double void_x		= VOID_POS[void_id][0];
					double void_y		= VOID_POS[void_id][1];
					double particle_x	= xp[pid];
					double particle_y	= yp[pid];
					int void_id_left	= 0;
					int void_id_right	= 0;
					
					if( (void_id != NO_NEIGHBOUR) && (void_id != i) && (void_y > particle_y) )
					{
						if( k == 0)
						{						
							void_id_left	= NUM_VOIDS_PARTICLE[pid] - 1;
							void_id_right	= k + 1;
						}
						else if (k == NUM_VOIDS_PARTICLE[pid] - 1)
						{
							void_id_left	= k - 1;
							void_id_right	= 0;							
						}
						else
						{
							void_id_left	= k - 1;
							void_id_right	= k + 1;
						}
						
						if(PARTICLE_BOUNDS_THE_VOID[pid][void_id_left] == i)
						{
							VOIDP[i][countP][countV]			=	void_id;
							LEFT_NEIGH_VOIDP[i][countP][countV]	=	1;
							countV++;
						}

						if(PARTICLE_BOUNDS_THE_VOID[pid][void_id_right] == i)
						{
							VOIDP[i][countP][countV]			=	void_id;
							RIGHT_NEIGH_VOIDP[i][countP][countV]=	1;
							countV++;
						}
					}
				}
				countP++;				
			}
		}
	}
	printf("End of Particle Void Location Function\n");
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
//												SUPPORT FUNCTIONS
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
bool isneighbour(int pid, int nid)
{
	bool neighbour=false;
	for (int j = 0;j < MAX_NEIGHBOUR; j++)
	{  
		if (CNUM[pid][j] == nid)
		{
			neighbour=true;
		}
	}
	return neighbour;
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
bool pnpoly(int npol, double x_poly[], double y_poly[], double x, double y)
{
	int i,j,k;	
	int c = 0;
	for (i = 0, j = npol-1; i < npol; j = i++)
	{
		if ( ( ((y_poly[i] <= y) && (y_poly[j] > y)) || ((y_poly[j] <= y) && (y_poly[i] > y)) ) )
		{
			if ( x < (x_poly[j] - x_poly[i]) * (y - y_poly[i]) / (y_poly[j] - y_poly[i]) + x_poly[i] )
			{
				 c	=	!c;		
			}
		}
	}
	return c;
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
int pip(int n_par, int par_pol[])
{	
	//checking whether any particle exists inside a polygon
	int i,j,k;
	int c = -1;
	int test = 0;
	double x_poly[n_par],y_poly[n_par];
	
	for(i = 0;i < n_par; i++)
	{
		x_poly[i] = xp[par_pol[i]];
		y_poly[i] = yp[par_pol[i]];		
	}	
	
	for(j = 0;j < NUM; j++)
	{	
		// Identifies particle j which is not part of the polygon formed by the new void
		bool no_coin=true;
		for(k = 0;k < n_par; k++)	//n_par = countp + 1
		{
			if(par_pol[k] == j)
			{
				no_coin=false;
			}
		}
		// says j is not part of the chain
		
		// Is that j inside the polygon? If yes, return the particle id 
		if(no_coin)
		{
			test = pnpoly(n_par, x_poly, y_poly, xp[j], yp[j]);			
			if(test == 1)
			{
				c = j;				
			}
		}		
	}
	return c;
}

// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
bool no_duplicate(double dummyx, double dummyy,int count)
{
   bool no_dup	= true;
   double tol	= 1e-5;
   
	if (count > 0)
	{
		for (int i = 0; i < count; i++)
		{  
			if (fabs(dummyx - VOID_POS[i][0]) < tol && fabs(dummyy - VOID_POS[i][1]) < tol)
			{
				no_dup = false;
			}        	
		}
	}
    return (no_dup);
}
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
int Num_touches (int pid, int count_p, int par_void[])
{
	int touches	=	0;
	
	for (int i = 0;i < MAX_NEIGHBOUR; i++)
	{
		int nid = CNUM[pid][i];							  // neighbour of particle inside void loop
		if (nid != NO_NEIGHBOUR)
		{
			for(int k = 0;k < count_p; k++)
			{
				if ( nid == par_void[k] )
				{			  
					touches++;
				}			
			}
		}
	}
	return touches;
}
// =================================================================================================================
/*
int Lowest_particle_of_void(int void_id)
{	
	int num_i;
	double y_min = 1000000.0;
		
	int Num_particles = int(VOID_POS[void_id][5]);
	for(int k = 0;k < Num_particles; k++)
	{
		int pid = VOID_ENCIRCLED_BY[void_id][k];
		if(y_min > yp[pid]) 
		{
			y_min	=	yp[pid];
			num_i	=	pid;
		}		
	}
	return(num_i);
}
*/