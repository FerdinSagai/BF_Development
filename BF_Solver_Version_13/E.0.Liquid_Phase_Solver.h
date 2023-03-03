int Lowest_particle_of_void(int);
void Read_Module(void);
void liquid_sources(void);
void liquid_phase_velocity(void);
void liquid_flow(void);
bool void_duplic(int, int, int);
bool multi_par_void(int, int, int);
bool three_par_void(int, int, int);
void print_output_dl(void);
void liquid_position_transient(void);
void save_liquid_position(int);
void print_droplets_and_rivulets(int, int, int, double);
void print_rivulet_breakage(int, double, double);
void Liquid_Flowfield_Output(void);

double Liquid_Volume_Determination(int, int, int, int, int);
double Calculate_Detatchment_Angle (int, int, int, double, double, double);
double Rivulet_Centroid_Angle(int, double);
double Rivulet_Leading_Edge_Angle(int, double, double, double);
double Rivulet_Trailing_Edge_Angle(int, double, double, double);
double Rivulet_Breakage_Calculation(int, int, int, double, double, double, double, double, double, double);
// =============================================================================================================================================================================
void Liquid_Phase_Solver()
{  
	printf("Liquid Flow Solver.........\n");	
	Read_Module();
	print_output_dl();
	printf("Liquid Flow Module.........Initiated\n");	

  	for(ITER_LIQUID = 1;ITER_LIQUID <= 300; ITER_LIQUID++)
	{
		printf("\n**********Iteration=%d************\n",ITER_LIQUID);
		liquid_sources();
		liquid_position_transient(); 
		
		liquid_flow();
		liquid_phase_velocity();
		Liquid_Flowfield_Output();
		
		Previous_ITER_LIQUID = ITER_LIQUID;
	}
	liquid_update_count++;
	save_liquid_position(liquid_update_count); 		
	printf("Liquid Flow Module.........Completed\n");	
	print_output_dl(); 
}
// =============================================================================================================================================================================
void Read_Module()
{
	FILE *fp;
	fp = fopen("2.5.Void_Properties.int", "r");
		fscanf(fp,"%d\n",&NUM_OF_VOIDS);
		for (int i = 0; i < NUM_OF_VOIDS; i++)
		{
			fscanf(fp,"%lf %lf %le %lf %lf %lf\n",&VOID_POS[i][0],&VOID_POS[i][1],&VOID_POS[i][2],&VOID_POS[i][3],&VOID_POS[i][4],&VOID_POS[i][5] );
			fscanf(fp,"%d\n",&NUM_PARTICLE_VOIDS[i]);
			int p = 0;
			while(p < NUM_PARTICLE_VOIDS[i])
			{
				fscanf(fp,"%d\n",&VOID_ENCIRCLED_BY[i][p]);
				fscanf(fp,"%d\n",&VOIDP[i][p][0]);	
				fscanf(fp,"%d\n",&VOIDP[i][p][1]);
				fscanf(fp,"%d\n",&LEFT_NEIGH_VOIDP[i][p][0]);	
				fscanf(fp,"%d\n",&LEFT_NEIGH_VOIDP[i][p][1]);	
				fscanf(fp,"%d\n",&RIGHT_NEIGH_VOIDP[i][p][0]);
				fscanf(fp,"%d\n",&RIGHT_NEIGH_VOIDP[i][p][1]);
				p++;
			}	
		}		
	fclose(fp);
}
// =============================================================================================================================================================================
void liquid_sources()	
{
	int i, j, k;
	for(int i = 0; i < NUM_OF_VOIDS; i++)
	{
		int vid	=	i;
		for(int j = 0;j < VOID_SURR_NUM; j++)
		{
			int pid								= VOID_ENCIRCLED_BY[vid][j];		// Provides particle Id of particle encircling void i at position j
			LIQUID_ON_PARTICLE_OLD[pid][vid]	= 0.0;								// Volume of liquid on particle k which is an encircling memeber of void i is initialized
			LIQVELX_OLD[pid][vid]				= 0.0;								// X velocity of liquid on particle k which is an encircling memeber of void i is initialized
			LIQVELY_OLD[pid][vid]				= 0.0;								// Y velocity of liquid on particle k which is an encircling memeber of void i is initialized
		}
	}
	/* This loop determines the voids that would be filled by the rotameter and fills them with 80% liquid*/
	
	for(int i = 0;i < NUM_OF_VOIDS; i++)
	{	
		int vid	=	i;
		double void_x	= VOID_POS[vid][0];										// The X coordinate of the centroid of void i
		double void_y 	= VOID_POS[vid][1];										// The Y coordinate of the centroid of void i
		
		for(int r = 0; r < NUM_OF_ROTAMETERS; r++)
		{	
			int rotameter_id = r;
			ROTAMETER_X	= SOURCE_OF_LIQUID[rotameter_id][0];									// The X coordinate of the rotameter r
			ROTAMETER_Y	= SOURCE_OF_LIQUID[rotameter_id][1];									// The Y coordinate of the rotameter r
		
			double stream_Thickness = ROTAMETER_OPENING * 0.5; 
			
			if( ( void_x > (ROTAMETER_X - stream_Thickness)) && ( void_x < (ROTAMETER_X + stream_Thickness)) )
			{
				if( (void_y < ROTAMETER_Y) && (void_y < (ROTAMETER_Y - 3 * d[0])) )
				{					
					int Lowest_particle_id							=	Lowest_particle_of_void(vid);		// Gives the particle id of the lowest particle surrounding void i
					LIQ_IN_VOID_OLD[vid]								=	0.8 * VOID_POS[vid][2];			// 80% of the void volume (VOID_POS[i][2]) is filled with liquid
					LIQUID_ON_PARTICLE_OLD[Lowest_particle_id][vid]	=	LIQ_IN_VOID_OLD[vid];				// The liquid is placed on the lowermost particle
					LIQVELX_OLD[Lowest_particle_id][vid]				=	0.0;		
					LIQVELY_OLD[Lowest_particle_id][vid]				=	-(ROTAMETER_FLOWRATE/60000.0)*pow((PI*ROTAMETER_OPENING*ROTAMETER_OPENING/4.0),-1);			
					
					if(void_x >= xp[Lowest_particle_id])							// If the void is larger on the right then liquid will move to the right along the lowest particle
					{
						GOING_RIGHT_OLD[Lowest_particle_id][vid] = 1.0;
					}
					else if(void_x <  xp[Lowest_particle_id])							// If the void is larger on the left then liquid will move to the right along the lowest particle
					{
						GOING_RIGHT_OLD[Lowest_particle_id][vid] = 0.0;
					}
				}
			}
		}
	}
}
// =============================================================================================================================================================================
void liquid_position_transient()
{
	FILE *gp;
	// Prints the Consolidated results
	gp=fopen("7.2.Transient_Liquid Position.dat","a");
	
	if (ITER_LIQUID == 1 )
	{
		fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Size\"\n");
	}
	
	if(Previous_ITER_LIQUID != ITER_LIQUID)
	{
		fprintf(gp, "ZONE T = \"Liquid Position at %d \"  \n", ITER_LIQUID);
	}
	for(int i = 0; i < NUM_OF_VOIDS; i++)
	{
		int vid = i;
		for(int j =0; j < VOID_SURR_NUM; j++)	
		{
			int pid = VOID_ENCIRCLED_BY[i][j];
			fprintf( gp, "%e %e %e\n",LIQ_POS_X[pid][vid], LIQ_POS_Y[pid][vid], LIQ_IN_VOID[vid] );
		}
	}
	fclose(gp);	
}
// =============================================================================================================================================================================
void initialise_liquid_phase()
{
	init_liquid_phase	=	true;	
	for (int i = 0; i <= IMAX; i++)
	{
		for (int j = 0; j <= JMAX; j++)
		{ 
			ul[i][j]	=	0.00;
			vl[i][j]	=	0.00;
			epld[i][j]	=	0.00;
			Area_GL[i][j]	=	0.00;
			Area_SL[i][j]	=	0.00;
		}
	}
	
	liquid_update_count = 0;
	initial_liquid_volume	=	0.38*(1.732/4.0)*pow(dp,3.0)/5.0;
	CONTACT_ANGLE		=	CONTACT_ANGLE*PI/180.0;
	Previous_ITER_LIQUID=0;
	Rivulet_Update = 1;
	Droplet_Update = 1;
	Rivulet_Break_Update=1;
}
// =============================================================================================================================================================================
void liquid_flow()
{
	double liquid_mass;
	double liquidvelocity;
	double liquid_radius;
	double rivuletlength;
	double THETA_L;
	
	double free_fall_length;
	double vertical_dist;
	double FFtime;
	double TFtime;

	double NETA;
	double DELTA;
	double ZETTA;
	
	double L1;
	
	double a_gl,a_sl;
	double Cd_gl,Cd_sl;
	
	double phi1,phi2;
	double liquid_vol;
	
	double d1,dj;
	//bool DUP[NUM][NVOIDS];
	double total_liquid_volume = 0.0;
	double liquid_volume_compartment1 = 0.00; 		
	double liquid_volume_compartment2 = 0.00; 		
	double liquid_volume_compartment3 = 0.00; 		
	for(int i = 0; i < NUM_OF_VOIDS; i++)
	{															
		int vid_r	= i;																		// Receiver void i
		for(int j = 0; j < VOID_SURR_NUM; j++)
		{
			int pid_r = VOID_ENCIRCLED_BY[i][j];												// Receiving particle of void i
			if( pid_r != NO_NEIGHBOUR)
			{											
				// CONSOLIDATION OF DATA FOR VOID i
				double liq			=	0.0;
				double shift		=	0.0;
				double going_right	=	0.0;
				double liq_vel_x	=	0.0;
				double liq_vel_y	=	0.0;
				double liq_pos_x	=	0.0;
				double liq_pos_y	=	0.0;
				for(int k = 0; k < VOID_SURR_NUM; k++)
				{
					int pid_s = VOID_ENCIRCLED_BY[i][k];										// Source particle of void i
					if(multi_par_void(vid_r,pid_s,pid_r) || three_par_void(i,j,k))
					{
						for(int l = 0; l < 2; l++)							// Loop over all possible source voids
						{		
							int vid_s	= VOIDP[i][k][l];										// Source void of source particle of void i
							//==========DETERMINATION OF LIQUID VOLUME IN VOID i======================================================================================
							liquid_vol	=	Liquid_Volume_Determination(k, l, vid_s, vid_r, pid_s);
							liquid_mass	=	(rho_l * liquid_vol);
							//========================================================================================================================================					
							if(liquid_vol > 0.0)
							{	
								//==============DETERMINATION OF DRAG COEFFICIENT BETWEEN LIQUID AND PARTICLES========================================================
								liquidvelocity	=	pow((pow(LIQVELX[pid_s][vid_s],2.0)+pow(LIQVELY[pid_s][vid_s],2.0)),0.5);
								if(liquidvelocity > CUT_OFF_LIQUID_VELOCITY)										
								{
									Cd_sl	=	541.0/((rho_l*liquidvelocity*d[pid_s])/mu_l) + 33.0;
								}
								else
								{
									Cd_sl	=	33.0;
								}
								Cd_gl	=	4.4;
								//==============DETERMINATINING IF LIQUID FLOW IS IN RIVULET OR DROPLET FORM===========================================================
								double RMAX			=	0.5 * VOID_POS[vid_r][4];				// Maximum radius of the receiver void i
								double RMAX_TEST	=	pow( (0.75*liquid_vol/PI), (1.0/3.0) );		// Radius of the liquid flowing into the void i
								
								if(RMAX_TEST > RMAX) 							//! LIQUID enters as a RIVULET
								{	
									liquid_radius	=	RMAX;
									rivuletlength	=	abs((liquid_vol - 0.33*PI*pow(liquid_radius,3.0)*(2.0+3.0*cos(CONTACT_ANGLE) - pow(cos(CONTACT_ANGLE),3.0)))/(pow(liquid_radius,2.0)*(PI-CONTACT_ANGLE + cos(CONTACT_ANGLE)*sin(CONTACT_ANGLE))));
									if(ITER_LIQUID == 300 ) print_droplets_and_rivulets(1, pid_s, pid_r, liquid_radius);
								}
								else										// LIQUID enters as a DROPLET											
								{
									liquid_radius	=	RMAX_TEST;													
									rivuletlength	=	0;
									if(ITER_LIQUID == 300 ) print_droplets_and_rivulets(0, pid_s, pid_r, liquid_radius);
								}										
								
								a_sl					= PI*pow(liquid_radius*sin(CONTACT_ANGLE),2.0) + rivuletlength*2.0*(liquid_radius*sin(CONTACT_ANGLE));
								a_gl					= PI*liquid_radius*liquid_radius + rivuletlength*2.0*liquid_radius;
								
								double void_x		= VOID_POS[i][0];
								double void_y		= VOID_POS[i][1];
								int cell_x				= void_x/dx[1];
								int cell_y				= void_y/dy[1];
								Area_SL[cell_x][cell_y]	= a_sl;
								Area_GL[cell_x][cell_y]	= a_gl;
								CA						= CA + a_sl;								
								//==============DETERMINATINING THE ANGLE AT WHICH THE LIQUID DETACHES================================================================
								THETA_L	= Calculate_Detatchment_Angle (vid_r, pid_s, pid_r, liquid_vol, liquid_radius, liquidvelocity);
								//==============DETERMINATINING THE LOCATION ON THE RECEIVEING PARTICLE ONTO WHICH THE LIQUID DROPS===================================
								d1	=	fabs( (0.5 * d[pid_s]) + liquid_radius * cos(CONTACT_ANGLE) );	
								dj	=	fabs( (0.5 * d[pid_r]) + liquid_radius * cos(CONTACT_ANGLE) );  
								if(LEFT_NEIGH_VOIDP[i][k][l] == 1)
								{	
									d1 = d1;							
								}
								if(RIGHT_NEIGH_VOIDP[i][k][l] == 1)
								{												
									d1 = -d1;
								}
								double XLPOS	=	(xp[pid_s] + d1 * cos(THETA_L));
								double YLPOS	=	yp[pid_r] + ( sqrt( (pow(dj,2.0)) - pow((XLPOS - xp[pid_r]),2.0) ) );									
								//====================================================================================================================================	
								for(int q = 0; q < VOID_SURR_NUM; q++)
								{
									if(VOID_ENCIRCLED_BY[i][q]!=NO_NEIGHBOUR 
									&& q != j 
									&& q != k 
									&& !(yp[VOID_ENCIRCLED_BY[i][q]] > VOID_POS[i][1]))
									{
										if (XLPOS >= xp[VOID_ENCIRCLED_BY[i][q]] - (0.5 * d[pid_r]) 
										&&  XLPOS <= xp[VOID_ENCIRCLED_BY[i][q]] + (0.5 * d[pid_r])
										&&  yp[VOID_ENCIRCLED_BY[i][q]] > yp[pid_r] 
										&&  yp[VOID_ENCIRCLED_BY[i][q]] < yp[pid_s] )
										{										
											//DUP_PAR[pid_r][i]	=	true;
										} 														
									}
								}																		
								//==============DETERMINATINING THE VERTICAL AND FREEFALL DISTANCE AND EFFECTS========================================================
								if (XLPOS >= xp[pid_r] - (0.5 * d[pid_r]) 											// 
								&&  XLPOS <= xp[pid_r] + (0.5 * d[pid_r]) )
								//   && DUP_PAR[VOID_ENCIRCLED_BY[i][j]][i]!=true) // && DUP[VOID_ENCIRCLED_BY[i][k]][VOIDP[i][k][l]]==false
								{									
									vertical_dist	=	yp[pid_s] - (d1 * sin(THETA_L) ) - YLPOS ;
									FFtime			=	0;
									TFtime			=	fabs(vertical_dist/liquidvelocity);	
									phi1				=	0.5*Cd_gl*a_gl*rho_g;
									phi2				=	0.5*Cd_sl*a_sl*rho_l;
									
									if ( vertical_dist <= (0.5 * rivuletlength) )												// Free fall does not occurs
									{												
										free_fall_length	=	0.0;										
										liq_vel_x	=	0;	
										liq_vel_y	=	- pow( (rho_l*liquid_vol*gravity)/(phi1 + phi2) , 0.5);
									}
									else																						// Free fall occurs (Always true for droplet)
									{
										free_fall_length	=	fabs( vertical_dist - (0.5 * rivuletlength) );
							
										double LVX0	=	LIQVELX[pid_s][vid_s];
										double LVY0	=	LIQVELY[pid_s][vid_s];
							
										double void_x		= VOID_POS[i][0];
										double void_y		= VOID_POS[i][1];
										
										int cell_x = void_x/dx[1];
										int cell_y = void_y/dy[1];
										
										double u_gl	= u[cell_x][cell_y]	- ul[cell_x][cell_y];
										double v_gl	= v[cell_x][cell_y]	- vl[cell_x][cell_y];
										
										double Gas_drag_ForceX			=	phi1*pow(u_gl,2);
										double Gas_drag_ForceY			=	phi1*pow(v_gl,2);
										
										if(liquidvelocity > CUT_OFF_LIQUID_VELOCITY) 
										{	
											FFtime	=	fabs(free_fall_length/liquidvelocity);
											TFtime	=	fabs(vertical_dist/liquidvelocity);
										}
										else
										{
											FFtime	=	fabs(sqrt(2*free_fall_length/gravity));
											TFtime	=	fabs(sqrt(2*(vertical_dist)/gravity));
										}
										// Gas Effect needs to accounted for
										liq_vel_x=	LVX0	+ FFtime * (Gas_drag_ForceX/liquid_mass);		
										liq_vel_y=	LVY0	+ FFtime * (gravity - (Gas_drag_ForceY/liquid_mass) );
									}									
									//==============DETERMINATION OF DROPLET DYNAMICS/RIVULET BREAKAGE ON RECEIVING PARTICLE===========================
									if(rivuletlength > 0)														//! RIVULET..............
									{
										// Checking for overhangs
										if ( rivuletlength >  PI * 0.5 * d[pid_r] ) 							
										{
											L1	=	 PI * 0.5 * d[pid_r];
										}
										else
										{
											L1	=	rivuletlength;
										}								
												
										NETA	=	Rivulet_Centroid_Angle(pid_r,XLPOS);
										//===========================================================================================================
										int Side_flag = 0;
										if(XLPOS <= xp[pid_r]) 													//! RIVULET fallen towards left.....
										{
											Side_flag = - 1; 
										}
										else if(XLPOS > xp[pid_r]) 																	//! RIVULET fallen towards right....
										{
											Side_flag = 1; 										
										}
										
										if(L1 == rivuletlength)
										{
											if(XLPOS <= xp[pid_r]) 					//!  RIVULET GOES COMPLETELY TOWARDS LEFT.......
											{
												shift	=	0.0;
											}
											else if(XLPOS > xp[pid_r]) 				//!  RIVULET GOES COMPLETELY TOWARDS RIGHT.......
											{
												shift	=	1.0;									
											}												
										}	
										else 							//! LARGE RIVULET => BREAKAGE POSSIBLE
										{										
											shift = Rivulet_Breakage_Calculation(i, pid_r, Side_flag, NETA, liquid_radius, rivuletlength, L1, liquidvelocity, Cd_sl, Cd_gl);
										}
									}
									else if(rivuletlength <= 0)													//! DROPLET
									{
										if(XLPOS <= xp[pid_r]) 													//! DROPLET fallen towards left.....
										{
											shift	=	0.0;				
										}
										else																	//! DROPLET fallen towards right....
										{	
											shift	=	1.0;				
										}
									}	
							
									// CONSOLIDATION OF DATA
									liq 		+=	liquid_vol;
									going_right +=	(shift		*	liquid_vol);																		
									liq_vel_x 	+=	(liq_vel_x	*	liquid_vol);
									liq_vel_y 	+=	(liq_vel_y	*	liquid_vol);
									liq_pos_x 	+=	(XLPOS		*	liquid_vol);
									liq_pos_y 	+=	(YLPOS		*	liquid_vol);				
								}
							}								
						}
					}
				}
				
				if (liq > 0.0)
				{
					LIQUID_ON_PARTICLE_OLD[pid_r][vid_r]	=	liq;
					total_liquid_volume +=liq; 
					
					if	( (xp[pid_r]>0.0) && (xp[pid_r]<=0.02)
					&&(yp[pid_r]< 3*d[0]) )
					{
						liquid_volume_compartment1 +=liq; 		
					}
						
					if	( (xp[pid_r]>0.02) && (xp[pid_r]<0.04)
						&&(yp[pid_r]< 3*d[0]) )
					{
						liquid_volume_compartment2 +=liq; 		
					}
					
					if	( (xp[pid_r]>0.04) && (xp[pid_r]<0.06)
						&&(yp[pid_r]< 3*d[0]) )
					{
						liquid_volume_compartment3 +=liq; 		
					}						
				}			
				
				if(LIQUID_ON_PARTICLE_OLD[pid_r][vid_r] > 0.0)
				{								
					GOING_RIGHT_OLD[pid_r][vid_r]	=	going_right/LIQUID_ON_PARTICLE_OLD[pid_r][vid_r];
					LIQ_POS_X[pid_r][vid_r]			=	VOID_POS[vid_r][0];
					LIQ_POS_Y[pid_r][vid_r]			=	VOID_POS[vid_r][1];
					
					if(GOING_RIGHT_OLD[pid_r][vid_r] > 1.0 && GOING_RIGHT_OLD[pid_r][vid_r] != NO_NEIGHBOUR)
					{
						GOING_RIGHT_OLD[pid_r][vid_r]	=	1.0;									
					}
						
					if(GOING_RIGHT_OLD[pid_r][vid_r] < 0.0 && GOING_RIGHT_OLD[pid_r][vid_r] != NO_NEIGHBOUR) 
					{
						GOING_RIGHT_OLD[pid_r][vid_r] = 0.0;									
					}							
					LIQVELX_OLD[pid_r][vid_r]	=	liq_vel_x/LIQUID_ON_PARTICLE_OLD[pid_r][vid_r];
					LIQVELY_OLD[pid_r][vid_r]	=	liq_vel_y/LIQUID_ON_PARTICLE_OLD[pid_r][vid_r];
				}	
			}						
		}
		double liq_in_void = 0.0;
		for(int j = 0; j < VOID_SURR_NUM; j++)
		{
			int pid = VOID_ENCIRCLED_BY[vid_r][j];
			if(LIQUID_ON_PARTICLE_OLD[pid][vid_r] > 0.0 
			&& pid != NO_NEIGHBOUR 
			&& GOING_RIGHT_OLD[pid][vid_r] != NO_NEIGHBOUR)
			{
				liq_in_void	+=	LIQUID_ON_PARTICLE_OLD[pid][vid_r];
			}
		}
		LIQ_IN_VOID_OLD[vid_r] = liq_in_void;
	}	
	
// END OF LARGE LOOP
	for(int i = 0; i < NUM_OF_VOIDS; i++)
	{
		LIQ_IN_VOID[i] = LIQ_IN_VOID_OLD[i];
		double void_x		= VOID_POS[i][0];
		double void_y		= VOID_POS[i][1];

		int cell_x = void_x/dx[1];
		int cell_y = void_y/dy[1];

		epld[cell_x][cell_y] =  LIQ_IN_VOID[i];
		for(int j = 0; j < VOID_SURR_NUM; j++)
		{
			int pid = VOID_ENCIRCLED_BY[i][j];			
			if(GOING_RIGHT_OLD[pid][i] != NO_NEIGHBOUR && pid != NO_NEIGHBOUR)
			{
				GOING_RIGHT[pid][i]			=	GOING_RIGHT_OLD[pid][i];				
				LIQUID_ON_PARTICLE[pid][i]	=	LIQUID_ON_PARTICLE_OLD[pid][i];			
				LIQVELX[pid][i]				=	LIQVELX_OLD[pid][i];
				LIQVELY[pid][i]				=	LIQVELY_OLD[pid][i];
			}
		}
	}

	printf("contact_area = %e\n",CA);
	printf("Total liquid in Bed = %e\n",total_liquid_volume);
	printf("Compartmental Liquid Volume = %e %e %e\n",liquid_volume_compartment1, liquid_volume_compartment2, liquid_volume_compartment3);
}
// =============================================================================================================================================================================
void liquid_phase_velocity()
{
	for(int i = 0; i < NUM_OF_VOIDS; i++)
	{	
		int vid		= i;
		int void_x	= VOID_POS[vid][0];										// The X coordinate of the centroid of void i
		int void_y	= VOID_POS[vid][1];										// The Y coordinate of the centroid of void i
		for(int k = 0; k < VOID_SURR_NUM; k++)
		{
			int pid = VOID_ENCIRCLED_BY[vid][k];
			if( (pid != NO_NEIGHBOUR) && (LIQUID_ON_PARTICLE[pid][vid] != 0.00) )
			{
				double pid_x	= xp[pid];
				double pid_y	= yp[pid] ;
				
				int cell_x = pid_x/dx[1];
				int cell_y = pid_y/dy[1];
				
				ul[cell_x][cell_y]	= LIQVELX[pid][vid];
				vl[cell_x][cell_y]	= LIQVELY[pid][vid];
			}
		}
	}		
}
// =============================================================================================================================================================================
void Liquid_Flowfield_Output()
{
	int ijk = 10;
	FILE *gp;
	gp=fopen("7.2.Results_Steady_State_Liquid.dat","wt");
	fprintf(gp,"TITLE = Properties at Final state\n");
	fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Ul\",\"Vl\"\n");
	fprintf(gp, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",ijk, IMAX-1, JMAX-1);
	
	double yy = 0.5*dy[1];
	for (int j = 1;j <= JMAX-1; j++)
	{
		double xx  =0.5*dx[1];
		for (int i = 1;i <= IMAX-1;  i++)
		{			
			fprintf(gp,"%e %e %e %e\n",
			xx, 
			yy, 
			ul[i][j], 
			vl[i][j]);						
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}
	fclose(gp);
}
// =============================================================================================================================================================================
bool void_duplic(int VOIDI, int vi[], int count_void) 
{
	int i;
	bool dup = false;
	for(i=0; i<count_void; i++)
	{
		if(VOIDI==vi[i])
		{
			dup = true;
		}
			
	}
	return(dup);
}
//! ********************************************************************
bool multi_par_void(int vid,int pid_s,int pid_r)
{
	bool t	=	false;
	
	if( (pid_s != pid_r) 											// If Receiver and Source particle are not the same
	&& (pid_s != NO_NEIGHBOUR) 										// AND source particle exists
	&& (pid_r != NO_NEIGHBOUR) 										// AND Receiver particle exists
	&& (yp[pid_s] >= yp[pid_r]) 									// AND Source particle is higher than the receiver particle
	//&& !(yp[pid_s] > VOID_POS[vid][1]	&& yp[pid_r] > VOID_POS[vid][1])	// AND Source particle AND Receiver particle are lower than void centroid
	)
	{
		t=true;
	}
	return(t);
}
//! ******************************************************************
bool three_par_void(int i, int j, int k)
{
	bool t	=	false;
	int vid_r = i;
	int pid_r = VOID_ENCIRCLED_BY[i][j];
	int pid_s = VOID_ENCIRCLED_BY[i][k];
	
	double void_x = VOID_POS[vid_r][0];
	double void_y = VOID_POS[vid_r][1];
	
	if( (pid_s != pid_r) 						// If receiver and source particle are not the same
	&&  (pid_s != NO_NEIGHBOUR) 				// AND source particle is not void
	&&  (pid_r != NO_NEIGHBOUR)					// AND Receiver particle exists
	&&  fabs(VOID_POS[vid_r][5] - 3.0) < 1e-6		// Three particles in a void ? 
//	&&  (yp[pid_s] >= yp[pid_r]) 				// AND Source particle is higher than the receiver particle  // Error? Clarify
//	&&  (yp[pid_s] < yp[pid_r]) 				// AND Source particle is higher than the receiver particle  // Error? Clarify
	&&  (yp[pid_r] < void_y )			// AND Source particle is lower than void centroid
	&&  (yp[pid_s] < void_y )			// Receiver particle is lower than void centroid
	)
	{
		t	=	true;
/* 		if( void_x < VOID_POS[VOIDP[i][k][0]][0])
		{
			if( (1 - GOING_RIGHT[pid_s][VOIDP[i][k][0]]) > 0.0)
			{
				t	=	true;
			}
		}
		
		if( void_x  > VOID_POS[VOIDP[i][k][0]][0])
		{
			if( (GOING_RIGHT[pid_s][VOIDP[i][k][0]]) > 0.0)
			{
				t	=	true;
			}
		}		 */
	} 	
	return(t);								
}
// =============================================================================================================================================================================
void print_droplets_and_rivulets(int DRindicator, int pid_s, int pid_r, double radius)
{
	int i,j;
	FILE *gp1;
	FILE *gp2;

	// Prints the Consolidated results		
		if( DRindicator == 0)
		{
			gp1	=	fopen("7.3a.Droplets_Results.dat","a");	
				if (Droplet_Update == 1) 
				{
					fprintf(gp1, "VARIABLES=\"X (m)\",\"Y (m)\",\"Size\"\n");		
				}
				fprintf(gp1, "ZONE T = \"Droplets Position\"\n");
				fprintf(gp1,"%e %e %e\n",xp[pid_s],yp[pid_s],radius);
				fprintf(gp1,"%e %e %e\n",xp[pid_r],yp[pid_r],radius);
			fclose(gp1);	
			Droplet_Update++;
		}
		else if(DRindicator == 1)
		{
			gp2	=	fopen("7.3b.Rivulets_Results.dat","a");	
				if (Rivulet_Update == 1) 
				{
					fprintf(gp2, "VARIABLES=\"X (m)\",\"Y (m)\",\"Size\"\n");		
				}	
				fprintf(gp2, "ZONE T = \"Rivulets Position\"\n");
				fprintf(gp2,"%e %e %e\n",xp[pid_s],yp[pid_s],radius);
				fprintf(gp2,"%e %e %e\n",xp[pid_r],yp[pid_r],radius);
			fclose(gp2);
			Rivulet_Update++;
		}	
}
// =============================================================================================================================================================================
void print_rivulet_breakage(int vid_r, double l_rivulet, double radius)
{
	int i,j;
	FILE *gp1;
	
	double void_x = VOID_POS[vid_r][0];
	double void_y = VOID_POS[vid_r][1];

	gp1	=	fopen("7.3c.Rivulet_Breakage_Results.dat","a");	
		if (Rivulet_Break_Update == 1) 
		{
			fprintf(gp1, "VARIABLES=\"X (m)\",\"Y (m)\",\"Size\"\n");		
		}
		if(Previous_ITER_LIQUID != ITER_LIQUID)
		{
			fprintf(gp1, "ZONE T = \"Rivulet Break Position at %d \"  \n", ITER_LIQUID);
		}
		fprintf(gp1,"%e %e %e\n",void_x,void_y,l_rivulet);
		fprintf(gp1,"%e %e %e\n",void_x,void_y,radius);
	fclose(gp1);	
	Rivulet_Break_Update++;
}
// =============================================================================================================================================================================
void print_output_dl()
{
	int i,j;
	FILE *gp1;

	// Prints the Consolidated results
	gp1=fopen("7.0.Consolididated Results.dat","w");
	fprintf(gp1, "VARIABLES=\"X (m)\",\"Y (m)\",\"Size\"\n");
	fprintf(gp1, "ZONE T = \"Particle Position\"\n");

	for(int i = 0;i < NUM;i++)
   	{
		int pid = i;
		fprintf(gp1,"%lf %lf %lf\n",xp[pid],yp[pid],d[pid]);
   	}
	
	fprintf(gp1, "ZONE T = \"Void Position\"\n");
	for(int i = 0;i < NUM_OF_VOIDS; i++)
	{
		int vid = i;
		fprintf(gp1,"%e %e %e\n",VOID_POS[vid][0],VOID_POS[vid][1],VOID_POS[vid][2]);
	}
	
	fprintf(gp1, "ZONE T = \"Liquid Position\"\n");
	for(int i = 0;i < NUM_OF_VOIDS;i++)
	{
		int vid = i;
		for(int j = 0;j < VOID_SURR_NUM;j++)	
		{
			int pid = VOID_ENCIRCLED_BY[vid][j];
			fprintf(gp1,"%e %e %e\n",LIQ_POS_X[pid][vid], LIQ_POS_Y[pid][vid], LIQ_IN_VOID[vid]);
		}
	}
	fprintf(gp1, "ZONE T = \"Right Moving Liquid Position\"\n");
	for(int i = 0;i < NUM_OF_VOIDS;i++)
	{
		int vid = i;
		for(int j = 0;j < VOID_SURR_NUM;j++)	
		{
			int pid = VOID_ENCIRCLED_BY[vid][j];
			if(GOING_RIGHT[pid][vid] == 1)
			{
				fprintf(gp1,"%e %e %e\n",LIQ_POS_X[pid][vid], LIQ_POS_Y[pid][vid], LIQ_IN_VOID[vid]);
			}
		}
	}
	fprintf(gp1, "ZONE T = \"Left Moving Liquid Position\"\n");
	for(int i = 0;i < NUM_OF_VOIDS;i++)
	{
		int vid = i;
		for(int j = 0;j < VOID_SURR_NUM;j++)	
		{
			int pid = VOID_ENCIRCLED_BY[vid][j];
			if(GOING_RIGHT[pid][vid] == 0)
			{
				fprintf(gp1,"%e %e %e\n",LIQ_POS_X[pid][vid], LIQ_POS_Y[pid][vid], LIQ_IN_VOID[vid]);
			}
		}
	}
	fclose(gp1);
}

//=====================================================================================================================================================
void save_liquid_position(int update)
{
	FILE *gp;

	// Prints the Consolidated results
	gp=fopen("7.1.Liquid Position.dat","a");
	
	if (update == 1 )
	{
		fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Size\"\n");
	}
	fprintf(gp, "ZONE T = \"%d\" \n", update);
	for(int i = 0; i < NUM_OF_VOIDS; i++)
	{
		int vid = i;
		for(int j =0; j < VOID_SURR_NUM; j++)	
		{
			int pid = VOID_ENCIRCLED_BY[i][j];
			fprintf( gp, "%e %e %e\n",LIQ_POS_X[pid][vid], LIQ_POS_Y[pid][vid], LIQ_IN_VOID[vid] );
		}
	}
	fclose(gp);	
}
// =============================================================================================================================================================================
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
// =============================================================================================================================================================================
double Liquid_Volume_Determination(int k, int l, int vid_s, int vid_r, int pid_s)
{
	double liquid_vol = 0.0;
	if (vid_s != NO_NEIGHBOUR 															// Source void exists
	&&  VOID_POS[vid_r][1] <  VOID_POS[vid_s][1] 										// Source void is above the receiver void
	&&  LIQUID_ON_PARTICLE[pid_s][vid_s] > 0.0)											// Source void-particle has liquid
	{										
		if (LEFT_NEIGH_VOIDP[vid_r][k][l] == 1 												// The source void is to the left of the receiver void
		&& GOING_RIGHT[pid_s][vid_s] != NO_NEIGHBOUR )
		{
			liquid_vol = LIQUID_ON_PARTICLE[pid_s][vid_s] * GOING_RIGHT[pid_s][vid_s];												
		}
		
		if(RIGHT_NEIGH_VOIDP[vid_r][k][l]==1 												// The source void is to the right of the receiver void
		&& GOING_RIGHT[pid_s][vid_s] != NO_NEIGHBOUR )
		{
			liquid_vol = LIQUID_ON_PARTICLE[pid_s][vid_s] * (1 - GOING_RIGHT[pid_s][vid_s]);												
		}
	}
	return liquid_vol;
}
// =============================================================================================================================================================================
double Calculate_Detatchment_Angle(int void_id, int pid_s, int pid_r, double liquid_vol, double liquid_radius, double liquidvelocity)
{
	double valuetest 				=	0;
	double detatchment_angle 		=	0;
	double deltaANGLE				=	atan(liquid_radius*sin(CONTACT_ANGLE)/(0.5*d[pid_s]));		
	
	double liquid_mass				= liquid_vol * rho_l ;			
	double surface_tension_force	= 4 * PI * SURFACE_TENSION_LIQUID * pow(liquid_radius,2.0) * pow(sin(CONTACT_ANGLE),3.0) * cos(deltaANGLE/2);
	double centrifugal_force		= liquid_mass * pow(liquidvelocity,2.0) * liquid_radius * sin(CONTACT_ANGLE) / (liquid_radius * cos(CONTACT_ANGLE) + d[pid_s]/2.0);
	
	for (double thetaLtest=0;thetaLtest<=PI*2.0;thetaLtest=thetaLtest+PI/1000) // Changed PI/2.0 to PI*2.0
	{
		double gravitational_force		= liquid_mass * gravity * (d[pid_s] * sin(thetaLtest + deltaANGLE/2.0) + liquid_radius * cos(CONTACT_ANGLE) * cos(thetaLtest));
	
/* 		double void_x		= VOID_POS[void_id][0];
		double void_y		= VOID_POS[void_id][1];
		int cell_x 			= void_x/dx[1];
		int cell_y 			= void_y/dy[1];
		double u_g			= u[cell_x][cell_y];
		double v_g			= v[cell_x][cell_y];
		double u_mag		= sqrt( sqr(u_g) + sqr(v_g) );
		double gas_force	= 0.5 * Cd_gl * a_gl * rho_g * pow(u_mag,2);								 */
		
		valuetest = surface_tension_force - centrifugal_force - gravitational_force;
		if(valuetest <= 0) 
		{	
			detatchment_angle = thetaLtest;			
/* 			printf("%lf\n", thetaLtest*180/PI);	
			printf("%lf %lf %lf\n", surface_tension_force,centrifugal_force,gravitational_force);				 */
			goto RRR;
		}
	}
	RRR:
	if(detatchment_angle <= 0.0) 
	{
		detatchment_angle = 0.0;
	}
	
	if(detatchment_angle >= 90*PI/180.0) 
	{
		//detatchment_angle=90*PI/180.0;
	}
	
	if(liquidvelocity==0.0)
	{
		detatchment_angle=0.0;
	}
	
	if (fabs(VOID_POS[void_id][5] - 3.0) < 1.0e-6 
	&&  isneighbour(pid_r, pid_s)
	&&  yp[pid_r] < VOID_POS[void_id][1] 
	&&  yp[pid_s] < VOID_POS[void_id][1]
	&&  detatchment_angle > 0.0)
	{
		detatchment_angle=0.0;			
	}
	return detatchment_angle;
}
// =============================================================================================================================================================================
double Rivulet_Breakage_Calculation(int vid_r, int pid_r, int Side_flag, double NETA, double liquid_radius, double rivuletlength, double L1, double liquidvelocity, double Cd_sl, double Cd_gl)
{
	int BRK_DONE;
	double shift;
	double BETA_L_DASH;
	double L2;
	double ALPHA_L, BETA_L;
	double massR, massL;
	int cell_x, cell_y;
	double u_gl,v_gl;	
	
	double massEXTRAleft, massEXTRAright;
	double Vtmp1, Vtmp2;

	Vtmp1	=	(1 - CONTACT_ANGLE/PI) * PI * pow(liquid_radius,2.0) + sin(CONTACT_ANGLE) * cos(CONTACT_ANGLE) * pow(liquid_radius,2.0);
	Vtmp2	=	(1.0/6.0) * (PI*pow(liquid_radius,3.0)) * ( 2.0 + 3.0*cos(CONTACT_ANGLE) - pow(cos(CONTACT_ANGLE),3.0) );
	if( (0.5*rivuletlength) >= d[pid_r] * (0.50 * PI + (Side_flag*NETA) ) ) 													// Overhangs on the left
	{
		massEXTRAleft	=	(Vtmp1 * 0.5 * (rivuletlength - d[pid_r] * (0.5 * PI + (Side_flag*NETA) ) ) + Vtmp2) * rho_l;
	}
	else
	{
		massEXTRAleft	=	0;
	}

	if( (0.5*rivuletlength) >= d[pid_r] * (0.5 * PI - (Side_flag*NETA) ) ) 													// Overhangs on the right
	{
		massEXTRAright	=	(Vtmp1 * 0.5 * (rivuletlength - d[pid_r] * (0.5 * PI - (Side_flag*NETA) )) + Vtmp2) * rho_l;
	}
	else
	{
		massEXTRAright	=	0;
	}	

	double DELTA	=	Rivulet_Leading_Edge_Angle(pid_r, L1, NETA, liquid_radius);
	double ZETTA	=	Rivulet_Trailing_Edge_Angle(pid_r, L1, NETA, liquid_radius);
	
	double max_range	= (0.5*PI) + (Side_flag * ZETTA); 
	double d_zetta		= (-1 * Side_flag) * ZETTA/20.0; 
	for(BETA_L = (0.5*PI); BETA_L* Side_flag <= (max_range + d_zetta)* Side_flag; BETA_L = BETA_L - d_zetta)
	{	
		BRK_DONE		=	0;
		ALPHA_L		=	BETA_L - Side_flag * 0.5 * (ZETTA - DELTA + (0.5*PI));
		L2			=	2.0 * ( (0.5*d[pid_r]) + liquid_radius * cos(CONTACT_ANGLE)) * (max_range - BETA_L) *  Side_flag;
		
		massR		=	( Vtmp1 * L2 + Vtmp2) * rho_l;
		massL		=	( Vtmp1 * (L1 - L2) + Vtmp2) * rho_l;
		double AslR	=	2.0 * liquid_radius * L2 * sin(CONTACT_ANGLE);
		double AslL	=	2.0 * liquid_radius * (L1 - L2) * sin(CONTACT_ANGLE);

		double AglR	=	PI*liquid_radius*liquid_radius + 2.0 * liquid_radius * L2;
		double AglL	=	PI*liquid_radius*liquid_radius + 2.0 * liquid_radius * (L1 - L2);

		cell_x = (xp[pid_r] - d[pid_r])/dx[1];
		cell_y = (yp[pid_r])/dy[1];	
		u_gl	= u[cell_x][cell_y]	- ul[cell_x][cell_y];
		v_gl	= v[cell_x][cell_y]	- vl[cell_x][cell_y];
		double F_gas_drag_left	= 0.5*Cd_gl*rho_g*AglL * (u_gl*fabs(u_gl) + v_gl*fabs(v_gl));
		
		
		cell_x = (xp[pid_r] + d[pid_r])/dx[1];
		cell_y = (yp[pid_r])/dy[1];
		u_gl	= u[cell_x][cell_y]	- ul[cell_x][cell_y];
		v_gl	= v[cell_x][cell_y]	- vl[cell_x][cell_y];
		double F_gas_drag_right	= 0.5*Cd_gl*rho_g*AglR * (u_gl*fabs(u_gl) + v_gl*fabs(v_gl));
		
		double F_drag_left	= 0.5*Cd_sl*rho_l*AslL*pow(liquidvelocity,2.0);
		double F_drag_right	= 0.5*Cd_sl*rho_l*AslR*pow(liquidvelocity,2.0);	

		double F_Wt_left	= massL * gravity * sin(BETA_L - 0.5*PI);
		double F_Wt_right	= massR * gravity * cos(ALPHA_L);

		double F_ExWt_left	= massEXTRAleft * gravity;
		double F_ExWt_right	= massEXTRAright * gravity;

		double Force_left	=  F_Wt_left + F_drag_left  + F_ExWt_left + F_gas_drag_left;
		double Force_right	=  F_Wt_right- F_drag_right + F_ExWt_right + F_gas_drag_right;
		double Force_ST		= liquid_radius * SURFACE_TENSION_LIQUID * ( 2.0*PI*(1 - CONTACT_ANGLE/PI) + 2.0*sin(CONTACT_ANGLE) );	
		
		if( (Force_left + Force_right) >= Force_ST) 
		{	
			BRK_DONE	=	1;
			BETA_L_DASH	=	BETA_L;
			goto break_t;
		}
	}		
	break_t:	
	L2		=	2.0 * ( (0.5*d[pid_r]) + liquid_radius * cos(CONTACT_ANGLE)) * (max_range - BETA_L_DASH) *  Side_flag;
	if(BRK_DONE == 1 && VOID_POS[vid_r][0] >= xp[pid_r] - (d[pid_r]/4)) 
	{	
		
		shift	=	L2/L1;	
		print_rivulet_breakage(vid_r, L1, L2);
	}
	else
	{
		shift	=	0;
	}
	
	if ( BRK_DONE == 1 && VOID_POS[vid_r][0] <= xp[pid_r] + (d[pid_r]/4) ) 
	{	
		
		shift	=	(L1 - L2)/L1;
		print_rivulet_breakage(vid_r, L1, L1 - L2);
	}
	else
	{
		shift	=	1;																
	}		
	
	return shift;
}
// =============================================================================================================================================================================
double Rivulet_Centroid_Angle(int pid_r, double XLPOS)
{
	double NETA	=	asin( abs(XLPOS - xp[pid_r]) / (0.5 * d[pid_r]) );	
	return NETA;
}
// =============================================================================================================================================================================
double Rivulet_Leading_Edge_Angle(int pid_r, double L1, double NETA, double liquid_radius)
{
	// double DELTA	=	(0.5 * PI) - NETA - (0.5 * L1)/(0.5 * d[pid_r] + liquid_radius * cos(CONTACT_ANGLE));					// ORIGINAL
	double DELTA	=	(0.5 * PI) - NETA - 2 * asin( (0.25 * L1)/(0.5 * d[pid_r] + liquid_radius * cos(CONTACT_ANGLE)) );		// PROPOSED CORRECTION
	if(DELTA < 0) 
	{
		DELTA	=	0.00;
	}
	
	return DELTA;
}
// =============================================================================================================================================================================
double Rivulet_Trailing_Edge_Angle(int pid_r, double L1, double NETA, double liquid_radius)
{
	//double ZETTA	=	(0.5 * L1) / (0.5 * d[pid_r] + liquid_radius * cos(CONTACT_ANGLE)) - NETA;								// ORIGINAL
	double ZETTA	=	2 * asin( (0.25 * L1) / (0.5 * d[pid_r] + liquid_radius * cos(CONTACT_ANGLE)) ) - NETA;					// PROPOSED CORRECTION
	if(ZETTA >= (0.5 * PI) )
	{
		ZETTA	=	0.5 * PI;
	}

	if(ZETTA < 0) 
	{
		ZETTA	=	0.00;
	}
	return ZETTA;
}