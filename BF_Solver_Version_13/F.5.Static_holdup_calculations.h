 void calc_epfs_packed_bed(void);
 void calc_epfs_cohesive_bed(void);
//-----------Static holdup calculations-------------------------------------------			
void Static_holdup_correlations()
{
	double Area_SH;
	time_hsb	=	clock();
	calc_epfs_packed_bed();
	//calc_epfs_cohesive_bed();
	enforce_domain_physics();
	Area_SH = 0.0;
	for (int i = 0; i <= IMAX; i++)
	{
		for (int j = 0;j<=JMAX;j++)
		{	
			epfs_max	=	epfs_frac*(1.0-void_frac_solid(i,j));
			epfs_limit	=	0.95*epfs_max;
			if(epfs[i][j] >= epfs_limit)
			{				
				Area_SH = Area_SH + (dx[i]*dy[j]);
			}
		}
	}
	time_hsa	=	clock();
	cpu_time_static_holdup = ((double) (time_hsa - time_hsb)) / (nthreads*CLOCKS_PER_SEC);
	printf("\tStatic Holdup.......................:: CPU time = %lf seconds ||Area = %lf m2\n",cpu_time_static_holdup,Area_SH);	
}
/*------Choice of correlations for Static Holdup Computation-----------------*/
void calc_epfs_packed_bed() 										// Packed Bed
{
	int i, j;
	double ug_e, ug_w, vg_n, vg_s;
	double ug_P, vg_P;
	double uf_e, uf_w, vf_n, vf_s;
	double uf_P, vf_P;
	double mag_vg, mag_vf;
	double Dstar, Rep, Frf;
	double r1, r2, r3, r4, r5;
	double a1, a2, a3, a4, a5;

	FILE *sh;
	for (i=1;i<=IMAX-1;i++)
	{
		for (j=1;j<=JMAX-1;j++)
		{
			epfs_max		=	epfs_frac*(1.0-void_frac_solid(i,j));
			epfs_limit		=	0.95*epfs_max;
			if(epfs[i][j]>=epfs_limit)
			{
				epfs[i][j]=epfs_max;
			}
			else
			{
				// DEFINING SUPERFICIAL VELOCITY FOR GAS
				ug_e = u[i-1][j];
				ug_w = u[i][j];
				vg_n = v[i][j];
				vg_s = v[i][j-1];
				// DETERMINING THE VELOCITY AT CELL CENTER
				ug_P = 0.5*(ug_w + ug_e);
				vg_P = 0.5*(vg_n + vg_s);
				mag_vg = sqrt( sqr(ug_P) + sqr(vg_P) );

				// DEFINING INTERSTITIAL VELOCITY FOR FINES
				uf_e = uf[i-1][j];
				uf_w = uf[i][j];
				vf_n = vf[i][j];
				vf_s = vf[i][j-1];
				// DETERMINING THE VELOCITY AT CELL CENTER
				uf_P = 0.5*(uf_w + uf_e);
				vf_P = 0.5*(vf_n + vf_s);
				mag_vf = sqrt( sqr(uf_P) + sqr(vf_P) );

				// DETERMINING THE NON_DIMENSIONAL SIMILARITY PARAMETERS
				//double Dstar	= 2*void_frac(i,j)*dp/(3*(void_frac_solid(i,j)));
				Dstar	= void_frac(i,j)*dp / void_frac_solid(i,j);
				Rep	= rho_g*mag_vg*dp / (mu_g*void_frac_solid(i,j)); 
				Frf	= mag_vg / sqrt(gravity*Dstar);

				// Determination of ratios and their exponents
				r1	= dpf / dp;
				a1	= 0.568;

				r2	= Rep;
				a2	= -0.227;
				
				r3	= Frf * Frf;
				a3	= -1.165;

				r4	= (mag_vf*rho_f*epfd[i][j]) / (mag_vg*rho_g);
				//r4	= Gf / (mag_vg*rho_g);
				a4	= 0.172;

				r5	= rho_f / rho_g;	
				a5	= -0.270;

				double c5	= 63.028;
				double A1	= pow(r1,a1)*pow(r2,a2)*pow(r3,a3)*pow(r4,a4)*pow(r5,a5);

				// Determining Static holdup
				if( (void_frac_solid(i,j)<0.36) || (loading == 0.00) )
				{
					epfs[i][j]=0.0;
				}
				else
				{
					epfs[i][j]=(min(c5*A1,epfs_max));
				}
			}
		}
	}	
	// Assigning values for Left and Right Faces
	for (j=1;j<=JMAX-1;j++)
	{
		epfs[0][j]	= epfs[1][j];
		epfs[IMAX][j]	= epfs[IMAX-1][j];
	}
	// Assigning values for Top and Bottom Faces
	for (i=1;i<=IMAX-1;i++)
	{
		epfs[i][0]	= epfs[i][1];
		epfs[i][JMAX]	= epfs[i][JMAX-1];
	}
	// Assigning values for Corners
	epfs[0][0]		= epfs[1][1];
	epfs[IMAX][JMAX]	= epfs[IMAX-1][JMAX-1];
	epfs[0][JMAX]		= epfs[1][JMAX-1];
	epfs[IMAX][0]		= epfs[IMAX-1][1];
}
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
void calc_epfs_cohesive_bed() 									// Cohesive Bed
{
	int i, j;
	double ug_e, ug_w, vg_n, vg_s;
	double ug_P, vg_P;
	double uf_e, uf_w, vf_n, vf_s;
	double uf_P, vf_P;
	double mag_vg, mag_vf;
	double Dstar, Rep, Frf;
	double r1, r2, r3, r4, r5;
	double a1, a2, a3, a4, a5;

	FILE *sh;
	for (i=1;i<=IMAX-1;i++)
	{
		for (j=1;j<=JMAX-1;j++)
		{
			epfs_max		=	epfs_frac*(1.0-void_frac_solid(i,j));
			epfs_limit		=	0.95*epfs_max;
			if(epfs[i][j]>=epfs_limit)
			{
				epfs[i][j]=epfs_max;
			}
			else
			{
				// DEFINING SUPERFICIAL VELOCITY FOR GAS
				ug_e = u[i-1][j];
				ug_w = u[i][j];
				vg_n = v[i][j];
				vg_s = v[i][j-1];
				// DETERMINING THE VELOCITY AT CELL CENTER
				ug_P = 0.5*(ug_w + ug_e);
				vg_P = 0.5*(vg_n + vg_s);
				mag_vg = sqrt( sqr(ug_P) + sqr(vg_P) );

				// DEFINING INTERSTITIAL VELOCITY FOR FINES
				uf_e = uf[i-1][j];
				uf_w = uf[i][j];
				vf_n = vf[i][j];
				vf_s = vf[i][j-1];
				// DETERMINING THE VELOCITY AT CELL CENTER
				uf_P = 0.5*(uf_w + uf_e);
				vf_P = 0.5*(vf_n + vf_s);
				mag_vf = sqrt( sqr(uf_P) + sqr(vf_P) );

				// DETERMINING THE NON_DIMENSIONAL SIMILARITY PARAMETERS
				//Dstar	= 2*void_frac(i,j)*dp/(3*(void_frac_solid(i,j)));
				Dstar	= void_frac(i,j)*dp / void_frac_solid(i,j);
				Rep	= rho_g*mag_vg*dp / (mu_g*void_frac_solid(i,j)); 
				Frf	= mag_vg / sqrt(gravity*Dstar);

				// Determination of ratios and their exponents
				r1	= dpf / dp;
				a1	= 1.562;

				//r2	= void_frac(i,j)*Rep;
				r2	= Rep;
				a2	= -0.62;
				
				//r3	= void_frac(i,j)*Frf * void_frac(i,j)*Frf;
				r3	= Frf * Frf;
				a3	= -1.416;

				r4	= Gf / (mag_vg*rho_g);
				a4	= 0.578;

				r5	= rho_f / rho_g;	
				a5	= 1.056;

				double c5	= 0.0012;
				double A1	= pow(r1,a1)*pow(r2,a2)*pow(r3,a3)*pow(r4,a4)*pow(r5,a5);

				// Determining Static holdup
				if( (void_frac_solid(i,j)<0.36) || (loading == 0.00) )
				{
					epfs[i][j]=0.0;
				}
				else
				{
					epfs[i][j]=(min(c5*A1,epfs_max));
				}
			}
		}
	}	
	// Assigning values for Left and Right Faces
	for (j=1;j<=JMAX-1;j++)
	{
		epfs[0][j]		= epfs[1][j];
		epfs[IMAX][j]	= epfs[IMAX-1][j];
	}
	// Assigning values for Top and Bottom Faces
	for (i=1;i<=IMAX-1;i++)
	{
		epfs[i][0]		= epfs[i][1];
		epfs[i][JMAX]	= epfs[i][JMAX-1];
	}
	// Assigning values for Corners
	epfs[0][0]			= epfs[1][1];
	epfs[IMAX][JMAX]	= epfs[IMAX-1][JMAX-1];
	epfs[0][JMAX]		= epfs[1][JMAX-1];
	epfs[IMAX][0]		= epfs[IMAX-1][1];
}