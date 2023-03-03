double velocity_magnitude_gas(int, int);
double velocity_magnitude_solid(int, int);
double velocity_magnitude_liquid(int, int);
double velocity_magnitude_fines(int, int);
void Forcefield_Computation(void);
double Fg_p(int, int);
double Fg_l(int, int);
double Fg_f(int, int);
double Fp_l(int, int);
double Fp_f(int, int);
double Fl_f(int, int);
// ============================================================================================
double differential_velocity_magnitude_gas_solid(int i, int j)
{
	double u_e, u_w, v_n, v_s;
	double u_P, v_P;
	double mag_v;
	if ( (us[i-1][j] >0.00 ) || (us[i][j] >0.00 ) || (vs[i][j-1] >0.00 ) || (vs[i][j] >0.00 ) )
	{
		printf("SOLID Velocity detected %d %d %lf %lf\n", i, j, us[i][j], vs[i][j]);
	}

	u_e	= u[i-1][j]	-	us[i-1][j];
	u_w	= u[i][j]	-	us[i][j];
	v_n	= v[i][j]	-	vs[i][j];
	v_s	= v[i][j-1]	-	vs[i][j-1];
	
	u_P	= 0.5 * (u_w + u_e);
	v_P	= 0.5 * (v_n + v_s);
	mag_v	= sqrt( sqr(u_P) + sqr(v_P) );
	
	return(mag_v);	
}

double differential_velocity_magnitude_gas_liquid(int i, int j)
{
	double u_e, u_w, v_n, v_s;
	double u_P, v_P;
	double mag_v;

	u_e	= u[i-1][j]	-	ul[i-1][j];
	u_w	= u[i][j]	-	ul[i][j];
	v_n	= v[i][j]	-	vl[i][j];
	v_s	= v[i][j-1]	-	vl[i][j-1];
	
	u_P	= 0.5 * (u_w + u_e);
	v_P	= 0.5 * (v_n + v_s);
	mag_v	= sqrt( sqr(u_P) + sqr(v_P) );
	
	return(mag_v);	
}

double differential_velocity_magnitude_gas_fines(int i, int j)
{
	double u_e, u_w, v_n, v_s;
	double u_P, v_P;
	double mag_v;

	u_e	= u[i-1][j]	-	uf[i-1][j];
	u_w	= u[i][j]	-	uf[i][j];
	v_n	= v[i][j]	-	vf[i][j];
	v_s	= v[i][j-1]	-	vf[i][j-1];
	
	u_P	= 0.5 * (u_w + u_e);
	v_P	= 0.5 * (v_n + v_s);
	mag_v	= sqrt( sqr(u_P) + sqr(v_P) );
	
	return(mag_v);	
}

double differential_velocity_magnitude_solid_liquid(int i, int j)
{
	double u_e, u_w, v_n, v_s;
	double u_P, v_P;
	double mag_v;

	u_e	= us[i-1][j]	-	ul[i-1][j];
	u_w	= us[i][j]	-	ul[i][j];
	v_n	= vs[i][j]	-	vl[i][j];
	v_s	= vs[i][j-1]	-	vl[i][j-1];
	
	u_P	= 0.5 * (u_w + u_e);
	v_P	= 0.5 * (v_n + v_s);
	mag_v	= sqrt( sqr(u_P) + sqr(v_P) );
	
	return(mag_v);	
}

double differential_velocity_magnitude_solid_fines(int i, int j)
{
	double u_e, u_w, v_n, v_s;
	double u_P, v_P;
	double mag_v;

	u_e	= us[i-1][j]	-	uf[i-1][j];
	u_w	= us[i][j]	-	uf[i][j];
	v_n	= vs[i][j]	-	vf[i][j];
	v_s	= vs[i][j-1]	-	vf[i][j-1];
	
	u_P	= 0.5 * (u_w + u_e);
	v_P	= 0.5 * (v_n + v_s);
	mag_v	= sqrt( sqr(u_P) + sqr(v_P) );
	
	return(mag_v);	
}

double differential_velocity_magnitude_liquid_fines(int i, int j)
{
	double u_e, u_w, v_n, v_s;
	double u_P, v_P;
	double mag_v;

	u_e	= ul[i-1][j]	-	uf[i-1][j];
	u_w	= ul[i][j]	-	uf[i][j];
	v_n	= vl[i][j]	-	vf[i][j];
	v_s	= vl[i][j-1]	-	vf[i][j-1];
	
	u_P	= 0.5 * (u_w + u_e);
	v_P	= 0.5 * (v_n + v_s);
	mag_v	= sqrt( sqr(u_P) + sqr(v_P) );
	
	return(mag_v);	
}
//===============================================================================
void Forcefield_Computation()
{
	int i, j;
	double mag_vg;
	double mag_vs;
	double mag_vf;
	double mag_vl;
	
	max_Fgp = 0.00;
	max_Fgl = 0.00;
	max_Fgf = 0.00;
	max_Fpl = 0.00;
	max_Fpf = 0.00;
	max_Flf = 0.00;
	
	for (i = 0; i < IMAX; i++)
	{
		for (j = 0; j < JMAX; j++)
		{
			Forcefield_GP[i][j]	=	0.00;
			Forcefield_GL[i][j]	=	0.00;
			Forcefield_GF[i][j]	=	0.00;
			Forcefield_PL[i][j]	=	0.00;
			Forcefield_PF[i][j]	=	0.00;
			Forcefield_LF[i][j]	=	0.00;
		}
	}
	
	for (i = 0; i < IMAX; i++)
	{
		for (j = 0; j < JMAX; j++)
		{
			double mag_Vgp	= differential_velocity_magnitude_gas_solid(i,j);
			double mag_Vgl	= differential_velocity_magnitude_gas_liquid(i,j);
			double mag_Vgf	= differential_velocity_magnitude_gas_fines(i,j);
			double mag_Vpl	= differential_velocity_magnitude_solid_liquid(i,j);
			double mag_Vpf	= differential_velocity_magnitude_solid_fines(i,j);
			double mag_Vlf	= differential_velocity_magnitude_liquid_fines(i,j);
			
				
			Forcefield_GP[i][j]	=	Fg_p(i,j);
			Forcefield_GL[i][j]	=	Fg_l(i,j);
			Forcefield_GF[i][j]	=	Fg_f(i,j);
			Forcefield_PL[i][j]	=	Fp_l(i,j);
			Forcefield_PF[i][j]	=	Fp_f(i,j);
			Forcefield_LF[i][j]	=	Fl_f(i,j);
			
			
			max_Fgp = max(max_Fgp, Forcefield_GP[i][j] * mag_Vgp);
			max_Fgl = max(max_Fgl, Forcefield_GL[i][j] * mag_Vgl);
			max_Fgf = max(max_Fgf, Forcefield_GF[i][j] * mag_Vgf);
			max_Fpl = max(max_Fpl, Forcefield_PL[i][j]);
			max_Fpf = max(max_Fpf, Forcefield_PF[i][j] * mag_Vpf);
			max_Flf = max(max_Flf, Forcefield_LF[i][j]);
		}
	}	
}
//-------------------------------------------------------------------------------
//----------------Force between particle and gas---------------------------------
//-------------------------------------------------------------------------------
double Fg_p(int i, int j)
{
	double Rep ,Cd;
	double mag_v;
	double aA, bB;
	// CALCULATING INTERSTITIAL VELOCITY FOR GAS & SOLID
	mag_v	= differential_velocity_magnitude_gas_solid(i,j);
	// CALCULATING REYNOLDS NUMBER AND Cd
	Rep	= rho_g*mag_v*dp/mu_g; 
	if ( (Rep <= 1.0) && (Rep > 0) )
		Cd	=	(24.0/Rep);
	else if ( (Rep <= 1000) && (Rep > 1.0) )
		Cd	=	(24.0/Rep)*(1.0 + 0.15*pow(Rep,0.687));
	else
		Cd	=	0.44;
	// CALCULATING COEFFICIENTS
	if(void_frac(i,j) < 0.8)
	{
		aA=150.0 * (mu_g/sqr(dp)) * ( sqr(void_frac_solid(i,j))/cube(void_frac(i,j)) );
		bB=1.75 * (rho_g/dp) * ( (void_frac_solid(i,j))/cube(void_frac(i,j)) );
	}
	else
	{
		aA = 0.00; 
		bB = 0.75*(Cd*rho_g/dp)*pow(void_frac(i,j),-4.65)*void_frac_solid(i,j);
	}
	// OTHER EXCEPTION HANDLING
	if (vfrac[i][j]==0.0)
	{
		return (0);
	}
	else
	{	
		return (aA + bB*mag_v);
	}
}
//-------------------------------------------------------------------------------
//----------------Force between liquid and gas---------------------------------
//-------------------------------------------------------------------------------
double Fg_l(int i, int j)
{
	double mag_v	=	differential_velocity_magnitude_gas_liquid(i,j);
	double Cd_gl	=	4.40;
	double Force = 0.5 * Cd_gl * rho_l * Area_GL[i][j] * sqr(mag_v) ;
	return (Force);
}
//-------------------------------------------------------------------------------
//----------------Force between fine and gas-------------------------------------
//-------------------------------------------------------------------------------
double Fg_f(int i, int j)
{
	double Ref, Cd;
	double mag_v;
	double aA, bB;
	if(loading > 0.0)
	{
		mag_v	=	differential_velocity_magnitude_gas_fines(i,j);
	}
	else
	{
		mag_v	=	0.0;
	}
		
	Ref=void_frac(i,j)*rho_g*mag_v*dpf/mu_g; 

	if ((Ref<=1.0) && (Ref>0))
		Cd=(24.0/Ref);
	else if ((Ref<=1000) && (Ref>1.0))
		Cd=(24.0/Ref)*(1.0 + 0.15*pow(Ref,0.687));
	else if(Ref==0.0) 
		Cd=0.0;
	else
		Cd=0.44;

	bB	=	0.75 * (Cd*rho_g/dpf) * pow(void_frac(i,j),-4.65)*void_frac_fines_dynamic(i,j);
		
	if (loading==0.0)
	{
		return (0);
	}
	else
	{
		return (bB * mag_v);
	}
}
//-------------------------------------------------------------------------------
//----------------Force between particle and gas---------------------------------
//-------------------------------------------------------------------------------
double Fp_l(int i, int j)
{
	double Cd_sl;
	double mag_v	=	differential_velocity_magnitude_solid_liquid(i,j);
	double Rel = rho_l * mag_v * d[0]/mu_l;
	if(mag_v > 0.00)
	{
		Cd_sl	=	541.0/( mag_v + 33.0);
	}
	else
	{
		Cd_sl	=	33.0;
	}
	double Force = 0.5 * Cd_sl * rho_l * Area_SL[i][j] * sqr(mag_v) ;
	return (Force);
}
//-------------------------------------------------------------------------------
//----------------Force between fine and particle--------------------------------
//-------------------------------------------------------------------------------
double Fp_f(int i, int j)
{
	double D,fk,Fr,s;
	double mag_v;
	double aA, bB;

	mag_v	=	differential_velocity_magnitude_solid_fines(i,j);
	
	// CALCULATING EQUIVALENT DIAMETER
	D	=	(2.0/3.0)*dp*(1 - void_frac_solid(i,j)) / (void_frac_solid(i,j));
	if(void_frac_solid(i,j) == 0) D = 0.0000095; 

	// CALCULATING FROUDE NUMBER AND fk (Additional Pressure loss coefficient)
	Fr	=	(mag_v)/(sqrt(D*gravity)); 
	if (Fr == 0)
		fk	=	0;
	else
		fk	=	14.98/pow(Fr,1.33);

	// CALCULATING COEFFICIENTS
	bB	=	(fk/(2.0*D))*rho_f*epfd[i][j];

	// OTHER EXCEPTION HANDLING	
	if ((loading==0.0)||(epfs[i][j]==0.24))
	{
		return (0);
	}
	else
	{
		return (bB*mag_v);
	}
}
//-------------------------------------------------------------------------------
//----------------Force between particle and gas---------------------------------
//-------------------------------------------------------------------------------
double Fl_f(int i, int j)
{
	return (0);
}