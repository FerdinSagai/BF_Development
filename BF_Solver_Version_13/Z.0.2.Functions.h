#define PI acos(-1.0)
int sign(double x)
{ 
	int s;
	
	if (x == 0)
	{
		return (0); 	
	}
	else
	{
		s = x/fabs(x);
		return (s); 
	}
} 

double sqr(double x)
{ 
  return (x*x); 
} 
 
int sqri(int x) 
{ 
  return (x*x); 
} 

int mini(int x, int y) 
{ 
  if (x<y) return x; 
  else return y; 
 
}

int maxi(int x, int y) 
{ 
  if (x>y) return x; 
  else return y; 
 
} 

double min(double x, double y) 
{ 
  if (x<y) return x; 
  else return y; 
 
} 
 
double cube(double x) 
{ 
  return (x*x*x); 
} 
 
double max(double x, double y) 
{ 
  if (x>y) return x; 
  else return y; 
 
} 
 
double fiver(double x)
{
  return (cube(x)*sqr(x));
}
 
double interp1D(double L, double FL, double R, double FR)
{
  double gradient=(FR-FL)/(L+R);
  double F=FR-gradient*R;
  return (F);
}

//---------------------------Raceway radius determination-------------------------------
double decrease_corr(double W, double H, double v_b, double D_T, double dp, double phi_s, double rho_s, double eps_g)
{
  double k=4.2; //3.3612;
  double phi_w=(10.76/180.0)*PI;
  double mu_w=tan(phi_w);
  double d_eff=dp*phi_s;
  double rho_eff=eps_g*rho_g + (1.0-eps_g)*rho_s;

  double C1=k*pow(rho_g*sqr(v_b)*D_T/(rho_eff*gravity*d_eff*W),0.6);
  double C2=pow(D_T/H,-0.12);
  double C3=pow(mu_w,-0.24);

  double D_r=D_T*C1*C2*C3;

  return (D_r);
}
//---------------------------Raceway radius determination-------------------------------
double increase_corr(double W, double H, double v_b, double D_T, double dp, double phi_s, double rho_s, double eps_g)
{
  double k=164;
  double phi_w=(10.76/180.0)*PI;
  double mu_w=tan(phi_w);
  double d_eff=dp*phi_s;
  double rho_eff=eps_g*rho_g + (1.0-eps_g)*rho_s;

  double C1=pow(rho_g*sqr(v_b*D_T)/(rho_eff*gravity*d_eff*H*W),0.8);
  double C2=pow(mu_w,-0.25);

  double D_r=D_T*k*C1*C2;

  return (D_r);
}
/*---------------------------------------------------------------------*/
/*------------------Calculation of Void Fraction-----------------------*/
/*---------------------------------------------------------------------*/
// Void fraction due to PACKED BED
 double void_frac_solid(int i, int j)
{  
	double chi;
   
   
	chi=1 - vfrac[i][j];   
	if(chi == 0) return (0.3); 
	else if (chi == 1) return (0.85); 
	else return (chi); 
}

// Void fraction due to TOTAL FINES
double void_frac_fines(int i, int j)
{
	return (epfs[i][j] + epfd[i][j]);
}

// Void fraction
double void_frac(int i, int j)
{
	if(1.0-void_frac_solid(i,j)-void_frac_fines(i,j) <= 0) return (0.1);  
	else if (1.0-void_frac_solid(i,j)-void_frac_fines(i,j) >= 1) return (0.85); 
	else return (1.0-void_frac_solid(i,j)-void_frac_fines(i,j));
}

// Void fraction due to DYNAMIC FINES
double void_frac_fines_dynamic(int i, int j)
{
	return (epfd[i][j]);  
}
// Void fraction due to DYNAMIC FINES
double void_frac_liquid(int i, int j)
{
	return (epfd[i][j]);  
}
/*---------------------------------------------------------------------*/
/*--------------------Calculation of Viscosity-------------------------*/
/*---------------------------------------------------------------------*/
double mut(int i, int j)
{
   return (c_mu*rho_g*(sqr(ke[i][j])/de[i][j]));
} 

double mue(int i, int j)
{
   return (mu_g + (c_mu*rho_g*(sqr(ke[i][j])/de[i][j])));
} 
//---------------------------------------------------------------------
void enforce_domain_physics()
{
	int i, j;
	
	// Setting Void Fractions in inlet and outlet regions
	/*
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
	*/
		
	for(i = 0; i <= IMAX; i++)
	{
		for(j = 0;j <= JMAX; j++)
		{
			if(State[i][j] == -1)				// WALL OR IMPERMEABLE SOLID
			{	
				vfrac[i][j]	= 0.00;			
				u[i][j]		= 0.00;
				v[i][j]		= 0.00;
				f[i][j]		= 0.00;			
				g[i][j]		= 0.00;

				uf[i][j]	= 0.00;
				vf[i][j]	= 0.00;
				ff[i][j]	= 0.00;
				gf[i][j]	= 0.00;

				epfs[i][j]	= 0.0;
				epfd[i][j]	= 0.0;
				epfd1[i][j]	= 0.0;
				epfd2[i][j]	= 0.0;
			}
			else if(State[i][j] == 1)			// INLET REGIONS
			{  	
				vfrac[i][j]	= Raceway_Voidage;				
				u[i][j]		= vin;
				v[i][j]		= 0.00;
				f[i][j]		= u[i][j];
				g[i][j]		= v[i][j];

				uf[i][j]	= vin;
				vf[i][j]	= 0.00;
				ff[i][j]	= uf[i][j];
				gf[i][j]	= vf[i][j];

				epfs[i][j]	= 0.0;
				epfd[i][j]	= loading;
				epfd1[i][j]	= loading;
				epfd2[i][j]	= loading;
				
				ke[i][j]	= 0.009*sqr(vin);
				de[i][j]	= pow(0.009*sqr(vin),1.5)/lm;
				
				ke1[i][j]	= 0.009*sqr(vin);
				de1[i][j]	= pow(0.009*sqr(vin),1.5)/lm;
			}
			else	if(State[i][j] == 2)			// OUTLET REGIONS
			{
				vfrac[i][j]	= Raceway_Voidage;
				p[i][j]		= 0.00;
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