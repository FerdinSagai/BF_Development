//-------------------------------------------------------------------------------
//--------------------------Next run input files---------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//--------------------------Writing Report files---------------------------------
//-------------------------------------------------------------------------------
void report_data_GFS()
{
	printf("Reporting Writing Initiated ...\n");
	//-----------------------Initial Calculations-------------------------------------------
	FILE  *fp;
	fp=fopen("2.0.Reporting.rep","w");	
		fprintf(fp,"\t\t\t\tREPORTS\n");	
		fprintf(fp,"-------------------------------------------------------------------------\n");
		fprintf(fp,"DOMAIN\n");	
		fprintf(fp,"\tDomain Dimensions(Lx X Ly X Lz) ..............................%lf m X %lf X %lf m\n", Lx, Ly, Lz);
		fprintf(fp,"\tDomain Discretization(Nx X Ny X Nz) ..........................%d X %d X %d\n", IMAX, JMAX, KMAX);
		fprintf(fp,"\tDomain Discretization(dx X dy X dz) ..........................%lf X %lf X %lf\n", Lx/IMAX, Ly/JMAX,Lz/KMAX);
		fprintf(fp,"\tX and Y Wall velocity at the top boundary ....................%lf %lf\n",u_top,v_top);
		fprintf(fp,"\tX and Y Wall velocity at the bottom boundary .................%lf %lf\n",u_bottom,v_bottom);
		fprintf(fp,"\tX and Y Wall velocity at the left boundary ...................%lf %lf\n",u_left,v_left);
		fprintf(fp,"\tX and Y Wall velocity at the right boundary ..................%lf %lf\n",u_right,v_right);
		
		fprintf(fp,"\n");
		
		fprintf(fp,"GAS PROPERTIES\n");
		fprintf(fp,"\tDensity of the gas ...........................................%lf kg/m^3\n",rho_g);
		fprintf(fp,"\tViscosity of the gas .........................................%lf Pa s\n",mu_g);
		fprintf(fp,"\tFlowrate of the gas inflow via Tuyere ........................%lf lpm\n",gas_flowrate);
		fprintf(fp,"\tInflow Velocity of gas .......................................%lf m/s\n",vin);
		fprintf(fp,"\tSuperficial Velocity of gas ..................................%lf m/s\n",v_super);
		
		fprintf(fp,"\n");
		
		fprintf(fp,"SOLID DETAILS\n");
		fprintf(fp,"\tSize of the particle .........................................%lf m\n",dp);
		fprintf(fp,"\tSphericity of the particle ...................................%lf\n",phi_s);
		fprintf(fp,"\tDensity of the packing .......................................%lf kg/m^3 \n",rho_s);
		fprintf(fp,"\tCoefficient of friction ......................................%lf\n",coeff_friction);
		fprintf(fp,"\tCoefficient of restitution ...................................%lf\n",coeff_restitution);
		fprintf(fp,"\tNumber of particles ..........................................%d\n",NUM);
		fprintf(fp,"\tVolume of particle ...........................................%e m^3\n",PI*dp*dp);
		fprintf(fp,"\tVolume of Bed ................................................%e m^3\n",(NUM)*PI*dp*dp);
		fprintf(fp,"\tVolume of Domain .............................................%lf m^3\n",(Lx*Ly));
		fprintf(fp,"\tSolid Fraction ...............................................%lf\n",((NUM)*PI*dp*dp)/(Lx*Ly));
		fprintf(fp,"\tPorosity .....................................................%lf\n",1.00-((NUM)*PI*dp*dp)/(Lx*Ly));
/* 		fprintf(fp,"\tMass of Bed ..................................................%e kg\n",mass);
		fprintf(fp,"\tMoment of Inertia of Bed .....................................%e kgm^2\n",MMOI); */
		
		
		fprintf(fp,"\n");
		
		fprintf(fp,"LIQUID DETAILS\n");
		fprintf(fp,"\tFlowrate of the liquid inflow via Rotameter ..................%lf\n",ROTAMETER_FLOWRATE);
		fprintf(fp,"\tDensity of the liquid ........................................%lf kg/m^3\n",rho_l);
		fprintf(fp,"\tViscosity of the liquid ......................................%lf Pa s\n",mu_l);
		fprintf(fp,"\tContact angle between liquid and solid .......................%lf degree s\n",CONTACT_ANGLE);
		fprintf(fp,"\tSurface tension of the liquid ................................%lf N/m\n",SURFACE_TENSION_LIQUID);
		
		fprintf(fp,"\n");
		
		fprintf(fp,"FINES PROPERTIES\n");
		fprintf(fp,"\tLoading of the fines .........................................%lf \n",(Gf/(vin*rho_f)));
		fprintf(fp,"\tFlow rate of fines ...........................................%lf kg/m^2\n",Gf);
		fprintf(fp,"\tDiameter of the fine particle ................................%lf m\n",dpf);
		fprintf(fp,"\tSphericity of the fines ......................................%lf kg/m^3\n",phi_f);
		fprintf(fp,"\tDensity of the fines .........................................%lf kg/m^3\n",rho_f);
		fprintf(fp,"\tViscosity of the fines .......................................%lf Pa s\n",mu_f);

		fprintf(fp,"\n");

		
		fprintf(fp,"TURBULENCE MODEL CONSTANTS\n");
		fprintf(fp,"\tConstant C1 ..................................................%lf\n",c_1);
		fprintf(fp,"\tConstant C2 ..................................................%lf\n",c_2);
		fprintf(fp,"\tConstant c_mu ................................................%lf\n",c_mu);
		fprintf(fp,"\tConstant Sigma k .............................................%lf\n",sigma_k);
		fprintf(fp,"\tConstant Sigma e .............................................%lf\n",sigma_e);
		fprintf(fp,"\tMixing length assocaiated with turbulence in the raceway .....%lf m\n",lm);

		fprintf(fp,"\n");
		
		fprintf(fp,"TUYERE\n");	
		fprintf(fp,"Number of Tuyeres ..............................................%d \n",n_tuyeres);
		for (int i = 0; i < n_tuyeres; i++)
		{
			fprintf(fp,"\tTuyere Number ........................................%d \n",i+1);
			fprintf(fp,"\tProtrusion of the Tuyere into the bed ................%lf m\n",Tuyere_protrusion[i+1]);
			fprintf(fp,"\tHeight of the Tuyere from bottom of the vessel .......%lf m\n",Tuyere_height[i+1]);
			fprintf(fp,"\tSize of the Tuyere opening ...........................%lf m\n",Tuyere_opening[i+1]);			
		}

		fprintf(fp,"\n");		

		fprintf(fp,"ROTAMETERS DETAILS\n");
		fprintf(fp,"Number of Rotameters ...........................................%d \n",NUM_OF_ROTAMETERS);
		for (int i = 0; i < NUM_OF_ROTAMETERS; i++)
		{
			fprintf(fp,"\tRotameter Number .....................................%d \n",i+1);			
			fprintf(fp,"\tDistance of the Rotameter into the bed ...............%lf m\n",SOURCE_OF_LIQUID[i][0]);
			fprintf(fp,"\tHeight of the Rotameter from bottom of the vessel ....%lf m\n",SOURCE_OF_LIQUID[i][1]);
			fprintf(fp,"\tSize of the Rotameter opening ........................%lf m\n",ROTAMETER_OPENING);			
		}	
		
		fprintf(fp,"\tTuyere co-ordinates 1 =%lf %lf\n",x1_t, y1_t);
		fprintf(fp,"\tTuyere co-ordinates 2 =%lf %lf\n",x2_t, y2_t);
		fprintf(fp,"\tTuyere co-ordinates 3 =%lf %lf\n",x3_t, y3_t);
		fprintf(fp,"\tTuyere co-ordinates 4 =%lf %lf\n",x4_t, y4_t);
		fprintf(fp,"\tTuyere Cells :: ii01= %d ii02=%d jj01=%d jj02=%d\n",ii01, ii02, jj01, jj02);	
		fprintf(fp,"\n");

		fprintf(fp,"SIMILARITY NUMBERS\n");
		fprintf(fp,"\tThe Bed height is : %lf m\n",Bed_height);
		fprintf(fp,"\tThe flow Reynolds number is : %lf\n",Re);
		fprintf(fp,"\tThe CFL number in the X direction is : %lf\n",CFL_x);
		fprintf(fp,"\tThe CFL number in the Y direction is : %lf\n",CFL_y);
		fprintf(fp,"\n");

		fprintf(fp,"RACEWAY PARAMETERS\n");	
		fprintf(fp,"\tThe Raceway Isostress boundary is : %lf \n",iso_stress);
		fprintf(fp,"\tThe Raceway Center is : %lf %lf \n",Raceway_Center_X, Raceway_Center_Y);
		fprintf(fp,"\tThe Raceway Radius is : %lf m\n",raceway_radius);
		fprintf(fp,"\tThe Min Max Raceway X Range is : %lf %lf \n",Raceway_Center_X - raceway_radius, Raceway_Center_X + raceway_radius);
		fprintf(fp,"\tThe Min Max Raceway Y Range is : %lf %lf \n",Raceway_Center_Y - raceway_radius, Raceway_Center_Y + raceway_radius);
		fprintf(fp,"\n");		
		
		fprintf(fp,"COHESIVE BLOCKS\n");
		fprintf(fp,"\n");
		

		
		
		fprintf(fp,"CONSTANTS\n");
		fprintf(fp,"\tGravity = %lf m/s^2\n",gravity);
		fprintf(fp,"\tMark = %lf\n",mark);

		fprintf(fp,"\n");
		
		
		
		fprintf(fp,"ISO STRESS RACEWAY MODELLING\n");
		fprintf(fp,"\tThe Raceway Isostress boundary is : %lf \n",iso_stress);
		fprintf(fp,"\tThe Raceway Center is : %lf %lf \n",Raceway_Center_X, Raceway_Center_Y);
		fprintf(fp,"The Raceway Radius is : %lf \n",raceway_radius);
		fprintf(fp,"\tThe Min Max Raceway X Range is : %lf %lf \n",Raceway_Center_X - raceway_radius, Raceway_Center_X + raceway_radius);
		fprintf(fp,"\tThe Min Max Raceway Y Range is : %lf %lf \n",Raceway_Center_Y - raceway_radius, Raceway_Center_Y + raceway_radius);
		fprintf(fp,"\n");
		fprintf(fp,"\n");
		fprintf(fp,"DEM MODEL CONSTANTS\n");
		fprintf(fp,"\tSpring constant (particle-normal) = %lf\n",Kn);
		fprintf(fp,"\tSpring constant (particle-tangential) = %lf\n",Kt);
		//fprintf(fp,"\tDamping constant (particle-normal) = %lf\n",Cn);
		//fprintf(fp,"\tDamping constant (particle-tangential) = %lf\n",Ct);	
		fprintf(fp,"\tSpring constant (wall-normal) = %lf\n",KnWall);
		fprintf(fp,"\tSpring constant (wall-tangential) = %lf\n",KtWall);
		fprintf(fp,"\tDamping constant (wall-normal) = %lf\n",CnWall);
		fprintf(fp,"\tDamping constant (wall-tangential) = %lf\n",CtWall);
		fprintf(fp,"\n");
		
		fprintf(fp,"MOVING BEDS\n");
		fprintf(fp,"\tTime step = %lf s\n",STATICtime);
		fprintf(fp,"\tReal time for DEM movement = %lf  s\n",sim_time);
		fprintf(fp,"\tInsertion rate = %d\n",insertion_rate);
		fprintf(fp,"\tDischarge rate = %d\n",discharge_rate);
		fprintf(fp,"\tBed Bottom for Removal = %lf m\n",L_down);
		fprintf(fp,"\tLow Bottom Wall = %lf m\n",Y_bottom_wall);
		fprintf(fp,"\tLimiting Height for insertion= %lf m\n",limitH);
		fprintf(fp,"\n");
		
		fprintf(fp,"CONVERGENCE AND RELAXATION\n");
		fprintf(fp,"\tSpatial convergence limit = %lf\n",tol);
		fprintf(fp,"\tFine mass balance limit = %lf\n",mass_bal);
		fprintf(fp,"\tRelaxation parameter in computation of U = %lf\n",relaxu);
		fprintf(fp,"\tRelaxation parameter in computation of V = %lf\n",relaxv);
		fprintf(fp,"\tRelaxation parameter in computation of P = %lf\n",relaxp);
		fprintf(fp,"\tRelaxation parameter in computation of KE = %lf\n",relaxke);
		fprintf(fp,"\tRelaxation parameter in computation of DE = %lf\n",relaxde);
		fprintf(fp,"\tRelaxation parameter in computation of dynamic holdup = %lf\n",relaxepfd);
		fprintf(fp,"\n");
	fclose(fp);
	printf("Reporting Writing Initiated .................. Completed\n");
	printf("\n");
	printf("\n");
}
//-------------------------------------------------------------------------------
//--------------------------Prerun Checks-------------------------------------
//-------------------------------------------------------------------------------
void prerun_checks()
{
	printf("Pre-run Checks Initiated ...\n");
	//-------------------------------------------------------------------------------
	int i, j;
	FILE *gp;			
	
	gp=fopen("2.1.Prerun Checks.dat","w");
	fprintf(gp,"TITLE = States of the domain\n");
	fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"State\",\"Void Fraction\"\n");
	fprintf(gp, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n",Transient_output_count, IMAX, JMAX);
	double yy = 0.5*dy[1];
	for (int j = 1; j <= JMAX; j++)
	{
		double xx  = 0.5*dx[1];
		for (int i = 1; i <= IMAX; i++)
		{
			fprintf(gp,"%e %e %d %e\n",xx, yy,State[i][j],vfrac[i][j]);
			// UPDATE X POSITION
			xx	=	xx + 0.5*(dx[i]+dx[i+1]); 
		}
		// UPDATE Y POSITION
		yy	=	yy + 0.5*(dy[j]+dy[j+1]);
	}
	fclose(gp);
	//-------------------------------------------------------------------------------
	gp=fopen("Raceway_Map.inp","w");
	for (int j = JMAX; j >= 1; j--)
	{
		for (int i = 1; i <= IMAX; i++)
		{
			fprintf(gp,"%.2f ",vfrac[i][j]);
		}
		fprintf(gp,"\n");
	}
	fclose(gp);
	//-------------------------------------------------------------------------------
	printf("Pre-run Checks Initiated .................. Completed\n");
	printf("\n");
	printf("\n");
}