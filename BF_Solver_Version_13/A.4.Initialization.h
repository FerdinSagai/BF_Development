//---------------------------------------------------
void initialise()
{
	printf("Overall Initialization ...\n");
	/*------Local variables declaration--------------------------------------------------------*/
	tid						= omp_get_thread_num();
	nthreads				= omp_get_max_threads();

	status_gas_phase		=	false;	
	status_solid_phase		=	false;	
	status_liquid_phase		=	false;	
	status_fines_phase		=	false;	

	init_gas_phase			=	false;	
	init_solid_phase		=	false;
	init_liquid_phase		=	false;
	init_fines_phase		=	false;
	
	steady_state_gas		=	false;
	steady_state_solid		=	false;
	steady_state_liquid		=	false;
	steady_state_fines		=	false;
//-----------------------Initializaing iteration variables ----------------------------   	
	Residue					=	0.00;
	mass_fine				=	0.00;;
	epfs_frac				=	0.60;		// Should be input
	loading					=	0.00;
	epfs_total				=	0.00;
//-----------------------Tuyere--------------------------------------------------------
	rp						=	0.5*dp;
	x1_t					=	0.0; 
	x2_t					=	r;
	x3_t					=	0.0;
	x4_t					=	r;
	y1_t					=	h-0.5*D_T;
	y2_t					=	h-0.5*D_T;
	y3_t					=	h+0.5*D_T;
	y4_t					=	h+0.5*D_T;
//-----------------------Stress based raceway computations------------------------------	
	double H_B				= 	Ly -(h + 0.5*D_T);
	double rho_eff			= 	Domain_Voidage*rho_g + (1.0 - Domain_Voidage)*rho_s;
	double Fr				= 	vin*vin/(9.81*dp*phi_s);
	iso_stress				= 	4085 - 24.67*vin;
//-------Convergence and Status Monitor initializations--------------------------------------		
	max_Ug					= 	0.00;
	max_Vg					= 	0.00;
	max_KE					= 	0.00;
	max_de					= 	0.00;
	max_Uf					= 	0.00;
	max_Vf					= 	0.00;
	max_P					= 	0.00;
	min_Stress				= 	1000000.0;	
	max_epfs 				= 	0.00;
	max_epfd 				= 	0.00;
	Conv_U 					= 	0;
	Conv_V 					= 	0;
	Conv_P 					= 	0;
	Conv_KE 				= 	0;
	Conv_de 				= 	0;
	Conv_Uf 				= 	0;
	Conv_Vf 				= 	0;
	Conv_epfd 				= 	0;
//-------Output counter initializations--------------------------------------		
	Particle_position_count	=	1;
	Convergence_output_count=	1;
	Performance_scheme_count=	1; 
	Performance_method_count=	1;
	printf("Overall Initialization ...Complete\n");
}
//+=+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++