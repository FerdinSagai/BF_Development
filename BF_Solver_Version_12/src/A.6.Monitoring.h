//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void save_convergence_data()
{
	FILE *fp;
	fp = fopen("6.0.Convergence.dat","a");
	if(Convergence_output_count==1){

		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Residue\",\"Mass of fines IN\",\"Mass of fines OUT\",\"Mass balance of fines\"\n");
	}
		fprintf(fp,"%d %e %e %e %e\n",Convergence_output_count,Residue,mass_in_fine,mass_out_fine,mass_fine);
	fclose(fp);
	Convergence_output_count++;
	//========================================================================================================
	fp = fopen("6.1a.Performance_Scheme.dat","a");
	if(Performance_scheme_count==1)
	{
		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Velocity X\",\"Gas Velocity Y\",\"Gas Pressure\",\"Gas Turbulent KE\",\"Gas Turbulent DE\",\"Fines Velocity X\",\"Fines Velocity Y\",,\"Dynamic Holdup\"\n");
	}
		fprintf(fp,"%d %d %d %d %d %d %d %d %d\n",Performance_scheme_count,Conv_U,Conv_V,Conv_P,Conv_KE,Conv_de,Conv_Uf,Conv_Vf,Conv_epfd);
	fclose(fp);
	Performance_scheme_count++;
	//========================================================================================================
	fp = fopen("6.2b.Performance_Implementation.dat","a");
	if(Performance_method_count == 1)
	{
		fprintf(fp, "VARIABLES=\"Iteration\",\"Gas Velocity X\",\"Gas Velocity Y\",\"Gas Pressure\",\"Gas Turbulent KE\",\"Gas Turbulent DE\",\"Fines Velocity X\",\"Fines Velocity Y\",\"Dynamic Holdup\"\n");
	}
		fprintf(fp,"%d %e %e %e %e %e %e %e %e\n",Performance_method_count,cpu_time_U_gas,cpu_time_V_gas,cpu_time_Pressure,cpu_time_KE_gas,cpu_time_de_gas,cpu_time_U_fines,cpu_time_V_fines,cpu_time_dynamic_holdup);
	fclose(fp);
	Performance_method_count++;
	//========================================================================================================
	Conv_U = 0;
	Conv_V = 0;
	Conv_P = 0;
	Conv_KE = 0;
	Conv_de = 0;
	Conv_Uf = 0;
	Conv_Vf = 0;
	Conv_epfd =0;
}
