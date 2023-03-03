void Fines_Convergence_Check()
{
	mass_out_fine		= 0.0;
	mass_in_fine		= Gf*(D_T*Lz);
	Growth_Fines		= 0.00;
	for (int i = 1; i <= IMAX-1; i++)
	{
		mass_out_fine = mass_out_fine + (epfd[i][JMAX-1]*vf[i][JMAX-1]*(dx[i]*Lz)*rho_f);
	}
	mass_fine		= fabs(mass_in_fine - mass_out_fine);
	Growth_Fines	= mass_fine/mass_in_fine;
}