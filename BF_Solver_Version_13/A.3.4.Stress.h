//!************************* Solver *************************************
void stress_calc()
{
	int i,j,loop,iter,st,st1,fn,fn1;
	double temp1,temp2,temp3,temp4,y1,y2,x1,x2,sai1,sai2,sigma1,sigma2;
	double xw,yw,sigmaw,saiw,xwb,ywb,sigmawb,saiwb,sigma1b,sigma2b,sai1b,sai2b;
	


int l1,l2;
int n_tu,m_tu;
double xl,yl;


double saiwL, saiwR;

double eps, delta, rhos, gy, deltaw, sur_charge;

double sigma_xx[ni][nj],sigma_yy[ni][nj],tau_xy[ni][nj];
double sigma_alpha[ni][nj],tau_alpha[ni][nj];
double sigma_beta[ni][nj],tau_beta[ni][nj];
	xl= Lx;
	// m_tu= 102 ;//        !121
	yl=Ly;
	// n_tu=202  ;//        !241 
	//!***************** For Stress ********************

	l1 = 35;// !m 
	l2 =m+n+l1;

	double dx_MOC=xl/double(m);

	for(i=0;i<=m;i++)
	{
		for(j=m;j>=0;j--)
		{			
		xa[i][j] = xl-i*dx_MOC; // !0.02   !0.1  !0.02
		ya[i][j] = 0.0 ;
		// printf("\ni=%d\tj=%d\txa=%e\tya=%e",i,j,xa[i][j],ya[i][j]);
		}
	}

	rhos =(1 - Domain_Voidage) * rho_s;//  !1500
	gy = 9.8 ;
	sur_charge=0.0;
	delta = 28.213*(PI/180);// !24.55*(3.14/180);//  !30.65
	deltaw = 7.36*(PI/180) ;//!9.05*(3.14/180);// !10.76
	eps = (0.5*PI - delta)/2;

	for(i=0;i<=m;i++)
	{
		for(j=m;j>=0;j--)
		{
		sai[i][j] = 0.5*PI ;
		sigma[i][j] = (rhos*gy*ya[i][j])/(1+sin(delta)) ;
		// printf("\ni=%d\tj=%d\tsai=%e\tsigma=%e",i,j,sai[i][j],sigma[i][j]);
		}
	}
	saiwL = 0.5*(PI -(asin(sin(deltaw)/sin(delta)) - deltaw));// !78.42*(3.14/180)  !2.736/2 
	saiwR = 0.5*(PI +(asin(sin(deltaw)/sin(delta)) - deltaw));//  !101.58*(3.14/180)  !1.7729
	
	//!*************** Stress calculation *********************** 
	// This code represent the region "OBA"
	st=1; fn=m;
	for(iter=1;iter<=m;iter++) 
	{ 
		st1=m;
		for(i=st;i<=fn;i++) 
		{
			j = st1;
			//printf("\niter=%d\ti=%d\tj=%d",iter,i,j);
			x1 = xa[i][j-1];
			x2 = xa[i-1][j];
			y1 = ya[i][j-1];
			y2 = ya[i-1][j];

			sigma1 = sigma[i][j-1];
			sigma2 = sigma[i-1][j];
			sai1 = sai[i][j-1];
			sai2 = sai[i-1][j];
			sai1b = sai1;
			sai2b = sai2;
			sigma1b = sigma1;
			sigma2b = sigma2;
			//printf("\nx1=%e\tx2=%e\ty1=%e\ty2=%e\tsai1=%e\tsai2=%e\tsigma1=%e\tsigma2=%e",x1,x2,y1,y2,sai1,sai2,sigma1,sigma2);

			xa[i][j] = ((y2 - y1) - (x2*tan(sai2b + eps)) + (x1*tan(sai1b - eps)))/(tan(sai1b - eps) - tan(sai2b + eps));
			ya[i][j] = y1 + (xa[i][j] - x1)*tan(sai1b - eps);
			if(sur_charge==0.0)
			{
				sai[i][j] = 1.57;
				sigma[i][j] = (rhos*gy*ya[i][j])/(1+sin(delta));
			}
			else
			{
				temp1 = (sigma2-sigma1) ;
				temp2 = 2*tan(delta)*((sigma1b*sai1) + (sigma2b*sai2)) ;
				temp3 = rhos*gy*(y2 - y1 - (2*tan(delta)*xa[i][j]) + tan(delta)*(x2+x1));
				temp4 = (2*tan(delta)*(sigma1b+sigma2b));
				sai[i][j] = (temp1 + temp2 - temp3)/temp4 ;

				sigma[i][j] = sigma1 + 2*tan(delta)*sigma1b*(sai[i][j] - sai1)+ rhos*gy*(ya[i][j] - y1 - tan(delta)*(xa[i][j] - x1));
			}

			//  sai1b=(sai1+sai[i][j])/2;
			// sai2b=(sai2+sai[i][j])/2;
			// sigma1b = (sigma1+sigma[i][j])/2;
			// sigma2b = (sigma2+sigma[i][j])/2;

			//printf("\nxa=%e\tya=%e\tsigma=%e\tsai=%e",xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
			st1=st1-1;

			sigma_xx[i][j]=sigma[i][j]*(1.0+(sin(delta)*cos(2.0*sai[i][j])));
			sigma_yy[i][j]=sigma[i][j]*(1.0-(sin(delta)*cos(2.0*sai[i][j])));
			tau_xy[i][j]=sigma[i][j]*sin(delta)*sin(2.0*sai[i][j]);

			sigma_alpha[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_alpha[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
			sigma_beta[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_beta[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
		}//I loop ends
		st=st+1;
	}// Iter loop ends


	// fclose(fp); 
	//!*********alpha characteristics *****************
	// This code represent the region "OBD"
	st=m+1; fn=2*m;
	st1=1; fn1=m;
	for(j=st;j<=fn;j++) 
		{
		for(i=st1;i<=fn1;i++) 
		{
			//printf("\nRegion OBD\n \ni=%d\tj=%d",i,j);
			if(i==st1)
			{
				xa[i][j]=xl;
				x1=xa[i][j-1];
				y1=ya[i][j-1];

				sigma1=sigma[i][j-1];
				sai[i][j]= saiwL;
				sai1=(sai[i][j]+sai[i][j-1])/2;
				ya[i][j]=y1+(xa[i][j]-x1)*tan(sai1-eps);
				sigma[i][j]=sigma1+2*tan(delta)*sigma1*(sai[i][j]-sai1)+ rhos*gy*(ya[i][j] - y1 - tan(delta)*fabs(xa[i][j] - x1));
				//printf("\nix1=%d\tjx1=%d\tx1=%e\ty1=%e\txa=%e\tya=%e\tsigma=%e\tsai=%e",i-1,j,xa[i-1][j],ya[i-1][j],xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
			}
			else
			{
				x1 = xa[i][j-1];
				x2 = xa[i-1][j];
				y1 = ya[i][j-1];
				y2 = ya[i-1][j];

				sigma1 = sigma[i][j-1];
				sigma2 = sigma[i-1][j];
				sai1 = sai[i][j-1];
				sai2 = sai[i-1][j];
				sai1b = sai1;
				sai2b = sai2;
				sigma1b = sigma1;
				sigma2b = sigma2;
				//printf("\nx1=%e\tx2=%e\ty1=%e\ty2=%e\tsai1=%e\tsai2=%e\tsigma1=%e\tsigma2=%e",x1,x2,y1,y2,sai1,sai2,sigma1,sigma2);
				
				for(loop = 1;loop<=10;loop++) 
				{
					xa[i][j] = ((y2 - y1) - (x2*tan(sai2b + eps)) + (x1*tan(sai1b - eps)))/(tan(sai1b - eps) - tan(sai2b + eps));
					ya[i][j] = y1 + (xa[i][j] - x1)*tan(sai1b - eps);

					temp1 = (sigma2-sigma1) ;
					temp2 = 2*tan(delta)*((sigma1b*sai1) + (sigma2b*sai2)) ;
					temp3 = rhos*gy*((y2 - y1) - (2*tan(delta)*xa[i][j]) + tan(delta)*(x2+x1));
					temp4 = (2*tan(delta)*(sigma1b+sigma2b));
					sai[i][j] = (temp1 + temp2 - temp3)/temp4 ;

					sigma[i][j] = sigma1 + 2*tan(delta)*sigma1b*(sai[i][j] - sai1)+ rhos*gy*((ya[i][j] - y1) - tan(delta)*(xa[i][j] - x1));
					//printf("\nt1=%e\tt2=%e\ttt2=%e\tt3=%e\ttt3=%e",sigma1,2*tan(delta)*sigma1b,(sai[i][j] - sai1),(ya[i][j] - y1),tan(delta)*(xa[i][j] - x1));
					//printf("\nxa=%e\tya=%e\tsigma=%e\tsai=%e",xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
					sai1b=(sai1+sai[i][j])/2;
					sai2b=(sai2+sai[i][j])/2;
					sigma1b = (sigma1+sigma[i][j])/2;
					sigma2b = (sigma2+sigma[i][j])/2;
				}
			}
			sigma_xx[i][j]=sigma[i][j]*(1.0+(sin(delta)*cos(2.0*sai[i][j])));
			sigma_yy[i][j]=sigma[i][j]*(1.0-(sin(delta)*cos(2.0*sai[i][j])));
			tau_xy[i][j]=sigma[i][j]*sin(delta)*sin(2.0*sai[i][j]);

			sigma_alpha[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_alpha[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
			sigma_beta[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_beta[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
		}// J loop ends
		st1=st1+1; 
	}
	
	// This code represent the region "ABE"
	st=m+1; fn=2*m;
	st1=1; fn1=m;
	for(i=st;i<=fn;i++) 
	{
		for(j=st1;j<=fn1;j++) 
		{
			//printf("\n\nRegion ABE\ni=%d\tj=%d",i,j);
			if(j==st1)
			{
				xa[i][j]=0.0;
				x2=xa[i-1][j];
				y2=ya[i-1][j];

				sigma2=sigma[i-1][j];
				sai[i][j]=saiwR;
				sai2=(sai[i][j]+sai[i-1][j])/2;
				ya[i][j]=y2+(xa[i][j]-x2)*tan(sai2+eps);
				//printf("\nix1=%e\ty1=%e\txa-x1=%e\tya=%e\ttan=%e",x1,y1,(xa[i][j]-x1),ya[i][j] ,tan(sai1+eps));

				sigma[i][j]=sigma2-2*tan(delta)*sigma2*(sai[i][j]-sai2)+ rhos*gy*(ya[i][j] - y2 + tan(delta)*(xa[i][j] - x2));
				//printf("\nix1=%d\tjx1=%d\tx1=%e\ty1=%e\txa=%e\tya=%e\tsigma=%e\tsai=%e",i-1,j,xa[i-1][j],ya[i-1][j],xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
			}
			else
			{
				x1 = xa[i][j-1];
				x2 = xa[i-1][j];
				y1 = ya[i][j-1];
				y2 = ya[i-1][j];

				sigma1 = sigma[i][j-1];
				sigma2 = sigma[i-1][j];
				sai1 = sai[i][j-1];
				sai2 = sai[i-1][j];
				sai1b = sai1;
				sai2b = sai2;
				sigma1b = sigma1;
				sigma2b = sigma2;
				//printf("\nx1=%e\tx2=%e\ty1=%e\ty2=%e\tsai1=%e\tsai2=%e\tsigma1=%e\tsigma2=%e",x1,x2,y1,y2,sai1,sai2,sigma1,sigma2);
				for(loop = 1;loop<=10;loop++) 
				{  	
					xa[i][j] = ((y2 - y1) - (x2*tan(sai2b + eps)) + (x1*tan(sai1b - eps)))/(tan(sai1b - eps) - tan(sai2b + eps));
					ya[i][j] = y1 + (xa[i][j] - x1)*tan(sai1b - eps);

					temp1 = (sigma2-sigma1) ;
					temp2 = 2*tan(delta)*((sigma1b*sai1) + (sigma2b*sai2)) ;
					temp3 = rhos*gy*((y2 - y1) - (2*tan(delta)*xa[i][j]) + tan(delta)*(x2+x1));
					temp4 = (2*tan(delta)*(sigma1b+sigma2b));
					sai[i][j] = (temp1 + temp2 - temp3)/temp4 ;

					sigma[i][j] = sigma1 + 2*tan(delta)*sigma1b*(sai[i][j] - sai1)+ rhos*gy*(ya[i][j] - y1 - tan(delta)*(xa[i][j] - x1));
					//printf("\nxa=%e\tya=%e\tsigma=%e\tsai=%e",xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
					sai1b=(sai1+sai[i][j])/2;
					sai2b=(sai2+sai[i][j])/2;
					sigma1b = (sigma1+sigma[i][j])/2;
					sigma2b = (sigma2+sigma[i][j])/2;
				}
			}
			sigma_xx[i][j]=sigma[i][j]*(1.0+(sin(delta)*cos(2.0*sai[i][j])));
			sigma_yy[i][j]=sigma[i][j]*(1.0-(sin(delta)*cos(2.0*sai[i][j])));
			tau_xy[i][j]=sigma[i][j]*sin(delta)*sin(2.0*sai[i][j]);

			sigma_alpha[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_alpha[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
			sigma_beta[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_beta[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
		}// I loop ends
		st1=st1+1; 
	}

	// This code represent the region below "DBE"
	st=m+1; fn=6*m; st1=m+1; fn1=2*m+1;
	for(i=st;i<=fn;i++) 
	{
		for(j=st1;j<=fn1;j++) 
		{
			//printf("\n\nRegion below DBE\ni=%d\tj=%d",i,j);
			if(j==fn1)
			{
				xa[i][j]=xl;
				x1=xa[i][j-1];
				y1=ya[i][j-1];
				sigma1=sigma[i][j-1];
				sai[i][j]= saiwL;
				sai1=(sai[i][j]+sai[i][j-1])/2;
				ya[i][j]=y1+(xa[i][j]-x1)*tan(sai1-eps);
				sigma[i][j]=sigma1+2*tan(delta)*sigma1*(sai[i][j]-sai1)+ rhos*gy*(ya[i][j] - y1 - tan(delta)*fabs(xa[i][j] - x1));
				//printf("\n'Region j==fn1'\nix1=%d\tjx1=%d\tx1=%e\ty1=%e\txa=%e\tya=%e\tsigma=%e\tsai=%e",i-1,j,xa[i-1][j],ya[i-1][j],xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
			}
			else if(i>=(2*m+1)&&j==st1)
			{
				//printf("\n\nRegion below DBE\ni=%d\tj=%d",i,j);
				xa[i][j]=0.0;
				x2=xa[i-1][j];
				y2=ya[i-1][j];

				sigma2=sigma[i-1][j];
				sai[i][j]=saiwR;
				sai2=(sai[i][j]+sai[i-1][j])/2;
				ya[i][j]=y2+(xa[i][j]-x2)*tan(sai2+eps);
				//printf("\n'Region j==st1'\nix1=%e\ty1=%e\txa-x1=%e\tya=%e\ttan=%e",x1,y1,(xa[i][j]-x1),ya[i][j] ,tan(sai1+eps));

				sigma[i][j]=sigma2-2*tan(delta)*sigma2*(sai[i][j]-sai2)+ rhos*gy*(ya[i][j] - y2 + tan(delta)*(xa[i][j] - x2));
				//printf("\nix1=%d\tjx1=%d\tx1=%e\ty1=%e\txa=%e\tya=%e\tsigma=%e\tsai=%e",i-1,j,xa[i-1][j],ya[i-1][j],xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
			}
			else
			{
				x1 = xa[i][j-1];
				x2 = xa[i-1][j];
				y1 = ya[i][j-1];
				y2 = ya[i-1][j];

				sigma1 = sigma[i][j-1];
				sigma2 = sigma[i-1][j];
				sai1 = sai[i][j-1];
				sai2 = sai[i-1][j];
				sai1b = sai1;
				sai2b = sai2;
				sigma1b = sigma1;
				sigma2b = sigma2;
				//printf("\nx1=%e\tx2=%e\ty1=%e\ty2=%e\tsai1=%e\tsai2=%e\tsigma1=%e\tsigma2=%e",x1,x2,y1,y2,sai1,sai2,sigma1,sigma2);
				for(loop = 1;loop<=10;loop++) 
				{  	
					xa[i][j] = ((y2 - y1) - (x2*tan(sai2b + eps)) + (x1*tan(sai1b - eps)))/(tan(sai1b -eps) - tan(sai2b + eps));
					ya[i][j] = y2 + ( xa[i][j]-x2)*tan(sai2b +eps);

					temp1 = (sigma2-sigma1) ;
					temp2 = 2*tan(delta)*((sigma1b*sai1) + (sigma2b*sai2)) ;
					temp3 = rhos*gy*((y2 - y1) - (2*tan(delta)*xa[i][j]) + tan(delta)*(x2+x1));
					temp4 = (2*tan(delta)*(sigma1b+sigma2b));
					sai[i][j] = (temp1 + temp2 - temp3)/temp4 ;

					sigma[i][j] = sigma1 + 2*tan(delta)*sigma1b*(sai[i][j] - sai1)+ rhos*gy*(ya[i][j] - y1 - tan(delta)*(xa[i][j] - x1));
					//printf("\nxa=%e\tya=%e\tsigma=%e\tsai=%e",xa[i][j] ,ya[i][j] ,sigma[i][j],sai[i][j]);
					sai1b=(sai1+sai[i][j])/2;
					sai2b=(sai2+sai[i][j])/2;
					sigma1b = (sigma1+sigma[i][j])/2;
					sigma2b = (sigma2+sigma[i][j])/2;
					} 
			} 
			sigma_xx[i][j]=sigma[i][j]*(1.0+(sin(delta)*cos(2.0*sai[i][j])));
			sigma_yy[i][j]=sigma[i][j]*(1.0-(sin(delta)*cos(2.0*sai[i][j])));
			tau_xy[i][j]=sigma[i][j]*sin(delta)*sin(2.0*sai[i][j]);

			sigma_alpha[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_alpha[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
			sigma_beta[i][j]=0.5*(sigma_xx[i][j]+sigma_yy[i][j])+0.5*(sigma_xx[i][j]-sigma_yy[i][j])*cos(2.0*sai[i][j])-tau_xy[i][j]*sin(2.0*sai[i][j]);
			tau_beta[i][j]=0.5*(sigma_xx[i][j]-sigma_yy[i][j])*sin(2.0*sai[i][j])-tau_xy[i][j]*cos(2.0*sai[i][j]);
		}// J loop ends
		if(i<(2*m+1))
		{
			fn1=fn1+1; st1=st1;
		} 
		else
		{
			st1=st1+1; fn1=fn1;
		}
	}
}
//!************************** Interpolation ***********************
 //!*********************interpolation printout******************************************  
void interp_printout()
{
	int i,j,st,fn,st1,fn1,zn,iter;
	FILE *fp,*fp1;

	fp=fopen("Stress_Sai_MOC.m","w");

		fprintf(fp,"x=[");
		st=0; fn=m; st1=m; fn1=0; zn=0;
		for(iter = 0;iter<=3*(m+n);iter++)
		{
			i=st;
			{
				for(j = st1;j>=fn1;j--)
				{
					fprintf(fp,"\t%e",xa[i][j]);
					// printf("\niter=%d i=%d j=%d xa=%e",iter,i,j,xa[i][j]);
					i=i+1;
				}
			} // J loop ends

			if(iter==zn)//odd number
			{
				st=st+1; fn=fn; st1=st1; fn1=fn1+1;
			}
			else
			{
				st=st; fn=fn+1; st1=st1+1; fn1=fn1;
				zn=zn+2;
			}
			//zn=zn+2;
		}// Iter loop ends
		
		fprintf(fp,"];");
		fprintf(fp,"\ny=[");
		st=0; fn=m; st1=m; fn1=0; zn=0;
		
		for(iter = 0;iter<=3*(m+n);iter++)
		{
			i=st;
			{
				for(j = st1;j>=fn1;j--)
				{
					fprintf(fp,"\t%e",ya[i][j]);
					i=i+1;
				}
			} // J loop ends

			if(iter==zn)//odd number
			{
				st=st+1; fn=fn; st1=st1; fn1=fn1+1;
			}
			else
			{
				st=st; fn=fn+1; st1=st1+1; fn1=fn1;
				zn=zn+2;
			}
			//zn=zn+2;
		}// Iter loop ends
		
		fprintf(fp,"];");
		fprintf(fp,"\nsigma=[");
		st=0; fn=m; st1=m; fn1=0; zn=0;
		
		for(iter = 0;iter<=3*(m+n);iter++)
		{
			i=st;
			{
				for(j = st1;j>=fn1;j--)
				{
					fprintf(fp,"\t%e",sigma[i][j]);
					i=i+1;}
				} // J loop ends
				
				if(iter==zn)//odd number
				{
					st=st+1; fn=fn; st1=st1; fn1=fn1+1;
				}
				else
				{
					st=st; fn=fn+1; st1=st1+1; fn1=fn1;
					zn=zn+2;
				}
				//zn=zn+2;
		}// Iter loop ends
		
		fprintf(fp,"];");
		fprintf(fp,"\nsai=[");
		st=0; fn=m; st1=m; fn1=0; zn=0;
		
		for(iter = 0;iter<=3*(m+n);iter++)
		{
			i=st;
			{
				for(j = st1;j>=fn1;j--)
				{
					fprintf(fp,"\t%e",sai[i][j]);
					i=i+1;
				}
			} // J loop ends

			if(iter==zn)//odd number
			{
				st=st+1; fn=fn; st1=st1; fn1=fn1+1;
			}
			else
			{
				st=st; fn=fn+1; st1=st1+1; fn1=fn1;
			zn=zn+2;
			}
			//zn=zn+2;
		}// Iter loop ends
		fprintf(fp,"];");
		fprintf(fp,"\ndx=[");
		for(i = 1;i<IMAX;i++) 
		{
			fprintf(fp,"\t%e",x[i]);
		}

		fprintf(fp,"];");
		fprintf(fp,"\ndy=[");
		
		for(j = 1;j<JMAX;j++) 
		{
			fprintf(fp,"\t%e",y[j]);
		}
		fprintf(fp,"];");
		fprintf(fp,"\n[x1,y1]=meshgrid(dx,dy);");
		fprintf(fp,"\n[x1,y1,sigma_w]=griddata(x,%lf-y,sigma,x1,y1);",Ly);
		fprintf(fp,"\n[x1,y1,sai_w]=griddata(x,%lf-y,sai,x1,y1);",Ly);

	
		fprintf(fp,"\n cs = contour(x1,y1,sigma_w);");
		fprintf(fp,"\nclabel(cs)");

		fprintf(fp,"\nfid1=fopen('1.4.Stress_Distribution.inp','w+');");
		fprintf(fp,"\n for i=1:%d,",JMAX-1);
		fprintf(fp,"\n for j=1:%d,",IMAX-1);
		fprintf(fp,"\nfprintf(fid1, '%%e %%e \\n',sigma_w(i,j),sai_w(i,j));");

		fprintf(fp,"\nend");
		fprintf(fp,"\nend");
		fprintf(fp,"\nfclose(fid1);");
	fclose(fp) ;
}


