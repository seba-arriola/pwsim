#include "pwsim.h"

/* Falta descripcion de la funcion*/
int NumberOfPulsingSubs(int i)
{
    double pulsing_percent = 1.;
    double R_hypo2 = sqrt( pow((xhip - X[i]),2.) + pow((yhip - Y[i]),2.) + pow((zhip - Z[i]),2.) )/1000.;
    int Rmax = R_hypo2/((mean_length+mean_width)/2.); 
    int Rmin = Rmax - (pulsing_percent/2.)*N;

    if (Rmin<0)
    {
        Rmin=0;
    }
    int n=0;
    double r;
    for(int k=0;k<N;k++)
    {
        r = sqrt( pow((xhip - X[k]),2.) + pow((yhip - Y[k]),2.) + pow((zhip - Z[k]),2.) )/1000.;
        r = r/((mean_length+mean_width)/2);
        if (r>Rmin && r<=Rmax)
        {
            n=n+1;
        }
    }
    
    if (n==0)
        {
             n=1;
        }

    return n;
}


/* Function to calculate radiation patterns, seismic moment and corner
 * frequencies.
 * */
void calc()
{
	double to_rad = M_PI/180;

	Rpp_OH = sqrt(4./15.);
	Rpsv_OH = (1./2.)*sqrt(pow(sin(mean_rake*to_rad),2.)*( 14./15. +1./3.*pow(sin(2*mean_dip*to_rad),2.) ) + pow(cos(mean_rake*to_rad),2.)*(4./15.+  2./3.*pow(cos(mean_dip*to_rad),2.) ));
	Rpsh_OH = (1./2.)*sqrt(2./3.*pow(cos(mean_rake*to_rad),2.)*(1.+pow(sin(mean_dip*to_rad),2.)) + 1./3.*pow(sin(mean_rake*to_rad),2.)*(1.+pow(cos(2.*mean_dip*to_rad),2.)));
	M0 = pow(10. , (3./2.)*(Mw+10.7) );
	m0 = (double*)malloc(N*sizeof(double));
	
	int *NoOfActiveSubs = (int*)malloc(N*sizeof(int));

	for(int i=0;i<N;i++)
	{
		m0[i]=M0*slip[i]/sum_slip;
		NoOfActiveSubs[i]=NumberOfPulsingSubs(i);
	}


	f0s = 4.9 * pow(10.,6.) * beta * pow((dsigma/M0),(1./3.));
	f0p = 4.9 * pow(10.,6.) * alpha * pow((dsigma/M0),(1./3.));

	fcs =(double*)malloc(N*sizeof(double));
	fcp =(double*)malloc(N*sizeof(double));

	mergeSort(trup,0,N-1);

	int tpos;

	for(int k=0;k<N;k++)
	{
		//tpos = binarySearch(trup, 0, N-1, Tr[k])+1;
		tpos = NoOfActiveSubs[k];
		fcs[k] = 4.9*pow(10.,6.)*beta*pow( (dsigma/(tpos*M0/N)),(1./3.) );
		fcp[k] = 4.9*pow(10.,6.)*alpha*pow( (dsigma/(tpos*M0/N)),(1./3.) );
	}


	Fs = sample;
	T_S = sample*waveform_time;
	nfft = nextpow2(sample*waveform_time);
	nnyq = (nfft+2.)/2.;
	double fnyq = 1./(2.*(1./sample));


	f=(double*)malloc(nnyq*sizeof(double));

	for(int i=0;i<nnyq;i++)
	{
		f[i]=i*fnyq/(nnyq-1);
	}

	return;

}


/* Function to calcule the slip amplification from the arrival of Sv waves at surface.
 * Receives an incidence angle and two pointers for storing vertical and radial
 * amplification
 * */
void FS_SV(double theta, double *AR, double *AZ)
{
	double j0 = theta;
	double C = 1./sqrt(3.);
	double i0 = asin(sin(j0)/C);

	if (theta < asin(C))
	{
		double f1 = -( pow(cos(2.*j0),2.) - C*C*sin(2.*j0)*sin(2.*i0))/( pow(cos(2.*j0),2.) + C*C *sin(2.*j0)* sin(2.*i0));
		double f2 = 2*C * sin(2.*j0)*cos(2.*j0)/( pow(cos(2.*j0),2.) + C*C * sin(2.*j0)*sin(2.*i0));
		*AR = fabs(1. - f1 + f2*(sin(i0)/cos(j0)));
		*AZ = fabs(1. + f1 + f2*(cos(i0)/sin(j0)));
	}

	else if(theta==asin(C))
	{
		*AR=2.;
		*AZ=2.;
	}
	else
	{

		double ar = pow(cos(2.*j0),2.) * cos(j0);
		double br = sqrt( pow(sin(j0),2.) - C*C)*pow(sin(2.*j0),2.);
		double Rr = sqrt(ar*ar + br*br);

		double az = 2.*sqrt( pow(sin(j0),2.) - C*C )*sin(j0)*sin(2.*j0);
		double bz =  -(pow(cos(2*j0),2.));
		double Rz = sqrt(az*az + bz*bz);

		*AR = fabs((2. *  cos(2.*j0)/( (  pow(cos(2.*j0),4.) + 4.*( pow(sin(j0),2.) - C*C)*pow(sin(2.*j0),2.) * pow(sin(j0),2.))))* Rr/cos(j0));

		*AZ = fabs(2.*sqrt(pow(sin(j0),2.) - C*C)*sin(2.*j0))/( pow(cos(2.*j0),4.) + 4.*( pow(sin(j0),2.) - C*C)*pow(sin(2.*j0),2.)* pow(sin(j0),2.))*(Rz/sin(j0));
	}
	return;
}


/* Function to calcule the slip amplification from the arrival of P waves at surface.
 * Receives an incidence angle and two pointers for storing vertical and radial
 * amplification
 * */
void FS_P(double theta, double *AR, double *AZ)
{
	double j0 = theta;
	double k = sqrt(3.);
	double i0 = asin(sin(j0)/k);

	double ar = (sin(2.*j0)*sin(2.*i0)-pow(k*cos(2.*i0),2.))/(sin(2.*j0)*sin(2.*i0)+pow(k*cos(2.*i0),2.));
    double br = 2.*k*sin(2.*j0) * cos(2.*i0)/(sin(2.*j0)*sin(2.*i0)+  pow(k*cos(2.*i0),2.)    );

    *AR = fabs(1 + ar + br*cos(i0)/sin(j0));

    *AZ = fabs(1 - ar + br*sin(i0)/cos(j0));
	return;
}





/* Function to calcule the geometric and anelastic attenuation
 * */
double P(double f,double R,double Q0,double v)
{
	double Z=1.0;
	double Q = Q0*pow((f+0.00001),Qexp); //input parameter
	
	int aux_ind,j ; 
	
	for(j=0;j<=N_attpar;j++)
	{
		if(R>=R_att_aux_1[j] && R<R_att_aux_1[j+1])
		{
			aux_ind = j+1;
			break;	
		}
	}
	
	
	double *R_aux2=(double *)malloc( (aux_ind+1)*sizeof(double ));
	R_aux2[aux_ind]=R;
	R_aux2[0]=R_att[0];
	
	for(j=aux_ind;j>1;j--)
	{
		R_aux2[j-1] = R_att_aux_1[j-1];
		
	}
	
	for(j=aux_ind;j>0;j--)
	{
		Z = Z*pow(( R_aux2[j-1]   / R_aux2[j]     ),p_att[j-1]   );
	}
		

	return Z*exp(-M_PI*f*R/(Q*v));
	
	free(R_aux2);
}

/*
 * Term related to filtering high frequency at/near surface
 * */
double D(double f, double kappa_sta)
{
	return exp(-M_PI*kappa_sta*f);
}

/* Function for giving shape to the acceleration specter
 * */
void Ysp(double M0,double f[],double fc,double R,double C,double Q0,double v,int n_nyq, double kappa_sta,double Ysp_ret[], double Gamma)
{
	for(int i=0;i<n_nyq;i++)
	{
		Ysp_ret[i]=pow(10.0,-20.0)*C*M0*(1./(1. + pow((f[i]/fc),(Gamma))   )) * pow((2.*M_PI*f[i]),2.0)*P(f[i],R,Q0,v)*D(f[i],kappa_sta);
	}

	return;

}


int readTFfromVmodel(double **vs_sta, double **vp_sta,double **espesor_sta,char *filename)
{
	
	int L;
	int C,lee;
	L=countLines(filename);
	*espesor_sta=(double*)malloc(L*sizeof(double));
	*vp_sta=(double*)malloc(L*sizeof(double));
	*vs_sta=(double*)malloc(L*sizeof(double));
	double A1,A2, A3;

	FILE	*archivo = fopen(filename,"rt");
	if(archivo == NULL)
	{
		printf("Local velocity model file %s cannot be opened\n", filename);
		exit(EXIT_FAILURE);
	}
	C=0;

	//read local velocity model
	do
	{
		fscanf(archivo,"%lf  %lf  %lf", &A1,&A2,&A3);
		lee=feof(archivo);
		
		if (lee==1) 
			break;  //stops if cannot read
			
		(*espesor_sta)[C]=A1;
		(*vp_sta)[C]=A2;
		(*vs_sta)[C]=A3;
		C++;
	}while(1);
	fclose(archivo); 
	return L;
	
}

/* Function to generate the acceleration and integrated velocity waveforms, which are
 * written to disk.
 * Receives and station name, a sample frequency, the station coordinates and the
 * amount of points in the generated waveforms (this is function of the total time
 * and the Fs).
 * */
void generateWaveform(char *sta_name, char *sta_model, int stat, double x_est, double y_est, double z_est, double kappa_sta,double *R, double *Csv_r,double *Csv_z,double *Csh,double *Cp_r,double *Cp_z,double *phi,double *theta,double *dtotal,double *Rpp,double *Rpsv,double *Rpsh,double *i_e,double *tpa,double *tsa,double *offset_p,double *offset_s,double *ap,double *bp,double *as,double *bs,double complex *SFsv_r ,double complex *SFsh, double complex *SFsv_z, double complex *SFp_r, double complex *SFp_z,double complex *SN,double *Ysv_r,double *Ysh,double *Ysv_z,double *Yp_r,double *Yp_z,double *asv_r,double *asv_z,double *ash,double *ap_r,double *ap_z, double Gamma )
{

	int i;
	double R_distance;

	R_distance = sqrt( pow((xhip - x_est),2.) + pow((yhip - y_est),2.) + pow((zhip - z_est),2.) )/1000.;


	FN(T_S,Fs, R_distance,SN);


	//arrays for transfer functions
	double complex *FT_SV;
	double complex *FT_P;
	double complex *FT_SH;
	double complex **FT_P_SV;
	double *vs_sta, *vp_sta,*espesor_sta;

	double complex *FT_H;
	double complex *FT_V;


    double max_ts = 0.0;

    double min_offtp = 999999.0;
    double angle_tf=0.0;
    for(i=0;i<N;i++)
	{
		angles( x_est, X[i],  y_est, Y[i], z_est, Z[i], &theta[i],&phi[i],&dtotal[i],&tpa[i],&tsa[i],&i_e[i] );
		//printf("%.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %d \t %s\n",theta[i],tpa[i],tsa[i],dtotal[i],i_e[i],i,sta_name);
		if(max_ts<tsa[i])
			max_ts=tsa[i];
		Rpp[i] = cos(rake[i])*sin(dip[i])*pow(sin(i_e[i]),2.)*sin(2.*(phi[i]-strike[i])) - cos(rake[i])*cos(dip[i])*sin(2.*i_e[i])*cos(phi[i]-strike[i]) + sin(rake[i])*sin(2.*dip[i])*(pow(cos(i_e[i]),2.)-pow(sin(i_e[i]),2.)*pow(sin(phi[i]-strike[i]),2.)) + sin(rake[i])*cos(2.*dip[i])*sin(2.*i_e[i])*sin(phi[i]-strike[i]);
		Rpsv[i] = sin(rake[i])*cos(2.*dip[i])*cos(2.*i_e[i])*sin(phi[i]-strike[i]) -cos(rake[i])*cos(dip[i])*cos(2.*i_e[i])*cos(phi[i]-strike[i]) + 1./2.*cos(rake[i])*sin(dip[i])*sin(2.*i_e[i])*sin(2.*(phi[i]-strike[i])) - 1./2.*sin(rake[i])*sin(2.*dip[i])*sin(2.*i_e[i])*(1.+ pow(sin(phi[i]-strike[i]),2.) );
		Rpsh[i] = cos(rake[i])*cos(dip[i])*cos(i_e[i])*sin(phi[i]-strike[i]) + cos(rake[i])*sin(dip[i])*sin(i_e[i])*cos(2.*(phi[i]-strike[i])) + sin(rake[i])*cos(2.*dip[i])*cos(i_e[i])*cos(phi[i]-strike[i]) - 1./2.*sin(rake[i])*sin(2.*dip[i])*sin(i_e[i])*sin(2.*(phi[i]-strike[i]));
		angle_tf=angle_tf+theta[i]/N;
	}

	//if you need to apply TF for some stations, this functions must be d
	switch(applyTF)
	{
		case 1:
			switch(stat)
				case 1:
				{
					int L=readTFfromVmodel(&vs_sta, &vp_sta,&espesor_sta,sta_model);
					f[0]=f[1];
					FT_P_SV = SATF_P_SV(f,nnyq,vs_sta,vp_sta,espesor_sta,L,am_p,am_sv,angle_tf,rho_tf);
					f[0]=0.0;
					FT_P  = FT_P_SV[0];
					FT_SV = FT_P_SV[1];

					FT_SH = SATF_SH(f,nnyq,vs_sta,espesor_sta,L,am_sh,rho_tf);
				}
				case 2:
				{
					double complex ** Amps = interpolateFA(sta_model);
					FT_H=Amps[0];
					FT_V=Amps[1];
				}
	}

	//For long distance stations, it could be necessary to calculate max_ts
	int final_size = (int)(max_ts/(1./Fs)) + T_S;

    //final accelerations
    double FS_SV_r=1.0;
	double FS_SV_z=1.0;
	double FS_SH_t=1.0;
	double FS_P_r=1.0;
	double FS_P_z=1.0;

	double *total_ac_N;
	double *total_ac_E;
	double *total_ac_V;


	calacv(calcfs,stat, &total_ac_N,&total_ac_E,&total_ac_V, R , x_est, y_est, z_est,FS_SV_r, FS_SV_z, FS_SH_t, FS_P_r, FS_P_z, Csv_r,Rpsv,Rpsh,Rpp,Csv_z,Csh,Cp_r,Cp_z,offset_p,offset_s,tpa,tsa,&min_offtp,kappa_sta,Ysv_r,Ysh,Ysv_z,Yp_r,Yp_z,SFsv_r,SFsh,SFsv_z, SFp_r, SFp_z,SN,FT_SV,FT_P,FT_SH,FT_P_SV,FT_H,  FT_V, ap, bp, as, bs, theta,   asv_r, asv_z, ash, ap_r, ap_z, phi, final_size,Gamma, R_distance);



	pthread_mutex_lock(&mutex2);
	fprintf(fppp,"%s\t%.5f\n",sta_name,min_offtp);
	pthread_mutex_unlock(&mutex2);

	char output_acc[50];
	sprintf(output_acc,"output/synt_%s.dat",sta_name);


	FILE *fp = fopen(output_acc,"w");
	if(fp == NULL)
	{
		printf("Cannot write synthetic waveforms to file %s\n", output_acc);
		exit(EXIT_FAILURE);
	}

	
	for(i=0;i<final_size;i++)
	{
		fprintf(fp,"%.8f\t%.8f\t%.8f\n",total_ac_N[i],total_ac_E[i],total_ac_V[i]);
	}

	fclose(fp);



	


	free(total_ac_N);
	free(total_ac_E);
	free(total_ac_V);


	if(applyTF==1)
	{
		if(stat==1)
		{
			free(vs_sta);
			free(vp_sta);
			free(espesor_sta);
			free(FT_P);
			free(FT_SV);
			free(FT_SH);
		}
		if(stat==2)
		{
			free(FT_V);
			free(FT_H);
		}
	}

	malloc_trim(0);
	return;

}


//thread function to calculate waveforms in different cores
void CalculateWaveform(int s, int e)
{
	double *R 	   = (double*)malloc(N*sizeof(double));
	double *Csv_r  = (double*)malloc(N*sizeof(double));
	double *Csv_z  = (double*)malloc(N*sizeof(double));
	double *Csh    = (double*)malloc(N*sizeof(double));
	double *Cp_r   = (double*)malloc(N*sizeof(double));
	double *Cp_z   = (double*)malloc(N*sizeof(double));
	double *phi    = (double*)malloc(N*sizeof(double));
	double *theta  = (double*)malloc(N*sizeof(double));
	double *dtotal = (double*)malloc(N*sizeof(double));
	double *Rpp    = (double*)malloc(N*sizeof(double));
	double *Rpsv   = (double*)malloc(N*sizeof(double));
	double *Rpsh   = (double*)malloc(N*sizeof(double));
	double *i_e    = (double*)malloc(N*sizeof(double));
	double *tpa	   =(double*)malloc(N*sizeof(double));
    double *tsa    =(double*)malloc(N*sizeof(double));
	double *offset_p = (double*)malloc(N*sizeof(double));
	double *offset_s = (double*)malloc(N*sizeof(double));
	double *ap       =(double*)malloc(N*sizeof(double));
	double *bp       =(double*)malloc(N*sizeof(double));
	double *as       =(double*)malloc(N*sizeof(double));
	double *bs       =(double*)malloc(N*sizeof(double));
	double complex *SFsv_r = (double complex*)malloc(nfft*sizeof(double complex));
	double complex *SFsh   = (double complex*)malloc(nfft*sizeof(double complex));
	double complex *SFsv_z = (double complex*)malloc(nfft*sizeof(double complex));
	double complex *SFp_r  = (double complex*)malloc(nfft*sizeof(double complex));
	double complex *SFp_z  = (double complex*)malloc(nfft*sizeof(double complex));
	double complex *SN = (double complex *) malloc(sizeof(double complex ) * nextpow2(T_S));
	double *Ysv_r	= (double*)malloc(nnyq*sizeof(double));
	double *Ysh		= (double*)malloc(nnyq*sizeof(double));
	double *Ysv_z	= (double*)malloc(nnyq*sizeof(double));
	double *Yp_r	= (double*)malloc(nnyq*sizeof(double));
	double *Yp_z	= (double*)malloc(nnyq*sizeof(double));
	double *asv_r = (double*)malloc(T_S*sizeof(double));
	double *asv_z = (double*)malloc(T_S*sizeof(double));
	double *ash   = (double*)malloc(T_S*sizeof(double));
	double *ap_r  = (double*)malloc(T_S*sizeof(double));
	double *ap_z  = (double*)malloc(T_S*sizeof(double));

	for(int k=s;k<=e;k++)
	{
		generateWaveform(sta_names[k], sta_models[k], TFstat[k], sta_x[k], sta_y[k], sta_z[k],kappa[k],R,Csv_r,Csv_z,Csh ,Cp_r ,Cp_z ,phi ,theta ,dtotal ,Rpp ,Rpsv ,Rpsh ,i_e ,tpa ,tsa ,offset_p ,offset_s ,ap ,bp ,as ,bs,SFsv_r,SFsh,SFsv_z,SFp_r,SFp_z,SN, Ysv_r, Ysh, Ysv_z, Yp_r, Yp_z, asv_r, asv_z, ash, ap_r, ap_z, gamma_sta[k]  );
	}

	free(R);
	free(Csv_r);
	free(Csv_z);
	free(Csh);
	free(Cp_r);
	free(Cp_z);
	free(phi);
	free(theta);
	free(dtotal);
	free(Rpp);
	free(Rpsv);
	free(Rpsh);
	free(i_e);
	free(tpa);
	free(tsa);
	free(offset_p);
	free(offset_s);
	free(ap);
	free(bp);
	free(as);
	free(bs);
	free(SFsv_r);
	free(SFsh);
	free(SFsv_z);
	free(SFp_r);
	free(SFp_z);
	free(SN);
	free(Ysv_r);
	free(Ysh);
	free(Ysv_z);
	free(Yp_r);
	free(Yp_z);
	free(asv_r);
	free(asv_z);
	free(ash );
	free(ap_r);
	free(ap_z);

	return;

}
