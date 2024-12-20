#include "pwsim.h"


void calc_radpat(int wtd,double *Rpsv_factor,double *Rpsh_factor,double *Rpp_factor, int index,double *Rpsv,double *Rpsh,double *Rpp)
{			

			
	switch(wtd)
    {
    case 0: //radiation pattern is calculated for every subfault.
        *Rpsv_factor=Rpsv[index];
        *Rpsh_factor=Rpsh[index];
        *Rpp_factor=Rpp[index] ;
        break;

    case 1: //use average value
        *Rpsv_factor=Rpsv_OH;
        *Rpsh_factor=Rpsh_OH;
        *Rpp_factor=Rpp_OH;
        break;
    }

    return;
	
	
}


void TFapplied(int applyTF, int stat, int index, double complex *aSFsv_r, double complex *aSFsh, double complex *aSFsv_z, double complex *aSFp_r, double complex *aSFp_z,  double complex *aFT_SV,  double complex *aFT_SH,  double complex *aFT_P,  double complex *aFT_H,  double complex *aFT_V)
{
	double complex factor1=0.0;
	double complex factor2=0.0;
	double complex factor3=0.0;
	double complex factor4=0.0;
	double complex factor5=0.0;
	
	switch(applyTF)
	{
		case 0: //don't apply transfer functions
			return;
			
		case 1: //apply
			switch(stat)
			{
				case 1:
				factor1 =   creal(aFT_SV[index]) ;
				factor2 =   creal(aFT_SH[index]) ;
				factor3 =   creal(aFT_SV[index]) ;
				factor4 =   creal(aFT_P[index])  ;
				factor5 =   creal(aFT_P[index])  ;
				
				
				case 2:
				factor1 =   aFT_H[index] ;
				factor2 =   aFT_H[index] ;
				factor3 =   aFT_V[index] ;
				factor4 =   aFT_H[index] ;
				factor5 =   aFT_V[index] ;
				
			}
		
	}
	
	*aSFsv_r = *aSFsv_r*factor1 ;
	*aSFsh   = *aSFsh*factor2 ;
	*aSFsv_z = *aSFsv_z*factor3;
	*aSFp_r  = *aSFp_r*factor4 ;
	*aSFp_z  = *aSFp_z*factor5 ;
				
}

//Calculate ace or vel
void calacv(int usar_FS, int stat, double **total_ac_N,double **total_ac_E,double **total_ac_V,double *R,double x_est, double y_est, double z_est,double FS_SV_r, double FS_SV_z, double FS_SH_t, double FS_P_r,double FS_P_z, double *Csv_r,double *Rpsv,double *Rpsh,double *Rpp,double *Csv_z,double *Csh,double *Cp_r,double *Cp_z,double *offset_p,double *offset_s,double *tpa,double *tsa,double *min_offtp, double kappa_sta,double *Ysv_r,double *Ysh,double *Ysv_z,double *Yp_r,double *Yp_z,double complex *SFsv_r,double complex *SFsh, double complex *SFsv_z, double complex *SFp_r, double complex *SFp_z,double complex *SN,double complex *FT_SV,double complex *FT_P,double complex *FT_SH,double complex **FT_P_SV,double complex *FT_H,double complex *FT_V,double *ap,double *bp,double *as,double *bs,double *theta,double *asv_r,double *asv_z,double *ash,double *ap_r,double *ap_z, double *phi, int final_size, double Gamma, double R_aux_jb)
{
	
	//create variables for this local function
	double EPsvr, EPsh, EPsvz;
	double EPpr, EPpz;
	double Hp, Hs;
	double ar,ah,az;
    double aN,aE,aV;
	
	double AoVexp=2.0;

	double **asvr  = (double**)calloc(N,sizeof(double*));
	double **asvz  = (double**)calloc(N,sizeof(double*));
    double **ashh  = (double**)calloc(N,sizeof(double*));
    double **apr   = (double**)calloc(N,sizeof(double*));
    double **apz   = (double**)calloc(N,sizeof(double*));

	//memory allocation for acc 
	*total_ac_N=(double*)malloc(final_size*sizeof(double));
	*total_ac_E=(double*)malloc(final_size*sizeof(double));
	*total_ac_V=(double*)malloc(final_size*sizeof(double));
	
	
	int i, j;
	
	//initializes acc arrays to 0.0 
	for(j=0;j<final_size;j++)
	{
		(*total_ac_N)[j]=0.0;
		(*total_ac_E)[j]=0.0;
		(*total_ac_V)[j]=0.0;


	}

	double Rpsv_factor=0.0;
	double Rpsh_factor=0.0;
	double Rpp_factor=0.0;
	
	
	
	for(i=0;i<N;i++)
	{
		ap[i]=0.0;
		bp[i]=0.0;
		as[i]=0.0;
		bs[i]=0.0;
		
		//if moment is 0 skip this subfault
		if(m0[i]==0.0) continue;
		
		//use Joyner-Boore or cartesian distance
		R[i] = sqrt(  (X[i]-x_est)*(X[i]-x_est) + (Y[i]-y_est)*(Y[i]-y_est) + (Z[i]-z_est)*(Z[i]-z_est) )/1000.0;



		EPsvr=cos(theta[i]);
		EPsh=1.;
		EPsvz=sin(theta[i]);
		EPpr = -1.*sin(theta[i]);
		EPpz = cos(theta[i]);

		if(usar_FS>0)
		{
			FS_SV(theta[i],&FS_SV_r,&FS_SV_z);
			FS_SH_t=2.;
			FS_P(theta[i],&FS_P_r, &FS_P_z);
		}


		calc_radpat(radpat,&Rpsv_factor,&Rpsh_factor,&Rpp_factor,i,Rpsv,Rpsh,Rpp);
		
		Csv_r[i] = (Rpsv_factor* fabs(FS_SV_r) * EPsvr ) / (4.0 * M_PI * rho  * pow(beta,3.0) );
		Csv_z[i] = (Rpsv_factor * fabs(FS_SV_z) * EPsvz ) / (4.0 * M_PI * rho  * pow(beta,3.0) );
		Csh[i]   = (Rpsh_factor * fabs(FS_SH_t) * EPsh  ) / (4.0 * M_PI * rho  * pow(beta,3.0) );
		Cp_r[i]  = (Rpp_factor * fabs(FS_P_r) *EPpr)/(4.0 * M_PI * rho  * pow(alpha,3.0) );
		Cp_z[i]  = (Rpp_factor * fabs(FS_P_z) *EPpz)/(4.0 * M_PI * rho  * pow(alpha,3.0) );
		

		offset_p[i]=Tr[i]+tpa[i];
		offset_s[i]=Tr[i]+tsa[i];
		if(*min_offtp>offset_p[i])
			*min_offtp=offset_p[i];

		
		Ysp(m0[i],f,fcs[i],R[i],Csv_r[i],Qso,beta, nnyq, kappa_sta,Ysv_r,Gamma);
		Ysp(m0[i],f,fcs[i],R[i],Csh[i],  Qso,beta, nnyq, kappa_sta,Ysh,Gamma);
        Ysp(m0[i],f,fcs[i],R[i],Csv_z[i],Qso,beta, nnyq, kappa_sta,Ysv_z,Gamma);
		Ysp(m0[i],f,fcp[i],R[i],Cp_r[i], Qpo,alpha,nnyq, kappa_sta,Yp_r,Gamma);
        Ysp(m0[i],f,fcp[i],R[i],Cp_z[i], Qpo,alpha,nnyq, kappa_sta,Yp_z,Gamma);
		
		
		for(j=0;j<nnyq;j++)
		{

			ap[i]=ap[i]+pow((    pow(  f[j], AoVexp   )    /(1.+pow((f[j]/f0p),Gamma)))*D(f[j],kappa_sta),2.);
			bp[i]=bp[i]+pow((    pow(  f[j], AoVexp   )    /(1.+pow((f[j]/fcp[i]),Gamma)))*D(f[j],kappa_sta),2.);

			as[i]=as[i]+pow((    pow(  f[j] , AoVexp   )    /(1.+pow((f[j]/f0s),Gamma)))*D(f[j],kappa_sta),2.);
			bs[i]=bs[i]+pow((    pow(  f[j] , AoVexp   )    /(1.+pow((f[j]/fcs[i]),Gamma)))*D(f[j],kappa_sta),2.);


            SFsv_r[j] = SN[j] * Ysv_r[j];
            SFsv_z[j] = SN[j] * Ysv_z[j];
            SFsh[j]   = SN[j] * Ysh[j];
            SFp_r[j]  = SN[j] * Yp_r[j];
            SFp_z[j]  = SN[j] * Yp_z[j];
            
            TFapplied(applyTF,stat,j,&SFsv_r[j], &SFsh[j], &SFsv_z[j], &SFp_r[j], &SFp_z[j],  FT_SV,  FT_SH,  FT_P,  FT_H,  FT_V);			
			
			
			if(j<1) continue;
			
			SFsv_r[nfft-j] = conj( SFsv_r[j]);
			SFsv_z[nfft-j] = conj(SFsv_z[j]);
			SFsh[nfft-j]  =  conj(SFsh[j]);
			SFp_r[nfft-j] =  conj(SFp_r[j]);
			SFp_z[nfft-j] =  conj(SFp_z[j]);
			

		}

		// CG modifications in the Scaling Factor H, now is the same employed in EXSIM
		Hp = sqrt((N)*ap[i]/bp[i]);
		Hs = sqrt((N)*as[i]/bs[i]);

		ITF(SFsv_r,Fs,nfft,T_S,Hs,asv_r );
		ITF(SFsv_z,Fs,nfft,T_S,Hs,asv_z);
        ITF(SFsh,  Fs,nfft,T_S,Hs,ash);

        ITF(SFp_r,Fs,nfft,T_S,Hp,ap_r);
        ITF(SFp_z,Fs,nfft,T_S,Hp,ap_z);


        asvr[i]=offset(offset_s[i],asv_r,Fs,T_S,final_size);
		asvz[i]=offset(offset_s[i],asv_z,Fs,T_S,final_size);
		ashh[i]=offset(offset_s[i],ash,  Fs,T_S,final_size);

		apr[i]=offset(offset_p[i],ap_r,Fs,T_S,final_size);
        apz[i]=offset(offset_p[i],ap_z,Fs,T_S,final_size);


        for(j=0;j<final_size;j++)
        {

			ar = asvr[i][j] + apr[i][j];
			ah = ashh[i][j];
			az = asvz[i][j]+ apz[i][j];

			
			if(only_SH==1)
			{
				aN=ah;
				aE=0;
				aV=0;
				
			}
			else // if you need the 3-components (P, SV, SH in EW, NS, and UD)
			{
				aN = cos(phi[i])*ar - sin(phi[i])*ah;
				aE = sin(phi[i])*ar + cos(phi[i])*ah;
				aV = az;
			}
			
			



			// in m/s/s
			(*total_ac_N)[j]=(*total_ac_N)[j]+aN/100.0;
			(*total_ac_E)[j]=(*total_ac_E)[j]+aE/100.0;
			(*total_ac_V)[j]=(*total_ac_V)[j]+aV/100.0;


		}


	}


	for(i=0;i<N;i++)
	{
		free(asvr[i]);
		free(apr[i]);
		free(ashh[i]);
		free(asvz[i]);
		free(apz[i]);
	}
	free(asvr);
	free(asvz);
    free(ashh);
    free(apr);
    free(apz);

}
