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

    default: //unknown value: fall back to per-subfault pattern
        *Rpsv_factor=Rpsv[index];
        *Rpsh_factor=Rpsh[index];
        *Rpp_factor=Rpp[index] ;
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
            if(stat == 0)
            {
                factor1 = 1.0;
                factor2 = 1.0;
                factor3 = 1.0;
                factor4 = 1.0;
                factor5 = 1.0;
            }
            else if(stat == 1)
            {
                factor1 =   creal(aFT_SV[index]) ;
                factor2 =   creal(aFT_SH[index]) ;
                factor3 =   creal(aFT_SV[index]) ;
                factor4 =   creal(aFT_P[index])  ;
                factor5 =   creal(aFT_P[index])  ;
            }
            else if(stat == 2)
            {
                factor1 =   aFT_H[index] ;
                factor2 =   aFT_H[index] ;
                factor3 =   aFT_V[index] ;
                factor4 =   aFT_H[index] ;
                factor5 =   aFT_V[index] ;
            }
            break;

    }

    *aSFsv_r = *aSFsv_r*factor1 ;
    *aSFsh   = *aSFsh*factor2 ;
    *aSFsv_z = *aSFsv_z*factor3;
    *aSFp_r  = *aSFp_r*factor4 ;
    *aSFp_z  = *aSFp_z*factor5 ;

}

//Calculate acceleration waveforms for one station
void calacv(int usar_FS, double **total_ac_N,double **total_ac_E,double **total_ac_V, WaveContext *w)
{
    double FS_SV_r = 1.0;
    double FS_SV_z = 1.0;
    double FS_SH_t = 1.0;
    double FS_P_r  = 1.0;
    double FS_P_z  = 1.0;

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
    *total_ac_N=(double*)malloc(w->final_size*sizeof(double));
    *total_ac_E=(double*)malloc(w->final_size*sizeof(double));
    *total_ac_V=(double*)malloc(w->final_size*sizeof(double));


    int i, j;

    //initializes acc arrays to 0.0
    for(j=0;j<w->final_size;j++)
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
        w->ap[i]=0.0;
        w->bp[i]=0.0;
        w->as[i]=0.0;
        w->bs[i]=0.0;

        //if moment is 0 skip this subfault
        if(m0[i]==0.0) continue;

        //cartesian distance from subfault to station, in km
        w->R[i] = sqrt(  (X[i]-w->x_est)*(X[i]-w->x_est) + (Y[i]-w->y_est)*(Y[i]-w->y_est) + (Z[i]-w->z_est)*(Z[i]-w->z_est) )/M_PER_KM;



        EPsvr=cos(w->theta[i]);
        EPsh=1.;
        EPsvz=sin(w->theta[i]);
        EPpr = -1.*sin(w->theta[i]);
        EPpz = cos(w->theta[i]);

        if(usar_FS>0)
        {
            FS_SV(w->theta[i],&FS_SV_r,&FS_SV_z);
            FS_SH_t=2.;
            FS_P(w->theta[i],&FS_P_r, &FS_P_z);
        }


        calc_radpat(radpat,&Rpsv_factor,&Rpsh_factor,&Rpp_factor,i,w->Rpsv,w->Rpsh,w->Rpp);

        w->Csv_r[i] = (Rpsv_factor* fabs(FS_SV_r) * EPsvr ) / (4.0 * M_PI * rho  * pow(beta,3.0) );
        w->Csv_z[i] = (Rpsv_factor * fabs(FS_SV_z) * EPsvz ) / (4.0 * M_PI * rho  * pow(beta,3.0) );
        w->Csh[i]   = (Rpsh_factor * fabs(FS_SH_t) * EPsh  ) / (4.0 * M_PI * rho  * pow(beta,3.0) );
        w->Cp_r[i]  = (Rpp_factor * fabs(FS_P_r) *EPpr)/(4.0 * M_PI * rho  * pow(alpha,3.0) );
        w->Cp_z[i]  = (Rpp_factor * fabs(FS_P_z) *EPpz)/(4.0 * M_PI * rho  * pow(alpha,3.0) );


        w->offset_p[i]=Tr[i]+w->tpa[i];
        w->offset_s[i]=Tr[i]+w->tsa[i];
        if(w->min_offtp>w->offset_p[i])
            w->min_offtp=w->offset_p[i];


        Ysp(m0[i],f,fcs[i],w->R[i],w->Csv_r[i],Qso,beta, nnyq, w->kappa_sta,w->Ysv_r,w->Gamma);
        Ysp(m0[i],f,fcs[i],w->R[i],w->Csh[i],  Qso,beta, nnyq, w->kappa_sta,w->Ysh,w->Gamma);
        Ysp(m0[i],f,fcs[i],w->R[i],w->Csv_z[i],Qso,beta, nnyq, w->kappa_sta,w->Ysv_z,w->Gamma);
        Ysp(m0[i],f,fcp[i],w->R[i],w->Cp_r[i], Qpo,alpha,nnyq, w->kappa_sta,w->Yp_r,w->Gamma);
        Ysp(m0[i],f,fcp[i],w->R[i],w->Cp_z[i], Qpo,alpha,nnyq, w->kappa_sta,w->Yp_z,w->Gamma);


        for(j=0;j<nnyq;j++)
        {

            w->ap[i]=w->ap[i]+pow((    pow(  f[j], AoVexp   )    /(1.+pow((f[j]/f0p),w->Gamma)))*D(f[j],w->kappa_sta),2.);
            w->bp[i]=w->bp[i]+pow((    pow(  f[j], AoVexp   )    /(1.+pow((f[j]/fcp[i]),w->Gamma)))*D(f[j],w->kappa_sta),2.);

            w->as[i]=w->as[i]+pow((    pow(  f[j] , AoVexp   )    /(1.+pow((f[j]/f0s),w->Gamma)))*D(f[j],w->kappa_sta),2.);
            w->bs[i]=w->bs[i]+pow((    pow(  f[j] , AoVexp   )    /(1.+pow((f[j]/fcs[i]),w->Gamma)))*D(f[j],w->kappa_sta),2.);


            w->SFsv_r[j] = w->SN[j] * w->Ysv_r[j];
            w->SFsv_z[j] = w->SN[j] * w->Ysv_z[j];
            w->SFsh[j]   = w->SN[j] * w->Ysh[j];
            w->SFp_r[j]  = w->SN[j] * w->Yp_r[j];
            w->SFp_z[j]  = w->SN[j] * w->Yp_z[j];

            TFapplied(applyTF,w->stat,j,&w->SFsv_r[j], &w->SFsh[j], &w->SFsv_z[j], &w->SFp_r[j], &w->SFp_z[j],  w->FT_SV,  w->FT_SH,  w->FT_P,  w->FT_H,  w->FT_V);


            if(j<1) continue;

            w->SFsv_r[nfft-j] = conj( w->SFsv_r[j]);
            w->SFsv_z[nfft-j] = conj(w->SFsv_z[j]);
            w->SFsh[nfft-j]  =  conj(w->SFsh[j]);
            w->SFp_r[nfft-j] =  conj(w->SFp_r[j]);
            w->SFp_z[nfft-j] =  conj(w->SFp_z[j]);


        }

        // CG modifications in the Scaling Factor H, now is the same employed in EXSIM
        Hp = sqrt((N)*w->ap[i]/w->bp[i]);
        Hs = sqrt((N)*w->as[i]/w->bs[i]);

        ITF(w->SFsv_r,Fs,nfft,T_S,Hs,w->asv_r );
        ITF(w->SFsv_z,Fs,nfft,T_S,Hs,w->asv_z);
        ITF(w->SFsh,  Fs,nfft,T_S,Hs,w->ash);

        ITF(w->SFp_r,Fs,nfft,T_S,Hp,w->ap_r);
        ITF(w->SFp_z,Fs,nfft,T_S,Hp,w->ap_z);


        asvr[i]=offset(w->offset_s[i],w->asv_r,Fs,T_S,w->final_size);
        asvz[i]=offset(w->offset_s[i],w->asv_z,Fs,T_S,w->final_size);
        ashh[i]=offset(w->offset_s[i],w->ash,  Fs,T_S,w->final_size);

        apr[i]=offset(w->offset_p[i],w->ap_r,Fs,T_S,w->final_size);
        apz[i]=offset(w->offset_p[i],w->ap_z,Fs,T_S,w->final_size);


        for(j=0;j<w->final_size;j++)
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
                aN = cos(w->phi[i])*ar - sin(w->phi[i])*ah;
                aE = sin(w->phi[i])*ar + cos(w->phi[i])*ah;
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
