#include "pwsim.h"

#ifdef __GLIBC__
#include <malloc.h>
#endif

/* Number of subfaults that pulse simultaneously with subfault `i`. */
int NumberOfPulsingSubs(int i)
{
    double pulsing_percent = 1.;
    double R_hypo2 = sqrt( pow((xhip - X[i]),2.) + pow((yhip - Y[i]),2.) + pow((zhip - Z[i]),2.) )/M_PER_KM;
    if ((mean_length+mean_width) == 0.0)
        return 1;
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
        r = sqrt( pow((xhip - X[k]),2.) + pow((yhip - Y[k]),2.) + pow((zhip - Z[k]),2.) )/M_PER_KM;
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
    double to_rad = DEG2RAD;

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


    f0s = BRUNE_CONST * beta * pow((dsigma/M0),(1./3.));
    f0p = BRUNE_CONST * alpha * pow((dsigma/M0),(1./3.));

    fcs =(double*)malloc(N*sizeof(double));
    fcp =(double*)malloc(N*sizeof(double));

    mergeSort(trup,0,N-1);

    int tpos;

    for(int k=0;k<N;k++)
    {
        tpos = NoOfActiveSubs[k];
        fcs[k] = BRUNE_CONST*beta*pow( (dsigma/(tpos*M0/N)),(1./3.) );
        fcp[k] = BRUNE_CONST*alpha*pow( (dsigma/(tpos*M0/N)),(1./3.) );
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

    free(NoOfActiveSubs);

    return;

}


/* Function to calculate the slip amplification from the arrival of Sv waves at surface.
 * Receives an incidence angle and two pointers for storing vertical and radial
 * amplification
 * */
void FS_SV(double theta, double *AR, double *AZ)
{
    double j0 = theta;
    double C = 1.0/sqrt(3.0);  // S/P velocity ratio for a Poisson solid
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


/* Function to calculate the slip amplification from the arrival of P waves at surface.
 * Receives an incidence angle and two pointers for storing vertical and radial
 * amplification
 * */
void FS_P(double theta, double *AR, double *AZ)
{
    double j0 = theta;
    double k = sqrt(3.); // P/S velocity ratio for a Poisson solid
    double i0 = asin(sin(j0)/k);

    double ar = (sin(2.*j0)*sin(2.*i0)-pow(k*cos(2.*i0),2.))/(sin(2.*j0)*sin(2.*i0)+pow(k*cos(2.*i0),2.));
    double br = 2.*k*sin(2.*j0) * cos(2.*i0)/(sin(2.*j0)*sin(2.*i0)+  pow(k*cos(2.*i0),2.)    );

    *AR = fabs(1 + ar + br*cos(i0)/sin(j0));

    *AZ = fabs(1 - ar + br*sin(i0)/cos(j0));
    return;
}





/* Function to calculate the geometric and anelastic attenuation
 * */
double P(double f,double R,double Q0,double v)
{
    double Z=1.0;
    double Q = Q0*pow((f+0.00001),Qexp); //input parameter

    int aux_ind = N_attpar;
    int j ;

    for(j=0;j<N_attpar;j++)
    {
        if(R>=R_att_aux_1[j] && R<R_att_aux_1[j+1])
        {
            aux_ind = j+1;
            break;
        }
    }


    double R_aux2[N_attpar+1];
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

    double result = Z*exp(-M_PI*f*R/(Q*v));

    return result;
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


int readTFfromVmodel(double **vs_sta, double **vp_sta,double **espesor_sta,const char *filename)
{

    int L;
    int C;
    L=countLines(filename);
    *espesor_sta=(double*)malloc(L*sizeof(double));
    *vp_sta=(double*)malloc(L*sizeof(double));
    *vs_sta=(double*)malloc(L*sizeof(double));
    double A1,A2, A3;

    FILE    *archivo = fopen(filename,"rt");
    if(archivo == NULL)
    {
        printf("Local velocity model file %s cannot be opened\n", filename);
        exit(EXIT_FAILURE);
    }
    C=0;

    //read local velocity model
    while(fscanf(archivo,"%lf  %lf  %lf", &A1,&A2,&A3)==3)
    {
        (*espesor_sta)[C]=A1;
        (*vp_sta)[C]=A2;
        (*vs_sta)[C]=A3;
        C++;
    }
    fclose(archivo);

    if(C == 0)
    {
        printf("Local velocity model file %s has no valid layers\n", filename);
        exit(EXIT_FAILURE);
    }
    return C;

}

/* Function to generate the acceleration and integrated velocity waveforms, which are
 * written to disk.
 * Receives and station name, a sample frequency, the station coordinates and the
 * amount of points in the generated waveforms (this is function of the total time
 * and the Fs).
 * */
void generateWaveform(const char *sta_name, const char *sta_model, WaveContext *w)
{

    int i;
    double R_distance;

    R_distance = sqrt( pow((xhip - w->x_est),2.) + pow((yhip - w->y_est),2.) + pow((zhip - w->z_est),2.) )/M_PER_KM;


    FN(T_S,Fs, R_distance,w->SN);


    //arrays for transfer functions
    double complex *FT_SV = NULL;
    double complex *FT_P = NULL;
    double complex *FT_SH = NULL;
    double complex **FT_P_SV = NULL;
    double *vs_sta = NULL, *vp_sta = NULL, *espesor_sta = NULL;

    double complex *FT_H = NULL;
    double complex *FT_V = NULL;


    double max_ts = 0.0;

    w->min_offtp = LARGE_TIME;
    double angle_tf=0.0;
    for(i=0;i<N;i++)
    {
        angles( w->x_est, X[i],  w->y_est, Y[i], w->z_est, Z[i], &w->theta[i],&w->phi[i],&w->dtotal[i],&w->tpa[i],&w->tsa[i],&w->i_e[i] );
        if(max_ts<w->tsa[i])
            max_ts=w->tsa[i];
        w->Rpp[i] = cos(rake[i])*sin(dip[i])*pow(sin(w->i_e[i]),2.)*sin(2.*(w->phi[i]-strike[i])) - cos(rake[i])*cos(dip[i])*sin(2.*w->i_e[i])*cos(w->phi[i]-strike[i]) + sin(rake[i])*sin(2.*dip[i])*(pow(cos(w->i_e[i]),2.)-pow(sin(w->i_e[i]),2.)*pow(sin(w->phi[i]-strike[i]),2.)) + sin(rake[i])*cos(2.*dip[i])*sin(2.*w->i_e[i])*sin(w->phi[i]-strike[i]);
        w->Rpsv[i] = sin(rake[i])*cos(2.*dip[i])*cos(2.*w->i_e[i])*sin(w->phi[i]-strike[i]) -cos(rake[i])*cos(dip[i])*cos(2.*w->i_e[i])*cos(w->phi[i]-strike[i]) + 1./2.*cos(rake[i])*sin(dip[i])*sin(2.*w->i_e[i])*sin(2.*(w->phi[i]-strike[i])) - 1./2.*sin(rake[i])*sin(2.*dip[i])*sin(2.*w->i_e[i])*(1.+ pow(sin(w->phi[i]-strike[i]),2.) );
        w->Rpsh[i] = cos(rake[i])*cos(dip[i])*cos(w->i_e[i])*sin(w->phi[i]-strike[i]) + cos(rake[i])*sin(dip[i])*sin(w->i_e[i])*cos(2.*(w->phi[i]-strike[i])) + sin(rake[i])*cos(2.*dip[i])*cos(w->i_e[i])*cos(w->phi[i]-strike[i]) - 1./2.*sin(rake[i])*sin(2.*dip[i])*sin(w->i_e[i])*sin(2.*(w->phi[i]-strike[i]));
        angle_tf=angle_tf+w->theta[i]/N;
    }

    //if you need to apply TF for some stations, these functions must be done
    switch(applyTF)
    {
        case 1:
            if(w->stat == 1)
            {
                int L=readTFfromVmodel(&vs_sta, &vp_sta,&espesor_sta,sta_model);
                double *f_tf = (double*)malloc(nnyq*sizeof(double));
                memcpy(f_tf, f, nnyq*sizeof(double));
                f_tf[0]=f_tf[1];
                FT_P_SV = SATF_P_SV(f_tf,nnyq,vs_sta,vp_sta,espesor_sta,L,damping_p,damping_sv,angle_tf,rho_tf);
                free(f_tf);
                FT_P  = FT_P_SV[0];
                FT_SV = FT_P_SV[1];
                free(FT_P_SV);
                FT_P_SV = NULL;

                FT_SH = SATF_SH(f,nnyq,vs_sta,espesor_sta,L,damping_sh,rho_tf);
            }
            else if(w->stat == 2)
            {
                double complex ** Amps = interpolateFA(sta_model);
                FT_H=Amps[0];
                FT_V=Amps[1];
                free(Amps);
            }
            break;
    }

    //For long distance stations, it could be necessary to calculate max_ts
    int final_size = (int)(max_ts/(1./Fs)) + T_S;
    w->final_size = final_size;

    double *total_ac_N;
    double *total_ac_E;
    double *total_ac_V;

    double *total_ac_N_nofs = NULL;
    double *total_ac_E_nofs = NULL;
    double *total_ac_V_nofs = NULL;

    w->FT_SV = FT_SV;
    w->FT_P  = FT_P;
    w->FT_SH = FT_SH;
    w->FT_H  = FT_H;
    w->FT_V  = FT_V;

    calacv(calcfs, &total_ac_N,&total_ac_E,&total_ac_V, w);

    if(calcfs == 2)
    {
        calacv(0, &total_ac_N_nofs,&total_ac_E_nofs,&total_ac_V_nofs, w);
    }


    pthread_mutex_lock(&mutex2);
    fprintf(fppp,"%s\t%.5f\n",sta_name,w->min_offtp);
    pthread_mutex_unlock(&mutex2);

    char output_acc[512];
    snprintf(output_acc,sizeof(output_acc),"output/synt_%s.dat",sta_name);


    FILE *fp = fopen(output_acc,"w");
    if(fp == NULL)
    {
        printf("Cannot write synthetic waveforms to file %s\n", output_acc);
        exit(EXIT_FAILURE);
    }


    for(i=0;i<final_size;i++)
    {
        if(calcfs == 2)
        {
            fprintf(fp,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",
                total_ac_N[i],total_ac_E[i],total_ac_V[i],
                total_ac_N_nofs[i],total_ac_E_nofs[i],total_ac_V_nofs[i]);
        }
        else
        {
            fprintf(fp,"%.8f\t%.8f\t%.8f\n",total_ac_N[i],total_ac_E[i],total_ac_V[i]);
        }
    }

    fclose(fp);




    free(total_ac_N);
    free(total_ac_E);
    free(total_ac_V);
    if(calcfs == 2)
    {
        free(total_ac_N_nofs);
        free(total_ac_E_nofs);
        free(total_ac_V_nofs);
    }


    if(applyTF==1)
    {
        if(w->stat==1)
        {
            free(vs_sta);
            free(vp_sta);
            free(espesor_sta);
            free(FT_P);
            free(FT_SV);
            free(FT_SH);
        }
        if(w->stat==2)
        {
            free(FT_V);
            free(FT_H);
        }
    }

#ifdef __GLIBC__
    malloc_trim(0);
#endif
    return;

}


//thread function to calculate waveforms in different cores
void CalculateWaveform(int s, int e)
{
    WaveContext w;

    w.R        = (double*)malloc(N*sizeof(double));
    w.Csv_r    = (double*)malloc(N*sizeof(double));
    w.Csv_z    = (double*)malloc(N*sizeof(double));
    w.Csh      = (double*)malloc(N*sizeof(double));
    w.Cp_r     = (double*)malloc(N*sizeof(double));
    w.Cp_z     = (double*)malloc(N*sizeof(double));
    w.phi      = (double*)malloc(N*sizeof(double));
    w.theta    = (double*)malloc(N*sizeof(double));
    w.dtotal   = (double*)malloc(N*sizeof(double));
    w.Rpp      = (double*)malloc(N*sizeof(double));
    w.Rpsv     = (double*)malloc(N*sizeof(double));
    w.Rpsh     = (double*)malloc(N*sizeof(double));
    w.i_e      = (double*)malloc(N*sizeof(double));
    w.tpa      = (double*)malloc(N*sizeof(double));
    w.tsa      = (double*)malloc(N*sizeof(double));
    w.offset_p = (double*)malloc(N*sizeof(double));
    w.offset_s = (double*)malloc(N*sizeof(double));
    w.ap       = (double*)malloc(N*sizeof(double));
    w.bp       = (double*)malloc(N*sizeof(double));
    w.as       = (double*)malloc(N*sizeof(double));
    w.bs       = (double*)malloc(N*sizeof(double));
    w.SFsv_r   = (double complex*)malloc(nfft*sizeof(double complex));
    w.SFsh     = (double complex*)malloc(nfft*sizeof(double complex));
    w.SFsv_z   = (double complex*)malloc(nfft*sizeof(double complex));
    w.SFp_r    = (double complex*)malloc(nfft*sizeof(double complex));
    w.SFp_z    = (double complex*)malloc(nfft*sizeof(double complex));
    w.SN       = (double complex*)malloc(sizeof(double complex) * nextpow2(T_S));
    w.Ysv_r    = (double*)malloc(nnyq*sizeof(double));
    w.Ysh      = (double*)malloc(nnyq*sizeof(double));
    w.Ysv_z    = (double*)malloc(nnyq*sizeof(double));
    w.Yp_r     = (double*)malloc(nnyq*sizeof(double));
    w.Yp_z     = (double*)malloc(nnyq*sizeof(double));
    w.asv_r    = (double*)malloc(T_S*sizeof(double));
    w.asv_z    = (double*)malloc(T_S*sizeof(double));
    w.ash      = (double*)malloc(T_S*sizeof(double));
    w.ap_r     = (double*)malloc(T_S*sizeof(double));
    w.ap_z     = (double*)malloc(T_S*sizeof(double));
    w.FT_SV    = NULL;
    w.FT_P     = NULL;
    w.FT_SH    = NULL;
    w.FT_H     = NULL;
    w.FT_V     = NULL;

    for(int k=s;k<=e;k++)
    {
        w.stat     = TFstat[k];
        w.x_est    = sta_x[k];
        w.y_est    = sta_y[k];
        w.z_est    = sta_z[k];
        w.kappa_sta= kappa[k];
        w.Gamma    = gamma_sta[k];

        generateWaveform(sta_names[k], sta_models[k], &w);
    }

    free(w.R);
    free(w.Csv_r);
    free(w.Csv_z);
    free(w.Csh);
    free(w.Cp_r);
    free(w.Cp_z);
    free(w.phi);
    free(w.theta);
    free(w.dtotal);
    free(w.Rpp);
    free(w.Rpsv);
    free(w.Rpsh);
    free(w.i_e);
    free(w.tpa);
    free(w.tsa);
    free(w.offset_p);
    free(w.offset_s);
    free(w.ap);
    free(w.bp);
    free(w.as);
    free(w.bs);
    free(w.SFsv_r);
    free(w.SFsh);
    free(w.SFsv_z);
    free(w.SFp_r);
    free(w.SFp_z);
    free(w.SN);
    free(w.Ysv_r);
    free(w.Ysh);
    free(w.Ysv_z);
    free(w.Yp_r);
    free(w.Yp_z);
    free(w.asv_r);
    free(w.asv_z);
    free(w.ash );
    free(w.ap_r);
    free(w.ap_z);

    return;

}
