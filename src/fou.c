#include "pwsim.h"

/* Function to calcule the Fourier transform of a tappered white noise window
 * Receives the number of points in the time window, the sample frequency,
 * distance to hypocenter and moment magnitude
 * */
void FN(int TS,int Fs,double R_hypo,double complex out[])
{

	double dt = 1./Fs;

	double *t = (double*)malloc(TS*sizeof(double));
	double *taper = (double*)malloc(TS*sizeof(double));
	double *s = (double*)malloc(TS*sizeof(double));
	
	double e, n, ft;
	
	e = envelope_e;
	n = envelope_n;
	ft = envelope_ft;
	


	double b = -(e*log(n))/(1+e*(log(e)-1));
	double c = b/e;
	double a = pow((exp(1.)/e),b);

	double Tgm;

	//duration according to  Ghofrani (2013) and Yoshi (2013)
	if(Mw >= 8.9)
		Tgm = 97.67 + 0.128*R_hypo;
	else if(Mw >= 7.5 && Mw < 8.9)
		Tgm = 0.0015*pow(10.,(0.5*Mw))+0.51*pow(R_hypo,0.3);
	else
		Tgm = 0.0015*pow(10.,(0.5*Mw))+0.02*pow(R_hypo,1.04);


	int i;
	double ax;
	double r=0.05;

	int nfft = nextpow2(TS);

	double complex *in;
	fftw_plan pt;

	in = (double complex *) malloc(sizeof(double complex) * nfft);


	for(i=0;i<TS;i++)
	{
		t[i]=i*dt;

		//tukey window, tapper
		ax=i/(TS-1.);
		if(ax<r/2.)
			taper[i]= 0.5*( 1. + cos(  (2.*M_PI/r)*(ax-r/2.)    )  );
		else if(ax>=r/2. && ax<1.-r/2)
			taper[i]=1.;
		else
			taper[i]= 0.5*( 1. + cos(  (2.*M_PI/r)*(ax-1.+r/2.)    )  );

		//envol
		s[i] = ruido[i]*taper[i]*( a*pow((t[i]/(ft*Tgm)),b) * exp(-c*(t[i]/(ft*Tgm))) );
		in[i]=s[i];
	}
	pthread_mutex_lock(&mutex);
	pt = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	pthread_mutex_unlock(&mutex);
	fftw_execute(pt);



	double factor=0.0;


	for(i=0;i<nfft;i++)
	{
		out[i]=out[i]/Fs;
		factor=factor+ cabs(out[i])*cabs(out[i])/nfft;
	}

	factor=sqrt(factor);

	for(i=0;i<nfft;i++)
	{
		out[i]=out[i]/factor;
	}

	free(t);
	free(taper);
	free(s);
	free(in);
	pthread_mutex_lock(&mutex);
	fftw_destroy_plan(pt);
	pthread_mutex_unlock(&mutex);

	return;

}

/* Inverse Fourier transform
 * */
void ITF(double complex EF[],double Fs,int nfft,int nt,double H,double ifft[])
{
	fftw_plan pi;

	double complex *final_out=(double complex *)malloc(nfft*sizeof(double complex ));
	pthread_mutex_lock(&mutex);
	pi = fftw_plan_dft_1d(nfft, EF, final_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	pthread_mutex_unlock(&mutex);
	fftw_execute(pi);

	for(int i=0;i<nt;i++)
		ifft[i]=creal(final_out[i])*Fs*H/nfft;


	free(final_out);
	pthread_mutex_lock(&mutex);
	fftw_destroy_plan(pi);
	pthread_mutex_unlock(&mutex);
	return;
}
