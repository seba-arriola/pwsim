#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <malloc.h>


//pthread for mutual exclusion in fftw3 plans and writing P times
extern pthread_mutex_t mutex;
extern pthread_mutex_t mutex2;

//File for P waves times
extern FILE *fppp;

//Number of points in the finite fault
extern int N;

//number of layers in velocity model
extern int M;

//number of stations 
extern int V;

//number of subfaults with slip>0
extern int NR;

//calculate amplification factors due to the free surface
extern int calcfs;

//seismic moment and moment magnitude
extern double M0;
extern double Mw;

//Medium parameters
extern double alpha;
extern double beta;
extern double rho;




//Quality factors
extern double Qso,Qpo;
extern double Qexp;
extern int N_attpar;
//stress drop
extern double dsigma;

extern double *R_att;
extern double *R_att_aux_1;
extern double *p_att;

//hipocenter coordinates in meters, xhip=yhip=0.0
extern double xhip, yhip, zhip;

//corner frequency, according to Brune
extern double f0s, f0p;

//total time in seconds for waveforms
extern int waveform_time;

//Numbers of cores to use.
extern int n_cores;

extern int nfft;
extern int nnyq;
extern int Fs;
extern int T_S;

//Hypocenter geographical coordinates, used as reference 
extern double ALATO;
extern double ALNGO;

//parameters for filtering acceleration waveforms
extern int filter_acc;
extern double fmin_acc;
extern double fmax_acc;

//parameters for filtering velocity waveforms
extern int velocity;
extern int filter_vel;
extern double fmin_vel;
extern double fmax_vel;

//number of noise waveform to take an average
extern int N_simul;

//
extern int only_SH;

//apply transfer functions
extern int applyTF;
extern double rho_tf;
extern double am_p,am_sv,am_sh;
//input filenames
extern char finite_fault[512];
extern char vel_model[512];
extern char stations_file[512];
extern char envelope_file[512];
extern char attenuation_file[512];
extern double envelope_e;
extern double envelope_n;
extern double envelope_ft;

//Fault parameters
extern double *X;
extern double *Y;
extern double *Z;
extern double *slip;
extern double *strike;
extern double *dip;
extern double *rake;
extern double *length;
extern double *width;
extern double *Tr;
extern double *trup;
extern double *m0;
extern double mean_dip;
extern double mean_rake;
extern double mean_length;
extern double mean_width;
extern double sum_slip;

//corner frequencies for every subfault, P and S waves
extern double *fcs;
extern double *fcp;
extern double *f;

//Radiation pattern to use
extern int radpat;

// Average radiation patterns according to Onishi Horike
extern double Rpp_OH;
extern double Rpsv_OH;
extern double Rpsh_OH;

//global velocity model parameters
extern double *prof;
extern double *vp;
extern double *vs;

//stations parameters
extern int sample;
extern double *ruido;
extern int *TFstat;
extern char **sta_names;
extern char **sta_models;
extern double *sta_x;
extern double *sta_y;
extern double *sta_z;
extern double *kappa;
extern double *gamma_sta;

//seed for randoms
extern int seed;

//thread struct
typedef struct
{
	int starting;
	int ending;
	pthread_t pid;
} Args;

#define BUFFER_SIZE 1000

int isEmpty(const char *str);

int removeEmptyLines(char *path);

/* function that counts lines in a text file
 * Receives the name of the text file
 * */
int countLines(char *filename);


/* Gauss-Kruger's method to convert lat-lon to x-y
 * Receives two pointers x and y to store the x-y calculated coordinates  
 * and the longitude lon_x and latitude lat_y to be converted
 * */
void geoToXY(double *x, double *y, double lon_x, double lat_y);


/* Gauss-Kruger's method to convert x-y to lat-lon
 * Receives two pointers ALNGDG and ALATDG to store the lon-lat 
 * calculated coordinates and the longitude lon_x and latitude lat_y
 * to be converted
 * */
void xyToGEO(double *ALNGDG, double *ALATDG, double x, double y);


/* Function to generate additive white Gaussian Noise samples with zero mean and a standard deviation of 1. 
 * This code was developed by Cagri Tanriover and published in https://www.embeddedrelated.com/
*/
double AWGN_generator();


/* Initializes subfaults, global velocity model and stations parameters 
 * from files.
 * */
void readInputs();


void readParams(int argcasd, char *argvasd);



/* Function to remove linear trends.
 * It receives an array with its size and return a pointer to a new 
 * array with the linear trend removed.
 * */
double *detrend(double a[],int a_size);



/* Receives an integer and return the higher next pow of 2.
 * */
int nextpow2(int x);



/*////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
 * A standard implementation of mergesort algorithm and iterative
 * binary search. The three function was mainly taken from 
 * geeksforgeeks.org
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////*/

void merge(double arr[], int l, int m, int r);

 
void mergeSort(double arr[], int l, int r);


int binarySearch(double arr[], int l, int r, double x);

/*/////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////*/

/* Falta descripcion de la funcion*/
int NumberOfPulsingSubs(int i);

/* Function to calculate radiation patterns, seismic moment and corner
 * frequencies.
 * */
void calc();

/* Function to calculate the angle of incidence and the angle of coordinate's rotation 
 * Receives the station coordinates and subfault coordinates, and pointers to 
 * variables which stores the calculated angles.
 * */
void angles(double xst, double xso, double yst, double yso, double zst,  double zso, double *theta, double *phi, double *d_total,double *tp, double *ts, double *i_e );

/* Function to calcule the slip amplification from the arrival of Sv waves at surface.
 * Receives an incidence angle and two pointers for storing vertical and radial
 * amplification
 * */
void FS_SV(double theta, double *AR, double *AZ);

/* Function to calcule the slip amplification from the arrival of P waves at surface.
 * Receives an incidence angle and two pointers for storing vertical and radial
 * amplification
 * */
void FS_P(double theta, double *AR, double *AZ);


/* Function to calcule the Fourier transform of a tappered white noise window
 * Receives the number of points in the time window, the sample frequency,
 * distance to hypocenter and moment magnitude
 * */
void FN(int TS,int Fs,double R_hypo,double complex out[]);

/* Function to calcule the geometric and anelastic attenuation
 * */
double P(double f,double R,double Q0,double v);

/* 
 * Term related to filtering high frequency at/near surface
 * */
double D(double f, double kappa_sta);

/* Function for giving shape to the acceleration specter
 * */
void Ysp(double M0,double f[],double fc,double R,double C,double Q0,double v,int n_nyq, double kappa_sta,double Ysp_ret[], double Gamma);

/* Inverse Fourier transform
 * */
void ITF(double complex EF[],double Fs,int nfft,int nt,double H,double ifft[]);

/* Function that moves acceleration waveform in a time t0 adding padding zeroes.
 * Receives a time t0, an acceleration array, the sample frequency, the size
 * of the acceleration array, and the final size of the array with zeroes.
 * */
double *offset(double t0,double a[],double Fs,int a_size,int final_size);



/*////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
 * A standard implementation to find the inverse of a complex square matrix
 * using Gauss-Jordan elimination algorithm. 
 * The base of this code was taken from codesansar.com and modified for 
 * using arbitrary complex matrix size.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////*/
double complex** inverter(double complex**a, int n);
/*/////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////*/


/* Auxiliary function for P and SV transfer functions. It returns a 
 * 2x2 complex matrix.
 * */
double complex**halfspace(double cs,double cp,double ds,double dp,double rho_hs,double k,double w,double fact);

/* Auxiliary function for P and SV transfer functions. It returns a 
 * 4x4 complex matrix.
 * */
double complex **stiflayer(double k,double wc,double ccs,double ccp,double h,double rho_sti,double ds,double dp);


/* Auxiliary function for P and SV transfer functions. It returns a 
 * complex matrix which size depends on the local velocity model.
 * */
double complex **globalmatrix(double k, double w, double cs[],double cp[],double h[],int model_size, double rho_satf,double beta1, double beta2);

/* Smooth a single vector value depending on its index, using a moving average.
 * The index of the value and the order determinates the amount of neighbors
 * to use.
 * */
double complex smooth(double complex *A,int A_size, int ind, double complex *cum_sum, int a_order);

/* Calculates the P and SV transfer functions
 * */
double complex ** SATF_P_SV(double *freq,int freq_size,double *vs_sta,double *vp_sta,double *espesor_sta,int model_size,double am_p, double am_sv,double angle,double rho_satf);

/* Calculates the SH transfer function
 * */
double complex * SATF_SH(double *freq,int freq_size,double *vs_sta,double *espesor_sta,int model_size,double am_sh, double rho_satf);

/* Function to generate the acceleration and integrated velocity waveforms, which are
 * written to disk.
 * Receives and station name, a sample frequency, the station coordinates and the 
 * amount of points in the generated waveforms (this is function of the total time
 * and the Fs).
 * */
void generateWaveform(char *sta_name, char *sta_model, int stat, double x_est, double y_est, double z_est, double kappa_sta,double *R, double *Csv_r,double *Csv_z,double *Csh,double *Cp_r,double *Cp_z,double *phi,double *theta,double *dtotal,double *Rpp,double *Rpsv,double *Rpsh,double *i_e,double *tpa,double *tsa,double *offset_p,double *offset_s,double *ap,double *bp,double *as,double *bs,double complex *SFsv_r ,double complex *SFsh, double complex *SFsv_z, double complex *SFp_r, double complex *SFp_z,double complex *SN,double *Ysv_r,double *Ysh,double *Ysv_z,double *Yp_r,double *Yp_z,double *asv_r,double *asv_z,double *ash,double *ap_r,double *ap_z, double Gamma );


//thread function to calculate waveforms in different cores
void CalculateWaveform(int s, int e);

//thread function to calculate waveforms in different cores
void *generateWaveform_thread(void *ptr);

//interpolate transfers functions amplitudes provided manually by user
double complex ** interpolateFA(char *filename);

//Calculate ace or vel
void calacv(int usar_FS, int stat, double **total_ac_N,double **total_ac_E,double **total_ac_V,double *R,double x_est, double y_est, double z_est,double FS_SV_r, double FS_SV_z, double FS_SH_t, double FS_P_r,double FS_P_z, double *Csv_r,double *Rpsv,double *Rpsh,double *Rpp,double *Csv_z,double *Csh,double *Cp_r,double *Cp_z,double *offset_p,double *offset_s,double *tpa,double *tsa,double *min_offtp, double kappa_sta,double *Ysv_r,double *Ysh,double *Ysv_z,double *Yp_r,double *Yp_z,double complex *SFsv_r,double complex *SFsh, double complex *SFsv_z, double complex *SFp_r, double complex *SFp_z,double complex *SN,double complex *FT_SV,double complex *FT_P,double complex *FT_SH,double complex **FT_P_SV,double complex *FT_H,double complex *FT_V,double *ap,double *bp,double *as,double *bs,double *theta, double *asv_r,double *asv_z,double *ash,double *ap_r,double *ap_z,double *phi, int final_size, double Gamma, double R_aux_jb);
